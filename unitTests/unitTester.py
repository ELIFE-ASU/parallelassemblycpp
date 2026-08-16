#!/usr/bin/env python3
"""Build AssemblyCpp and run its regression manifest."""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "AssemblyCpp"
DEFAULT_MANIFEST = TEST_DIRECTORY / "regression_cases.tsv"
DEFAULT_PATHWAY_MANIFEST = TEST_DIRECTORY / "pathway_cases.tsv"
MANIFEST_HEADER = ("molecule", "expected_assembly_index")
PATHWAY_MANIFEST_HEADER = ("molecule", "expected_pathway")
PATHWAY_KEYS = {"file_graph", "remnant", "duplicates", "removed_edges"}
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")


class TestConfigurationError(RuntimeError):
    """Raised when the manifest, a fixture, or a build setting is invalid."""


@dataclass(frozen=True)
class TestCase:
    name: str
    source: Path
    expected: int
    expected_pathway: Path | None = None


@dataclass(frozen=True)
class TestResult:
    case: TestCase
    actual: int | None
    duration_seconds: float
    failure: str | None = None
    error: str | None = None

    @property
    def status(self) -> str:
        if self.error is not None:
            return "ERROR"
        if self.failure is not None or self.actual != self.case.expected:
            return "FAIL"
        return "PASS"


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def resolve_test_path(path: Path) -> Path:
    if path.is_absolute() or path.exists():
        return path.resolve()

    test_relative = TEST_DIRECTORY / path
    if test_relative.exists():
        return test_relative.resolve()

    return path.resolve()


def resolve_fixture(name: str, fixture_directory: Path) -> Path:
    requested = Path(name)
    candidate = requested if requested.is_absolute() else fixture_directory / requested

    if candidate.is_file():
        return candidate.resolve()

    mol_candidate = candidate.with_suffix(".mol")
    if mol_candidate.is_file():
        return mol_candidate.resolve()

    raise TestConfigurationError(
        f"fixture {name!r} was not found as {candidate} or {mol_candidate}"
    )


def load_manifest(path: Path) -> tuple[Path, list[TestCase]]:
    manifest = resolve_test_path(path)
    seen: dict[str, int] = {}
    cases: list[TestCase] = []

    try:
        with manifest.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != MANIFEST_HEADER:
                raise TestConfigurationError(
                    f"invalid header in {manifest}: expected {MANIFEST_HEADER}, got {header}"
                )

            for row in reader:
                line_number = reader.line_num
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != 2:
                    raise TestConfigurationError(
                        f"invalid row in {manifest}:{line_number}: expected 2 columns"
                    )

                name, expected_text = (value.strip() for value in row)
                if not name:
                    raise TestConfigurationError(
                        f"empty fixture name in {manifest}:{line_number}"
                    )
                if name in seen:
                    raise TestConfigurationError(
                        f"duplicate fixture {name!r} in {manifest}:{line_number}; "
                        f"first declared on line {seen[name]}"
                    )

                try:
                    expected = int(expected_text)
                except ValueError as error:
                    raise TestConfigurationError(
                        f"invalid assembly index in {manifest}:{line_number}: "
                        f"{expected_text!r}"
                    ) from error

                seen[name] = line_number
                cases.append(
                    TestCase(
                        name=name,
                        source=resolve_fixture(name, manifest.parent),
                        expected=expected,
                    )
                )
    except OSError as error:
        raise TestConfigurationError(f"cannot read {manifest}: {error}") from error

    if not cases:
        raise TestConfigurationError(f"manifest contains no test cases: {manifest}")

    return manifest, cases


def parse_pathway_document(path: Path) -> dict[str, object]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error

    if not isinstance(document, dict):
        raise ValueError(f"pathway document must be a JSON object: {path}")
    if set(document) != PATHWAY_KEYS:
        raise ValueError(
            f"invalid pathway keys in {path}: expected {sorted(PATHWAY_KEYS)}, "
            f"got {sorted(document)}"
        )
    for key in PATHWAY_KEYS:
        if not isinstance(document[key], list):
            raise ValueError(f"pathway field {key!r} must be an array in {path}")

    return document


def load_pathway_manifest(
    path: Path, cases: Sequence[TestCase]
) -> tuple[Path, list[TestCase]]:
    manifest = resolve_test_path(path)
    case_names = {case.name for case in cases}
    pathways: dict[str, Path] = {}

    try:
        with manifest.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != PATHWAY_MANIFEST_HEADER:
                raise TestConfigurationError(
                    f"invalid header in {manifest}: expected "
                    f"{PATHWAY_MANIFEST_HEADER}, got {header}"
                )

            for row in reader:
                line_number = reader.line_num
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != 2:
                    raise TestConfigurationError(
                        f"invalid row in {manifest}:{line_number}: expected 2 columns"
                    )

                name, relative_path = (value.strip() for value in row)
                if name in pathways:
                    raise TestConfigurationError(
                        f"duplicate pathway case {name!r} in "
                        f"{manifest}:{line_number}"
                    )
                if name not in case_names:
                    raise TestConfigurationError(
                        f"pathway case {name!r} is not in the regression manifest"
                    )

                expected_pathway = (manifest.parent / relative_path).resolve()
                try:
                    parse_pathway_document(expected_pathway)
                except ValueError as error:
                    raise TestConfigurationError(str(error)) from error
                pathways[name] = expected_pathway
    except OSError as error:
        raise TestConfigurationError(f"cannot read {manifest}: {error}") from error

    if not pathways:
        raise TestConfigurationError(
            f"pathway manifest contains no test cases: {manifest}"
        )

    return manifest, [
        TestCase(
            name=case.name,
            source=case.source,
            expected=case.expected,
            expected_pathway=pathways.get(case.name),
        )
        for case in cases
    ]


def audit_test_data(manifest: Path, cases: Sequence[TestCase], verbose: bool) -> None:
    fixture_directory = manifest.parent
    mol_fixtures = {path.resolve() for path in fixture_directory.glob("*.mol")}
    referenced_mol = {
        case.source for case in cases if case.source.suffix.lower() == ".mol"
    }
    fixture_only = sorted(mol_fixtures - referenced_mol)

    cases_by_hash: dict[str, list[TestCase]] = defaultdict(list)
    for case in cases:
        digest = hashlib.sha256(case.source.read_bytes()).hexdigest()
        cases_by_hash[digest].append(case)

    shared_content = [group for group in cases_by_hash.values() if len(group) > 1]
    conflicting_content = [
        group for group in shared_content if len({case.expected for case in group}) > 1
    ]
    if conflicting_content:
        details = "; ".join(
            ", ".join(f"{case.name}={case.expected}" for case in group)
            for group in conflicting_content
        )
        raise TestConfigurationError(
            f"byte-identical fixtures have conflicting expectations: {details}"
        )

    graph_cases = sum(case.source.suffix.lower() != ".mol" for case in cases)
    print(f"Manifest: {manifest}")
    print(f"Regression cases: {len(cases)} ({len(cases) - graph_cases} mol, {graph_cases} graph)")
    print(f"Molecule fixtures: {len(mol_fixtures)}")
    print(f"Fixture-only molecules: {len(fixture_only)}")
    print(f"Shared-content case groups: {len(shared_content)} (consistent expectations)")
    print(f"Pathway golden cases: {sum(case.expected_pathway is not None for case in cases)}")

    if verbose and fixture_only:
        print("Fixture-only molecule names:")
        for path in fixture_only:
            print(f"  {path.stem}")


def resolve_executable(path: Path) -> Path:
    if path.is_file():
        return path.resolve()

    located = shutil.which(str(path))
    if located is not None:
        return Path(located).resolve()

    raise TestConfigurationError(
        f"executable not found: {path}. Run again with --build or compile it first."
    )


def run_cli_command(
    executable: Path, arguments: Sequence[str], working_directory: Path
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            [str(executable), *arguments],
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise TestConfigurationError(
            f"CLI check timed out after {error.timeout:g} seconds: "
            f"{shlex.join([str(executable), *arguments])}"
        ) from error
    except OSError as error:
        raise TestConfigurationError(f"cannot run CLI check: {error}") from error


def require_cli(
    condition: bool,
    message: str,
    completed: subprocess.CompletedProcess[str] | None = None,
) -> None:
    if condition:
        return

    diagnostics = ""
    if completed is not None:
        diagnostics = format_process_diagnostics(completed)
    suffix = f"; {diagnostics}" if diagnostics else ""
    raise TestConfigurationError(f"CLI check failed: {message}{suffix}")


def read_first_line_assembly_index(path: Path) -> int | None:
    """Return an index only when the first output line is fully numeric."""
    if not path.is_file():
        return None
    lines = path.read_text().splitlines()
    if not lines:
        return None
    match = ASSEMBLY_INDEX_PATTERN.search(lines[0])
    if match is None or match.end() != len(lines[0]):
        return None
    return int(match.group(1))


def read_last_intermediate_index(path: Path) -> int | None:
    """Return the assembly index in the final well-formed intermediate row."""
    if not path.is_file():
        return None
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if not rows or len(rows[-1]) != 2 or not rows[-1][1].lstrip("-").isdigit():
        return None
    return int(rows[-1][1])


def run_cli_checks(executable: Path) -> int:
    """Exercise help, validation, aliases, input handling, and output flags."""
    scenarios = 0
    help_tokens = (
        "Usage:",
        "--runtime=<TICKS>",
        "--enum-max=<COUNT>",
        "--pathway=<0|1>",
        "--remove-hydrogens=<0|1>",
        "--compensate-disjoint=<0|1>",
        "--memory-report=<0|1>",
        "--write-intermediate-mas=<0|1>",
        "Outputs:",
        "Compatibility:",
    )

    with tempfile.TemporaryDirectory(prefix="assemblycpp-cli-") as directory:
        working_directory = Path(directory)

        for help_option in ("--help", "-h"):
            completed = run_cli_command(
                executable, [help_option], working_directory
            )
            require_cli(
                completed.returncode == 0,
                f"{help_option} should exit successfully",
                completed,
            )
            missing_tokens = [
                token for token in help_tokens if token not in completed.stdout
            ]
            require_cli(
                not missing_tokens,
                f"{help_option} output is missing {missing_tokens}",
                completed,
            )
            require_cli(
                "pathwayFolder" not in completed.stdout,
                f"{help_option} still documents the removed pathwayFolder option",
                completed,
            )
            scenarios += 1

        completed = run_cli_command(executable, [], working_directory)
        require_cli(
            completed.returncode == 2,
            "a missing INPUT should be a command-line error",
            completed,
        )
        require_cli(
            "missing required INPUT" in completed.stderr,
            "the missing-input error should explain what is required",
            completed,
        )
        scenarios += 1

        completed = run_cli_command(
            executable, ["missing-input-file"], working_directory
        )
        require_cli(
            completed.returncode == 1,
            "a missing input file should fail the calculation",
            completed,
        )
        require_cli(
            "input file not found" in completed.stderr,
            "a missing input file should produce a clear error",
            completed,
        )
        scenarios += 1

        compatibility_options = (
            "-runtime=1000000000",
            "--runTime=1000000000",
            "-runTime=1000000000",
            "-enumMax=1000000",
            "-removeHydrogens=0",
            "-compensateDisjoint=0",
            "-disjointCompensation=0",
            "-memTest=0",
            "-testMemory=0",
            "-writeIntermediateMAs=0",
        )
        for option in compatibility_options:
            completed = run_cli_command(
                executable, ["--help", option], working_directory
            )
            require_cli(
                completed.returncode == 0,
                f"compatibility option {option!r} should remain accepted",
                completed,
            )
            scenarios += 1

        invalid_cases = (
            (["input", "--does-not-exist=1"], "unknown option"),
            (["input", "--pathway"], "requires"),
            (["input", "--pathway="], "expects 0 or 1"),
            (["input", "--pathway=2"], "expects 0 or 1"),
            (["input", "--remove-hydrogens=yes"], "expects 0 or 1"),
            (["input", "--enum-max=0"], "expects a value from 1"),
            (["input", "--enum-max=12junk"], "non-negative whole number"),
            (["input", "--runtime=-1"], "non-negative whole number"),
            (["input", f"--runtime={'9' * 100}"], "non-negative whole number"),
            (
                ["input", "--pathway=0", "--pathway=1"],
                "specified more than once",
            ),
            (["first-input", "second-input"], "expected one INPUT"),
        )
        for arguments, error_text in invalid_cases:
            completed = run_cli_command(executable, arguments, working_directory)
            require_cli(
                completed.returncode == 2,
                f"invalid arguments {arguments!r} should be rejected",
                completed,
            )
            require_cli(
                error_text in completed.stderr,
                f"invalid arguments {arguments!r} should report {error_text!r}",
                completed,
            )
            scenarios += 1

        source = TEST_DIRECTORY / "butane.mol"
        input_path = working_directory / "input.mol"
        shutil.copy2(source, input_path)
        completed = run_cli_command(
            executable,
            [
                "--runtime=1000000000",
                "--enum-max=1000000",
                "--pathway=0",
                "--remove-hydrogens=0",
                "--compensate-disjoint=1",
                "--memory-report=0",
                "--write-intermediate-mas=1",
                input_path.name,
            ],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "canonical options before a .mol input should run successfully",
            completed,
        )
        require_cli(
            (working_directory / "inputOut").is_file(),
            "a .mol input should create output without .mol in the output name",
            completed,
        )
        require_cli(
            (working_directory / "inputIntermediateMAs").is_file(),
            "--write-intermediate-mas=1 should create its output",
            completed,
        )
        require_cli(
            not (working_directory / "inputPathway").exists(),
            "--pathway=0 should suppress pathway output",
            completed,
        )
        require_cli(
            not (working_directory / "memUsage").exists(),
            "--memory-report=0 should suppress memory output",
            completed,
        )
        scenarios += 1

        dash_input = working_directory / "-dash-input.mol"
        shutil.copy2(source, dash_input)
        completed = run_cli_command(
            executable,
            ["--pathway=0", "--", dash_input.name],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "-- should allow an input name beginning with a dash",
            completed,
        )
        require_cli(
            (working_directory / "-dash-inputOut").is_file(),
            "the dash-prefixed input should create the expected output",
            completed,
        )
        scenarios += 1

        cutoff_cases = (
            (
                "runtime-limit",
                TEST_DIRECTORY / "butane.mol",
                "--runtime=0",
                "status: runtime limit reached",
            ),
            (
                "enumeration-limit",
                TEST_DIRECTORY / "113.mol",
                "--enum-max=1",
                "status: enumeration limit reached",
            ),
        )
        for name, source_path, limit_option, expected_status in cutoff_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            shutil.copy2(source_path, case_directory / "input.mol")
            completed = run_cli_command(
                executable,
                ["input.mol", "--pathway=0", limit_option],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"{name} scenario should return its best result successfully",
                completed,
            )
            output_path = case_directory / "inputOut"
            index = read_first_line_assembly_index(output_path)
            require_cli(
                index is not None and index != 2_147_483_647,
                f"{name} scenario should put a finite assembly index on the first line",
                completed,
            )
            require_cli(
                expected_status in output_path.read_text().splitlines(),
                f"{name} scenario should record {expected_status!r}",
                completed,
            )
            scenarios += 1

        explicit_hydrogen_mol = "\n".join(
            (
                "Explicit hydrogens",
                "AssemblyCpp CLI test",
                "",
                "  6  5  0  0  0  0  0  0  0  0999 V2000",
                "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "    0.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "    1.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                "    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
                "    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "    3.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "  1  3  1  0  0  0  0",
                "  2  3  1  0  0  0  0",
                "  3  4  1  0  0  0  0",
                "  4  5  1  0  0  0  0",
                "  4  6  1  0  0  0  0",
                "M  END",
                "",
            )
        )
        hydrogen_cases = (
            ("hydrogens-default", [], ["C", "C"], 1),
            (
                "hydrogens-on",
                ["--remove-hydrogens=1"],
                ["C", "C"],
                1,
            ),
            (
                "hydrogens-off",
                ["--remove-hydrogens=0"],
                ["H", "H", "C", "C", "H", "H"],
                5,
            ),
        )
        for (
            name,
            hydrogen_options,
            expected_colours,
            expected_edge_count,
        ) in hydrogen_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input.mol").write_text(explicit_hydrogen_mol)
            completed = run_cli_command(
                executable,
                ["input.mol", "--pathway=1", *hydrogen_options],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"explicit-hydrogen scenario {name!r} should succeed",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"explicit-hydrogen scenario {name!r} omitted its pathway",
                completed,
            )
            pathway = parse_pathway_document(pathway_path)
            graph = pathway["file_graph"][0]
            require_cli(
                graph["VertexColours"] == expected_colours
                and len(graph["Edges"]) == expected_edge_count,
                f"explicit-hydrogen scenario {name!r} transformed the wrong "
                "atoms or bonds",
                completed,
            )
            scenarios += 1

        native_hydrogen_graph = "\n".join(
            (
                "native-hydrogens",
                "4",
                "1 3 2 3 3 4",
                "H H C C",
                "1 1 1",
                "",
            )
        )
        expected_native_graph = None
        for name, hydrogen_option in (
            ("native-hydrogens-on", "--remove-hydrogens=1"),
            ("native-hydrogens-off", "--remove-hydrogens=0"),
        ):
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input").write_text(native_hydrogen_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=1", hydrogen_option],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"native-graph hydrogen scenario {name!r} should succeed",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"native-graph hydrogen scenario {name!r} omitted its pathway",
                completed,
            )
            graph = parse_pathway_document(pathway_path)["file_graph"][0]
            require_cli(
                graph["VertexColours"] == ["H", "H", "C", "C"]
                and len(graph["Edges"]) == 3,
                "--remove-hydrogens should not transform native graph inputs",
                completed,
            )
            if expected_native_graph is None:
                expected_native_graph = graph
            else:
                require_cli(
                    graph == expected_native_graph,
                    "native graph output should be identical with hydrogen "
                    "removal on or off",
                    completed,
                )
            scenarios += 1

        all_hydrogen_mol = "\n".join(
            (
                "Hydrogen",
                "AssemblyCpp CLI test",
                "",
                "  2  1  0  0  0  0  0  0  0  0999 V2000",
                "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",
                "  1  2  1  0  0  0  0",
                "M  END",
                "",
            )
        )
        empty_graph_results: list[tuple[int, int]] = []
        for name, compensation_option in (
            ("empty-compensation-off", "--compensate-disjoint=0"),
            ("empty-compensation-on", "--compensate-disjoint=1"),
        ):
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input.mol").write_text(all_hydrogen_mol)
            completed = run_cli_command(
                executable,
                [
                    "input.mol",
                    "--pathway=0",
                    "--remove-hydrogens=1",
                    compensation_option,
                    "--write-intermediate-mas=1",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"empty-graph compensation scenario {name!r} should succeed",
                completed,
            )
            final_index = read_first_line_assembly_index(case_directory / "inputOut")
            intermediate_index = read_last_intermediate_index(
                case_directory / "inputIntermediateMAs"
            )
            require_cli(
                final_index is not None and intermediate_index == final_index,
                f"empty-graph compensation scenario {name!r} should report "
                "consistent indices",
                completed,
            )
            empty_graph_results.append((final_index, intermediate_index))
            scenarios += 1
        require_cli(
            empty_graph_results[0] == empty_graph_results[1],
            "disjoint compensation should not change an empty processed graph",
        )

        disconnected_graph = "\n".join(
            (
                "disconnected",
                "4",
                "1 2 3 4",
                "C C C C",
                "1 1",
                "",
            )
        )
        compensation_cases = (
            (
                "canonical-off",
                "--compensate-disjoint=0",
                "--write-intermediate-mas=1",
                1,
            ),
            (
                "canonical-on",
                "--compensate-disjoint=1",
                "--write-intermediate-mas=1",
                0,
            ),
            (
                "legacy-parser",
                "-compensateDisjoint=1",
                "-writeIntermediateMAs=1",
                0,
            ),
            (
                "legacy-docs",
                "-disjointCompensation=1",
                "--write-intermediate-mas=1",
                0,
            ),
        )
        for (
            name,
            compensation_option,
            intermediate_option,
            expected_index,
        ) in compensation_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input").write_text(disconnected_graph)
            completed = run_cli_command(
                executable,
                [
                    "input",
                    "-pathway=0",
                    compensation_option,
                    intermediate_option,
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"disjoint-compensation scenario {name!r} should succeed",
                completed,
            )

            output_path = case_directory / "inputOut"
            intermediate_path = case_directory / "inputIntermediateMAs"
            require_cli(
                output_path.is_file() and intermediate_path.is_file(),
                f"disjoint-compensation scenario {name!r} omitted an output file",
                completed,
            )
            output_match = ASSEMBLY_INDEX_PATTERN.search(output_path.read_text())
            require_cli(
                output_match is not None and int(output_match.group(1)) == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong final index",
                completed,
            )
            last_intermediate = read_last_intermediate_index(intermediate_path)
            require_cli(
                last_intermediate == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong intermediate index",
                completed,
            )
            require_cli(
                not (case_directory / "inputPathway").exists(),
                "legacy -pathway=0 should suppress pathway output",
                completed,
            )
            scenarios += 1

        for enum_limit, expect_limit_status in ((1, True), (2, False)):
            case_directory = working_directory / f"enum-boundary-{enum_limit}"
            case_directory.mkdir()
            (case_directory / "input").write_text(disconnected_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=0", f"--enum-max={enum_limit}"],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"enum boundary {enum_limit} should return a best result",
                completed,
            )
            output_path = case_directory / "inputOut"
            status_present = (
                output_path.is_file()
                and "status: enumeration limit reached"
                in output_path.read_text().splitlines()
            )
            require_cli(
                status_present == expect_limit_status,
                f"enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            scenarios += 1

        output_failure_cases = [
            (
                "pathway-output-failure",
                "inputPathway",
                ["--pathway=1", "--write-intermediate-mas=0", "--memory-report=0"],
            ),
            (
                "intermediate-output-failure",
                "inputIntermediateMAs",
                ["--pathway=0", "--write-intermediate-mas=1", "--memory-report=0"],
            ),
        ]
        if sys.platform.startswith("linux"):
            output_failure_cases.append(
                (
                    "memory-output-failure",
                    "memUsage",
                    [
                        "--pathway=0",
                        "--write-intermediate-mas=0",
                        "--memory-report=1",
                    ],
                )
            )

        for name, target_name, output_options in output_failure_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            shutil.copy2(source, case_directory / "input.mol")
            (case_directory / target_name).mkdir()
            completed = run_cli_command(
                executable,
                ["input.mol", *output_options],
                case_directory,
            )
            require_cli(
                completed.returncode != 0,
                f"{name} should fail when {target_name!r} cannot be opened",
                completed,
            )
            require_cli(
                f"could not open output file '{target_name}'" in completed.stderr,
                f"{name} should identify the output it could not open",
                completed,
            )
            scenarios += 1

        disabled_output_directory = working_directory / "disabled-output-sentinels"
        disabled_output_directory.mkdir()
        shutil.copy2(source, disabled_output_directory / "input.mol")
        sentinels = {
            "inputPathway": "existing pathway sentinel\n",
            "inputIntermediateMAs": "existing intermediate sentinel\n",
            "memUsage": "existing memory sentinel\n",
        }
        for filename, content in sentinels.items():
            (disabled_output_directory / filename).write_text(content)
        completed = run_cli_command(
            executable,
            [
                "input.mol",
                "--pathway=0",
                "--write-intermediate-mas=0",
                "--memory-report=0",
            ],
            disabled_output_directory,
        )
        require_cli(
            completed.returncode == 0,
            "disabled output flags should not obstruct a successful calculation",
            completed,
        )
        for filename, content in sentinels.items():
            require_cli(
                (disabled_output_directory / filename).read_text() == content,
                f"disabled output flag unexpectedly overwrote {filename!r}",
                completed,
            )
        scenarios += 1

    for memory_option, expected_on_linux in (
        (None, False),
        ("--memory-report=0", False),
        ("--memory-report=1", True),
        ("-memTest=1", True),
        ("-testMemory=1", True),
    ):
        with tempfile.TemporaryDirectory(
            prefix="assemblycpp-memory-"
        ) as directory:
            working_directory = Path(directory)
            shutil.copy2(TEST_DIRECTORY / "butane.mol", working_directory / "input.mol")
            arguments = ["input.mol", "--pathway=0"]
            if memory_option is not None:
                arguments.append(memory_option)
            completed = run_cli_command(executable, arguments, working_directory)
            require_cli(
                completed.returncode == 0,
                f"memory scenario {memory_option or 'default'} should succeed",
                completed,
            )

            memory_path = working_directory / "memUsage"
            should_exist = expected_on_linux and sys.platform.startswith("linux")
            require_cli(
                memory_path.exists() == should_exist,
                f"memory scenario {memory_option or 'default'} created an unexpected report",
                completed,
            )
            if should_exist:
                require_cli(
                    re.fullmatch(r"VmPeak:\s+\d+\s+kB\n?", memory_path.read_text())
                    is not None,
                    "the Linux memory report should contain a VmPeak value in kB",
                    completed,
                )
            scenarios += 1

    return scenarios


def build_executable(executable: Path, compiler: str) -> Path:
    executable = executable.resolve()
    executable.parent.mkdir(parents=True, exist_ok=True)

    compiler_command = shlex.split(compiler)
    if not compiler_command:
        raise TestConfigurationError("the compiler command is empty")
    if shutil.which(compiler_command[0]) is None:
        raise TestConfigurationError(f"compiler not found: {compiler_command[0]}")

    command = [
        *compiler_command,
        str(REPOSITORY_ROOT / "v5" / "main.cpp"),
        "-std=c++23",
        "-O3",
        "-mpopcnt",
        "-march=x86-64-v3",
    ]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        command.append(f"-I{Path(conda_prefix) / 'include'}")
    command.extend(["-o", str(executable)])

    print(f"Building: {shlex.join(command)}", flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise TestConfigurationError(
            f"build failed with exit code {completed.returncode}"
        )

    return executable


def format_process_diagnostics(completed: subprocess.CompletedProcess[str]) -> str:
    diagnostics = []
    if completed.returncode != 0:
        diagnostics.append(f"process exited with code {completed.returncode}")
    if completed.stdout.strip():
        diagnostics.append(f"stdout: {completed.stdout.strip()}")
    if completed.stderr.strip():
        diagnostics.append(f"stderr: {completed.stderr.strip()}")
    return "; ".join(diagnostics)


def compare_pathway_output(actual_path: Path, expected_path: Path) -> str | None:
    if not actual_path.is_file():
        return f"expected pathway output was not created: {actual_path}"

    try:
        actual = parse_pathway_document(actual_path)
    except ValueError as error:
        return str(error)
    expected = parse_pathway_document(expected_path)
    if actual == expected:
        return None

    expected_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
    actual_lines = json.dumps(actual, indent=2, sort_keys=True).splitlines()
    difference = "\n".join(
        difflib.unified_diff(
            expected_lines,
            actual_lines,
            fromfile=str(expected_path),
            tofile=actual_path.name,
            lineterm="",
        )
    )
    return f"pathway output differs from its golden file:\n{difference}"


def run_test_case(executable: Path, case: TestCase, timeout: float) -> TestResult:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "-", Path(case.name).name)
    started = time.perf_counter()

    try:
        with tempfile.TemporaryDirectory(prefix=f"assemblycpp-{safe_name}-") as directory:
            working_directory = Path(directory)
            copied_input = working_directory / case.source.name
            shutil.copy2(case.source, copied_input)
            input_argument = (
                copied_input.with_suffix("")
                if copied_input.suffix.lower() == ".mol"
                else copied_input
            )

            completed = subprocess.run(
                [
                    str(executable),
                    str(input_argument),
                    f"--pathway={int(case.expected_pathway is not None)}",
                ],
                cwd=working_directory,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
            duration = time.perf_counter() - started
            output_path = Path(f"{input_argument}Out")

            if completed.returncode != 0:
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=format_process_diagnostics(completed),
                )
            if not output_path.is_file():
                diagnostics = format_process_diagnostics(completed)
                suffix = f"; {diagnostics}" if diagnostics else ""
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=f"expected output file was not created: {output_path}{suffix}",
                )

            output = output_path.read_text()
            match = ASSEMBLY_INDEX_PATTERN.search(output)
            if match is None:
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=f"assembly index was not found in {output_path.name}",
                )

            actual_index = int(match.group(1))
            pathway_failure = None
            if case.expected_pathway is not None:
                pathway_failure = compare_pathway_output(
                    Path(f"{input_argument}Pathway"), case.expected_pathway
                )

            return TestResult(
                case=case,
                actual=actual_index,
                duration_seconds=duration,
                failure=pathway_failure,
            )
    except subprocess.TimeoutExpired:
        return TestResult(
            case=case,
            actual=None,
            duration_seconds=time.perf_counter() - started,
            error=f"timed out after {timeout:g} seconds",
        )
    except OSError as error:
        return TestResult(
            case=case,
            actual=None,
            duration_seconds=time.perf_counter() - started,
            error=str(error),
        )


def run_test_cases(
    executable: Path,
    cases: Sequence[TestCase],
    timeout: float,
    jobs: int,
    on_result: Callable[[TestResult], None] | None = None,
) -> list[TestResult]:
    results: list[TestResult] = []

    def record(result: TestResult) -> None:
        results.append(result)
        if on_result is not None:
            on_result(result)

    if jobs == 1:
        for case in cases:
            record(run_test_case(executable, case, timeout))
        return results

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for result in executor.map(
            lambda case: run_test_case(executable, case, timeout), cases
        ):
            record(result)

    return results


def print_result_header() -> None:
    print()
    print(
        f"{'Molecule':<32} {'Expected':>8} {'Actual':>8} "
        f"{'Time (s)':>10}  Status"
    )
    print("-" * 78)


def print_result(result: TestResult) -> None:
    actual = "-" if result.actual is None else str(result.actual)
    print(
        f"{result.case.name:<32.32} {result.case.expected:>8} "
        f"{actual:>8} {result.duration_seconds:>10.3f}  {result.status}"
    )
    if result.error:
        print(f"  {result.error}")
    if result.failure:
        print(f"  {result.failure}")


def print_summary(results: Sequence[TestResult]) -> None:
    passed = sum(result.status == "PASS" for result in results)
    failed = sum(result.status == "FAIL" for result in results)
    errors = sum(result.status == "ERROR" for result in results)
    duration = sum(result.duration_seconds for result in results)

    print(
        f"Summary: {passed} passed, {failed} failed, {errors} errors, "
        f"{len(results)} total ({duration:.3f}s cumulative)"
    )


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build AssemblyCpp and run its regression manifest."
    )
    parser.add_argument(
        "executable",
        nargs="?",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=f"AssemblyCpp executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"tab-separated regression manifest (default: {DEFAULT_MANIFEST.name})",
    )
    parser.add_argument(
        "--pathway-manifest",
        type=Path,
        default=DEFAULT_PATHWAY_MANIFEST,
        help=(
            "tab-separated pathway golden manifest "
            f"(default: {DEFAULT_PATHWAY_MANIFEST.name})"
        ),
    )
    parser.add_argument(
        "--audit",
        action="store_true",
        help="validate the manifest and report fixture coverage without running tests",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="compile v5/main.cpp for x86-64-v3 before running the tests",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler command used with --build (default: CXX or c++)",
    )
    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=1,
        help="number of test processes to run concurrently (default: 1)",
    )
    parser.add_argument(
        "--limit",
        type=positive_int,
        help="run only the first N cases",
    )
    parser.add_argument(
        "--pathways-only",
        action="store_true",
        help="run only cases with expected pathway golden files",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=300.0,
        help="per-case timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print every passing test or list fixture-only molecules during an audit",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)

    try:
        manifest, all_cases = load_manifest(arguments.manifest)
        _, all_cases = load_pathway_manifest(
            arguments.pathway_manifest, all_cases
        )
        if arguments.audit:
            audit_test_data(manifest, all_cases, arguments.verbose)
            return 0

        cases = (
            [case for case in all_cases if case.expected_pathway is not None]
            if arguments.pathways_only
            else all_cases
        )
        cases = cases[: arguments.limit]
        executable_path = arguments.executable
        if arguments.build:
            executable_path = build_executable(executable_path, arguments.compiler)
        executable = resolve_executable(executable_path)
        cli_scenarios = run_cli_checks(executable)
    except TestConfigurationError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"CLI checks: {cli_scenarios} passed", flush=True)
    print(
        f"Running {len(cases)} tests with {arguments.jobs} worker(s); "
        f"timeout={arguments.timeout:g}s",
        flush=True,
    )
    if arguments.verbose:
        print_result_header()
    results = run_test_cases(
        executable=executable,
        cases=cases,
        timeout=arguments.timeout,
        jobs=arguments.jobs,
        on_result=print_result if arguments.verbose else None,
    )

    failures = [result for result in results if result.status != "PASS"]
    if failures and not arguments.verbose:
        print_result_header()
        for result in failures:
            print_result(result)
    print_summary(results)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
