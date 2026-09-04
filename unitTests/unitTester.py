"""Run AssemblyCpp CLI and regression checks."""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "AssemblyCpp"
DEFAULT_MANIFEST = TEST_DIRECTORY / "regression_cases.tsv"
DEFAULT_PATHWAY_MANIFEST = TEST_DIRECTORY / "pathway_cases.tsv"
MANIFEST_HEADER = ("molecule", "expected_assembly_index")
PATHWAY_MANIFEST_HEADER = ("molecule", "expected_pathway")
PATHWAY_KEYS = {"file_graph", "remnant", "duplicates", "removed_edges"}
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
MOLFILE_SUFFIXES = (".mol", ".sdf")


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
    if not math.isfinite(parsed) or parsed <= 0:
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


def has_molfile_suffix(path: Path) -> bool:
    return path.name[-4:].lower() in MOLFILE_SUFFIXES


def molecule_output_path(path: Path, output_suffix: str) -> Path:
    input_name = path.name[:-4] if has_molfile_suffix(path) else path.name
    return path.parent / f"{input_name}{output_suffix}"


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
                    f"invalid header in {manifest}: expected {MANIFEST_HEADER}, "
                    f"got {header}"
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
        # A wrong JSON shape is malformed file content, not caller misuse.
        raise ValueError(  # noqa: TRY004
            f"pathway document must be a JSON object: {path}"
        )
    if set(document) != PATHWAY_KEYS:
        raise ValueError(
            f"invalid pathway keys in {path}: expected {sorted(PATHWAY_KEYS)}, "
            f"got {sorted(document)}"
        )
    for key in PATHWAY_KEYS:
        if not isinstance(document[key], list):
            # Keep every pathway schema violation under the ValueError contract.
            raise ValueError(  # noqa: TRY004
                f"pathway field {key!r} must be an array in {path}"
            )

    return document


def load_pathway_manifest(path: Path, cases: Sequence[TestCase]) -> list[TestCase]:
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
                        f"duplicate pathway case {name!r} in {manifest}:{line_number}"
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

    return [
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
    mol_fixtures = {
        path.resolve()
        for path in fixture_directory.iterdir()
        if path.is_file() and has_molfile_suffix(path)
    }
    referenced_mol = {case.source for case in cases if has_molfile_suffix(case.source)}
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

    graph_cases = sum(not has_molfile_suffix(case.source) for case in cases)
    print(f"Manifest: {manifest}")
    print(
        f"Regression cases: {len(cases)} "
        f"({len(cases) - graph_cases} MOL/SDF, {graph_cases} graph)"
    )
    print(f"Molecule fixtures: {len(mol_fixtures)}")
    print(f"Fixture-only molecules: {len(fixture_only)}")
    print(
        "Shared-content case groups:",
        len(shared_content),
        "(consistent expectations)",
    )
    print(
        "Pathway golden cases: "
        f"{sum(case.expected_pathway is not None for case in cases)}"
    )

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


def make_mask_capacity_graph(component_sizes: Sequence[int]) -> str:
    """Build a bounded-work native graph whose final atom appears in edge zero."""
    atom_colours: list[str] = []
    component_edges: list[list[tuple[int, int]]] = []
    first_atom = 1
    for component, size in enumerate(component_sizes):
        if size < 3:
            raise ValueError("capacity-test components must be cycles of size >= 3")
        vertices = list(range(first_atom, first_atom + size))
        atom_colours.extend([f"C{component}"] * size)
        component_edges.append(
            list(zip(vertices, vertices[1:] + vertices[:1], strict=True))
        )
        first_atom += size

    # Visit the last component first and start at its wraparound edge. This
    # makes the first connected-subgraph seed include the highest atom index,
    # so a 513-atom case necessarily exercises AtomMask word eight as well as
    # EdgeMask word eight.
    last_edges = component_edges[-1]
    edges = [last_edges[-1], *last_edges[:-1]]
    for cycle_edges in component_edges[:-1]:
        edges.extend(cycle_edges)

    return "\n".join(
        (
            "mask-capacity-boundary",
            str(len(atom_colours)),
            " ".join(f"{first} {second}" for first, second in edges),
            " ".join(atom_colours),
            " ".join("1" for _ in edges),
            "",
        )
    )


def run_cli_checks(executable: Path) -> int:
    """Exercise help, validation, aliases, input handling, and output flags."""
    scenarios = 0
    help_tokens = (
        "Usage:",
        "--runtime=<TICKS>",
        "--enum-max=<COUNT>",
        "--pathway=<0|1>",
        "--run-strings=<0|1>",
        "--accept-palindromes=<0|1>",
        "--parallel=<auto|on|off>",
        (
            "Select parallel search automatically, require it, or disable it. "
            "Default: off."
        ),
        "--threads=<auto|N>",
        "--remove-hydrogens=<0|1>",
        "--verbose=<0|1>",
        "--compensate-disjoint=<0|1>",
        "--memory-report=<0|1>",
        "--write-intermediate-mas=<0|1>",
        "Outputs:",
        "Legacy options:",
    )
    telemetry_supported: bool | None = None

    with tempfile.TemporaryDirectory(prefix="assemblycpp-cli-") as directory:
        working_directory = Path(directory)

        for help_option in ("--help", "-h"):
            completed = run_cli_command(executable, [help_option], working_directory)
            require_cli(
                completed.returncode == 0,
                f"{help_option} should exit successfully",
                completed,
            )
            help_has_telemetry = "--telemetry=<0|1>" in completed.stdout
            if telemetry_supported is None:
                telemetry_supported = help_has_telemetry
            require_cli(
                help_has_telemetry == telemetry_supported,
                "help aliases disagree about telemetry support",
                completed,
            )
            expected_help_tokens = help_tokens + (
                ("--telemetry=<0|1>",) if telemetry_supported else ()
            )
            missing_tokens = [
                token for token in expected_help_tokens if token not in completed.stdout
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
            "INPUT is required" in completed.stderr,
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
            "-runStrings=0",
            "--runStrings=0",
            "-acceptPalindromes=0",
            "-palindrome=0",
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

        valid_execution_options = (
            "--parallel=auto",
            "--parallel=on",
            "--parallel=off",
            "--threads=auto",
            "--threads=1",
            "--threads=2147483647",
        )
        for option in valid_execution_options:
            completed = run_cli_command(
                executable, ["--help", option], working_directory
            )
            require_cli(
                completed.returncode == 0,
                f"execution option {option!r} should be accepted",
                completed,
            )
            scenarios += 1

        invalid_cases = [
            (["input", "--does-not-exist=1"], "unknown option"),
            (["input", "--pathway"], "requires a value"),
            (["input", "--pathway="], "expected 0 or 1"),
            (["input", "--pathway=2"], "expected 0 or 1"),
            (["input", "--run-strings=2"], "expected 0 or 1"),
            (["input", "--accept-palindromes=yes"], "expected 0 or 1"),
            (["input", "--parallel"], "requires a value"),
            (["input", "--parallel="], "expected auto, on, or off"),
            (["input", "--parallel=ON"], "expected auto, on, or off"),
            (["input", "--parallel=1"], "expected auto, on, or off"),
            (["input", "--threads"], "requires a value"),
            (["input", "--threads="], "expected a non-negative integer"),
            (["input", "--threads=0"], "expected auto or an integer from 1"),
            (["input", "--threads=-1"], "expected a non-negative integer"),
            (["input", "--threads=2junk"], "expected a non-negative integer"),
            (
                ["input", "--threads=2147483648"],
                "expected auto or an integer from 1",
            ),
            (["input", "--remove-hydrogens=yes"], "expected 0 or 1"),
            (["input", "--verbose=2"], "expected 0 or 1"),
            (["input", "--enum-max=0"], "expected an integer from 1"),
            (["input", "--enum-max=12junk"], "non-negative integer"),
            (["input", "--runtime=-1"], "non-negative integer"),
            (["input", f"--runtime={'9' * 100}"], "non-negative integer"),
            (
                ["input", "--pathway=0", "--pathway=1"],
                "may be specified only once",
            ),
            (
                ["input", "--run-strings=0", "-runStrings=1"],
                "may be specified only once",
            ),
            (
                ["input", "--parallel=auto", "--parallel=off"],
                "may be specified only once",
            ),
            (
                ["input", "--threads=auto", "--threads=2"],
                "may be specified only once",
            ),
            (["first-input", "second-input"], "expected one INPUT"),
        ]
        if telemetry_supported:
            invalid_cases.append((["input", "--telemetry=2"], "expected 0 or 1"))
        else:
            invalid_cases.append((["input", "--telemetry=1"], "unknown option"))
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

        string_directory = working_directory / "string-assembly"
        string_directory.mkdir()
        string_input = string_directory / "strings.txt"
        string_cases = tuple(
            symbol_count * "0" + symbol_count * "1" + symbol_count * "2"
            for symbol_count in (5, 10, 15, 20, 25)
        )
        string_input.write_text("\n".join(string_cases) + "\n")
        completed = run_cli_command(
            executable,
            [
                string_input.name,
                "--run-strings=1",
                "--accept-palindromes=0",
                "--pathway=1",
            ],
            string_directory,
        )
        require_cli(
            completed.returncode == 0,
            "string assembly should process a line-oriented input",
            completed,
        )
        string_output = string_directory / "strings.txtOut"
        require_cli(
            string_output.is_file(),
            "string assembly did not create INPUTOut",
            completed,
        )
        string_output_text = string_output.read_text()
        string_indices = [
            int(match.group(1))
            for match in ASSEMBLY_INDEX_PATTERN.finditer(string_output_text)
        ]
        require_cli(
            string_indices == [11, 14, 17, 17, 20],
            "string assembly disagrees with the upstream reference corpus: "
            f"{string_indices}",
            completed,
        )
        require_cli(
            string_output_text.count("time elapsed:") == len(string_cases),
            "string assembly should record timing for every input line",
            completed,
        )
        for line_index, value in enumerate(string_cases):
            pathway_path = string_directory / f"strings.txt_{line_index}_Pathway"
            require_cli(
                pathway_path.is_file(),
                f"string line {line_index} did not create its pathway",
                completed,
            )
            pathway = json.loads(pathway_path.read_text())
            require_cli(
                pathway["file_graph"][0]["Fragments"] == [value]
                and isinstance(pathway["duplicates"], list),
                f"string line {line_index} pathway has the wrong structure",
                completed,
            )
        scenarios += 1

        no_pathway_directory = working_directory / "string-no-pathway"
        no_pathway_directory.mkdir()
        no_pathway_input = no_pathway_directory / "input"
        no_pathway_input.write_text("abab\n")
        completed = run_cli_command(
            executable,
            ["input", "-runStrings=1", "--pathway=0"],
            no_pathway_directory,
        )
        require_cli(
            completed.returncode == 0,
            "the legacy string flag should run successfully",
            completed,
        )
        require_cli(
            read_first_line_assembly_index(no_pathway_directory / "inputOut") == 2,
            "legacy string mode calculated the wrong index for 'abab'",
            completed,
        )
        require_cli(
            not (no_pathway_directory / "input_0_Pathway").exists(),
            "--pathway=0 should suppress string pathway output",
            completed,
        )
        scenarios += 1

        reversal_directory = working_directory / "string-reversal"
        reversal_directory.mkdir()
        (reversal_directory / "input").write_text("abcxcba\n")
        reversal_indices: list[int | None] = []
        for reversal_option in ("--accept-palindromes=0", "-acceptPalindromes=1"):
            completed = run_cli_command(
                executable,
                [
                    "input",
                    "--run-strings=1",
                    reversal_option,
                    "--pathway=0",
                ],
                reversal_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"string reversal option {reversal_option!r} should succeed",
                completed,
            )
            reversal_indices.append(
                read_first_line_assembly_index(reversal_directory / "inputOut")
            )
        require_cli(
            reversal_indices == [6, 4],
            "--accept-palindromes did not enable reversal equivalence: "
            f"{reversal_indices}",
        )
        scenarios += 1

        unique_string_directory = working_directory / "string-unique-pathway"
        unique_string_directory.mkdir()
        (unique_string_directory / "input").write_text("abcdef\n")
        completed = run_cli_command(
            executable,
            ["input", "--run-strings=1", "--pathway=1"],
            unique_string_directory,
        )
        require_cli(
            completed.returncode == 0,
            "a unique string should still produce a valid pathway document",
            completed,
        )
        unique_pathway = json.loads(
            (unique_string_directory / "input_0_Pathway").read_text()
        )
        require_cli(
            unique_pathway["remnant"][0] == {"Fragments": ["abcdef"], "Positions": [0]}
            and unique_pathway["duplicates"] == [],
            "a no-copy pathway should preserve the whole string as its remnant",
            completed,
        )
        scenarios += 1

        source = TEST_DIRECTORY / "butane.mol"
        input_path = working_directory / "input.mol"
        shutil.copy2(source, input_path)
        canonical_options = [
            "--runtime=1000000000",
            "--enum-max=1000000",
            "--pathway=0",
            "--parallel=off",
            "--threads=2",
            "--remove-hydrogens=0",
            "--verbose=0",
            "--compensate-disjoint=1",
            "--memory-report=0",
            "--write-intermediate-mas=1",
        ]
        if telemetry_supported:
            canonical_options.append("--telemetry=0")
        completed = run_cli_command(
            executable,
            [*canonical_options, input_path.name],
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
        require_cli(
            not (working_directory / "inputTelemetry.json").exists(),
            "--telemetry=0 should suppress telemetry output",
            completed,
        )
        require_cli(
            "Graph:" not in completed.stdout and "  Atom " not in completed.stdout,
            "quiet mode should suppress the parsed graph dump",
            completed,
        )
        scenarios += 1

        extension_directory = working_directory / "molfile-extension-coverage"
        extension_directory.mkdir()
        for case_index, extension in enumerate(
            (".MOL", ".mOl", ".sdf", ".SDF", ".sDf"),
            start=1,
        ):
            case_directory = extension_directory / f"case-{case_index}"
            case_directory.mkdir()
            input_name = f"input{extension}"
            extension_input_path = case_directory / input_name
            if extension.lower() == ".sdf":
                extension_input_path.write_text(
                    source.read_text()
                    + "$$$$\n"
                    + (TEST_DIRECTORY / "alanine.mol").read_text()
                    + "$$$$\n"
                )
            else:
                shutil.copy2(source, extension_input_path)
            completed = run_cli_command(
                executable,
                [
                    input_name,
                    "--pathway=0",
                    "--parallel=off",
                    "--verbose=1",
                    "--remove-hydrogens=0",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"a {extension} input should run successfully",
                completed,
            )
            require_cli(
                f"Input: {input_name}\n" in completed.stdout
                and "Molfile: 4 atoms, 3 bonds" in completed.stdout,
                f"a {extension} input should use MOL parsing",
                completed,
            )
            require_cli(
                (case_directory / "inputOut").is_file(),
                f"a {extension} input should omit its suffix from output names",
                completed,
            )
            require_cli(
                read_first_line_assembly_index(case_directory / "inputOut") == 2,
                f"a {extension} input should use its first V2000 structure",
                completed,
            )
            require_cli(
                not (case_directory / f"{input_name}Out").exists(),
                f"a {extension} input should not retain its suffix in output names",
                completed,
            )
            scenarios += 1

        native_suffix_directory = extension_directory / "native-final-suffix"
        native_suffix_directory.mkdir()
        native_suffix_name = "input.sdf.txt"
        shutil.copy2(
            TEST_DIRECTORY / "graphio_test",
            native_suffix_directory / native_suffix_name,
        )
        completed = run_cli_command(
            executable,
            [native_suffix_name, "--pathway=0", "--verbose=1"],
            native_suffix_directory,
        )
        require_cli(
            completed.returncode == 0,
            "a native graph with non-final .sdf text should run successfully",
            completed,
        )
        require_cli(
            "Molfile:" not in completed.stdout and "Graph:" in completed.stdout,
            "only a final MOL/SDF suffix should select MOL parsing",
            completed,
        )
        require_cli(
            (native_suffix_directory / f"{native_suffix_name}Out").is_file(),
            "a native graph should retain its full filename in output names",
            completed,
        )
        scenarios += 1

        format_parity_directory = working_directory / "input-format-parity"
        mol_directory = format_parity_directory / "mol"
        native_directory = format_parity_directory / "native"
        mol_directory.mkdir(parents=True)
        native_directory.mkdir(parents=True)
        shutil.copy2(TEST_DIRECTORY / "alanine.mol", mol_directory / "input.mol")
        (native_directory / "input").write_text(
            "alanine native graph\n6\n1 2 2 3 3 4 3 5 2 6\nN C C O O C\n1 1 2 1 1\n"
        )

        equivalent_results: list[tuple[int | None, dict[str, object]]] = []
        for format_name, input_name, case_directory in (
            ("MOL", "input.mol", mol_directory),
            ("native graph", "input", native_directory),
        ):
            completed = run_cli_command(
                executable,
                [
                    input_name,
                    "--pathway=1",
                    "--parallel=off",
                    "--remove-hydrogens=0",
                    "--memory-report=0",
                    "--write-intermediate-mas=0",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the equivalent {format_name} input should run successfully",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"the equivalent {format_name} input omitted its pathway",
                completed,
            )
            assembly_index = read_first_line_assembly_index(case_directory / "inputOut")
            require_cli(
                assembly_index is not None,
                f"the equivalent {format_name} input omitted its assembly index",
                completed,
            )
            equivalent_results.append(
                (
                    assembly_index,
                    parse_pathway_document(pathway_path),
                )
            )

        require_cli(
            equivalent_results[0][0] == equivalent_results[1][0],
            "equivalent MOL and native graph inputs should produce the same "
            "assembly index",
        )
        require_cli(
            equivalent_results[0][1] == equivalent_results[1][1],
            "equivalent MOL and native graph inputs should produce the same pathway",
        )
        scenarios += 1

        precedence_directory = working_directory / "native-input-precedence"
        precedence_directory.mkdir()
        shutil.copy2(
            TEST_DIRECTORY / "graphio_test",
            precedence_directory / "input",
        )
        shutil.copy2(
            TEST_DIRECTORY / "tridecane.mol",
            precedence_directory / "input.mol",
        )
        completed = run_cli_command(
            executable,
            ["input", "--pathway=0", "--verbose=1"],
            precedence_directory,
        )
        require_cli(
            completed.returncode == 0,
            "an exact native input with a .mol sibling should run successfully",
            completed,
        )
        require_cli(
            read_first_line_assembly_index(precedence_directory / "inputOut") == 5,
            "an exact native input should take precedence over its .mol sibling",
            completed,
        )
        require_cli(
            "Input: input\n" in completed.stdout
            and "Input: input.mol\n" not in completed.stdout,
            "verbose output should identify the exact native input",
            completed,
        )
        scenarios += 1

        completed = run_cli_command(
            executable,
            ["--pathway=0", "--verbose=1", input_path.name],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "--verbose=1 should run successfully",
            completed,
        )
        require_cli(
            "Input: input.mol" in completed.stdout
            and "Molfile: 4 atoms, 3 bonds" in completed.stdout
            and "Graph: 4 atoms, 3 bonds" in completed.stdout
            and "  Atom 1 (C):" in completed.stdout,
            "--verbose=1 should print the input summary and parsed graph",
            completed,
        )
        scenarios += 1

        parity_indices: list[int | None] = []
        for pathway_enabled in (False, True):
            mode = int(pathway_enabled)
            case_directory = working_directory / f"pathway-parity-{mode}"
            case_directory.mkdir()
            shutil.copy2(
                TEST_DIRECTORY / "ketoconazole.mol",
                case_directory / "input.mol",
            )
            parity_options = [
                "input.mol",
                f"--pathway={mode}",
                "--memory-report=0",
                "--write-intermediate-mas=0",
            ]
            if telemetry_supported:
                parity_options.append("--telemetry=0")
            completed = run_cli_command(
                executable,
                parity_options,
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"ketoconazole pathway={mode} parity run should succeed",
                completed,
            )
            parity_indices.append(
                read_first_line_assembly_index(case_directory / "inputOut")
            )
            require_cli(
                (case_directory / "inputPathway").is_file() == pathway_enabled,
                f"ketoconazole pathway={mode} created unexpected pathway output",
                completed,
            )
        require_cli(
            parity_indices == [22, 22],
            "ketoconazole should have assembly index 22 in both pathway modes",
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
            # Cyclosporin exercises two active words in both mask domains.
            (
                "wide-runtime-limit",
                TEST_DIRECTORY / "cyclosporin.mol",
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

        explicit_hydrogen_mol = (
            "Explicit hydrogens\n"
            "AssemblyCpp CLI test\n"
            "\n"
            "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    0.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    3.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  3  1  0  0  0  0\n"
            "  2  3  1  0  0  0  0\n"
            "  3  4  1  0  0  0  0\n"
            "  4  5  1  0  0  0  0\n"
            "  4  6  1  0  0  0  0\n"
            "M  END\n"
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

        native_hydrogen_graph = "native-hydrogens\n4\n1 3 2 3 3 4\nH H C C\n1 1 1\n"
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

        all_hydrogen_mol = (
            "Hydrogen\n"
            "AssemblyCpp CLI test\n"
            "\n"
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "M  END\n"
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

        disconnected_graph = "disconnected\n4\n1 2 3 4\nC C C C\n1 1\n"
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
                output_match is not None
                and int(output_match.group(1)) == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong "
                "final index",
                completed,
            )
            last_intermediate = read_last_intermediate_index(intermediate_path)
            require_cli(
                last_intermediate == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong "
                "intermediate index",
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

        connected_graph = "connected\n3\n1 2 2 3\nC C C\n1 1\n"
        for enum_limit, expect_limit_status in ((2, True), (3, False)):
            case_directory = working_directory / f"enum-connected-boundary-{enum_limit}"
            case_directory.mkdir()
            (case_directory / "input").write_text(connected_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=0", f"--enum-max={enum_limit}"],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"connected enum boundary {enum_limit} should return a best result",
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
                f"connected enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            scenarios += 1

        for enum_limit, expected_index, expect_limit_status in (
            (86, 12, True),
            (87, 8, False),
        ):
            case_directory = working_directory / f"enum-deep-boundary-{enum_limit}"
            case_directory.mkdir()
            shutil.copy2(
                TEST_DIRECTORY / "113.mol",
                case_directory / "input.mol",
            )
            deep_options = [
                "input.mol",
                "--pathway=0",
                f"--enum-max={enum_limit}",
            ]
            if telemetry_supported:
                deep_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                deep_options,
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"deep enum boundary {enum_limit} should return a best result",
                completed,
            )
            output_path = case_directory / "inputOut"
            output_lines = output_path.read_text().splitlines()
            require_cli(
                read_first_line_assembly_index(output_path) == expected_index,
                f"deep enum boundary {enum_limit} returned the wrong index",
                completed,
            )
            require_cli(
                ("status: enumeration limit reached" in output_lines)
                == expect_limit_status,
                f"deep enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            if telemetry_supported:
                phases = json.loads(
                    (case_directory / "inputTelemetry.json").read_text()
                )["memory"]["phases"]
                expected_later_activations = 0 if expect_limit_status else 1
                require_cli(
                    phases["dag_conversion"]["activations"]
                    == expected_later_activations
                    and phases["assembly_search"]["activations"]
                    == expected_later_activations,
                    f"deep enum boundary {enum_limit} reported incorrect phase "
                    "activity",
                    completed,
                )
            scenarios += 1

        # Exercise scalar, wide, adaptive, later-word, and pre-fragment
        # equivalence-quotient paths.
        for edge_count, expected_index, active_words, cache_outcome in (
            (64, 6, 1, "scalar-lookups"),
            (65, 7, 2, "equivalence-quotient"),
            (127, 10, 2, "wide-hits"),
            (128, 7, 2, "adaptive-fallback"),
            (129, 8, 3, "adaptive-fallback"),
        ):
            wide_graph = "\n".join(
                (
                    "wide-path",
                    str(edge_count + 1),
                    " ".join(
                        f"{vertex} {vertex + 1}" for vertex in range(1, edge_count + 1)
                    ),
                    " ".join("C" for _ in range(edge_count + 1)),
                    " ".join("1" for _ in range(edge_count)),
                    "",
                )
            )
            wide_directory = working_directory / f"wide-initial-dag-{edge_count}"
            wide_directory.mkdir()
            (wide_directory / "input").write_text(wide_graph)
            wide_options = ["input", "--pathway=0"]
            if telemetry_supported:
                wide_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                wide_options,
                wide_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the {edge_count}-edge initial DAG scenario should succeed",
                completed,
            )
            require_cli(
                read_first_line_assembly_index(wide_directory / "inputOut")
                == expected_index,
                f"the {edge_count}-edge initial DAG scenario returned the wrong index",
                completed,
            )
            if not telemetry_supported:
                scenarios += 1
                continue
            telemetry = json.loads((wide_directory / "inputTelemetry.json").read_text())
            counters = telemetry["counters"]
            graph = telemetry["processed_graph"]
            residual = telemetry["caches"]["residual_decomposition"]
            canonical = telemetry["caches"]["canonical_mask"]
            require_cli(
                graph["edges"] == edge_count
                and graph["active_mask_words"] == active_words,
                f"the {edge_count}-edge telemetry reported the wrong mask width",
                completed,
            )
            require_cli(
                residual["eligible_for_processed_graph"] is True,
                f"the {edge_count}-edge telemetry reported the wrong cache eligibility",
                completed,
            )
            require_cli(
                counters["retained_mask_attempts"]
                == counters["retained_masks"]
                + counters["duplicate_mask_attempts"]
                + counters["rejected_masks"],
                f"the {edge_count}-edge retained-mask counters are inconsistent",
                completed,
            )
            require_cli(
                counters["canonicalisation_calls"]
                == canonical["hits"] + canonical["misses"],
                f"the {edge_count}-edge canonical counters are inconsistent",
                completed,
            )
            require_cli(
                residual["lookups"] == residual["hits"] + residual["misses"],
                f"the {edge_count}-edge residual-cache counters are inconsistent",
                completed,
            )
            require_cli(
                counters["retained_masks"] > 0
                and counters["matching_visits"] > 0
                and counters["canonicalisation_calls"] > 0,
                f"the {edge_count}-edge telemetry did not exercise search counters",
                completed,
            )
            require_cli(
                residual["requests"] > 0
                and residual["eligible_requests"] == residual["requests"]
                and residual["small_molecule_bypasses"] == 0
                and residual["wide_molecule_bypasses"] == 0
                and residual["eligible_requests"]
                == residual["small_residual_bypasses"]
                + residual["first_occurrence_bypasses"]
                + residual["runtime_disabled_bypasses"]
                + residual["lookups"],
                f"the {edge_count}-edge case did not exercise an eligible cache path",
                completed,
            )
            if cache_outcome == "scalar-lookups":
                require_cli(
                    residual["lookups"] > 0
                    and residual["admissions"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not exercise scalar caching",
                    completed,
                )
            elif cache_outcome == "wide-hits":
                require_cli(
                    residual["lookups"] > 0
                    and residual["hits"] > 0
                    and residual["admissions"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not exercise wide cache hits",
                    completed,
                )
            elif cache_outcome == "equivalence-quotient":
                require_cli(
                    counters["matching_visits"] < 5000
                    and residual["first_occurrence_bypasses"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not reduce equivalent "
                    "matchings before adaptive fallback",
                    completed,
                )
            else:
                require_cli(
                    residual["first_occurrence_bypasses"] > 0
                    and residual["runtime_disabled_bypasses"] > 0,
                    f"the {edge_count}-edge case did not exercise adaptive fallback",
                    completed,
                )
            phases = telemetry["memory"]["phases"]
            require_cli(
                set(phases)
                == {
                    "input_setup",
                    "initial_enumeration",
                    "dag_conversion",
                    "assembly_search",
                    "output",
                },
                f"the {edge_count}-edge telemetry omitted a search phase",
                completed,
            )
            if sys.platform.startswith("linux"):
                require_cli(
                    all(
                        phase["peak_rss_kib"] is None
                        or (
                            phase["peak_rss_kib"] >= phase["start_rss_kib"]
                            and phase["peak_rss_kib"] >= phase["end_rss_kib"]
                        )
                        for phase in phases.values()
                    ),
                    f"the {edge_count}-edge phase RSS peaks are inconsistent",
                    completed,
                )
            scenarios += 1

        # Cross the former fixed 512-bit cap in both mask domains without
        # introducing a large connected-subgraph search. Each component is a
        # small cycle with its own atom colour; the enumeration cap is reached
        # only after every one-edge EdgeMask has been retained and the first
        # high-index AtomMask seed has expanded once.
        for atom_count, edge_count, component_sizes in (
            (512, 512, [4] * 128),
            (513, 513, [4] * 127 + [5]),
        ):
            capacity_directory = (
                working_directory / f"dynamic-mask-capacity-{atom_count}a-{edge_count}e"
            )
            capacity_directory.mkdir()
            capacity_graph = make_mask_capacity_graph(component_sizes)
            (capacity_directory / "input").write_text(capacity_graph)
            capacity_options = [
                "input",
                "--pathway=0",
                f"--enum-max={edge_count + 1}",
            ]
            if telemetry_supported:
                capacity_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                capacity_options,
                capacity_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "should succeed",
                completed,
            )
            output_path = capacity_directory / "inputOut"
            output_lines = output_path.read_text().splitlines()
            require_cli(
                read_first_line_assembly_index(output_path) == edge_count - 1,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "returned the wrong bounded-search index",
                completed,
            )
            require_cli(
                "status: enumeration limit reached" in output_lines,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "did not reach its deterministic enumeration boundary",
                completed,
            )
            if telemetry_supported:
                telemetry = json.loads(
                    (capacity_directory / "inputTelemetry.json").read_text()
                )
                graph = telemetry["processed_graph"]
                counters = telemetry["counters"]
                require_cli(
                    graph
                    == {
                        "atoms": atom_count,
                        "edges": edge_count,
                        "active_mask_words": (edge_count + 63) // 64,
                    },
                    f"the {atom_count}-atom/{edge_count}-edge capacity telemetry "
                    "reported the wrong dynamic mask dimensions",
                    completed,
                )
                require_cli(
                    counters["retained_masks"] == edge_count + 1
                    and counters["rejected_masks"] == 1,
                    f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                    "did not traverse every singleton mask and the first child",
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
        if telemetry_supported:
            output_failure_cases.append(
                (
                    "telemetry-output-failure",
                    "inputTelemetry.json",
                    ["--pathway=0", "--telemetry=1", "--memory-report=0"],
                )
            )
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
            "inputTelemetry.json": "existing telemetry sentinel\n",
        }
        for filename, content in sentinels.items():
            (disabled_output_directory / filename).write_text(content)
        disabled_options = [
            "input.mol",
            "--pathway=0",
            "--write-intermediate-mas=0",
            "--memory-report=0",
        ]
        if telemetry_supported:
            disabled_options.append("--telemetry=0")
        completed = run_cli_command(
            executable,
            disabled_options,
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
        with tempfile.TemporaryDirectory(prefix="assemblycpp-memory-") as directory:
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
                f"memory scenario {memory_option or 'default'} created an "
                "unexpected report",
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


def compiler_command(compiler: str) -> list[str]:
    compiler_command = shlex.split(compiler)
    if not compiler_command:
        raise TestConfigurationError("the compiler command is empty")
    if shutil.which(compiler_command[0]) is None:
        raise TestConfigurationError(f"compiler not found: {compiler_command[0]}")
    return compiler_command


def run_cpp_unit_test(
    compiler: str,
    source_stem: str,
    label: str,
    *,
    stack_limit_bytes: int | None = None,
) -> None:
    """Build and run one standalone C++ unit-test executable."""
    command_prefix = compiler_command(compiler)
    prefix = f"assemblycpp-{label.replace(' ', '-')}-tests-"
    with tempfile.TemporaryDirectory(prefix=prefix) as directory:
        test_executable = Path(directory) / source_stem
        command = [
            *command_prefix,
            str(TEST_DIRECTORY / f"{source_stem}.cpp"),
            "-std=c++20",
            "-O2",
            "-mpopcnt",
            "-march=x86-64-v3",
            "-o",
            str(test_executable),
        ]
        print(f"Building {label} tests: {shlex.join(command)}", flush=True)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise TestConfigurationError(
                f"{label} test build failed with exit code {completed.returncode}"
            )

        run_options: dict[str, object] = {}
        if stack_limit_bytes is not None and sys.platform.startswith("linux"):
            # resource is unavailable on Windows, so keep this import platform-local.
            import resource  # noqa: PLC0415

            def limit_stack() -> None:
                resource.setrlimit(
                    resource.RLIMIT_STACK,
                    (stack_limit_bytes, stack_limit_bytes),
                )

            run_options["preexec_fn"] = limit_stack
        completed = subprocess.run([str(test_executable)], check=False, **run_options)
        if completed.returncode != 0:
            raise TestConfigurationError(
                f"{label} tests failed with exit code {completed.returncode}"
            )
        print(f"{label.capitalize()} tests: passed", flush=True)


def run_mask_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(compiler, "activeWordMaskTester", "mask")


def run_tree_canon_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(
        compiler,
        "treeCanonTester",
        "tree canon",
        stack_limit_bytes=64 * 1024,
    )


def run_cyclic_canon_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(compiler, "cyclicCanonTester", "cyclic canon")


def build_executable(executable: Path, compiler: str) -> Path:
    executable = executable.resolve()
    executable.parent.mkdir(parents=True, exist_ok=True)

    run_mask_unit_tests(compiler)
    run_tree_canon_unit_tests(compiler)
    run_cyclic_canon_unit_tests(compiler)
    command_prefix = compiler_command(compiler)

    command = [
        *command_prefix,
        str(REPOSITORY_ROOT / "v5" / "main.cpp"),
        "-std=c++20",
        "-O3",
        "-mpopcnt",
        "-march=x86-64-v3",
        "-DASSEMBLY_ENABLE_TELEMETRY",
    ]
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
        with tempfile.TemporaryDirectory(
            prefix=f"assemblycpp-{safe_name}-"
        ) as directory:
            working_directory = Path(directory)
            copied_input = working_directory / case.source.name
            shutil.copy2(case.source, copied_input)
            input_argument = (
                copied_input.with_suffix("")
                if copied_input.name.endswith(".mol")
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
            output_path = molecule_output_path(copied_input, "Out")

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
                missing_output_error = (
                    f"expected output file was not created: {output_path}{suffix}"
                )
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=missing_output_error,
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
                    molecule_output_path(copied_input, "Pathway"),
                    case.expected_pathway,
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
    print(f"{'Molecule':<32} {'Expected':>8} {'Actual':>8} {'Time (s)':>10}  Status")
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
        description="Run AssemblyCpp CLI and regression checks."
    )
    parser.add_argument(
        "executable",
        nargs="?",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=(
            "AssemblyCpp executable "
            f"(default: {DEFAULT_EXECUTABLE.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=(
            "regression manifest (TSV; default: "
            f"{DEFAULT_MANIFEST.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--pathway-manifest",
        type=Path,
        default=DEFAULT_PATHWAY_MANIFEST,
        help=(
            "pathway manifest (TSV; default: "
            f"{DEFAULT_PATHWAY_MANIFEST.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--audit",
        action="store_true",
        help="check manifests and fixture coverage without running tests",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="build x86-64-v3 test executables before testing",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler used by --build (uses $CXX, then c++)",
    )
    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=1,
        help="parallel test processes (default: 1)",
    )
    parser.add_argument(
        "--limit",
        type=positive_int,
        help="run the first N regression cases",
    )
    parser.add_argument(
        "--pathways-only",
        action="store_true",
        help="select regression cases with pathway golden files",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=300.0,
        help="timeout per case in seconds (default: 300)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="list passing checks; in audit mode, list fixture-only molecules",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)

    try:
        manifest, all_cases = load_manifest(arguments.manifest)
        all_cases = load_pathway_manifest(arguments.pathway_manifest, all_cases)
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
