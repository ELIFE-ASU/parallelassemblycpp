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
                    f"-pathway={int(case.expected_pathway is not None)}",
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
        help="compile v5/main.cpp before running the tests",
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
    except TestConfigurationError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

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
