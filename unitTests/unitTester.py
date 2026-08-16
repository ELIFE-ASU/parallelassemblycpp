#!/usr/bin/env python3
"""Build AssemblyCpp and run its molecule regression batteries."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "AssemblyCpp"
DEFAULT_MOLECULES = TEST_DIRECTORY / "batteryTest2"
DEFAULT_EXPECTED = TEST_DIRECTORY / "batteryTest2Base"
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")


class TestConfigurationError(RuntimeError):
    """Raised when a test input or build setting is invalid."""


@dataclass(frozen=True)
class TestCase:
    name: str
    source: Path
    expected: int


@dataclass(frozen=True)
class TestResult:
    case: TestCase
    actual: int | None
    duration_seconds: float
    error: str | None = None

    @property
    def status(self) -> str:
        if self.error is not None:
            return "ERROR"
        if self.actual == self.case.expected:
            return "PASS"
        return "FAIL"


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


def read_nonempty_lines(path: Path) -> list[str]:
    try:
        return [line.strip() for line in path.read_text().splitlines() if line.strip()]
    except OSError as error:
        raise TestConfigurationError(f"cannot read {path}: {error}") from error


def resolve_data_file(path: Path) -> Path:
    if path.is_absolute() or path.exists():
        return path.resolve()

    test_relative = TEST_DIRECTORY / path
    if test_relative.exists():
        return test_relative.resolve()

    return path.resolve()


def resolve_molecule(name: str, data_directory: Path) -> Path:
    requested = Path(name)
    candidate = requested if requested.is_absolute() else data_directory / requested

    if candidate.is_file():
        return candidate.resolve()

    mol_candidate = candidate.with_suffix(".mol")
    if mol_candidate.is_file():
        return mol_candidate.resolve()

    raise TestConfigurationError(
        f"molecule {name!r} was not found as {candidate} or {mol_candidate}"
    )


def load_test_cases(
    molecule_list: Path, expected_results: Path, limit: int | None
) -> list[TestCase]:
    molecule_list = resolve_data_file(molecule_list)
    expected_results = resolve_data_file(expected_results)
    molecule_names = read_nonempty_lines(molecule_list)
    if not molecule_names:
        raise TestConfigurationError(f"molecule list is empty: {molecule_list}")

    expected_values: list[int] = []
    for line_number, value in enumerate(
        read_nonempty_lines(expected_results), start=1
    ):
        try:
            expected_values.append(int(value))
        except ValueError as error:
            raise TestConfigurationError(
                f"invalid assembly index in {expected_results}:{line_number}: {value!r}"
            ) from error

    if len(molecule_names) != len(expected_values):
        raise TestConfigurationError(
            "molecule and expected-result files have different lengths: "
            f"{len(molecule_names)} != {len(expected_values)}"
        )

    selected = zip(molecule_names, expected_values)
    if limit is not None:
        selected = list(selected)[:limit]

    return [
        TestCase(
            name=name,
            source=resolve_molecule(name, molecule_list.parent),
            expected=expected,
        )
        for name, expected in selected
    ]


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


def run_test_case(executable: Path, case: TestCase, timeout: float) -> TestResult:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "-", Path(case.name).name)
    started = time.perf_counter()

    try:
        with tempfile.TemporaryDirectory(prefix=f"assemblycpp-{safe_name}-") as directory:
            working_directory = Path(directory)
            copied_input = working_directory / case.source.name
            shutil.copy2(case.source, copied_input)

            if copied_input.suffix.lower() == ".mol":
                input_argument = copied_input.with_suffix("")
            else:
                input_argument = copied_input

            completed = subprocess.run(
                [str(executable), str(input_argument), "-pathway=0"],
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
                    error=f"assembly index was not found in {output_path.name}: {output.strip()}",
                )

            return TestResult(
                case=case,
                actual=int(match.group(1)),
                duration_seconds=duration,
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
            lambda case: run_test_case(executable, case, timeout),
            cases,
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
        f"{actual:>8} {result.duration_seconds:>10.3f}  {result.status}",
        flush=True,
    )
    if result.error:
        print(f"  {result.error}", flush=True)


def print_summary(results: Sequence[TestResult]) -> None:
    passed = sum(result.status == "PASS" for result in results)
    failed = sum(result.status == "FAIL" for result in results)
    errors = sum(result.status == "ERROR" for result in results)
    duration = sum(result.duration_seconds for result in results)

    print()
    print(
        f"Summary: {passed} passed, {failed} failed, {errors} errors, "
        f"{len(results)} total ({duration:.3f}s cumulative)"
    )


def print_results(results: Sequence[TestResult]) -> None:
    print_result_header()
    for result in results:
        print_result(result)
    print_summary(results)


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build AssemblyCpp and run a molecule regression battery."
    )
    parser.add_argument(
        "executable",
        nargs="?",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=f"AssemblyCpp executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "molecule_list",
        nargs="?",
        type=Path,
        default=DEFAULT_MOLECULES,
        help=f"molecule list (default: {DEFAULT_MOLECULES.name})",
    )
    parser.add_argument(
        "expected_results",
        nargs="?",
        type=Path,
        default=DEFAULT_EXPECTED,
        help=f"expected assembly indices (default: {DEFAULT_EXPECTED.name})",
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
        "--timeout",
        type=positive_float,
        default=3600.0,
        help="per-case timeout in seconds (default: 3600)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)

    try:
        executable_path = arguments.executable
        if arguments.build:
            executable_path = build_executable(executable_path, arguments.compiler)
        executable = resolve_executable(executable_path)
        cases = load_test_cases(
            arguments.molecule_list,
            arguments.expected_results,
            arguments.limit,
        )
    except TestConfigurationError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(
        f"Running {len(cases)} tests with {arguments.jobs} worker(s); "
        f"timeout={arguments.timeout:g}s",
        flush=True,
    )
    print_result_header()
    results = run_test_cases(
        executable=executable,
        cases=cases,
        timeout=arguments.timeout,
        jobs=arguments.jobs,
        on_result=print_result,
    )
    print_summary(results)
    return 0 if all(result.status == "PASS" for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
