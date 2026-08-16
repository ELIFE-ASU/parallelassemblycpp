#!/usr/bin/env python3
"""Run a small, repeatable AssemblyCpp speed benchmark."""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "AssemblyCpp"
DEFAULT_INPUT = REPOSITORY_ROOT / "unitTests" / "ketoconazole.mol"
DEFAULT_EXPECTED_ASSEMBLY_INDEX = 22
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
CLOCK_TICKS_PATTERN = re.compile(r"^time elapsed:\s*(\d+)\s*$", re.MULTILINE)


class BenchmarkError(RuntimeError):
    """Raised when the benchmark cannot produce a valid measurement."""


@dataclass(frozen=True)
class Measurement:
    wall_seconds: float
    clock_ticks: int


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return parsed


def non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def resolve_file(path: Path, description: str) -> Path:
    candidates = [path]
    if not path.is_absolute():
        candidates.append(REPOSITORY_ROOT / path)

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    raise BenchmarkError(f"{description} not found: {path}")


def resolve_executable(path: Path) -> Path:
    if path.is_file():
        return path.resolve()

    located = shutil.which(str(path))
    if located is not None:
        return Path(located).resolve()

    raise BenchmarkError(
        f"executable not found: {path}. Run again with --build or compile it first."
    )


def build_executable(executable: Path, compiler: str) -> Path:
    executable = executable.resolve()
    executable.parent.mkdir(parents=True, exist_ok=True)

    compiler_command = shlex.split(compiler)
    if not compiler_command:
        raise BenchmarkError("the compiler command is empty")
    if shutil.which(compiler_command[0]) is None:
        raise BenchmarkError(f"compiler not found: {compiler_command[0]}")

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
    try:
        completed = subprocess.run(command, check=False)
    except OSError as error:
        raise BenchmarkError(f"could not run compiler: {error}") from error
    if completed.returncode != 0:
        raise BenchmarkError(f"build failed with exit code {completed.returncode}")

    return executable


def parse_measurement(
    output_path: Path,
    completed: subprocess.CompletedProcess[str],
    wall_seconds: float,
    expected_assembly_index: int | None,
) -> Measurement:
    if completed.returncode != 0:
        diagnostics = completed.stderr.strip() or completed.stdout.strip()
        suffix = f": {diagnostics}" if diagnostics else ""
        raise BenchmarkError(
            f"AssemblyCpp exited with code {completed.returncode}{suffix}"
        )
    if not output_path.is_file():
        raise BenchmarkError(f"AssemblyCpp did not create {output_path.name}")

    try:
        output = output_path.read_text(encoding="utf-8")
    except OSError as error:
        raise BenchmarkError(f"could not read {output_path}: {error}") from error

    assembly_index_match = ASSEMBLY_INDEX_PATTERN.search(output)
    if assembly_index_match is None:
        raise BenchmarkError(f"assembly index not found in {output_path.name}")
    assembly_index = int(assembly_index_match.group(1))
    if (
        expected_assembly_index is not None
        and assembly_index != expected_assembly_index
    ):
        raise BenchmarkError(
            f"expected assembly index {expected_assembly_index}, got {assembly_index}"
        )

    clock_ticks_match = CLOCK_TICKS_PATTERN.search(output)
    if clock_ticks_match is None:
        raise BenchmarkError(
            f"elapsed clock ticks not found in {output_path.name}; "
            "rebuild the executable with --build"
        )

    return Measurement(
        wall_seconds=wall_seconds,
        clock_ticks=int(clock_ticks_match.group(1)),
    )


def run_once(
    executable: Path,
    input_name: str,
    output_path: Path,
    working_directory: Path,
    expected_assembly_index: int | None,
    timeout: float,
) -> Measurement:
    command = [
        str(executable),
        "--pathway=0",
        "--memory-report=0",
        "--write-intermediate-mas=0",
        "--",
        input_name,
    ]
    output_path.unlink(missing_ok=True)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise BenchmarkError(
            f"AssemblyCpp timed out after {error.timeout:g} seconds"
        ) from error
    except OSError as error:
        raise BenchmarkError(f"could not run AssemblyCpp: {error}") from error
    wall_seconds = time.perf_counter() - started

    return parse_measurement(
        output_path,
        completed,
        wall_seconds,
        expected_assembly_index,
    )


def run_benchmark(
    executable: Path,
    source: Path,
    expected_assembly_index: int | None,
    runs: int,
    warmup: int,
    timeout: float,
) -> list[Measurement]:
    with tempfile.TemporaryDirectory(prefix="assemblycpp-benchmark-") as directory:
        working_directory = Path(directory)
        input_path = working_directory / source.name
        shutil.copy2(source, input_path)

        output_name = f"{input_path.name.removesuffix('.mol')}Out"
        output_path = working_directory / output_name

        for current in range(1, warmup + 1):
            print(f"Warm-up {current}/{warmup}...", flush=True)
            run_once(
                executable,
                input_path.name,
                output_path,
                working_directory,
                expected_assembly_index,
                timeout,
            )

        measurements = []
        for current in range(1, runs + 1):
            measurement = run_once(
                executable,
                input_path.name,
                output_path,
                working_directory,
                expected_assembly_index,
                timeout,
            )
            measurements.append(measurement)
            print(
                f"Run {current}/{runs}: {measurement.wall_seconds:.6f} s wall, "
                f"{measurement.clock_ticks} clock ticks",
                flush=True,
            )

    return measurements


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run serial AssemblyCpp calculations and summarize their speed. "
            "The default input is unitTests/ketoconazole.mol."
        )
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=f"AssemblyCpp executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"molfile or native graph to benchmark (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--expected",
        type=int,
        help=(
            "expected assembly index; defaults to 22 for the default input and "
            "is unchecked for a custom input"
        ),
    )
    parser.add_argument(
        "--runs",
        type=positive_int,
        default=5,
        help="number of measured runs (default: 5)",
    )
    parser.add_argument(
        "--warmup",
        type=non_negative_int,
        default=1,
        help="number of unmeasured warm-up runs (default: 1)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=300.0,
        help="per-run timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="compile v5/main.cpp before benchmarking",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler command used with --build (default: CXX or c++)",
    )
    return parser


def print_summary(measurements: Sequence[Measurement]) -> None:
    wall_times = [measurement.wall_seconds for measurement in measurements]
    clock_ticks = [measurement.clock_ticks for measurement in measurements]

    print("\nSummary")
    print(f"  wall min:       {min(wall_times):.6f} s")
    print(f"  wall median:    {statistics.median(wall_times):.6f} s")
    print(f"  wall mean:      {statistics.mean(wall_times):.6f} s")
    print(f"  clock min:      {min(clock_ticks)} ticks")
    print(f"  clock median:   {statistics.median(clock_ticks):g} ticks")
    print(f"  clock mean:     {statistics.mean(clock_ticks):.1f} ticks")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)

    try:
        source = resolve_file(arguments.input, "benchmark input")
        executable_path = arguments.executable
        if arguments.build:
            executable_path = build_executable(
                executable_path, arguments.compiler
            )
        executable = resolve_executable(executable_path)
        expected_assembly_index = arguments.expected
        if expected_assembly_index is None and source == DEFAULT_INPUT.resolve():
            expected_assembly_index = DEFAULT_EXPECTED_ASSEMBLY_INDEX

        print(f"Executable: {executable}")
        print(f"Input: {source}")
        if expected_assembly_index is None:
            print("Assembly index: checked for presence only")
        else:
            print(f"Expected assembly index: {expected_assembly_index}")
        print(
            f"Runs: {arguments.runs} measured, {arguments.warmup} warm-up",
            flush=True,
        )

        measurements = run_benchmark(
            executable=executable,
            source=source,
            expected_assembly_index=expected_assembly_index,
            runs=arguments.runs,
            warmup=arguments.warmup,
            timeout=arguments.timeout,
        )
    except (BenchmarkError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted", file=sys.stderr)
        return 130

    print_summary(measurements)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
