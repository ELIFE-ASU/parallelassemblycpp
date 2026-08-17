"""Run repeatable single-input or corpus AssemblyCpp speed benchmarks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_DIRECTORY = Path(__file__).resolve().parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "AssemblyCpp"
DEFAULT_INPUT = REPOSITORY_ROOT / "unitTests" / "ketoconazole.mol"
DEFAULT_EXPECTED_ASSEMBLY_INDEX = 22
DEFAULT_MANIFEST = BENCHMARK_DIRECTORY / "cases.tsv"
MANIFEST_HEADER = (
    "name",
    "input",
    "expected_assembly_index",
    "expectation",
    "suites",
    "workload",
)
KNOWN_SUITES = ("quick", "full", "profile", "scaling")
KNOWN_EXPECTATIONS = ("reviewed", "provisional")
CASE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
CLOCK_TICKS_PATTERN = re.compile(r"^time elapsed:\s*(\d+)\s*$", re.MULTILINE)


class BenchmarkError(RuntimeError):
    """Raised when the benchmark cannot produce a valid measurement."""


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    source: Path
    expected_assembly_index: int | None
    expectation: str
    suites: tuple[str, ...]
    workload: str


@dataclass(frozen=True)
class PreparedCase:
    case: BenchmarkCase
    input_name: str
    output_path: Path
    working_directory: Path


@dataclass(frozen=True)
class Measurement:
    wall_seconds: float
    clock_ticks: int
    assembly_index: int


@dataclass(frozen=True)
class MetricSummary:
    minimum: float
    median: float
    mean: float
    mad: float
    p95: float


@dataclass(frozen=True)
class CaseResult:
    case: BenchmarkCase
    measurements: tuple[Measurement, ...]
    baseline_measurements: tuple[Measurement, ...] = ()


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


def resolve_manifest_input(path: Path, manifest: Path) -> Path:
    candidates = [path]
    if not path.is_absolute():
        candidates = [manifest.parent / path, REPOSITORY_ROOT / path, path]

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    raise BenchmarkError(f"benchmark input not found in {manifest.name}: {path}")


def resolve_executable(path: Path) -> Path:
    if path.is_file():
        return path.resolve()

    located = shutil.which(str(path))
    if located is not None:
        return Path(located).resolve()

    raise BenchmarkError(
        f"executable not found: {path}. Run again with --build or compile it first."
    )


def ensure_distinct_executables(candidate: Path, baseline: Path) -> None:
    candidate = candidate.resolve()
    baseline = baseline.resolve()
    same_file = candidate == baseline
    if not same_file and candidate.exists():
        try:
            same_file = candidate.samefile(baseline)
        except OSError:
            pass
    if same_file:
        raise BenchmarkError(
            "candidate and baseline executables resolve to the same file"
        )


def executable_metadata(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        size = path.stat().st_size
    except OSError as error:
        raise BenchmarkError(f"could not fingerprint {path}: {error}") from error
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "size_bytes": size,
    }


def verify_executable_unchanged(
    path: Path, expected: dict[str, object], role: str
) -> None:
    if executable_metadata(path) != expected:
        raise BenchmarkError(f"{role} executable changed during the benchmark")


def cpu_description() -> str:
    processor = platform.processor().strip()
    if processor:
        return processor
    try:
        with Path("/proc/cpuinfo").open(encoding="utf-8") as stream:
            for line in stream:
                key, separator, value = line.partition(":")
                if (
                    separator
                    and key.strip() in {"model name", "Processor"}
                    and (description := value.strip())
                ):
                    return description
    except OSError:
        pass
    return platform.machine() or "unknown"


def load_manifest(path: Path) -> tuple[Path, list[BenchmarkCase]]:
    manifest = resolve_file(path, "benchmark manifest")
    cases: list[BenchmarkCase] = []
    names: set[str] = set()

    try:
        with manifest.open(encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != MANIFEST_HEADER:
                raise BenchmarkError(
                    f"invalid header in {manifest}: expected {MANIFEST_HEADER}, "
                    f"got {header}"
                )

            for line_number, row in enumerate(reader, start=2):
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != len(MANIFEST_HEADER):
                    raise BenchmarkError(
                        f"invalid row in {manifest}:{line_number}: expected "
                        f"{len(MANIFEST_HEADER)} columns"
                    )

                (
                    name,
                    input_text,
                    expected_text,
                    expectation,
                    suites_text,
                    workload,
                ) = (value.strip() for value in row)
                if not CASE_NAME_PATTERN.fullmatch(name):
                    raise BenchmarkError(
                        f"invalid benchmark name in {manifest}:{line_number}: {name!r}"
                    )
                if name in names:
                    raise BenchmarkError(
                        f"duplicate benchmark name in {manifest}:{line_number}: "
                        f"{name!r}"
                    )
                names.add(name)

                try:
                    expected = int(expected_text)
                except ValueError as error:
                    raise BenchmarkError(
                        f"invalid assembly index in {manifest}:{line_number}: "
                        f"{expected_text!r}"
                    ) from error

                if expectation not in KNOWN_EXPECTATIONS:
                    raise BenchmarkError(
                        f"invalid expectation in {manifest}:{line_number}: "
                        f"{expectation!r}"
                    )

                suites = tuple(
                    suite.strip() for suite in suites_text.split(",") if suite.strip()
                )
                if not suites or len(set(suites)) != len(suites):
                    raise BenchmarkError(
                        f"invalid suites in {manifest}:{line_number}: {suites_text!r}"
                    )
                unknown_suites = sorted(set(suites) - set(KNOWN_SUITES))
                if unknown_suites:
                    raise BenchmarkError(
                        f"unknown suites in {manifest}:{line_number}: "
                        f"{', '.join(unknown_suites)}"
                    )
                if not workload:
                    raise BenchmarkError(f"empty workload in {manifest}:{line_number}")

                source = resolve_manifest_input(Path(input_text), manifest)
                cases.append(
                    BenchmarkCase(
                        name=name,
                        source=source,
                        expected_assembly_index=expected,
                        expectation=expectation,
                        suites=suites,
                        workload=workload,
                    )
                )
    except OSError as error:
        raise BenchmarkError(f"could not read {manifest}: {error}") from error

    if not cases:
        raise BenchmarkError(f"benchmark manifest is empty: {manifest}")
    return manifest, cases


def select_cases(
    cases: Sequence[BenchmarkCase],
    suite: str | None,
    requested_names: Sequence[str],
) -> list[BenchmarkCase]:
    selected = [case for case in cases if suite is None or suite in case.suites]
    if requested_names:
        requested = set(requested_names)
        known = {case.name for case in cases}
        unknown = sorted(requested - known)
        if unknown:
            raise BenchmarkError(f"unknown benchmark case(s): {', '.join(unknown)}")
        selected = [case for case in selected if case.name in requested]
        missing = sorted(requested - {case.name for case in selected})
        if missing:
            raise BenchmarkError(
                f"case(s) not in the {suite!r} suite: {', '.join(missing)}"
            )

    if not selected:
        label = f" suite {suite!r}" if suite is not None else ""
        raise BenchmarkError(f"no benchmark cases selected from{label}")
    return selected


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
        "-mpopcnt",
        "-march=x86-64-v3",
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
        assembly_index=assembly_index,
    )


def prepare_cases(
    cases: Sequence[BenchmarkCase], working_root: Path
) -> list[PreparedCase]:
    prepared: list[PreparedCase] = []
    for index, case in enumerate(cases):
        working_directory = working_root / f"{index:03d}-{case.name}"
        working_directory.mkdir()
        input_path = working_directory / case.source.name
        shutil.copy2(case.source, input_path)
        output_name = f"{input_path.name.removesuffix('.mol')}Out"
        prepared.append(
            PreparedCase(
                case=case,
                input_name=input_path.name,
                output_path=working_directory / output_name,
                working_directory=working_directory,
            )
        )
    return prepared


def run_once(
    executable: Path,
    prepared: PreparedCase,
    timeout: float,
) -> Measurement:
    command = [
        str(executable),
        "--pathway=0",
        "--memory-report=0",
        "--write-intermediate-mas=0",
        "--",
        prepared.input_name,
    ]
    prepared.output_path.unlink(missing_ok=True)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=prepared.working_directory,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise BenchmarkError(
            f"{prepared.case.name} timed out after {error.timeout:g} seconds"
        ) from error
    except OSError as error:
        raise BenchmarkError(
            f"could not run AssemblyCpp for {prepared.case.name}: {error}"
        ) from error
    wall_seconds = time.perf_counter() - started

    try:
        return parse_measurement(
            prepared.output_path,
            completed,
            wall_seconds,
            prepared.case.expected_assembly_index,
        )
    except BenchmarkError as error:
        raise BenchmarkError(f"{prepared.case.name}: {error}") from error


def rotate_cases(cases: Sequence[PreparedCase], round_index: int) -> list[PreparedCase]:
    if not cases:
        return []
    offset = round_index % len(cases)
    return [*cases[offset:], *cases[:offset]]


def run_benchmarks(
    executable: Path,
    baseline_executable: Path | None,
    cases: Sequence[BenchmarkCase],
    runs: int,
    warmup: int,
    timeout: float,
) -> list[CaseResult]:
    with tempfile.TemporaryDirectory(prefix="assemblycpp-benchmark-") as directory:
        prepared_cases = prepare_cases(cases, Path(directory))
        measurements: dict[str, list[Measurement]] = {case.name: [] for case in cases}
        baseline_measurements: dict[str, list[Measurement]] = {
            case.name: [] for case in cases
        }
        case_count = len(prepared_cases)

        def executable_order(round_index: int) -> list[tuple[str, Path]]:
            if baseline_executable is None:
                return [("candidate", executable)]
            if round_index % 2 == 0:
                return [
                    ("baseline", baseline_executable),
                    ("candidate", executable),
                ]
            return [
                ("candidate", executable),
                ("baseline", baseline_executable),
            ]

        def validate_unchecked_pair(
            prepared: PreparedCase, paired: dict[str, Measurement]
        ) -> None:
            if (
                baseline_executable is not None
                and prepared.case.expected_assembly_index is None
                and paired["candidate"].assembly_index
                != paired["baseline"].assembly_index
            ):
                raise BenchmarkError(
                    f"{prepared.case.name}: candidate assembly index "
                    f"{paired['candidate'].assembly_index} does not match "
                    "baseline assembly index "
                    f"{paired['baseline'].assembly_index}"
                )

        for warmup_round in range(warmup):
            for prepared in rotate_cases(prepared_cases, warmup_round):
                paired_warmups: dict[str, Measurement] = {}
                for role, current_executable in executable_order(warmup_round):
                    role_suffix = f" {role}" if baseline_executable is not None else ""
                    print(
                        f"Warm-up {warmup_round + 1}/{warmup} "
                        f"[{prepared.case.name}{role_suffix}]...",
                        flush=True,
                    )
                    try:
                        paired_warmups[role] = run_once(
                            current_executable, prepared, timeout
                        )
                    except BenchmarkError as error:
                        raise BenchmarkError(
                            f"warm-up {warmup_round + 1}/{warmup} {role}: {error}"
                        ) from error
                validate_unchecked_pair(prepared, paired_warmups)

        for run_index in range(runs):
            ordered = rotate_cases(prepared_cases, warmup + run_index)
            for case_index, prepared in enumerate(ordered, start=1):
                paired_measurements: dict[str, Measurement] = {}
                for role, current_executable in executable_order(run_index):
                    try:
                        measurement = run_once(current_executable, prepared, timeout)
                    except BenchmarkError as error:
                        raise BenchmarkError(
                            f"round {run_index + 1}/{runs} {role}: {error}"
                        ) from error
                    paired_measurements[role] = measurement
                    if role == "candidate":
                        measurements[prepared.case.name].append(measurement)
                    else:
                        baseline_measurements[prepared.case.name].append(measurement)

                    if case_count == 1:
                        prefix = f"Run {run_index + 1}/{runs}"
                    else:
                        prefix = (
                            f"Round {run_index + 1}/{runs} "
                            f"[{case_index}/{case_count}] {prepared.case.name}"
                        )
                    if baseline_executable is not None:
                        prefix += f" {role}"
                    print(
                        f"{prefix}: {measurement.wall_seconds:.6f} s wall, "
                        f"{measurement.clock_ticks} clock ticks, "
                        f"AI {measurement.assembly_index}",
                        flush=True,
                    )

                validate_unchecked_pair(prepared, paired_measurements)

    return [
        CaseResult(
            case=case,
            measurements=tuple(measurements[case.name]),
            baseline_measurements=tuple(baseline_measurements[case.name]),
        )
        for case in cases
    ]


def percentile(values: Sequence[float], percentile_value: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot calculate a percentile of no values")
    position = (len(ordered) - 1) * percentile_value
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    fraction = position - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


def summarize(values: Sequence[float]) -> MetricSummary:
    if not values:
        raise ValueError("cannot summarize no values")
    median = statistics.median(values)
    return MetricSummary(
        minimum=min(values),
        median=median,
        mean=statistics.mean(values),
        mad=statistics.median(abs(value - median) for value in values),
        p95=percentile(values, 0.95),
    )


def result_summaries(
    result: CaseResult, *, baseline: bool = False
) -> tuple[MetricSummary, MetricSummary]:
    measurements = result.baseline_measurements if baseline else result.measurements
    return (
        summarize([measurement.wall_seconds for measurement in measurements]),
        summarize([float(measurement.clock_ticks) for measurement in measurements]),
    )


def aggregate_round_summaries(
    results: Sequence[CaseResult], *, baseline: bool = False
) -> tuple[MetricSummary, MetricSummary]:
    selected = [
        result.baseline_measurements if baseline else result.measurements
        for result in results
    ]
    run_count = len(selected[0])
    wall_totals = [
        sum(measurements[index].wall_seconds for measurements in selected)
        for index in range(run_count)
    ]
    clock_totals = [
        float(sum(measurements[index].clock_ticks for measurements in selected))
        for index in range(run_count)
    ]
    return summarize(wall_totals), summarize(clock_totals)


def paired_speedup_summary(
    candidate: Sequence[Measurement],
    baseline: Sequence[Measurement],
    attribute: str,
) -> MetricSummary | None:
    ratios = []
    for candidate_value, baseline_value in zip(candidate, baseline, strict=True):
        candidate_metric = float(getattr(candidate_value, attribute))
        baseline_metric = float(getattr(baseline_value, attribute))
        if candidate_metric <= 0 or baseline_metric <= 0:
            return None
        ratios.append(baseline_metric / candidate_metric)
    return summarize(ratios)


def aggregate_paired_speedup_summary(
    results: Sequence[CaseResult], attribute: str
) -> MetricSummary | None:
    run_count = len(results[0].measurements)
    ratios = []
    for index in range(run_count):
        candidate_total = sum(
            float(getattr(result.measurements[index], attribute)) for result in results
        )
        baseline_total = sum(
            float(getattr(result.baseline_measurements[index], attribute))
            for result in results
        )
        if candidate_total <= 0 or baseline_total <= 0:
            return None
        ratios.append(baseline_total / candidate_total)
    return summarize(ratios)


def geometric_mean(values: Sequence[float]) -> float:
    if not values or any(value <= 0 for value in values):
        raise ValueError("geometric mean requires positive values")
    return math.exp(statistics.mean(math.log(value) for value in values))


def equal_weight_paired_speedup_summary(
    results: Sequence[CaseResult], attribute: str
) -> MetricSummary | None:
    run_count = len(results[0].measurements)
    round_geometric_means = []
    for index in range(run_count):
        ratios = []
        for result in results:
            candidate = float(getattr(result.measurements[index], attribute))
            baseline = float(getattr(result.baseline_measurements[index], attribute))
            if candidate <= 0 or baseline <= 0:
                return None
            ratios.append(baseline / candidate)
        round_geometric_means.append(geometric_mean(ratios))
    return summarize(round_geometric_means)


def print_summary(results: Sequence[CaseResult]) -> None:
    if results[0].baseline_measurements:
        print_comparison_summary(results)
        return

    print("\nSummary")
    if len(results) == 1:
        wall, clock = result_summaries(results[0])
        print(f"  wall min:       {wall.minimum:.6f} s")
        print(f"  wall median:    {wall.median:.6f} s")
        print(f"  wall mean:      {wall.mean:.6f} s")
        print(f"  clock min:      {clock.minimum:g} ticks")
        print(f"  clock median:   {clock.median:g} ticks")
        print(f"  clock mean:     {clock.mean:.1f} ticks")
        print(f"  wall MAD:       {wall.mad:.6f} s")
        print(f"  wall p95:       {wall.p95:.6f} s")
        print(f"  clock MAD:      {clock.mad:g} ticks")
        print(f"  clock p95:      {clock.p95:g} ticks")
        return

    name_width = max(12, min(24, max(len(result.case.name) for result in results)))
    workload_width = max(
        20, min(32, max(len(result.case.workload) for result in results))
    )
    print(
        f"  {'Case':<{name_width}} {'Workload':<{workload_width}} "
        f"{'Wall median':>12} {'MAD':>10} {'p95':>12} "
        f"{'Clock median':>14}"
    )
    print("  " + "-" * (name_width + workload_width + 53))
    for result in results:
        wall, clock = result_summaries(result)
        print(
            f"  {result.case.name:<{name_width}.{name_width}} "
            f"{result.case.workload:<{workload_width}.{workload_width}} "
            f"{wall.median:>11.6f}s {wall.mad:>9.6f}s "
            f"{wall.p95:>11.6f}s {clock.median:>14g}"
        )

    wall, clock = aggregate_round_summaries(results)
    print("\nSuite round totals")
    print(f"  wall median:    {wall.median:.6f} s")
    print(f"  wall MAD:       {wall.mad:.6f} s")
    print(f"  wall p95:       {wall.p95:.6f} s")
    print(f"  clock median:   {clock.median:g} ticks")


def print_comparison_summary(results: Sequence[CaseResult]) -> None:
    print("\nComparison summary (speedup > 1 means candidate is faster)")
    name_width = max(12, min(24, max(len(result.case.name) for result in results)))
    print(
        f"  {'Case':<{name_width}} {'Baseline':>12} {'Candidate':>12} "
        f"{'Wall speedup':>14} {'Clock speedup':>15}"
    )
    print("  " + "-" * (name_width + 59))
    for result in results:
        candidate_wall, _ = result_summaries(result)
        baseline_wall, _ = result_summaries(result, baseline=True)
        wall_speedup = paired_speedup_summary(
            result.measurements,
            result.baseline_measurements,
            "wall_seconds",
        )
        clock_speedup = paired_speedup_summary(
            result.measurements,
            result.baseline_measurements,
            "clock_ticks",
        )
        assert wall_speedup is not None
        clock_text = "n/a"
        if clock_speedup is not None:
            clock_text = f"{clock_speedup.median:.4f}x"
        print(
            f"  {result.case.name:<{name_width}.{name_width}} "
            f"{baseline_wall.median:>11.6f}s "
            f"{candidate_wall.median:>11.6f}s "
            f"{wall_speedup.median:>13.4f}x {clock_text:>15}"
        )

    round_wall = aggregate_paired_speedup_summary(results, "wall_seconds")
    round_clock = aggregate_paired_speedup_summary(results, "clock_ticks")
    equal_wall = equal_weight_paired_speedup_summary(results, "wall_seconds")
    equal_clock = equal_weight_paired_speedup_summary(results, "clock_ticks")
    assert round_wall is not None
    assert equal_wall is not None
    print("\nCorpus comparison")
    print(
        f"  primary paired round-total wall:  {round_wall.median:.4f}x "
        f"(MAD {round_wall.mad:.4f}x)"
    )
    if round_clock is not None:
        print(
            f"  primary paired round-total clock: {round_clock.median:.4f}x "
            f"(MAD {round_clock.mad:.4f}x)"
        )
    print(
        f"  descriptive equal-weight wall:    {equal_wall.median:.4f}x "
        f"(MAD {equal_wall.mad:.4f}x)"
    )
    if equal_clock is not None:
        print(
            f"  descriptive equal-weight clock:   {equal_clock.median:.4f}x "
            f"(MAD {equal_clock.mad:.4f}x)"
        )


def write_json_report(
    path: Path,
    candidate_metadata: dict[str, object],
    baseline_metadata: dict[str, object] | None,
    manifest: Path | None,
    suite: str | None,
    runs: int,
    warmup: int,
    timeout: float,
    results: Sequence[CaseResult],
) -> None:
    def measurement_report(measurements: tuple[Measurement, ...]) -> dict[str, object]:
        wall = summarize([measurement.wall_seconds for measurement in measurements])
        clock = summarize(
            [float(measurement.clock_ticks) for measurement in measurements]
        )
        return {
            "measurements": [
                {"round": index, **asdict(value)}
                for index, value in enumerate(measurements, start=1)
            ],
            "wall_seconds": asdict(wall),
            "clock_ticks": asdict(clock),
        }

    cases = []
    for result in results:
        case_report: dict[str, object] = {
            "name": result.case.name,
            "input": str(result.case.source),
            "expected_assembly_index": result.case.expected_assembly_index,
            "expectation": result.case.expectation,
            "workload": result.case.workload,
            "suites": list(result.case.suites),
            "candidate": measurement_report(result.measurements),
            "baseline": None,
            "comparison": None,
        }
        if result.baseline_measurements:
            wall_speedup = paired_speedup_summary(
                result.measurements,
                result.baseline_measurements,
                "wall_seconds",
            )
            clock_speedup = paired_speedup_summary(
                result.measurements,
                result.baseline_measurements,
                "clock_ticks",
            )
            case_report["baseline"] = measurement_report(result.baseline_measurements)
            case_report["comparison"] = {
                "paired_wall_speedup": (
                    None if wall_speedup is None else asdict(wall_speedup)
                ),
                "paired_clock_speedup": (
                    None if clock_speedup is None else asdict(clock_speedup)
                ),
            }
        cases.append(case_report)

    aggregate_wall, aggregate_clock = aggregate_round_summaries(results)
    baseline_aggregate = None
    comparison = None
    if baseline_metadata is not None:
        baseline_wall, baseline_clock = aggregate_round_summaries(
            results, baseline=True
        )
        wall_speedup = aggregate_paired_speedup_summary(results, "wall_seconds")
        clock_speedup = aggregate_paired_speedup_summary(results, "clock_ticks")
        equal_wall_speedup = equal_weight_paired_speedup_summary(
            results, "wall_seconds"
        )
        equal_clock_speedup = equal_weight_paired_speedup_summary(
            results, "clock_ticks"
        )
        baseline_aggregate = {
            "suite_round_wall_seconds": asdict(baseline_wall),
            "suite_round_clock_ticks": asdict(baseline_clock),
        }
        comparison = {
            "primary_metric": "paired_round_wall_speedup",
            "paired_round_wall_speedup": (
                None if wall_speedup is None else asdict(wall_speedup)
            ),
            "paired_round_clock_speedup": (
                None if clock_speedup is None else asdict(clock_speedup)
            ),
            "equal_weight_wall_geomean_speedup": (
                None if equal_wall_speedup is None else asdict(equal_wall_speedup)
            ),
            "equal_weight_clock_geomean_speedup": (
                None if equal_clock_speedup is None else asdict(equal_clock_speedup)
            ),
        }
    report = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "platform": {
            "description": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "cpu_model": cpu_description(),
            "python": platform.python_version(),
        },
        "executables": {
            "candidate": candidate_metadata,
            "baseline": baseline_metadata,
        },
        "manifest": None if manifest is None else str(manifest),
        "suite": suite,
        "runs": runs,
        "warmup": warmup,
        "timeout_seconds": timeout,
        "schedule": {
            "case_order": "rotated by one position after each round",
            "comparison_order": (
                None
                if baseline_metadata is None
                else "baseline/candidate on odd rounds, candidate/baseline on even rounds"
            ),
            "measurement_array_order": "round order, starting at round 1",
        },
        "cases": cases,
        "candidate_aggregate": {
            "suite_round_wall_seconds": asdict(aggregate_wall),
            "suite_round_clock_ticks": asdict(aggregate_clock),
        },
        "baseline_aggregate": baseline_aggregate,
        "comparison": comparison,
    }
    try:
        with path.open("w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    except OSError as error:
        raise BenchmarkError(f"could not write JSON report {path}: {error}") from error


def print_case_list(cases: Sequence[BenchmarkCase]) -> None:
    print(f"{'Case':<24} {'Suites':<20} {'Expected':>8} {'Status':<12}  Workload")
    print("-" * 95)
    for case in cases:
        print(
            f"{case.name:<24.24} {','.join(case.suites):<20.20} "
            f"{case.expected_assembly_index:>8} "
            f"{case.expectation:<12.12}  {case.workload}"
        )
        print(f"  {case.source}")


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run serial AssemblyCpp calculations and summarize their speed. "
            "Without a corpus option, the default input remains "
            "unitTests/ketoconazole.mol."
        )
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=f"AssemblyCpp executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "--baseline-executable",
        type=Path,
        help=(
            "baseline executable for paired A/B measurements; --executable is "
            "the candidate"
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        help=f"single molfile or native graph (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--expected",
        type=int,
        help=(
            "expected assembly index for single-input mode; defaults to 22 for "
            "the default input and is unchecked for a custom input"
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help=f"benchmark corpus manifest (default with corpus options: {DEFAULT_MANIFEST})",
    )
    parser.add_argument(
        "--suite",
        choices=KNOWN_SUITES,
        help="run one named suite from the benchmark manifest",
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="NAME",
        help="run only this manifest case; may be repeated",
    )
    parser.add_argument(
        "--list-cases",
        action="store_true",
        help="list selected manifest cases and exit without building or running",
    )
    parser.add_argument(
        "--runs",
        type=positive_int,
        help="number of measured rounds per case (default: 5; paired: 6)",
    )
    parser.add_argument(
        "--warmup",
        type=non_negative_int,
        default=1,
        help="number of unmeasured rounds per case (default: 1)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=300.0,
        help="per-calculation timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="write measurements and summaries to this JSON file",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="compile v5/main.cpp for x86-64-v3 before benchmarking",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler command used with --build (default: CXX or c++)",
    )
    return parser


def resolve_requested_cases(
    arguments: argparse.Namespace,
) -> tuple[Path | None, list[BenchmarkCase]]:
    corpus_requested = bool(
        arguments.manifest or arguments.suite or arguments.case or arguments.list_cases
    )
    if arguments.input is not None and corpus_requested:
        raise BenchmarkError(
            "--input cannot be combined with manifest, suite, case, or list options"
        )
    if arguments.expected is not None and corpus_requested:
        raise BenchmarkError("--expected is available only in single-input mode")

    if corpus_requested:
        manifest_path = arguments.manifest or DEFAULT_MANIFEST
        manifest, cases = load_manifest(manifest_path)
        return manifest, select_cases(cases, arguments.suite, arguments.case)

    source = resolve_file(arguments.input or DEFAULT_INPUT, "benchmark input")
    expected = arguments.expected
    if expected is None and source == DEFAULT_INPUT.resolve():
        expected = DEFAULT_EXPECTED_ASSEMBLY_INDEX
    if (
        source == DEFAULT_INPUT.resolve()
        and expected == DEFAULT_EXPECTED_ASSEMBLY_INDEX
    ):
        expectation = "reviewed"
    elif expected is not None:
        expectation = "supplied"
    else:
        expectation = "unchecked"
    return None, [
        BenchmarkCase(
            name=source.name.removesuffix(".mol"),
            source=source,
            expected_assembly_index=expected,
            expectation=expectation,
            suites=(),
            workload="custom input" if arguments.input else "default",
        )
    ]


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)
    runs = arguments.runs
    if runs is None:
        runs = 6 if arguments.baseline_executable is not None else 5

    try:
        manifest, cases = resolve_requested_cases(arguments)
        if arguments.list_cases:
            print_case_list(cases)
            return 0

        baseline_executable = None
        if arguments.baseline_executable is not None:
            baseline_executable = resolve_executable(arguments.baseline_executable)

        executable_path = arguments.executable
        if arguments.build:
            if baseline_executable is not None:
                ensure_distinct_executables(
                    executable_path.resolve(), baseline_executable
                )
            executable_path = build_executable(executable_path, arguments.compiler)
        executable = resolve_executable(executable_path)
        if baseline_executable is not None:
            ensure_distinct_executables(executable, baseline_executable)
            if runs % 2 != 0:
                print(
                    "Warning: an odd number of paired rounds gives one "
                    "executable more first-run positions",
                    file=sys.stderr,
                )
        candidate_metadata = executable_metadata(executable)
        baseline_metadata = (
            None
            if baseline_executable is None
            else executable_metadata(baseline_executable)
        )

        if baseline_executable is None:
            print(f"Executable: {executable}")
        else:
            print(f"Candidate executable: {executable}")
            print(f"Baseline executable: {baseline_executable}")
        if manifest is None:
            print(f"Input: {cases[0].source}")
            if cases[0].expected_assembly_index is None:
                print("Assembly index: checked for presence only")
            else:
                print(f"Expected assembly index: {cases[0].expected_assembly_index}")
        else:
            print(f"Manifest: {manifest}")
            print(f"Suite: {arguments.suite or 'all selected cases'}")
            print(f"Cases: {len(cases)}")
            provisional = [
                case.name for case in cases if case.expectation == "provisional"
            ]
            if provisional:
                print(f"Provisional expectations: {', '.join(provisional)}")
        if len(cases) == 1:
            print(
                f"Runs: {runs} measured, {arguments.warmup} warm-up",
                flush=True,
            )
        else:
            print(
                f"Rounds: {runs} measured, {arguments.warmup} warm-up per case",
                flush=True,
            )

        results = run_benchmarks(
            executable=executable,
            baseline_executable=baseline_executable,
            cases=cases,
            runs=runs,
            warmup=arguments.warmup,
            timeout=arguments.timeout,
        )
        verify_executable_unchanged(executable, candidate_metadata, "candidate")
        if baseline_executable is not None:
            assert baseline_metadata is not None
            verify_executable_unchanged(
                baseline_executable, baseline_metadata, "baseline"
            )
        print_summary(results)
        if arguments.json_output is not None:
            write_json_report(
                path=arguments.json_output,
                candidate_metadata=candidate_metadata,
                baseline_metadata=baseline_metadata,
                manifest=manifest,
                suite=arguments.suite,
                runs=runs,
                warmup=arguments.warmup,
                timeout=arguments.timeout,
                results=results,
            )
            print(f"JSON report: {arguments.json_output}")
    except (BenchmarkError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted", file=sys.stderr)
        return 130

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
