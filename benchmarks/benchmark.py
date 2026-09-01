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
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import NoReturn, TypeGuard, cast

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
ENVIRONMENT_KEY_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
SEARCH_TELEMETRY_PHASES = frozenset(
    (
        "input_setup",
        "initial_enumeration",
        "dag_conversion",
        "assembly_search",
        "output",
    )
)
PARALLEL_TELEMETRY_COUNTERS = frozenset(
    (
        "retained_mask_attempts",
        "retained_masks",
        "duplicate_mask_attempts",
        "rejected_masks",
        "matching_visits",
        "canonicalisation_calls",
        "canonicalisation_mask_cache_hits",
        "canonicalisation_mask_cache_misses",
        "canonical_class_insertions",
        "canonical_class_reuses",
        "vf2_calls",
        "vf2_matches",
        "residual_decomposition_requests",
        "residual_cache_eligible_requests",
        "residual_cache_small_molecule_bypasses",
        "residual_cache_wide_molecule_bypasses",
        "residual_cache_small_residual_bypasses",
        "residual_cache_first_occurrence_bypasses",
        "residual_cache_runtime_disabled_bypasses",
        "residual_cache_lookups",
        "residual_cache_hits",
        "residual_cache_misses",
        "residual_cache_admissions",
        "assembly_cache_lookups",
        "assembly_cache_hits",
        "assembly_cache_misses",
        "assembly_cache_pruned_hits",
        "assembly_cache_updated_hits",
        "pair_bound_cache_lookups",
        "pair_bound_cache_hits",
        "pair_bound_cache_misses",
    )
)
ADAPTIVE_SPLITTING_POLICY = {
    "minimum_queued_tasks_per_worker": 8,
    "target_queued_tasks_per_worker": 16,
    "maximum_queued_tasks_per_worker": 32,
    "maximum_depth": 4,
    "warm_start": "largest_duplicate_first",
}
PARALLEL_SCHEDULER_SUM_FIELDS = (
    "depth_two_tasks_spawned",
    "depth_two_tasks_executed",
    "deeper_tasks_spawned",
    "deeper_tasks_executed",
    "task_steal_attempts",
    "task_steals",
    "local_task_executions",
    "scheduler_idle_waits",
    "scheduler_idle_nanoseconds",
    "deep_refill_activations",
    "proactive_tail_refills",
    "warm_start_branches",
)


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
    telemetry_path: Path
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
    telemetry: dict[str, object] | None = None


@dataclass(frozen=True)
class ExecutionConfig:
    """Launcher prefix and explicit environment overrides for one executable."""

    launcher: tuple[str, ...] = ()
    environment: tuple[tuple[str, str], ...] = ()


def positive_int(value: str) -> int:
    parsed_value = int(value)
    if parsed_value < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return parsed_value


def non_negative_int(value: str) -> int:
    parsed_value = int(value)
    if parsed_value < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return parsed_value


def positive_float(value: str) -> float:
    parsed_value = float(value)
    if not math.isfinite(parsed_value) or parsed_value <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed_value


def environment_assignment(value: str) -> tuple[str, str]:
    """Parse one conventional KEY=VALUE environment assignment."""
    key, separator, setting = value.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("must use KEY=VALUE syntax")
    if not ENVIRONMENT_KEY_PATTERN.fullmatch(key):
        raise argparse.ArgumentTypeError(
            "environment key must match [A-Za-z_][A-Za-z0-9_]*"
        )
    if "\x00" in setting:
        raise argparse.ArgumentTypeError("environment value must not contain NUL")
    return key, setting


def launcher_prefix(value: str) -> tuple[str, ...]:
    """Shell-split a launcher prefix without invoking a shell."""
    try:
        launcher = tuple(shlex.split(value))
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"invalid shell quoting: {error}") from error
    if not launcher:
        raise argparse.ArgumentTypeError("launcher prefix must not be empty")
    return launcher


def create_execution_config(
    launcher: tuple[str, ...] | None,
    environment: Sequence[tuple[str, str]],
    role: str,
) -> ExecutionConfig | None:
    """Validate and normalize one optional command execution configuration."""
    seen: set[str] = set()
    for key, _ in environment:
        if key in seen:
            raise BenchmarkError(f"duplicate {role} environment setting: {key}")
        seen.add(key)

    environment_values = dict(environment)
    launcher_path = environment_values.get("PATH", os.environ.get("PATH"))
    if launcher and shutil.which(launcher[0], path=launcher_path) is None:
        raise BenchmarkError(f"{role} launcher not found: {launcher[0]}")

    config = ExecutionConfig(launcher or (), tuple(environment))
    return None if not config.launcher and not config.environment else config


def execution_config_metadata(config: ExecutionConfig | None) -> dict[str, object]:
    """Return the exact explicit execution configuration used for a role."""
    effective = config or ExecutionConfig()
    return {
        "launcher": list(effective.launcher),
        "environment": dict(effective.environment),
    }


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


def execution_configs_alias(
    candidate: ExecutionConfig | None,
    baseline: ExecutionConfig | None,
) -> bool:
    """Return whether two configurations launch with identical semantics."""
    candidate = candidate or ExecutionConfig()
    baseline = baseline or ExecutionConfig()
    return candidate.launcher == baseline.launcher and dict(
        candidate.environment
    ) == dict(baseline.environment)


def ensure_distinct_executables(
    candidate: Path,
    baseline: Path,
    candidate_execution: ExecutionConfig | None = None,
    baseline_execution: ExecutionConfig | None = None,
) -> None:
    if paths_alias(candidate, baseline) and execution_configs_alias(
        candidate_execution, baseline_execution
    ):
        raise BenchmarkError(
            "candidate and baseline executables resolve to the same file and use "
            "the same execution configuration"
        )


def ensure_distinct_execution_identities(
    candidate_metadata: dict[str, object],
    baseline_metadata: dict[str, object],
    candidate_execution: ExecutionConfig | None,
    baseline_execution: ExecutionConfig | None,
) -> None:
    """Reject byte-identical binaries only when their configurations also alias."""
    if candidate_metadata.get("sha256") == baseline_metadata.get(
        "sha256"
    ) and execution_configs_alias(candidate_execution, baseline_execution):
        raise BenchmarkError(
            "candidate and baseline use the same binary fingerprint and execution "
            "configuration"
        )


def paths_alias(first: Path, second: Path) -> bool:
    """Return whether two paths resolve to, or hardlink, the same file."""
    first = first.resolve()
    second = second.resolve()
    if first == second:
        return True
    try:
        return first.samefile(second)
    except OSError:
        return False


def ensure_json_output_is_distinct(
    output: Path,
    candidate: Path,
    baseline: Path | None,
    telemetry: Path | None,
    manifest: Path | None,
    cases: Sequence[BenchmarkCase],
) -> None:
    """Refuse to overwrite an executable or benchmark corpus source."""
    protected_paths: list[tuple[str, Path]] = [
        ("candidate executable", candidate),
    ]
    if baseline is not None:
        protected_paths.append(("baseline executable", baseline))
    if telemetry is not None:
        protected_paths.append(("telemetry executable", telemetry))
    if manifest is not None:
        protected_paths.append(("benchmark manifest", manifest))
    protected_paths.extend(
        (f"benchmark input for {case.name!r}", case.source) for case in cases
    )

    for description, protected_path in protected_paths:
        if paths_alias(output, protected_path):
            raise BenchmarkError(
                f"--json-output resolves to or aliases the {description}: "
                f"{protected_path}"
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


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file used by the benchmark corpus."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise BenchmarkError(f"could not fingerprint {path}: {error}") from error
    return digest.hexdigest()


def benchmark_corpus_metadata(
    manifest: Path | None,
    cases: Sequence[BenchmarkCase],
) -> dict[str, object]:
    """Fingerprint the manifest and selected inputs represented by a report."""
    return {
        "manifest": (
            None
            if manifest is None
            else {"path": str(manifest), "sha256": file_sha256(manifest)}
        ),
        "inputs": [
            {
                "name": case.name,
                "path": str(case.source),
                "sha256": file_sha256(case.source),
            }
            for case in cases
        ],
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
    manifest_path = resolve_file(path, "benchmark manifest")
    cases: list[BenchmarkCase] = []
    names: set[str] = set()

    try:
        with manifest_path.open(encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != MANIFEST_HEADER:
                raise BenchmarkError(
                    f"invalid header in {manifest_path}: expected {MANIFEST_HEADER}, "
                    f"got {header}"
                )

            for line_number, row in enumerate(reader, start=2):
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != len(MANIFEST_HEADER):
                    raise BenchmarkError(
                        f"invalid row in {manifest_path}:{line_number}: expected "
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
                        f"invalid benchmark name in {manifest_path}:{line_number}: "
                        f"{name!r}"
                    )
                if name in names:
                    raise BenchmarkError(
                        f"duplicate benchmark name in {manifest_path}:{line_number}: "
                        f"{name!r}"
                    )
                names.add(name)

                try:
                    expected = int(expected_text)
                except ValueError as error:
                    raise BenchmarkError(
                        f"invalid assembly index in {manifest_path}:{line_number}: "
                        f"{expected_text!r}"
                    ) from error

                if expectation not in KNOWN_EXPECTATIONS:
                    raise BenchmarkError(
                        f"invalid expectation in {manifest_path}:{line_number}: "
                        f"{expectation!r}"
                    )

                suites = tuple(
                    suite.strip() for suite in suites_text.split(",") if suite.strip()
                )
                if not suites or len(set(suites)) != len(suites):
                    raise BenchmarkError(
                        f"invalid suites in {manifest_path}:{line_number}: "
                        f"{suites_text!r}"
                    )
                unknown_suites = sorted(set(suites) - set(KNOWN_SUITES))
                if unknown_suites:
                    raise BenchmarkError(
                        f"unknown suites in {manifest_path}:{line_number}: "
                        f"{', '.join(unknown_suites)}"
                    )
                if not workload:
                    raise BenchmarkError(
                        f"empty workload in {manifest_path}:{line_number}"
                    )

                source = resolve_manifest_input(Path(input_text), manifest_path)
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
        raise BenchmarkError(f"could not read {manifest_path}: {error}") from error

    if not cases:
        raise BenchmarkError(f"benchmark manifest is empty: {manifest_path}")
    return manifest_path, cases


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


def build_executable(
    executable: Path,
    compiler: str,
    *,
    telemetry: bool = False,
) -> Path:
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
        "-std=c++20",
        "-O3",
        "-mpopcnt",
        "-march=x86-64-v3",
    ]
    if telemetry:
        command.append("-DASSEMBLY_ENABLE_TELEMETRY")
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
        telemetry_name = f"{input_path.name.removesuffix('.mol')}Telemetry.json"
        prepared.append(
            PreparedCase(
                case=case,
                input_name=input_path.name,
                output_path=working_directory / output_name,
                telemetry_path=working_directory / telemetry_name,
                working_directory=working_directory,
            )
        )
    return prepared


def terminate_command(process: subprocess.Popen[str]) -> None:
    """Terminate a configured launch, including its POSIX process group."""
    killed_group = False
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
            killed_group = True
        except ProcessLookupError:
            killed_group = True
        except OSError:
            pass
    if not killed_group:
        process.kill()


def run_command(
    command: Sequence[str],
    working_directory: Path,
    timeout: float,
    environment: dict[str, str] | None,
) -> subprocess.CompletedProcess[str]:
    """Run a command and clean up its process group on timeout or interruption."""
    popen_arguments: dict[str, object] = {
        "cwd": working_directory,
        "env": environment,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        popen_arguments["start_new_session"] = True

    process = subprocess.Popen(command, **popen_arguments)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as error:
        terminate_command(process)
        stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(
            command,
            error.timeout,
            output=stdout,
            stderr=stderr,
        ) from error
    except BaseException:
        # start_new_session deliberately keeps launcher workers out of the
        # driver's foreground process group, so an interrupt must stop them
        # explicitly before it propagates to the caller.
        terminate_command(process)
        process.communicate()
        raise

    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def run_once(
    executable: Path,
    prepared: PreparedCase,
    timeout: float,
    *,
    execution: ExecutionConfig | None = None,
) -> Measurement:
    command = [
        *(execution.launcher if execution is not None else ()),
        str(executable),
        "--pathway=0",
        "--memory-report=0",
        "--write-intermediate-mas=0",
        "--",
        prepared.input_name,
    ]
    environment = None
    if execution is not None and execution.environment:
        environment = os.environ.copy()
        environment.update(execution.environment)
    prepared.output_path.unlink(missing_ok=True)
    started = time.perf_counter()
    try:
        if execution is None:
            completed = subprocess.run(
                command,
                cwd=prepared.working_directory,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        else:
            completed = run_command(
                command,
                prepared.working_directory,
                timeout,
                environment,
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


def parse_search_telemetry(path: Path) -> dict[str, object]:
    def is_nonnegative_integer(value: object) -> TypeGuard[int]:
        return type(value) is int and value >= 0

    def invalid_parallel(detail: str) -> NoReturn:
        raise BenchmarkError(f"invalid parallel telemetry in {path.name}: {detail}")

    def parallel_counters(
        value: object,
        context: str,
    ) -> dict[str, object]:
        if (
            not isinstance(value, dict)
            or not PARALLEL_TELEMETRY_COUNTERS.issubset(value)
            or any(
                not is_nonnegative_integer(value.get(name))
                for name in PARALLEL_TELEMETRY_COUNTERS
            )
        ):
            invalid_parallel(f"invalid {context} counters")
        if value["retained_mask_attempts"] != (
            value["retained_masks"]
            + value["duplicate_mask_attempts"]
            + value["rejected_masks"]
        ):
            invalid_parallel(f"inconsistent {context} retained-mask counters")
        if value["vf2_matches"] > value["vf2_calls"]:
            invalid_parallel(f"inconsistent {context} VF2 counters")
        if value["canonicalisation_calls"] != (
            value["canonicalisation_mask_cache_hits"]
            + value["canonicalisation_mask_cache_misses"]
        ):
            invalid_parallel(f"inconsistent {context} canonical counters")
        if value["canonicalisation_mask_cache_misses"] != (
            value["canonical_class_insertions"] + value["canonical_class_reuses"]
        ):
            invalid_parallel(f"inconsistent {context} canonical-class counters")
        if value["residual_cache_lookups"] != (
            value["residual_cache_hits"] + value["residual_cache_misses"]
        ):
            invalid_parallel(f"inconsistent {context} residual-cache counters")
        if value["residual_decomposition_requests"] != (
            value["residual_cache_eligible_requests"]
            + value["residual_cache_small_molecule_bypasses"]
            + value["residual_cache_wide_molecule_bypasses"]
        ):
            invalid_parallel(f"inconsistent {context} residual request counters")
        if value["residual_cache_eligible_requests"] != (
            value["residual_cache_small_residual_bypasses"]
            + value["residual_cache_first_occurrence_bypasses"]
            + value["residual_cache_runtime_disabled_bypasses"]
            + value["residual_cache_lookups"]
        ):
            invalid_parallel(f"inconsistent {context} residual path counters")
        if value["residual_cache_admissions"] > value["residual_cache_misses"]:
            invalid_parallel(f"inconsistent {context} residual admissions")
        if value["assembly_cache_lookups"] != (
            value["assembly_cache_hits"] + value["assembly_cache_misses"]
        ):
            invalid_parallel(f"inconsistent {context} assembly-cache counters")
        if value["assembly_cache_hits"] != (
            value["assembly_cache_pruned_hits"] + value["assembly_cache_updated_hits"]
        ):
            invalid_parallel(f"inconsistent {context} assembly-cache hits")
        if value["pair_bound_cache_lookups"] != (
            value["pair_bound_cache_hits"] + value["pair_bound_cache_misses"]
        ):
            invalid_parallel(f"inconsistent {context} pair-bound counters")
        return value

    def shared_assembly_cache_counters(
        value: object,
        rank_count: int,
    ) -> dict[str, object]:
        names = (
            "table_count",
            "hits",
            "misses",
            "collision_chain_steps",
            "allocated_bytes",
            "lock_acquisitions",
            "lock_waits",
            "lock_wait_nanoseconds",
        )
        if not isinstance(value, dict) or any(
            not is_nonnegative_integer(value.get(name)) for name in names
        ):
            invalid_parallel("invalid shared assembly-cache counters")
        if value["table_count"] > rank_count:
            invalid_parallel("shared assembly-cache table count exceeds ranks")
        if value["lock_acquisitions"] != value["hits"] + value["misses"]:
            invalid_parallel("inconsistent shared assembly-cache lookups")
        if value["lock_waits"] > value["lock_acquisitions"]:
            invalid_parallel("shared assembly-cache waits exceed acquisitions")
        if (value["misses"] == 0) != (value["allocated_bytes"] == 0):
            invalid_parallel("inconsistent shared assembly-cache allocation")
        if value["table_count"] == 0 and any(
            value[name] != 0 for name in names if name != "table_count"
        ):
            invalid_parallel("disabled shared assembly-cache has activity")
        return value

    def validate_rate(
        container: dict[str, object],
        name: str,
        numerator: int,
        denominator: int,
    ) -> None:
        value = container.get(name)
        if denominator == 0:
            if value is not None:
                raise BenchmarkError(f"invalid cache rate in {path.name}")
            return
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
            or value > 1
            or not math.isclose(
                float(value),
                numerator / denominator,
                rel_tol=5e-6,
                abs_tol=1e-12,
            )
        ):
            raise BenchmarkError(f"invalid cache rate in {path.name}")

    if not path.is_file():
        raise BenchmarkError(f"AssemblyCpp did not create {path.name}")
    try:
        telemetry = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BenchmarkError(f"could not parse {path.name}: {error}") from error
    if (
        not isinstance(telemetry, dict)
        or type(telemetry.get("schema_version")) is not int
        or telemetry["schema_version"] != 1
    ):
        raise BenchmarkError(f"unsupported search telemetry schema in {path.name}")

    required_counters = {
        "retained_mask_attempts",
        "retained_masks",
        "duplicate_mask_attempts",
        "rejected_masks",
        "matching_visits",
        "canonicalisation_calls",
        "vf2_calls",
        "vf2_matches",
    }
    counters = telemetry.get("counters")
    if not isinstance(counters, dict) or not required_counters <= counters.keys():
        raise BenchmarkError(f"missing search counters in {path.name}")
    if any(not is_nonnegative_integer(counters[name]) for name in required_counters):
        raise BenchmarkError(f"invalid search counter in {path.name}")
    if counters["retained_mask_attempts"] != sum(
        counters[name]
        for name in (
            "retained_masks",
            "duplicate_mask_attempts",
            "rejected_masks",
        )
    ):
        raise BenchmarkError(f"inconsistent retained-mask counters in {path.name}")
    if counters["vf2_matches"] > counters["vf2_calls"]:
        raise BenchmarkError(f"inconsistent VF2 counters in {path.name}")

    processed_graph = telemetry.get("processed_graph")
    if not isinstance(processed_graph, dict):
        raise BenchmarkError(f"missing processed graph telemetry in {path.name}")
    for name in ("atoms", "edges", "active_mask_words"):
        value = processed_graph.get(name)
        if not is_nonnegative_integer(value):
            raise BenchmarkError(f"invalid processed graph telemetry in {path.name}")
    expected_active_words = (processed_graph["edges"] + 63) // 64
    if processed_graph["active_mask_words"] != expected_active_words:
        raise BenchmarkError(f"inconsistent processed mask width in {path.name}")

    caches = telemetry.get("caches")
    required_caches = {
        "canonical_mask",
        "canonical_class",
        "residual_decomposition",
        "assembly_state",
        "pair_bound",
    }
    if not isinstance(caches, dict) or not required_caches <= caches.keys():
        raise BenchmarkError(f"missing cache telemetry in {path.name}")
    canonical = caches["canonical_mask"]
    canonical_class = caches["canonical_class"]
    residual = caches["residual_decomposition"]
    assembly_cache = caches["assembly_state"]
    pair_bound_cache = caches["pair_bound"]
    if not all(
        isinstance(cache, dict)
        for cache in (
            canonical,
            canonical_class,
            residual,
            assembly_cache,
            pair_bound_cache,
        )
    ):
        raise BenchmarkError(f"invalid cache telemetry in {path.name}")
    residual_counter_names = (
        "requests",
        "eligible_requests",
        "small_molecule_bypasses",
        "wide_molecule_bypasses",
        "small_residual_bypasses",
        "first_occurrence_bypasses",
        "runtime_disabled_bypasses",
        "lookups",
        "hits",
        "misses",
        "admissions",
    )
    for cache, names in (
        (canonical, ("hits", "misses")),
        (canonical_class, ("insertions", "reuses")),
        (residual, residual_counter_names),
        (
            assembly_cache,
            ("lookups", "hits", "misses", "pruned_hits", "updated_hits"),
        ),
        (pair_bound_cache, ("lookups", "hits", "misses")),
    ):
        if any(not is_nonnegative_integer(cache.get(name)) for name in names):
            raise BenchmarkError(f"invalid cache counter in {path.name}")
    if counters["canonicalisation_calls"] != (canonical["hits"] + canonical["misses"]):
        raise BenchmarkError(f"inconsistent canonical cache counters in {path.name}")
    if residual["lookups"] != (residual["hits"] + residual["misses"]):
        raise BenchmarkError(f"inconsistent residual cache counters in {path.name}")
    if type(residual.get("eligible_for_processed_graph")) is not bool:
        raise BenchmarkError(f"invalid residual cache eligibility in {path.name}")
    expected_residual_cache_eligibility = processed_graph["edges"] >= 31
    if residual["eligible_for_processed_graph"] != expected_residual_cache_eligibility:
        raise BenchmarkError(f"inconsistent residual cache eligibility in {path.name}")
    if residual["requests"] != (
        residual["eligible_requests"]
        + residual["small_molecule_bypasses"]
        + residual["wide_molecule_bypasses"]
    ):
        raise BenchmarkError(f"inconsistent residual request counters in {path.name}")
    if residual["eligible_for_processed_graph"]:
        if (
            residual["eligible_requests"] != residual["requests"]
            or residual["small_molecule_bypasses"] != 0
            or residual["wide_molecule_bypasses"] != 0
        ):
            raise BenchmarkError(f"inconsistent residual eligibility in {path.name}")
    elif (
        residual["eligible_requests"] != 0
        or residual["small_molecule_bypasses"] + residual["wide_molecule_bypasses"]
        != residual["requests"]
    ):
        raise BenchmarkError(f"inconsistent residual bypasses in {path.name}")
    if residual["eligible_requests"] != (
        residual["small_residual_bypasses"]
        + residual["first_occurrence_bypasses"]
        + residual["runtime_disabled_bypasses"]
        + residual["lookups"]
    ):
        raise BenchmarkError(f"inconsistent residual path counters in {path.name}")
    if residual["admissions"] > residual["misses"]:
        raise BenchmarkError(f"inconsistent residual admissions in {path.name}")
    if canonical["misses"] != (
        canonical_class["insertions"] + canonical_class["reuses"]
    ):
        raise BenchmarkError(f"inconsistent canonical class counters in {path.name}")
    for label, cache in (
        ("assembly", assembly_cache),
        ("pair-bound", pair_bound_cache),
    ):
        if cache["lookups"] != cache["hits"] + cache["misses"]:
            raise BenchmarkError(f"inconsistent {label} cache counters in {path.name}")
    if (
        assembly_cache["pruned_hits"] + assembly_cache["updated_hits"]
        != assembly_cache["hits"]
    ):
        raise BenchmarkError(f"inconsistent assembly hit counters in {path.name}")

    validate_rate(
        canonical,
        "hit_rate",
        canonical["hits"],
        canonical["hits"] + canonical["misses"],
    )
    validate_rate(
        canonical_class,
        "reuse_rate",
        canonical_class["reuses"],
        canonical_class["insertions"] + canonical_class["reuses"],
    )
    validate_rate(
        residual,
        "lookup_hit_rate",
        residual["hits"],
        residual["hits"] + residual["misses"],
    )
    validate_rate(
        residual,
        "request_hit_rate",
        residual["hits"],
        residual["requests"],
    )
    validate_rate(
        assembly_cache,
        "hit_rate",
        assembly_cache["hits"],
        assembly_cache["hits"] + assembly_cache["misses"],
    )
    validate_rate(
        pair_bound_cache,
        "hit_rate",
        pair_bound_cache["hits"],
        pair_bound_cache["hits"] + pair_bound_cache["misses"],
    )

    memory = telemetry.get("memory")
    if not isinstance(memory, dict):
        raise BenchmarkError(f"missing phase memory telemetry in {path.name}")
    phases = memory.get("phases")
    required_phases = {
        "input_setup",
        "initial_enumeration",
        "dag_conversion",
        "assembly_search",
        "output",
    }
    if not isinstance(phases, dict) or set(phases) != required_phases:
        raise BenchmarkError(f"missing phase memory telemetry in {path.name}")
    if (
        not isinstance(memory.get("method"), str)
        or type(memory.get("phase_peaks_are_absolute_not_additive")) is not bool
        or not memory["phase_peaks_are_absolute_not_additive"]
        or type(memory.get("phase_peaks_complete")) is not bool
    ):
        raise BenchmarkError(f"invalid phase memory metadata in {path.name}")

    memory_value_names = (
        "start_rss_kib",
        "peak_rss_kib",
        "end_rss_kib",
        "start_virtual_kib",
        "end_virtual_kib",
    )
    for phase_name in required_phases:
        phase = phases[phase_name]
        if not isinstance(phase, dict):
            raise BenchmarkError(f"invalid phase telemetry in {path.name}")
        if not is_nonnegative_integer(phase.get("clock_ticks")) or not (
            is_nonnegative_integer(phase.get("activations"))
        ):
            raise BenchmarkError(f"invalid phase counter in {path.name}")
        if "wall_nanoseconds" in phase and not is_nonnegative_integer(
            phase["wall_nanoseconds"]
        ):
            raise BenchmarkError(f"invalid phase wall time in {path.name}")
        if any(
            value is not None and not is_nonnegative_integer(value)
            for value in (phase.get(name) for name in memory_value_names)
        ):
            raise BenchmarkError(f"invalid phase memory value in {path.name}")
        start = phase.get("start_rss_kib")
        peak = phase.get("peak_rss_kib")
        end = phase.get("end_rss_kib")
        if peak is not None and (
            not is_nonnegative_integer(start)
            or not is_nonnegative_integer(end)
            or peak < start
            or peak < end
        ):
            raise BenchmarkError(f"invalid phase RSS peak in {path.name}")
        if phase["activations"] == 0 and any(
            phase.get(name) is not None for name in memory_value_names
        ):
            raise BenchmarkError(f"inactive phase has memory data in {path.name}")

    overall_peak = memory.get("overall_peak_resident_kib")
    process_virtual_peak = memory.get("process_peak_virtual_kib")
    if process_virtual_peak is not None and not is_nonnegative_integer(
        process_virtual_peak
    ):
        raise BenchmarkError(f"invalid process memory telemetry in {path.name}")
    if memory["phase_peaks_complete"]:
        activated_peaks = [
            phase["peak_rss_kib"]
            for phase in phases.values()
            if phase["activations"] > 0
        ]
        if (
            not activated_peaks
            or any(not is_nonnegative_integer(peak) for peak in activated_peaks)
            or not is_nonnegative_integer(overall_peak)
            or overall_peak != max(activated_peaks)
        ):
            raise BenchmarkError(f"invalid overall RSS peak in {path.name}")
    elif overall_peak is not None:
        raise BenchmarkError(f"partial phase memory has an overall peak in {path.name}")

    if "parallel" not in telemetry:
        return telemetry
    parallel = telemetry["parallel"]
    if not isinstance(parallel, dict):
        invalid_parallel("expected an object")

    if (
        memory["method"] != "disabled_parallel"
        or memory["phase_peaks_complete"]
        or overall_peak is not None
        or process_virtual_peak is not None
        or any(
            phase.get(name) is not None
            for phase in phases.values()
            for name in memory_value_names
        )
    ):
        invalid_parallel("parallel phase memory must be disabled")

    mode = parallel.get("mode")
    scope = parallel.get("aggregation_scope")
    if mode not in {"openmp", "mpi", "hybrid"}:
        invalid_parallel("invalid mode")
    if scope not in {"process", "all_mpi_ranks"}:
        invalid_parallel("invalid aggregation scope")
    if parallel.get("busy_timing_method") != "elapsed_minus_scheduler_idle_time":
        invalid_parallel("invalid busy timing method")
    if parallel.get("elapsed_timing_method") != "parallel_region_steady_clock":
        invalid_parallel("invalid elapsed timing method")
    if parallel.get("enabled") is not True:
        invalid_parallel("parallel telemetry is not enabled")

    rank_count = parallel.get("rank_count")
    worker_count = parallel.get("worker_count")
    if not is_nonnegative_integer(rank_count) or rank_count == 0:
        invalid_parallel("invalid rank count")
    if not is_nonnegative_integer(worker_count) or worker_count == 0:
        invalid_parallel("invalid worker count")
    if mode == "openmp":
        if scope != "process" or rank_count != 1:
            invalid_parallel("OpenMP telemetry must describe one process")
    elif scope != "all_mpi_ranks":
        invalid_parallel("MPI telemetry must aggregate all ranks")

    threads_per_rank = parallel.get("local_threads_per_rank")
    if (
        not isinstance(threads_per_rank, list)
        or len(threads_per_rank) != rank_count
        or any(
            not is_nonnegative_integer(threads) or threads == 0
            for threads in threads_per_rank
        )
        or sum(threads_per_rank) != worker_count
    ):
        invalid_parallel("inconsistent local thread counts")
    if mode == "mpi" and any(threads != 1 for threads in threads_per_rank):
        invalid_parallel("MPI telemetry must report one thread per rank")
    uniform_threads = threads_per_rank[0] if len(set(threads_per_rank)) == 1 else None
    local_threads = parallel.get("local_threads")
    if (
        local_threads is not None and not is_nonnegative_integer(local_threads)
    ) or local_threads != uniform_threads:
        invalid_parallel("inconsistent uniform local thread count")

    has_static_shards = "shard_ownership" in parallel
    has_dynamic_leases = "branch_scheduler" in parallel
    if has_static_shards == has_dynamic_leases:
        invalid_parallel("expected exactly one branch scheduling description")
    lease_size: int | None = None
    scheduler_complete = True
    uses_distributed_root_queue = False
    uses_legacy_rank_partition = False
    if has_static_shards:
        shard_ownership = parallel["shard_ownership"]
        if not isinstance(shard_ownership, dict):
            invalid_parallel("missing shard ownership")
        if (
            shard_ownership.get("strategy") != "root_branch_ordinal_modulo_worker_count"
            or type(shard_ownership.get("complete")) is not bool
            or not shard_ownership["complete"]
            or shard_ownership.get("shard_count") != worker_count
        ):
            invalid_parallel("inconsistent shard ownership")
    else:
        branch_scheduler = parallel["branch_scheduler"]
        if not isinstance(branch_scheduler, dict):
            invalid_parallel("missing branch scheduler")
        strategy = branch_scheduler.get("strategy")
        lease_size = branch_scheduler.get("lease_size")
        scheduler_complete = branch_scheduler.get("complete")
        adaptive_splitting = branch_scheduler.get("adaptive_splitting")
        uses_distributed_root_queue = strategy == "distributed_global_root_queue"
        uses_legacy_rank_partition = (
            strategy == "dynamic_leases_with_static_mpi_rank_partition"
        )
        scheduling_metadata_valid = False
        if uses_distributed_root_queue:
            root_queue = branch_scheduler.get("root_queue")
            scheduling_metadata_valid = (
                isinstance(root_queue, dict)
                and is_nonnegative_integer(root_queue.get("participant_count"))
                and root_queue.get("participant_count") == rank_count
            )
        elif uses_legacy_rank_partition:
            scheduling_metadata_valid = (
                is_nonnegative_integer(branch_scheduler.get("rank_partition_count"))
                and branch_scheduler.get("rank_partition_count") == rank_count
            )
        if (
            not scheduling_metadata_valid
            or not is_nonnegative_integer(lease_size)
            or lease_size == 0
            or type(scheduler_complete) is not bool
            or not isinstance(adaptive_splitting, dict)
            or any(
                (
                    not is_nonnegative_integer(adaptive_splitting.get(name))
                    or adaptive_splitting.get(name) != expected
                )
                if isinstance(expected, int)
                else adaptive_splitting.get(name) != expected
                for name, expected in ADAPTIVE_SPLITTING_POLICY.items()
            )
        ):
            invalid_parallel("inconsistent branch scheduler")
    if type(parallel.get("branch_scan_complete")) is not bool:
        invalid_parallel("invalid branch scan status")

    aggregate = parallel.get("aggregate")
    if not isinstance(aggregate, dict):
        invalid_parallel("missing aggregate")
    aggregate_counters = parallel_counters(
        aggregate.get("counters"),
        "aggregate",
    )
    # Schema v1 is additive and historical documents predate this object.
    if "shared_assembly_cache" in aggregate:
        shared_assembly_cache_counters(
            aggregate["shared_assembly_cache"],
            rank_count,
        )
    aggregate_integer_names = (
        "branch_assignments",
        "elapsed_nanoseconds",
        "worker_elapsed_nanoseconds",
        "worker_busy_nanoseconds",
        "task_queue_high_watermark",
        "maximum_task_depth_executed",
        *PARALLEL_SCHEDULER_SUM_FIELDS,
    )
    if any(
        not is_nonnegative_integer(aggregate.get(name))
        for name in aggregate_integer_names
    ):
        invalid_parallel("invalid aggregate measurement")
    aggregate_branch_candidates = aggregate.get("branch_candidates")
    if aggregate_branch_candidates is not None and not is_nonnegative_integer(
        aggregate_branch_candidates
    ):
        invalid_parallel("invalid aggregate branch count")
    if aggregate["worker_busy_nanoseconds"] > aggregate["worker_elapsed_nanoseconds"]:
        invalid_parallel("aggregate busy time exceeds elapsed worker time")
    aggregate_branch_leases = aggregate.get("branch_leases")
    if has_dynamic_leases and not is_nonnegative_integer(aggregate_branch_leases):
        invalid_parallel("invalid aggregate branch lease count")

    workers = parallel.get("workers")
    if not isinstance(workers, list) or len(workers) != worker_count:
        invalid_parallel("worker count does not match worker records")
    worker_counters = []
    worker_branch_candidates = []
    rank_local_ids: list[set[int]] = [set() for _ in range(rank_count)]
    global_ids = set()
    shard_ids = set()
    rank_branch_assignments = [0] * rank_count
    total_branch_leases = 0
    total_branch_assignments = 0
    total_elapsed = 0
    total_busy = 0
    maximum_elapsed = 0
    scheduler_sums = dict.fromkeys(PARALLEL_SCHEDULER_SUM_FIELDS, 0)
    maximum_task_queue_high_watermark = 0
    maximum_task_depth_executed = 0
    rank_offsets = []
    offset = 0
    for threads in threads_per_rank:
        rank_offsets.append(offset)
        offset += threads

    for worker_index, worker in enumerate(workers):
        if not isinstance(worker, dict):
            invalid_parallel(f"worker {worker_index} is not an object")
        rank = worker.get("rank")
        local_worker_index = worker.get("local_worker_index")
        global_worker_index = worker.get("global_worker_index")
        if (
            not is_nonnegative_integer(rank)
            or rank >= rank_count
            or not is_nonnegative_integer(local_worker_index)
            or local_worker_index >= threads_per_rank[rank]
            or not is_nonnegative_integer(global_worker_index)
            or global_worker_index >= worker_count
            or global_worker_index != rank_offsets[rank] + local_worker_index
        ):
            invalid_parallel(f"invalid worker identity at record {worker_index}")
        rank_local_ids[rank].add(local_worker_index)
        global_ids.add(global_worker_index)

        branch_leases = 0
        if has_static_shards:
            shard = worker.get("shard")
            if (
                not isinstance(shard, dict)
                or shard.get("index") != global_worker_index
                or shard.get("count") != worker_count
            ):
                invalid_parallel(f"invalid worker shard at record {worker_index}")
            shard_ids.add(shard["index"])
        else:
            branch_leases = worker.get("branch_leases")
            if uses_distributed_root_queue:
                root_queue = worker.get("root_queue")
                if (
                    not isinstance(root_queue, dict)
                    or not is_nonnegative_integer(root_queue.get("participant_rank"))
                    or root_queue.get("participant_rank") != rank
                    or not is_nonnegative_integer(root_queue.get("participant_count"))
                    or root_queue.get("participant_count") != rank_count
                ):
                    invalid_parallel(
                        "invalid worker root queue participant at record "
                        f"{worker_index}"
                    )
            else:
                rank_partition = worker.get("rank_partition")
                if (
                    not isinstance(rank_partition, dict)
                    or not is_nonnegative_integer(rank_partition.get("index"))
                    or rank_partition.get("index") != rank
                    or not is_nonnegative_integer(rank_partition.get("count"))
                    or rank_partition.get("count") != rank_count
                ):
                    invalid_parallel(
                        f"invalid worker rank partition at record {worker_index}"
                    )
            if not is_nonnegative_integer(branch_leases):
                invalid_parallel(f"invalid worker lease count at record {worker_index}")

        branch_candidates = worker.get("branch_candidates")
        branch_assignments = worker.get("branch_assignments")
        elapsed = worker.get("elapsed_nanoseconds")
        busy = worker.get("busy_nanoseconds")
        worker_scheduler_values = {
            name: worker.get(name) for name in PARALLEL_SCHEDULER_SUM_FIELDS
        }
        task_queue_high_watermark = worker.get("task_queue_high_watermark")
        task_depth = worker.get("maximum_task_depth_executed")
        if any(
            not is_nonnegative_integer(value)
            for value in (
                branch_candidates,
                branch_assignments,
                elapsed,
                busy,
                task_queue_high_watermark,
                task_depth,
                *worker_scheduler_values.values(),
            )
        ):
            invalid_parallel(f"invalid worker measurement at record {worker_index}")
        if branch_assignments > branch_candidates:
            invalid_parallel(
                f"worker assignments exceed branches at record {worker_index}"
            )
        if has_dynamic_leases and (
            (branch_assignments == 0 and branch_leases != 0)
            or (
                branch_assignments > 0
                and (
                    branch_leases == 0
                    or branch_leases > branch_assignments
                    or branch_assignments > branch_leases * lease_size
                )
            )
        ):
            invalid_parallel(f"inconsistent worker leases at record {worker_index}")
        if busy > elapsed:
            invalid_parallel(
                f"worker busy time exceeds elapsed time at record {worker_index}"
            )
        task_steal_attempts = worker_scheduler_values["task_steal_attempts"]
        task_steals = worker_scheduler_values["task_steals"]
        local_task_executions = worker_scheduler_values["local_task_executions"]
        transferred_tasks_executed = (
            worker_scheduler_values["depth_two_tasks_executed"]
            + worker_scheduler_values["deeper_tasks_executed"]
        )
        if task_steals > task_steal_attempts:
            invalid_parallel(
                f"worker task steals exceed attempts at record {worker_index}"
            )
        if local_task_executions + task_steals != transferred_tasks_executed:
            invalid_parallel(
                "worker task executions do not match transferred tasks at "
                f"record {worker_index}"
            )
        if not (
            (transferred_tasks_executed == 0 and task_depth == 0)
            or (
                transferred_tasks_executed > 0
                and 2 <= task_depth <= ADAPTIVE_SPLITTING_POLICY["maximum_depth"]
            )
        ):
            invalid_parallel(
                f"invalid worker maximum task depth at record {worker_index}"
            )
        if worker_scheduler_values["warm_start_branches"] > 1:
            invalid_parallel(
                f"invalid worker warm start count at record {worker_index}"
            )
        if busy != elapsed - min(
            worker_scheduler_values["scheduler_idle_nanoseconds"], elapsed
        ):
            invalid_parallel(f"inconsistent worker busy time at record {worker_index}")
        worker_branch_candidates.append(branch_candidates)
        rank_branch_assignments[rank] += branch_assignments
        total_branch_leases += branch_leases
        total_branch_assignments += branch_assignments
        total_elapsed += elapsed
        total_busy += busy
        maximum_elapsed = max(maximum_elapsed, elapsed)
        for name, value in worker_scheduler_values.items():
            scheduler_sums[name] += value
        maximum_task_queue_high_watermark = max(
            maximum_task_queue_high_watermark,
            task_queue_high_watermark,
        )
        maximum_task_depth_executed = max(
            maximum_task_depth_executed,
            task_depth,
        )

        graph = worker.get("processed_graph")
        if not isinstance(graph, dict):
            invalid_parallel(f"missing worker graph at record {worker_index}")
        if (
            any(
                not is_nonnegative_integer(graph.get(name))
                for name in ("atoms", "edges", "active_mask_words")
            )
            or type(graph.get("residual_cache_eligible")) is not bool
        ):
            invalid_parallel(f"invalid worker graph at record {worker_index}")
        if (
            graph["atoms"] != processed_graph["atoms"]
            or graph["edges"] != processed_graph["edges"]
            or graph["active_mask_words"] != processed_graph["active_mask_words"]
            or graph["residual_cache_eligible"]
            != residual["eligible_for_processed_graph"]
        ):
            invalid_parallel(f"inconsistent worker graph at record {worker_index}")

        worker_phases = worker.get("phases")
        if (
            not isinstance(worker_phases, dict)
            or set(worker_phases) != SEARCH_TELEMETRY_PHASES
        ):
            invalid_parallel(f"invalid worker phases at record {worker_index}")
        for phase in worker_phases.values():
            if not isinstance(phase, dict) or any(
                not is_nonnegative_integer(phase.get(name))
                for name in ("wall_nanoseconds", "activations")
            ):
                invalid_parallel(f"invalid worker phase at record {worker_index}")
            if phase["activations"] == 0 and phase["wall_nanoseconds"] != 0:
                invalid_parallel(
                    f"inactive worker phase has time at record {worker_index}"
                )

        worker_counters.append(
            parallel_counters(worker.get("counters"), f"worker {worker_index}")
        )

    if global_ids != set(range(worker_count)) or (
        has_static_shards and shard_ids != set(range(worker_count))
    ):
        invalid_parallel("worker or shard indexes are not contiguous")
    if any(
        rank_local_ids[rank] != set(range(threads_per_rank[rank]))
        for rank in range(rank_count)
    ):
        invalid_parallel("local worker indexes are not contiguous")
    if aggregate["branch_assignments"] != total_branch_assignments:
        invalid_parallel("aggregate branch assignments do not match workers")
    if has_dynamic_leases and aggregate_branch_leases != total_branch_leases:
        invalid_parallel("aggregate branch leases do not match workers")
    if aggregate["worker_elapsed_nanoseconds"] != total_elapsed:
        invalid_parallel("aggregate worker elapsed time does not match workers")
    if aggregate["worker_busy_nanoseconds"] != total_busy:
        invalid_parallel("aggregate worker busy time does not match workers")
    if aggregate["elapsed_nanoseconds"] < maximum_elapsed:
        invalid_parallel("aggregate elapsed time is shorter than a worker")
    for name, worker_sum in scheduler_sums.items():
        if aggregate[name] != worker_sum:
            invalid_parallel(f"aggregate scheduler field {name} does not match workers")
    if aggregate["task_queue_high_watermark"] != maximum_task_queue_high_watermark:
        invalid_parallel("aggregate task queue high-water mark does not match workers")
    if aggregate["maximum_task_depth_executed"] != maximum_task_depth_executed:
        invalid_parallel("aggregate maximum task depth does not match workers")
    aggregate_depth_two_spawned = aggregate["depth_two_tasks_spawned"]
    aggregate_depth_two_executed = aggregate["depth_two_tasks_executed"]
    aggregate_deeper_spawned = aggregate["deeper_tasks_spawned"]
    aggregate_deeper_executed = aggregate["deeper_tasks_executed"]
    aggregate_transferred_tasks_executed = (
        aggregate_depth_two_executed + aggregate_deeper_executed
    )
    if aggregate["task_steals"] > aggregate["task_steal_attempts"]:
        invalid_parallel("aggregate task steals exceed attempts")
    if (
        aggregate["local_task_executions"] + aggregate["task_steals"]
        != aggregate_transferred_tasks_executed
    ):
        invalid_parallel("aggregate task executions do not match transferred tasks")
    if aggregate_depth_two_spawned != aggregate_depth_two_executed:
        invalid_parallel("spawned depth-two tasks were not each executed once")
    if (
        aggregate_depth_two_spawned + aggregate_deeper_spawned
        != aggregate_transferred_tasks_executed
    ):
        invalid_parallel("spawned transferred tasks were not each executed once")
    aggregate_maximum_task_depth = aggregate["maximum_task_depth_executed"]
    if not (
        (
            aggregate_transferred_tasks_executed == 0
            and aggregate_maximum_task_depth == 0
        )
        or (
            aggregate_transferred_tasks_executed > 0
            and 2
            <= aggregate_maximum_task_depth
            <= ADAPTIVE_SPLITTING_POLICY["maximum_depth"]
        )
    ):
        invalid_parallel("invalid aggregate maximum task depth")
    for name in PARALLEL_TELEMETRY_COUNTERS:
        if aggregate_counters[name] != sum(
            counters[name] for counters in worker_counters
        ):
            invalid_parallel(f"aggregate counter {name} does not match workers")

    candidates_agree = len(set(worker_branch_candidates)) == 1
    expected_branch_candidates = (
        worker_branch_candidates[0] if candidates_agree else None
    )
    if aggregate_branch_candidates != expected_branch_candidates:
        invalid_parallel("aggregate branch count does not match workers")
    if parallel["branch_scan_complete"] and (
        not candidates_agree
        or not scheduler_complete
        or total_branch_assignments != expected_branch_candidates
    ):
        invalid_parallel("complete branch scan has incomplete assignments")
    if uses_legacy_rank_partition and parallel["branch_scan_complete"]:
        expected_branch_candidates = cast(int, expected_branch_candidates)
        for rank, assignments in enumerate(rank_branch_assignments):
            expected_assignments = (
                0
                if expected_branch_candidates <= rank
                else 1 + (expected_branch_candidates - 1 - rank) // rank_count
            )
            if assignments != expected_assignments:
                invalid_parallel(
                    f"rank {rank} branch assignments do not match partition"
                )
    legacy_counters = {
        "retained_mask_attempts": counters["retained_mask_attempts"],
        "retained_masks": counters["retained_masks"],
        "duplicate_mask_attempts": counters["duplicate_mask_attempts"],
        "rejected_masks": counters["rejected_masks"],
        "matching_visits": counters["matching_visits"],
        "canonicalisation_calls": counters["canonicalisation_calls"],
        "canonicalisation_mask_cache_hits": canonical["hits"],
        "canonicalisation_mask_cache_misses": canonical["misses"],
        "canonical_class_insertions": canonical_class["insertions"],
        "canonical_class_reuses": canonical_class["reuses"],
        "vf2_calls": counters["vf2_calls"],
        "vf2_matches": counters["vf2_matches"],
        "residual_decomposition_requests": residual["requests"],
        "residual_cache_eligible_requests": residual["eligible_requests"],
        "residual_cache_small_molecule_bypasses": residual["small_molecule_bypasses"],
        "residual_cache_wide_molecule_bypasses": residual["wide_molecule_bypasses"],
        "residual_cache_small_residual_bypasses": residual["small_residual_bypasses"],
        "residual_cache_first_occurrence_bypasses": residual[
            "first_occurrence_bypasses"
        ],
        "residual_cache_runtime_disabled_bypasses": residual[
            "runtime_disabled_bypasses"
        ],
        "residual_cache_lookups": residual["lookups"],
        "residual_cache_hits": residual["hits"],
        "residual_cache_misses": residual["misses"],
        "residual_cache_admissions": residual["admissions"],
        "assembly_cache_lookups": assembly_cache["lookups"],
        "assembly_cache_hits": assembly_cache["hits"],
        "assembly_cache_misses": assembly_cache["misses"],
        "assembly_cache_pruned_hits": assembly_cache["pruned_hits"],
        "assembly_cache_updated_hits": assembly_cache["updated_hits"],
        "pair_bound_cache_lookups": pair_bound_cache["lookups"],
        "pair_bound_cache_hits": pair_bound_cache["hits"],
        "pair_bound_cache_misses": pair_bound_cache["misses"],
    }
    if any(
        aggregate_counters[name] != legacy_counters[name]
        for name in PARALLEL_TELEMETRY_COUNTERS
    ):
        invalid_parallel("aggregate counters do not match legacy telemetry")
    return telemetry


def run_telemetry_once(
    executable: Path,
    prepared: PreparedCase,
    timeout: float,
    *,
    execution: ExecutionConfig | None = None,
) -> dict[str, object]:
    command = [
        *(execution.launcher if execution is not None else ()),
        str(executable),
        "--pathway=0",
        "--memory-report=0",
        "--write-intermediate-mas=0",
        "--telemetry=1",
        "--",
        prepared.input_name,
    ]
    environment = None
    if execution is not None and execution.environment:
        environment = os.environ.copy()
        environment.update(execution.environment)
    prepared.output_path.unlink(missing_ok=True)
    prepared.telemetry_path.unlink(missing_ok=True)
    try:
        if execution is None:
            completed = subprocess.run(
                command,
                cwd=prepared.working_directory,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        else:
            completed = run_command(
                command,
                prepared.working_directory,
                timeout,
                environment,
            )
    except subprocess.TimeoutExpired as error:
        raise BenchmarkError(
            f"{prepared.case.name} telemetry run timed out after "
            f"{error.timeout:g} seconds"
        ) from error
    except OSError as error:
        raise BenchmarkError(
            f"could not run AssemblyCpp telemetry for {prepared.case.name}: {error}"
        ) from error
    try:
        parse_measurement(
            prepared.output_path,
            completed,
            0.0,
            prepared.case.expected_assembly_index,
        )
        return parse_search_telemetry(prepared.telemetry_path)
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
    telemetry_executable: Path | None = None,
    candidate_execution: ExecutionConfig | None = None,
    baseline_execution: ExecutionConfig | None = None,
) -> list[CaseResult]:
    with tempfile.TemporaryDirectory(prefix="assemblycpp-benchmark-") as directory:
        prepared_cases = prepare_cases(cases, Path(directory))
        measurements: dict[str, list[Measurement]] = {case.name: [] for case in cases}
        baseline_measurements: dict[str, list[Measurement]] = {
            case.name: [] for case in cases
        }
        telemetry: dict[str, dict[str, object] | None] = {
            case.name: None for case in cases
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

        def execute(
            role: str,
            current_executable: Path,
            prepared: PreparedCase,
        ) -> Measurement:
            execution = (
                candidate_execution if role == "candidate" else baseline_execution
            )
            # Preserve the historical three-argument call for default callers and
            # monkeypatched PGO/benchmark test doubles.
            if execution is None:
                return run_once(current_executable, prepared, timeout)
            return run_once(
                current_executable,
                prepared,
                timeout,
                execution=execution,
            )

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
                        paired_warmups[role] = execute(
                            role, current_executable, prepared
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
                        measurement = execute(role, current_executable, prepared)
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
                        f"assembly index {measurement.assembly_index}",
                        flush=True,
                    )

                validate_unchecked_pair(prepared, paired_measurements)

        if telemetry_executable is not None:
            for index, prepared in enumerate(prepared_cases, start=1):
                print(
                    f"Telemetry [{index}/{case_count}] {prepared.case.name}...",
                    flush=True,
                )
                if candidate_execution is None:
                    telemetry[prepared.case.name] = run_telemetry_once(
                        telemetry_executable,
                        prepared,
                        timeout,
                    )
                else:
                    telemetry[prepared.case.name] = run_telemetry_once(
                        telemetry_executable,
                        prepared,
                        timeout,
                        execution=candidate_execution,
                    )

    return [
        CaseResult(
            case=case,
            measurements=tuple(measurements[case.name]),
            baseline_measurements=tuple(baseline_measurements[case.name]),
            telemetry=telemetry[case.name],
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
        20, min(44, max(len(result.case.workload) for result in results))
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
    print("\nComparison (speedup > 1.0 = candidate faster)")
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
        wall_speedup = cast(MetricSummary, wall_speedup)
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
    round_wall = cast(MetricSummary, round_wall)
    equal_wall = cast(MetricSummary, equal_wall)
    print("\nCorpus comparison")
    print(
        f"  paired suite wall (primary):       {round_wall.median:.4f}x "
        f"(MAD {round_wall.mad:.4f}x)"
    )
    if round_clock is not None:
        print(
            f"  paired suite clock (primary):     {round_clock.median:.4f}x "
            f"(MAD {round_clock.mad:.4f}x)"
        )
    print(
        f"  equal-weight wall (descriptive):    {equal_wall.median:.4f}x "
        f"(MAD {equal_wall.mad:.4f}x)"
    )
    if equal_clock is not None:
        print(
            f"  equal-weight clock (descriptive): {equal_clock.median:.4f}x "
            f"(MAD {equal_clock.mad:.4f}x)"
        )


def print_telemetry_summary(results: Sequence[CaseResult]) -> None:
    if not any(result.telemetry is not None for result in results):
        return
    print("\nUntimed search telemetry")
    name_width = max(12, min(24, max(len(result.case.name) for result in results)))
    print(
        f"  {'Case':<{name_width}} {'Masks':>10} {'Matches':>12} "
        f"{'Canon':>12} {'VF2':>8} {'Canon hit':>10} "
        f"{'Residual hit':>13} {'Peak RSS':>11}"
    )
    print("  " + "-" * (name_width + 82))
    for result in results:
        telemetry = result.telemetry
        if telemetry is None:
            continue
        counters = telemetry["counters"]
        caches = telemetry["caches"]
        memory = telemetry["memory"]
        counters = cast(dict[str, object], counters)
        caches = cast(dict[str, object], caches)
        memory = cast(dict[str, object], memory)
        canonical = caches["canonical_mask"]
        residual = caches["residual_decomposition"]
        canonical = cast(dict[str, object], canonical)
        residual = cast(dict[str, object], residual)

        def rate_text(value: object) -> str:
            return "n/a" if value is None else f"{float(value):.1%}"

        peak = memory.get("overall_peak_resident_kib")
        peak_text = "n/a" if peak is None else f"{int(peak):,} KiB"
        print(
            f"  {result.case.name:<{name_width}.{name_width}} "
            f"{int(counters['retained_masks']):>10,} "
            f"{int(counters['matching_visits']):>12,} "
            f"{int(counters['canonicalisation_calls']):>12,} "
            f"{int(counters['vf2_calls']):>8,} "
            f"{rate_text(canonical['hit_rate']):>10} "
            f"{rate_text(residual['lookup_hit_rate']):>13} "
            f"{peak_text:>11}"
        )

    parallel_results = [
        result
        for result in results
        if result.telemetry is not None
        and isinstance(result.telemetry.get("parallel"), dict)
    ]
    if not parallel_results:
        return

    print("\nParallel worker telemetry")
    print(
        f"  {'Case':<{name_width}} {'Topology':>12} {'Branches':>12} "
        f"{'Critical':>12} {'Imbalance':>10} {'Coverage':>10}"
    )
    print("  " + "-" * (name_width + 71))
    for result in parallel_results:
        telemetry = cast(dict[str, object], result.telemetry)
        parallel = cast(dict[str, object], telemetry["parallel"])
        aggregate = parallel["aggregate"]
        workers = parallel["workers"]
        aggregate = cast(dict[str, object], aggregate)
        workers = cast(list[object], workers)
        elapsed_values = [
            int(worker["elapsed_nanoseconds"])
            for worker in workers
            if isinstance(worker, dict)
        ]
        mean_elapsed = statistics.fmean(elapsed_values)
        imbalance = max(elapsed_values) / mean_elapsed if mean_elapsed else 1.0
        worker_elapsed = int(aggregate["worker_elapsed_nanoseconds"])
        worker_busy = int(aggregate["worker_busy_nanoseconds"])
        coverage = worker_busy / worker_elapsed if worker_elapsed else 0.0
        topology = f"{parallel['mode']}:{parallel['worker_count']}w"
        branch_candidates = aggregate["branch_candidates"]
        branch_text = (
            "n/a" if branch_candidates is None else f"{int(branch_candidates):,}"
        )
        critical_ms = int(aggregate["elapsed_nanoseconds"]) / 1_000_000
        print(
            f"  {result.case.name:<{name_width}.{name_width}} "
            f"{topology:>12} {branch_text:>12} {critical_ms:>9.3f} ms "
            f"{imbalance:>9.3f}x {coverage:>9.1%}"
        )


def write_json_report(
    path: Path,
    candidate_metadata: dict[str, object],
    baseline_metadata: dict[str, object] | None,
    telemetry_metadata: dict[str, object] | None,
    corpus_metadata: dict[str, object],
    manifest: Path | None,
    suite: str | None,
    runs: int,
    warmup: int,
    timeout: float,
    results: Sequence[CaseResult],
    candidate_execution: ExecutionConfig | None = None,
    baseline_execution: ExecutionConfig | None = None,
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
        candidate_report = measurement_report(result.measurements)
        candidate_report["telemetry"] = result.telemetry
        case_report: dict[str, object] = {
            "name": result.case.name,
            "input": str(result.case.source),
            "expected_assembly_index": result.case.expected_assembly_index,
            "expectation": result.case.expectation,
            "workload": result.case.workload,
            "suites": list(result.case.suites),
            "candidate": candidate_report,
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
        "schema_version": 2,
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
            "telemetry": telemetry_metadata,
        },
        "execution": {
            "candidate": execution_config_metadata(candidate_execution),
            "baseline": (
                None
                if baseline_metadata is None
                else execution_config_metadata(baseline_execution)
            ),
            "telemetry": (
                None
                if telemetry_metadata is None
                else execution_config_metadata(candidate_execution)
            ),
        },
        "corpus": corpus_metadata,
        "manifest": None if manifest is None else str(manifest),
        "suite": suite,
        "runs": runs,
        "warmup": warmup,
        "timeout_seconds": timeout,
        "telemetry": {
            "enabled": any(result.telemetry is not None for result in results),
            "collection": (
                "one untimed separate-instrumented-executable run per case"
                if any(result.telemetry is not None for result in results)
                else None
            ),
            "excluded_from_timing_aggregates": True,
        },
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
    print(f"{'Case':<24} {'Suites':<20} {'Index':>8} {'Expectation':<12}  Workload")
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
        description="Benchmark one input or a manifest suite."
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=(
            "AssemblyCpp executable "
            f"(default: {DEFAULT_EXECUTABLE.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--baseline-executable",
        type=Path,
        help=("baseline executable for paired comparisons (candidate: --executable)"),
    )
    parser.add_argument(
        "--candidate-launcher",
        type=launcher_prefix,
        metavar="COMMAND",
        help=(
            "command prefix for the candidate executable (for example, 'mpirun -n 2')"
        ),
    )
    parser.add_argument(
        "--baseline-launcher",
        type=launcher_prefix,
        metavar="COMMAND",
        help="command prefix for the baseline executable",
    )
    parser.add_argument(
        "--candidate-env",
        "--candidate-environment",
        dest="candidate_environment",
        type=environment_assignment,
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="candidate environment override (repeatable)",
    )
    parser.add_argument(
        "--baseline-env",
        "--baseline-environment",
        dest="baseline_environment",
        type=environment_assignment,
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="baseline environment override (repeatable)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help=(
            "single molfile or native graph file "
            f"(default: {DEFAULT_INPUT.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--expected",
        type=int,
        help=(
            "expected assembly index; defaults to 22 for the default input, "
            "otherwise unchecked"
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help=(
            "benchmark corpus manifest "
            "(default with corpus options: "
            f"{DEFAULT_MANIFEST.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--suite",
        choices=KNOWN_SUITES,
        help="select a suite from the benchmark manifest",
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="NAME",
        help="select a manifest case (repeatable)",
    )
    parser.add_argument(
        "--list-cases",
        action="store_true",
        help="list selected cases without building or running",
    )
    parser.add_argument(
        "--runs",
        type=positive_int,
        help="measured rounds (default: 5; paired: 6)",
    )
    parser.add_argument(
        "--warmup",
        type=non_negative_int,
        default=1,
        help="warm-up rounds (default: 1)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=600.0,
        help="per-calculation timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="write a JSON report",
    )
    parser.add_argument(
        "--telemetry",
        action="store_true",
        help="run each case once more with untimed telemetry",
    )
    parser.add_argument(
        "--telemetry-executable",
        type=Path,
        help=("instrumented executable for telemetry; required unless --build is used"),
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="compile v5/main.cpp for x86-64-v3 before benchmarking",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler used by --build (uses $CXX, then c++)",
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
        resolved_manifest, cases = load_manifest(manifest_path)
        return resolved_manifest, select_cases(
            cases,
            arguments.suite,
            arguments.case,
        )

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
        manifest_path, cases = resolve_requested_cases(arguments)
        if arguments.list_cases:
            print_case_list(cases)
            return 0

        baseline_executable = None
        if arguments.baseline_executable is not None:
            baseline_executable = resolve_executable(arguments.baseline_executable)

        candidate_execution = create_execution_config(
            arguments.candidate_launcher,
            arguments.candidate_environment,
            "candidate",
        )
        if baseline_executable is None and (
            arguments.baseline_launcher or arguments.baseline_environment
        ):
            raise BenchmarkError(
                "--baseline-launcher and --baseline-env require --baseline-executable"
            )
        baseline_execution = create_execution_config(
            arguments.baseline_launcher,
            arguments.baseline_environment,
            "baseline",
        )

        if arguments.telemetry_executable is not None and not arguments.telemetry:
            raise BenchmarkError("--telemetry-executable requires --telemetry")

        telemetry_path = None
        if arguments.telemetry:
            telemetry_path = arguments.telemetry_executable
            if telemetry_path is None:
                if not arguments.build:
                    raise BenchmarkError(
                        "--telemetry requires --telemetry-executable unless "
                        "--build is used"
                    )
                telemetry_path = arguments.executable.resolve().with_name(
                    arguments.executable.resolve().name + "Telemetry"
                )
            if telemetry_path.resolve() == arguments.executable.resolve():
                raise BenchmarkError(
                    "telemetry and timed candidate executables must be distinct"
                )
            ensure_distinct_executables(
                telemetry_path.resolve(),
                arguments.executable.resolve(),
            )
            if (
                baseline_executable is not None
                and telemetry_path.resolve() == baseline_executable.resolve()
            ):
                raise BenchmarkError(
                    "telemetry and baseline executables must be distinct"
                )
            if baseline_executable is not None:
                ensure_distinct_executables(
                    telemetry_path.resolve(),
                    baseline_executable,
                )

        executable_path = arguments.executable
        if arguments.build:
            if baseline_executable is not None:
                ensure_distinct_executables(
                    executable_path.resolve(),
                    baseline_executable,
                    candidate_execution,
                    baseline_execution,
                )
            executable_path = build_executable(
                executable_path,
                arguments.compiler,
            )
        executable = resolve_executable(executable_path)
        telemetry_executable = None
        if arguments.telemetry:
            telemetry_path = cast(Path, telemetry_path)
            if arguments.build:
                telemetry_path = build_executable(
                    telemetry_path,
                    arguments.compiler,
                    telemetry=True,
                )
            telemetry_executable = resolve_executable(telemetry_path)
        if baseline_executable is not None:
            ensure_distinct_executables(
                executable,
                baseline_executable,
                candidate_execution,
                baseline_execution,
            )
            if runs % 2 != 0:
                print(
                    "Warning: an odd number of paired rounds gives one "
                    "executable more first-run positions",
                    file=sys.stderr,
                )
        if telemetry_executable is not None:
            ensure_distinct_executables(executable, telemetry_executable)
            if baseline_executable is not None:
                ensure_distinct_executables(
                    baseline_executable,
                    telemetry_executable,
                )
        if arguments.json_output is not None:
            ensure_json_output_is_distinct(
                output=arguments.json_output,
                candidate=executable,
                baseline=baseline_executable,
                telemetry=telemetry_executable,
                manifest=manifest_path,
                cases=cases,
            )
        candidate_metadata = executable_metadata(executable)
        baseline_metadata = (
            None
            if baseline_executable is None
            else executable_metadata(baseline_executable)
        )
        if baseline_metadata is not None:
            ensure_distinct_execution_identities(
                candidate_metadata,
                baseline_metadata,
                candidate_execution,
                baseline_execution,
            )
        telemetry_metadata = (
            None
            if telemetry_executable is None
            else executable_metadata(telemetry_executable)
        )
        corpus_metadata = (
            None
            if arguments.json_output is None
            else benchmark_corpus_metadata(manifest_path, cases)
        )

        if baseline_executable is None:
            print(f"Executable: {executable}")
        else:
            print(f"Candidate executable: {executable}")
            print(f"Baseline executable: {baseline_executable}")
        if telemetry_executable is not None:
            print(f"Telemetry executable: {telemetry_executable}")
        if manifest_path is None:
            print(f"Input: {cases[0].source}")
            if cases[0].expected_assembly_index is None:
                print("Assembly index: present but not validated")
            else:
                print(f"Expected assembly index: {cases[0].expected_assembly_index}")
        else:
            print(f"Manifest: {manifest_path}")
            print(f"Suite: {arguments.suite or 'all'}")
            print(f"Cases: {len(cases)}")
            provisional = [
                case.name for case in cases if case.expectation == "provisional"
            ]
            if provisional:
                print(f"Provisional assembly indices: {', '.join(provisional)}")
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
        if arguments.telemetry:
            print(
                "Telemetry: one separate, untimed run per case",
                flush=True,
            )

        results = run_benchmarks(
            executable=executable,
            baseline_executable=baseline_executable,
            telemetry_executable=telemetry_executable,
            candidate_execution=candidate_execution,
            baseline_execution=baseline_execution,
            cases=cases,
            runs=runs,
            warmup=arguments.warmup,
            timeout=arguments.timeout,
        )
        verify_executable_unchanged(executable, candidate_metadata, "candidate")
        if baseline_executable is not None:
            baseline_metadata = cast(dict[str, object], baseline_metadata)
            verify_executable_unchanged(
                baseline_executable, baseline_metadata, "baseline"
            )
        if telemetry_executable is not None:
            telemetry_metadata = cast(dict[str, object], telemetry_metadata)
            verify_executable_unchanged(
                telemetry_executable,
                telemetry_metadata,
                "telemetry",
            )
        if (
            corpus_metadata is not None
            and benchmark_corpus_metadata(manifest_path, cases) != corpus_metadata
        ):
            raise BenchmarkError(
                "benchmark manifest or selected input changed during measurement"
            )
        print_summary(results)
        print_telemetry_summary(results)
        if arguments.json_output is not None:
            corpus_metadata = cast(dict[str, object], corpus_metadata)
            write_json_report(
                path=arguments.json_output,
                candidate_metadata=candidate_metadata,
                baseline_metadata=baseline_metadata,
                telemetry_metadata=telemetry_metadata,
                corpus_metadata=corpus_metadata,
                manifest=manifest_path,
                suite=arguments.suite,
                runs=runs,
                warmup=arguments.warmup,
                timeout=arguments.timeout,
                results=results,
                candidate_execution=candidate_execution,
                baseline_execution=baseline_execution,
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
