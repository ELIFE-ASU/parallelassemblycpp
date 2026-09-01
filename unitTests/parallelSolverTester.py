"""Check serial/parallel solver parity and parallel telemetry reduction."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent
DEFAULT_MANIFEST = REPOSITORY_ROOT / "benchmarks" / "cases.tsv"
SKIP_RETURN_CODE = 77
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
DEFAULT_BRANCH_LEASE_SIZE = 4
ROOT_ENUMERATION_PHASES = ("initial_enumeration", "dag_conversion")
ADAPTIVE_SPLITTING_POLICY = {
    "minimum_queued_tasks_per_worker": 8,
    "target_queued_tasks_per_worker": 16,
    "maximum_queued_tasks_per_worker": 32,
    "maximum_depth": 4,
    "warm_start": "largest_duplicate_first",
}
SCHEDULER_SUM_FIELDS = (
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


@dataclass(frozen=True)
class CaseSpec:
    name: str
    edges: int
    active_mask_words: int


@dataclass(frozen=True)
class SolverCase:
    name: str
    source: Path
    expected_index: int
    edges: int
    active_mask_words: int


@dataclass(frozen=True)
class ParallelTopology:
    mode: str
    threads_per_rank: tuple[int, ...]

    @property
    def rank_count(self) -> int:
        return len(self.threads_per_rank)

    @property
    def worker_count(self) -> int:
        return sum(self.threads_per_rank)

    @property
    def aggregation_scope(self) -> str:
        return "process" if self.mode == "openmp" else "all_mpi_ranks"

    @property
    def description(self) -> str:
        threads = self.threads_per_rank[0]
        return f"{self.mode} {self.rank_count}x{threads} ({self.worker_count} workers)"


@dataclass(frozen=True)
class ParallelTarget:
    label: str
    executable: Path
    topology: ParallelTopology
    mpiexec: Path | None = None
    numproc_flag: str = "-n"
    branch_lease_size: int | None = None


PARITY_CASES = (
    CaseSpec("mask-boundary-path-063b", 63, 1),
    CaseSpec("mask-boundary-path-064b", 64, 1),
    CaseSpec("mask-boundary-path-065b", 65, 2),
    CaseSpec("mask-boundary-path-127b", 127, 2),
    CaseSpec("mask-boundary-path-128b", 128, 2),
    CaseSpec("mask-boundary-path-129b", 129, 3),
    # Four disconnected amino-acid components exercise component preprocessing;
    # default hydrogen removal leaves 31 processed edges from 32 input bonds.
    CaseSpec("amino-acid-scale-04c", 31, 1),
)
TELEMETRY_CASE_NAMES = {
    "mask-boundary-path-129b",
    "amino-acid-scale-04c",
}
SPARSE_ADAPTIVE_TELEMETRY_CASE = SolverCase(
    name="sparse-1253",
    source=TEST_DIRECTORY / "1253.mol",
    expected_index=9,
    # Hydrogen removal leaves twelve processed bonds in this sparse fixture.
    edges=12,
    active_mask_words=1,
)
LATE_REFILL_ADAPTIVE_TELEMETRY_CASE = SolverCase(
    name="late-refill-amino-acid-scale-09c",
    source=(
        REPOSITORY_ROOT
        / "benchmarks"
        / "inputs"
        / "scaling"
        / "amino_acid_scaling_09c_077a.mol"
    ),
    expected_index=25,
    # Hydrogen removal leaves 66 processed bonds from 68 input bonds.
    edges=66,
    active_mask_words=2,
)
PATHWAY_PARITY_NAME = "ketoconazole-pathway"
PATHWAY_PARITY_SOURCE = TEST_DIRECTORY / "ketoconazole.mol"
PATHWAY_PARITY_EXPECTED_INDEX = 22
PATHWAY_PARITY_EXPECTED = TEST_DIRECTORY / "expected_pathways" / "ketoconazole.json"
MPI_TOPOLOGY = ParallelTopology("mpi", (1, 1))
HYBRID_TOPOLOGY = ParallelTopology("hybrid", (2, 2))


class TestFailureError(RuntimeError):
    """Raised when a solver result or telemetry invariant is incorrect."""


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


def resolve_path(path: Path) -> Path:
    if path.is_absolute():
        return path
    return (Path.cwd() / path).resolve()


def resolve_command_path(path: Path) -> Path:
    command = str(path)
    if path.is_absolute() or path.parent != Path():
        return resolve_path(path)
    discovered = shutil.which(command)
    return Path(discovered) if discovered is not None else resolve_path(path)


def load_cases(manifest_argument: Path) -> list[SolverCase]:
    manifest = resolve_path(manifest_argument)
    wanted = {spec.name: spec for spec in PARITY_CASES}
    found: dict[str, SolverCase] = {}
    try:
        with manifest.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            required = {"name", "input", "expected_assembly_index"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise TestFailureError(
                    f"invalid benchmark manifest header in {manifest}; "
                    f"required columns are {sorted(required)}"
                )
            for row in reader:
                name = row["name"].strip()
                spec = wanted.get(name)
                if spec is None:
                    continue
                if name in found:
                    raise TestFailureError(f"duplicate case {name!r} in {manifest}")
                source = (manifest.parent / row["input"].strip()).resolve()
                try:
                    expected_index = int(row["expected_assembly_index"])
                except ValueError as error:
                    raise TestFailureError(
                        f"invalid expected index for {name!r} in {manifest}"
                    ) from error
                found[name] = SolverCase(
                    name=name,
                    source=source,
                    expected_index=expected_index,
                    edges=spec.edges,
                    active_mask_words=spec.active_mask_words,
                )
    except OSError as error:
        raise TestFailureError(
            f"cannot read benchmark manifest {manifest}: {error}"
        ) from error

    missing = [spec.name for spec in PARITY_CASES if spec.name not in found]
    if missing:
        raise TestFailureError(
            f"benchmark manifest is missing cases: {', '.join(missing)}"
        )
    cases = [found[spec.name] for spec in PARITY_CASES]
    absent_fixtures = [str(case.source) for case in cases if not case.source.is_file()]
    if absent_fixtures:
        raise TestFailureError(f"missing fixture(s): {', '.join(absent_fixtures)}")
    return cases


def executable_status(paths: Sequence[tuple[str, Path | None]]) -> int | None:
    missing = [
        f"{name}={path}"
        for name, path in paths
        if path is not None and not path.is_file()
    ]
    if missing:
        print(f"SKIP: requested executable is absent: {', '.join(missing)}")
        return SKIP_RETURN_CODE
    unusable = [
        f"{name}={path}"
        for name, path in paths
        if path is not None and not os.access(path, os.X_OK)
    ]
    if unusable:
        raise TestFailureError(
            f"solver target is not executable: {', '.join(unusable)}"
        )
    return None


def output_stem(input_path: Path) -> str:
    return input_path.name.removesuffix(".mol")


def format_completed(completed: subprocess.CompletedProcess[str]) -> str:
    return (
        f"return code: {completed.returncode}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )


def terminate_command(process: subprocess.Popen[str]) -> None:
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
    arguments: Sequence[str],
    working_directory: Path,
    environment: Mapping[str, str],
    timeout: float,
) -> subprocess.CompletedProcess[str]:
    popen_arguments: dict[str, object] = {
        "cwd": working_directory,
        "env": environment,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        popen_arguments["start_new_session"] = True
    process = subprocess.Popen(arguments, **popen_arguments)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as error:
        terminate_command(process)
        stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(
            arguments,
            error.timeout,
            output=stdout,
            stderr=stderr,
        ) from error
    except BaseException:
        terminate_command(process)
        process.communicate()
        raise
    return subprocess.CompletedProcess(arguments, process.returncode, stdout, stderr)


def run_solver(
    executable: Path,
    case: SolverCase,
    timeout: float,
    *,
    topology: ParallelTopology | None = None,
    mpiexec: Path | None = None,
    numproc_flag: str = "-n",
    branch_lease_size: int | None = None,
    telemetry: bool = False,
) -> tuple[int, dict[str, Any] | None]:
    environment = os.environ.copy()
    if topology is not None:
        # A default-mode test must not inherit a developer's local override.
        environment.pop("ASSEMBLYCPP_BRANCH_LEASE_SIZE", None)
        if branch_lease_size is not None:
            require(branch_lease_size > 0, "branch lease size must be positive")
            environment["ASSEMBLYCPP_BRANCH_LEASE_SIZE"] = str(branch_lease_size)
    if topology is not None and topology.mode in {"openmp", "hybrid"}:
        require(
            len(set(topology.threads_per_rank)) == 1,
            f"{topology.mode}: test topology must use uniform thread counts",
        )
        worker_text = str(topology.threads_per_rank[0])
        environment.update(
            {
                "OMP_NUM_THREADS": worker_text,
                "OMP_THREAD_LIMIT": worker_text,
                "OMP_DYNAMIC": "FALSE",
            }
        )

    with tempfile.TemporaryDirectory(prefix="assemblycpp-parallel-") as temporary:
        working_directory = Path(temporary)
        input_name = "input.mol" if case.source.suffix == ".mol" else "input"
        input_path = working_directory / input_name
        shutil.copy2(case.source, input_path)
        solver_arguments = [str(executable), input_name, "--pathway=0"]
        if topology is not None:
            require(
                len(set(topology.threads_per_rank)) == 1,
                f"{topology.mode}: test topology must use uniform thread counts",
            )
            solver_arguments.extend(
                [
                    (
                        "--parallel=on"
                        if topology.worker_count > 1
                        else "--parallel=off"
                    ),
                    f"--threads={topology.threads_per_rank[0]}",
                ]
            )
        if telemetry:
            solver_arguments.append("--telemetry=1")
        if topology is not None and topology.rank_count > 1:
            require(mpiexec is not None, f"{topology.mode}: mpiexec is required")
            arguments = [
                str(mpiexec),
                numproc_flag,
                str(topology.rank_count),
                *solver_arguments,
            ]
        else:
            arguments = solver_arguments
        try:
            completed = run_command(
                arguments,
                working_directory,
                environment,
                timeout,
            )
        except subprocess.TimeoutExpired as error:
            configuration = topology.description if topology is not None else "serial"
            raise TestFailureError(
                f"{case.name}: {executable.name} exceeded {timeout:g}s with "
                f"{configuration} configuration"
            ) from error
        if completed.returncode != 0:
            raise TestFailureError(
                f"{case.name}: {executable.name} failed\n{format_completed(completed)}"
            )

        stem = output_stem(input_path)
        output_path = working_directory / f"{stem}Out"
        try:
            output_text = output_path.read_text(encoding="utf-8")
        except OSError as error:
            raise TestFailureError(
                f"{case.name}: cannot read solver output {output_path}: {error}\n"
                f"{format_completed(completed)}"
            ) from error
        match = ASSEMBLY_INDEX_PATTERN.search(output_text)
        if match is None:
            raise TestFailureError(
                f"{case.name}: assembly index is absent from {output_path}\n"
                f"output file:\n{output_text}\n{format_completed(completed)}"
            )
        index = int(match.group(1))

        document: dict[str, Any] | None = None
        if telemetry:
            telemetry_path = working_directory / f"{stem}Telemetry.json"
            try:
                parsed = json.loads(telemetry_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                raise TestFailureError(
                    f"{case.name}: cannot read telemetry {telemetry_path}: {error}"
                ) from error
            if not isinstance(parsed, dict):
                raise TestFailureError(f"{case.name}: telemetry root must be an object")
            document = parsed
        return index, document


def require(condition: bool, message: str) -> None:
    if not condition:
        raise TestFailureError(message)


def require_nonnegative_integer(value: object, path: str) -> int:
    require(
        isinstance(value, int) and not isinstance(value, bool) and value >= 0,
        f"{path} must be a non-negative integer, got {value!r}",
    )
    return value


def require_mapping(value: object, path: str) -> Mapping[str, Any]:
    require(isinstance(value, dict), f"{path} must be an object")
    return value


def validate_counter_sum(
    aggregate: object,
    worker_values: Sequence[object],
    path: str,
) -> None:
    aggregate_map = require_mapping(aggregate, path)
    require(bool(aggregate_map), f"{path} must not be empty")
    worker_maps = [
        require_mapping(value, f"worker[{index}].counters")
        for index, value in enumerate(worker_values)
    ]
    aggregate_keys = set(aggregate_map)
    for index, worker_map in enumerate(worker_maps):
        require(
            set(worker_map) == aggregate_keys,
            f"worker[{index}].counters keys differ from {path}",
        )
    for key, aggregate_value in aggregate_map.items():
        child_path = f"{path}.{key}"
        children = [worker_map[key] for worker_map in worker_maps]
        if isinstance(aggregate_value, dict):
            validate_counter_sum(aggregate_value, children, child_path)
            continue
        aggregate_integer = require_nonnegative_integer(aggregate_value, child_path)
        worker_integers = [
            require_nonnegative_integer(value, f"worker[{index}].counters.{key}")
            for index, value in enumerate(children)
        ]
        require(
            aggregate_integer == sum(worker_integers),
            f"{child_path}={aggregate_integer} does not equal worker sum "
            f"{sum(worker_integers)}",
        )


def validate_shared_assembly_cache(
    value: object,
    rank_count: int,
    path: str,
) -> dict[str, int]:
    counters = require_mapping(value, path)
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
    parsed = {
        name: require_nonnegative_integer(counters.get(name), f"{path}.{name}")
        for name in names
    }
    require(
        parsed["table_count"] <= rank_count,
        f"{path}.table_count exceeds rank count",
    )
    require(
        parsed["lock_acquisitions"] == parsed["hits"] + parsed["misses"],
        f"{path} lookup counters are inconsistent",
    )
    require(
        parsed["lock_waits"] <= parsed["lock_acquisitions"],
        f"{path}.lock_waits exceeds acquisitions",
    )
    require(
        (parsed["misses"] == 0) == (parsed["allocated_bytes"] == 0),
        f"{path} allocation counters are inconsistent",
    )
    if parsed["table_count"] == 0:
        require(
            all(parsed[name] == 0 for name in names if name != "table_count"),
            f"{path} reports activity without a table",
        )
    return parsed


def validate_root_enumeration_phases(
    document: Mapping[str, Any],
    workers: Sequence[Mapping[str, Any]],
    worker_ranks: Sequence[int],
    topology: ParallelTopology,
    prefix: str,
) -> None:
    """Require process-owned root setup to run once on each MPI rank."""
    memory = require_mapping(document.get("memory"), f"{prefix}: memory")
    aggregate_phases = require_mapping(memory.get("phases"), f"{prefix}: memory.phases")
    for phase_name in ROOT_ENUMERATION_PHASES:
        activations_per_rank = [0 for _ in range(topology.rank_count)]
        worker_activation_total = 0
        for index, worker in enumerate(workers):
            worker_path = f"{prefix}: parallel.workers[{index}]"
            phases = require_mapping(worker.get("phases"), f"{worker_path}.phases")
            phase = require_mapping(
                phases.get(phase_name), f"{worker_path}.phases.{phase_name}"
            )
            activations = require_nonnegative_integer(
                phase.get("activations"),
                f"{worker_path}.phases.{phase_name}.activations",
            )
            wall_nanoseconds = require_nonnegative_integer(
                phase.get("wall_nanoseconds"),
                f"{worker_path}.phases.{phase_name}.wall_nanoseconds",
            )
            require(
                activations != 0 or wall_nanoseconds == 0,
                f"{worker_path}.phases.{phase_name} has elapsed time without "
                "an activation",
            )
            activations_per_rank[worker_ranks[index]] += activations
            worker_activation_total += activations

        for rank, activations in enumerate(activations_per_rank):
            require(
                activations == 1,
                f"{prefix}: {phase_name} ran {activations} times on rank {rank}; "
                "root enumeration must run once per rank, not once per worker",
            )

        aggregate_phase = require_mapping(
            aggregate_phases.get(phase_name),
            f"{prefix}: memory.phases.{phase_name}",
        )
        aggregate_activations = require_nonnegative_integer(
            aggregate_phase.get("activations"),
            f"{prefix}: memory.phases.{phase_name}.activations",
        )
        require(
            aggregate_activations == worker_activation_total,
            f"{prefix}: memory.phases.{phase_name}.activations="
            f"{aggregate_activations} does not equal worker sum "
            f"{worker_activation_total}",
        )


def validate_parallel_telemetry(
    document: Mapping[str, Any],
    case: SolverCase,
    topology: ParallelTopology,
    expected_lease_size: int | None = None,
    *,
    require_adaptive_splitting: bool = False,
    require_depth_two_tasks: bool = False,
) -> None:
    requested_workers = topology.worker_count
    prefix = f"{case.name} ({topology.description})"
    parallel = require_mapping(document.get("parallel"), f"{prefix}: parallel")
    require(parallel.get("enabled") is True, f"{prefix}: parallel.enabled must be true")
    require(
        parallel.get("mode") == topology.mode,
        f"{prefix}: parallel.mode must be {topology.mode!r}",
    )
    require(
        parallel.get("aggregation_scope") == topology.aggregation_scope,
        f"{prefix}: parallel.aggregation_scope must be {topology.aggregation_scope!r}",
    )
    require(
        parallel.get("rank_count") == topology.rank_count,
        f"{prefix}: parallel.rank_count must be {topology.rank_count}",
    )
    require(
        parallel.get("worker_count") == requested_workers,
        f"{prefix}: parallel.worker_count must be {requested_workers}",
    )
    require(
        parallel.get("elapsed_timing_method") == "parallel_region_steady_clock",
        f"{prefix}: parallel.elapsed_timing_method is incorrect",
    )
    require(
        parallel.get("busy_timing_method") == "elapsed_minus_scheduler_idle_time",
        f"{prefix}: parallel.busy_timing_method is incorrect",
    )
    require(
        parallel.get("branch_scan_complete") is True,
        f"{prefix}: parallel.branch_scan_complete must be true",
    )
    branch_scheduler = require_mapping(
        parallel.get("branch_scheduler"), f"{prefix}: parallel.branch_scheduler"
    )
    require(
        "shard_ownership" not in parallel,
        f"{prefix}: legacy static shard ownership must be absent",
    )
    require(
        branch_scheduler.get("strategy") == "distributed_global_root_queue",
        f"{prefix}: parallel.branch_scheduler.strategy is incorrect",
    )
    lease_size = require_nonnegative_integer(
        branch_scheduler.get("lease_size"),
        f"{prefix}: parallel.branch_scheduler.lease_size",
    )
    require(lease_size > 0, f"{prefix}: branch lease size must be positive")
    if expected_lease_size is not None:
        require(
            lease_size == expected_lease_size,
            f"{prefix}: branch lease size {lease_size} does not match requested "
            f"size {expected_lease_size}",
        )
    elif topology.mode in {"mpi", "hybrid"}:
        require(
            lease_size == 1,
            f"{prefix}: adaptive distributed branch lease size must be 1",
        )
    root_queue = require_mapping(
        branch_scheduler.get("root_queue"),
        f"{prefix}: parallel.branch_scheduler.root_queue",
    )
    require(
        root_queue.get("participant_count") == topology.rank_count,
        f"{prefix}: root queue participant count must be {topology.rank_count}",
    )
    require(
        branch_scheduler.get("complete") is True,
        f"{prefix}: parallel.branch_scheduler.complete must be true",
    )
    adaptive_splitting: Mapping[str, Any] | None = None
    if "adaptive_splitting" in branch_scheduler:
        adaptive_splitting = require_mapping(
            branch_scheduler.get("adaptive_splitting"),
            f"{prefix}: parallel.branch_scheduler.adaptive_splitting",
        )
        for name, expected in ADAPTIVE_SPLITTING_POLICY.items():
            path = f"{prefix}: parallel.branch_scheduler.adaptive_splitting.{name}"
            if isinstance(expected, int):
                actual: object = require_nonnegative_integer(
                    adaptive_splitting.get(name), path
                )
            else:
                actual = adaptive_splitting.get(name)
            require(
                actual == expected,
                f"{path} must be {expected!r}, got {actual!r}",
            )
    require(
        not require_adaptive_splitting or adaptive_splitting is not None,
        f"{prefix}: adaptive splitting telemetry is absent",
    )
    require(
        not require_depth_two_tasks or adaptive_splitting is not None,
        f"{prefix}: cannot require depth-two tasks without adaptive telemetry",
    )

    require(
        len(set(topology.threads_per_rank)) == 1,
        f"{prefix}: expected test topology must have uniform local threads",
    )
    local_threads = topology.threads_per_rank[0]
    require(
        parallel.get("local_threads") == local_threads,
        f"{prefix}: parallel.local_threads must be {local_threads}",
    )
    require(
        parallel.get("local_threads_per_rank") == list(topology.threads_per_rank),
        f"{prefix}: parallel.local_threads_per_rank must be "
        f"{list(topology.threads_per_rank)}",
    )

    workers_value = parallel.get("workers")
    require(
        isinstance(workers_value, list), f"{prefix}: parallel.workers must be an array"
    )
    workers: list[Mapping[str, Any]] = [
        require_mapping(worker, f"{prefix}: parallel.workers[{index}]")
        for index, worker in enumerate(workers_value)
    ]
    require(
        len(workers) == requested_workers,
        f"{prefix}: expected {requested_workers} worker records, got {len(workers)}",
    )

    global_indices: list[int] = []
    worker_ranks: list[int] = []
    rank_local_indices = [set() for _ in range(topology.rank_count)]
    rank_offsets: list[int] = []
    offset = 0
    for threads in topology.threads_per_rank:
        rank_offsets.append(offset)
        offset += threads
    for index, worker in enumerate(workers):
        worker_path = f"{prefix}: parallel.workers[{index}]"
        rank = require_nonnegative_integer(worker.get("rank"), f"{worker_path}.rank")
        require(rank < topology.rank_count, f"{worker_path}.rank is out of range")
        local_index = require_nonnegative_integer(
            worker.get("local_worker_index"), f"{worker_path}.local_worker_index"
        )
        require(
            local_index < topology.threads_per_rank[rank],
            f"{worker_path}.local_worker_index is out of range",
        )
        global_index = require_nonnegative_integer(
            worker.get("global_worker_index"), f"{worker_path}.global_worker_index"
        )
        require(
            global_index == rank_offsets[rank] + local_index,
            f"{worker_path}: rank/local/global index mismatch",
        )
        rank_local_indices[rank].add(local_index)
        worker_ranks.append(rank)
        global_indices.append(global_index)
        require(
            "rank_partition" not in worker,
            f"{worker_path}: legacy static rank partition must be absent",
        )
        root_queue = require_mapping(
            worker.get("root_queue"), f"{worker_path}.root_queue"
        )
        require(
            root_queue.get("participant_rank") == rank,
            f"{worker_path}: root queue participant rank must equal rank",
        )
        require(
            root_queue.get("participant_count") == topology.rank_count,
            f"{worker_path}.root_queue.participant_count must be {topology.rank_count}",
        )

        processed = require_mapping(
            worker.get("processed_graph"), f"{worker_path}.processed_graph"
        )
        require(
            processed.get("edges") == case.edges,
            f"{worker_path}.processed_graph.edges must be {case.edges}",
        )
        require(
            processed.get("active_mask_words") == case.active_mask_words,
            f"{worker_path}.processed_graph.active_mask_words must be {case.active_mask_words}",
        )
    expected_indices = list(range(requested_workers))
    require(
        sorted(global_indices) == expected_indices,
        f"{prefix}: global worker indices are not contiguous",
    )
    for rank, local_indices in enumerate(rank_local_indices):
        require(
            local_indices == set(range(topology.threads_per_rank[rank])),
            f"{prefix}: rank {rank} local worker indices are not contiguous",
        )

    validate_root_enumeration_phases(
        document,
        workers,
        worker_ranks,
        topology,
        prefix,
    )

    aggregate = require_mapping(
        parallel.get("aggregate"), f"{prefix}: parallel.aggregate"
    )
    # This exact nesting is intentional: reductions must never be inferred from
    # the process-level telemetry counters.
    validate_counter_sum(
        aggregate.get("counters"),
        [worker.get("counters") for worker in workers],
        f"{prefix}: parallel.aggregate.counters",
    )
    shared_assembly_cache = validate_shared_assembly_cache(
        aggregate.get("shared_assembly_cache"),
        topology.rank_count,
        f"{prefix}: parallel.aggregate.shared_assembly_cache",
    )
    if case.name == "amino-acid-scale-04c" and local_threads > 1:
        require(
            shared_assembly_cache["table_count"] == topology.rank_count,
            f"{prefix}: expected one active shared assembly cache per rank",
        )
        require(
            shared_assembly_cache["hits"] + shared_assembly_cache["misses"] > 0,
            f"{prefix}: shared assembly cache emitted no lookup activity",
        )

    worker_branch_candidates: list[int] = []
    branch_leases = 0
    branch_assignments = 0
    scheduler_sums = dict.fromkeys(SCHEDULER_SUM_FIELDS, 0)
    maximum_task_queue_high_watermark = 0
    maximum_task_depth_executed = 0
    warm_starts_per_rank = [0 for _ in range(topology.rank_count)]
    worker_elapsed = 0
    worker_busy = 0
    maximum_worker_elapsed = 0
    for index, worker in enumerate(workers):
        worker_path = f"{prefix}: parallel.workers[{index}]"
        candidates = require_nonnegative_integer(
            worker.get("branch_candidates"), f"{worker_path}.branch_candidates"
        )
        assignments = require_nonnegative_integer(
            worker.get("branch_assignments"), f"{worker_path}.branch_assignments"
        )
        leases = require_nonnegative_integer(
            worker.get("branch_leases"), f"{worker_path}.branch_leases"
        )
        require(
            assignments <= candidates, f"{worker_path}: assignments exceed candidates"
        )
        require(
            (assignments == 0 and leases == 0)
            or (assignments > 0 and 0 < leases <= assignments <= leases * lease_size),
            f"{worker_path}: utilized lease count is inconsistent with assignments",
        )
        worker_branch_candidates.append(candidates)
        branch_leases += leases
        branch_assignments += assignments
        scheduler_values = {
            name: require_nonnegative_integer(worker.get(name), f"{worker_path}.{name}")
            for name in SCHEDULER_SUM_FIELDS
        }
        task_queue_high_watermark = require_nonnegative_integer(
            worker.get("task_queue_high_watermark"),
            f"{worker_path}.task_queue_high_watermark",
        )
        task_depth = require_nonnegative_integer(
            worker.get("maximum_task_depth_executed"),
            f"{worker_path}.maximum_task_depth_executed",
        )
        require(
            scheduler_values["task_steals"] <= scheduler_values["task_steal_attempts"],
            f"{worker_path}: successful task steals exceed attempts",
        )
        executed_tasks = (
            scheduler_values["depth_two_tasks_executed"]
            + scheduler_values["deeper_tasks_executed"]
        )
        require(
            scheduler_values["local_task_executions"] + scheduler_values["task_steals"]
            == executed_tasks,
            f"{worker_path}: local executions plus steals do not equal "
            "executed transferred tasks",
        )
        require(
            (executed_tasks == 0 and task_depth == 0)
            or (
                executed_tasks > 0
                and 2 <= task_depth <= ADAPTIVE_SPLITTING_POLICY["maximum_depth"]
            ),
            f"{worker_path}: maximum executed task depth is inconsistent with "
            "transferred work",
        )
        warm_starts = scheduler_values["warm_start_branches"]
        require(
            warm_starts <= 1,
            f"{worker_path}: worker warmed more than one root branch",
        )
        for name, value in scheduler_values.items():
            scheduler_sums[name] += value
        maximum_task_queue_high_watermark = max(
            maximum_task_queue_high_watermark,
            task_queue_high_watermark,
        )
        maximum_task_depth_executed = max(
            maximum_task_depth_executed,
            task_depth,
        )
        warm_starts_per_rank[worker_ranks[index]] += warm_starts
        elapsed = require_nonnegative_integer(
            worker.get("elapsed_nanoseconds"), f"{worker_path}.elapsed_nanoseconds"
        )
        busy = require_nonnegative_integer(
            worker.get("busy_nanoseconds"), f"{worker_path}.busy_nanoseconds"
        )
        require(busy <= elapsed, f"{worker_path}: busy time exceeds elapsed time")
        require(
            busy
            == elapsed - min(scheduler_values["scheduler_idle_nanoseconds"], elapsed),
            f"{worker_path}: busy time does not equal elapsed time minus "
            "scheduler idle time",
        )
        worker_elapsed += elapsed
        worker_busy += busy
        maximum_worker_elapsed = max(maximum_worker_elapsed, elapsed)

    aggregate_candidates = require_nonnegative_integer(
        aggregate.get("branch_candidates"),
        f"{prefix}: parallel.aggregate.branch_candidates",
    )
    aggregate_assignments = require_nonnegative_integer(
        aggregate.get("branch_assignments"),
        f"{prefix}: parallel.aggregate.branch_assignments",
    )
    aggregate_leases = require_nonnegative_integer(
        aggregate.get("branch_leases"),
        f"{prefix}: parallel.aggregate.branch_leases",
    )
    require(
        len(set(worker_branch_candidates)) == 1
        and aggregate_candidates == worker_branch_candidates[0],
        f"{prefix}: workers disagree on the scanned branch count",
    )
    require(
        aggregate_assignments == branch_assignments,
        f"{prefix}: aggregate branch assignment total is incorrect",
    )
    require(
        aggregate_leases == branch_leases,
        f"{prefix}: aggregate utilized lease total is incorrect",
    )
    require(aggregate_assignments > 0, f"{prefix}: no branch was assigned")
    require(
        aggregate_assignments == aggregate_candidates,
        f"{prefix}: complete branch scan did not assign every candidate",
    )
    aggregate_scheduler_values = {
        name: require_nonnegative_integer(
            aggregate.get(name), f"{prefix}: parallel.aggregate.{name}"
        )
        for name in SCHEDULER_SUM_FIELDS
    }
    aggregate_task_queue_high_watermark = require_nonnegative_integer(
        aggregate.get("task_queue_high_watermark"),
        f"{prefix}: parallel.aggregate.task_queue_high_watermark",
    )
    aggregate_maximum_task_depth = require_nonnegative_integer(
        aggregate.get("maximum_task_depth_executed"),
        f"{prefix}: parallel.aggregate.maximum_task_depth_executed",
    )
    for name, worker_sum in scheduler_sums.items():
        require(
            aggregate_scheduler_values[name] == worker_sum,
            f"{prefix}: aggregate {name} is not the worker sum",
        )
    require(
        aggregate_task_queue_high_watermark == maximum_task_queue_high_watermark,
        f"{prefix}: aggregate task queue high-water mark is not the worker maximum",
    )
    require(
        aggregate_maximum_task_depth == maximum_task_depth_executed,
        f"{prefix}: aggregate maximum task depth is not the worker maximum",
    )
    aggregate_depth_two_spawned = aggregate_scheduler_values["depth_two_tasks_spawned"]
    aggregate_depth_two_executed = aggregate_scheduler_values[
        "depth_two_tasks_executed"
    ]
    aggregate_deeper_spawned = aggregate_scheduler_values["deeper_tasks_spawned"]
    aggregate_deeper_executed = aggregate_scheduler_values["deeper_tasks_executed"]
    aggregate_executed_tasks = aggregate_depth_two_executed + aggregate_deeper_executed
    require(
        aggregate_scheduler_values["task_steals"]
        <= aggregate_scheduler_values["task_steal_attempts"],
        f"{prefix}: aggregate successful task steals exceed attempts",
    )
    require(
        aggregate_scheduler_values["local_task_executions"]
        + aggregate_scheduler_values["task_steals"]
        == aggregate_executed_tasks,
        f"{prefix}: aggregate local executions plus steals do not equal "
        "executed transferred tasks",
    )
    require(
        aggregate_depth_two_spawned == aggregate_depth_two_executed,
        f"{prefix}: spawned depth-two tasks were not each executed once",
    )
    require(
        aggregate_depth_two_spawned + aggregate_deeper_spawned
        == aggregate_executed_tasks,
        f"{prefix}: spawned transferred tasks were not each executed once",
    )
    require(
        (aggregate_executed_tasks == 0 and aggregate_maximum_task_depth == 0)
        or (
            aggregate_executed_tasks > 0
            and 2
            <= aggregate_maximum_task_depth
            <= ADAPTIVE_SPLITTING_POLICY["maximum_depth"]
        ),
        f"{prefix}: aggregate maximum executed task depth is inconsistent with "
        "transferred work",
    )
    if adaptive_splitting is not None:
        aggregate_warm_starts = aggregate_scheduler_values["warm_start_branches"]
        expected_warm_starts_per_rank = 1 if aggregate_candidates > 0 else 0
        for rank, warm_starts in enumerate(warm_starts_per_rank):
            require(
                warm_starts == expected_warm_starts_per_rank,
                f"{prefix}: rank {rank} warmed {warm_starts} root branches, "
                f"expected {expected_warm_starts_per_rank}",
            )
        require(
            aggregate_warm_starts
            == expected_warm_starts_per_rank * topology.rank_count,
            f"{prefix}: aggregate warm-start count is incorrect",
        )
        require(
            not require_depth_two_tasks or aggregate_depth_two_spawned > 0,
            f"{prefix}: adaptive frontier did not expose depth-two work",
        )
    aggregate_elapsed = require_nonnegative_integer(
        aggregate.get("elapsed_nanoseconds"),
        f"{prefix}: parallel.aggregate.elapsed_nanoseconds",
    )
    require(
        aggregate_elapsed >= maximum_worker_elapsed,
        f"{prefix}: aggregate elapsed time is shorter than a worker",
    )
    aggregate_worker_elapsed = require_nonnegative_integer(
        aggregate.get("worker_elapsed_nanoseconds"),
        f"{prefix}: parallel.aggregate.worker_elapsed_nanoseconds",
    )
    aggregate_worker_busy = require_nonnegative_integer(
        aggregate.get("worker_busy_nanoseconds"),
        f"{prefix}: parallel.aggregate.worker_busy_nanoseconds",
    )
    require(
        aggregate_worker_elapsed == worker_elapsed,
        f"{prefix}: aggregate worker elapsed time is not the worker sum",
    )
    require(
        aggregate_worker_busy == worker_busy,
        f"{prefix}: aggregate worker busy time is not the worker sum",
    )


def run_parity_suite(
    serial: Path,
    openmp: Path,
    cases: Sequence[SolverCase],
    repetitions: int,
    timeout: float,
) -> int:
    runs = 0
    for case in cases:
        observed: dict[str, list[int]] = {}
        configurations: tuple[tuple[str, Path, ParallelTopology | None], ...] = (
            ("serial", serial, None),
            ("openmp-1", openmp, ParallelTopology("openmp", (1,))),
            ("openmp-2", openmp, ParallelTopology("openmp", (2,))),
            ("openmp-4", openmp, ParallelTopology("openmp", (4,))),
        )
        for label, executable, topology in configurations:
            indices = [
                run_solver(
                    executable,
                    case,
                    timeout,
                    topology=topology,
                )[0]
                for _ in range(repetitions)
            ]
            runs += repetitions
            require(
                all(index == case.expected_index for index in indices),
                f"{case.name}: {label} returned {indices}, expected "
                f"{case.expected_index} on every run",
            )
            observed[label] = indices
        require(
            len({index for indices in observed.values() for index in indices}) == 1,
            f"{case.name}: serial/OpenMP index mismatch: {observed}",
        )
        print(f"PASS parity {case.name}: index {case.expected_index}")
    runs += run_invalid_lease_configuration_suite(openmp, cases[0], timeout)
    return runs


def run_pathway_parity_suite(target: ParallelTarget, timeout: float) -> int:
    topology = target.topology
    require(
        len(set(topology.threads_per_rank)) == 1,
        f"{topology.mode}: test topology must use uniform thread counts",
    )
    threads = topology.threads_per_rank[0]
    environment = os.environ.copy()
    environment.pop("ASSEMBLYCPP_BRANCH_LEASE_SIZE", None)
    if target.branch_lease_size is not None:
        require(target.branch_lease_size > 0, "branch lease size must be positive")
        environment["ASSEMBLYCPP_BRANCH_LEASE_SIZE"] = str(target.branch_lease_size)
    if topology.mode in {"openmp", "hybrid"}:
        thread_text = str(threads)
        environment.update(
            {
                "OMP_NUM_THREADS": thread_text,
                "OMP_THREAD_LIMIT": thread_text,
                "OMP_DYNAMIC": "FALSE",
            }
        )
    try:
        expected_pathway = json.loads(
            PATHWAY_PARITY_EXPECTED.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise TestFailureError(
            f"cannot read expected pathway {PATHWAY_PARITY_EXPECTED}: {error}"
        ) from error

    with tempfile.TemporaryDirectory(
        prefix="assemblycpp-parallel-pathway-"
    ) as temporary:
        working_directory = Path(temporary)
        input_name = "input.mol"
        shutil.copy2(PATHWAY_PARITY_SOURCE, working_directory / input_name)
        solver_arguments = [
            str(target.executable),
            input_name,
            "--parallel=on",
            f"--threads={threads}",
        ]
        if topology.rank_count > 1:
            require(
                target.mpiexec is not None,
                f"{topology.mode}: mpiexec is required",
            )
            arguments = [
                str(target.mpiexec),
                target.numproc_flag,
                str(topology.rank_count),
                *solver_arguments,
            ]
        else:
            arguments = solver_arguments
        try:
            completed = run_command(arguments, working_directory, environment, timeout)
        except subprocess.TimeoutExpired as error:
            raise TestFailureError(
                f"{PATHWAY_PARITY_NAME}: {target.label} pathway run exceeded "
                f"{timeout:g}s"
            ) from error
        if completed.returncode != 0:
            raise TestFailureError(
                f"{PATHWAY_PARITY_NAME}: {target.executable.name} failed\n"
                f"{format_completed(completed)}"
            )

        output_path = working_directory / "inputOut"
        pathway_path = working_directory / "inputPathway"
        try:
            output_text = output_path.read_text(encoding="utf-8")
            actual_pathway = json.loads(pathway_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise TestFailureError(
                f"{PATHWAY_PARITY_NAME}: cannot read parallel pathway output: "
                f"{error}\n"
                f"{format_completed(completed)}"
            ) from error

        match = ASSEMBLY_INDEX_PATTERN.search(output_text)
        require(
            match is not None,
            f"{PATHWAY_PARITY_NAME}: assembly index is absent from {output_path}",
        )
        index = int(match.group(1))
        require(
            index == PATHWAY_PARITY_EXPECTED_INDEX,
            f"{PATHWAY_PARITY_NAME}: index {index}, expected "
            f"{PATHWAY_PARITY_EXPECTED_INDEX}",
        )
        require(
            actual_pathway == expected_pathway,
            f"{PATHWAY_PARITY_NAME}: parallel pathway JSON differs from "
            f"{PATHWAY_PARITY_EXPECTED}",
        )

    print(f"PASS pathway {PATHWAY_PARITY_NAME}: {target.label} index and JSON parity")
    return 1


def run_openmp_execution_policy_suite(openmp: Path, timeout: float) -> int:
    """Check automatic fallback, forced parallelism, and hard blockers."""
    case = SPARSE_ADAPTIVE_TELEMETRY_CASE
    environment = os.environ.copy()
    environment.pop("ASSEMBLYCPP_BRANCH_LEASE_SIZE", None)
    environment.update(
        {
            "OMP_NUM_THREADS": "2",
            "OMP_THREAD_LIMIT": "2",
            "OMP_DYNAMIC": "FALSE",
        }
    )
    fallback_prefix = "parallel: serial fallback: "
    scenarios = (
        (
            "automatic-work-estimate",
            ("--parallel=auto", "--threads=2"),
            0,
            (
                fallback_prefix + "estimated work ",
                "root jobs x",
                "retained DAG nodes) is below",
            ),
            (),
            False,
        ),
        (
            "forced-parallel",
            ("--parallel=on", "--threads=2"),
            0,
            (),
            (fallback_prefix,),
            False,
        ),
        (
            "parallel-disabled",
            ("--parallel=off", "--threads=2"),
            0,
            (),
            (fallback_prefix,),
            False,
        ),
        (
            "automatic-finite-runtime",
            (
                "--parallel=auto",
                "--threads=2",
                "--runtime=1000000000",
            ),
            0,
            (fallback_prefix + "finite --runtime budgets require serial execution",),
            (),
            False,
        ),
        (
            "forced-finite-runtime",
            (
                "--parallel=on",
                "--threads=2",
                "--runtime=1000000000",
            ),
            1,
            (
                (
                    "error: --parallel=on cannot be honored: finite --runtime "
                    "budgets require serial execution"
                ),
            ),
            (),
            False,
        ),
        (
            "automatic-intermediate-output",
            (
                "--parallel=auto",
                "--threads=2",
                "--write-intermediate-mas=1",
            ),
            0,
            (fallback_prefix + "--write-intermediate-mas requires serial execution",),
            (),
            True,
        ),
        (
            "forced-intermediate-output",
            (
                "--parallel=on",
                "--threads=2",
                "--write-intermediate-mas=1",
            ),
            1,
            (
                (
                    "error: --parallel=on cannot be honored: "
                    "--write-intermediate-mas requires serial execution"
                ),
            ),
            (),
            False,
        ),
    )

    for (
        name,
        options,
        expected_returncode,
        required_stderr,
        forbidden_stderr,
        expect_intermediate,
    ) in scenarios:
        with tempfile.TemporaryDirectory(
            prefix=f"assemblycpp-policy-{name}-"
        ) as temporary:
            working_directory = Path(temporary)
            input_name = "input.mol"
            shutil.copy2(case.source, working_directory / input_name)
            arguments = [
                str(openmp),
                input_name,
                "--pathway=0",
                *options,
            ]
            try:
                completed = run_command(
                    arguments,
                    working_directory,
                    environment,
                    timeout,
                )
            except subprocess.TimeoutExpired as error:
                raise TestFailureError(
                    f"execution policy {name!r} exceeded {timeout:g}s"
                ) from error

            require(
                completed.returncode == expected_returncode,
                f"execution policy {name!r} returned {completed.returncode}, "
                f"expected {expected_returncode}\n{format_completed(completed)}",
            )
            for diagnostic in required_stderr:
                require(
                    diagnostic in completed.stderr,
                    f"execution policy {name!r} did not report {diagnostic!r}\n"
                    f"{format_completed(completed)}",
                )
            for diagnostic in forbidden_stderr:
                require(
                    diagnostic not in completed.stderr,
                    f"execution policy {name!r} unexpectedly reported "
                    f"{diagnostic!r}\n{format_completed(completed)}",
                )

            if expected_returncode == 0:
                output_path = working_directory / "inputOut"
                try:
                    output_text = output_path.read_text(encoding="utf-8")
                except OSError as error:
                    raise TestFailureError(
                        f"execution policy {name!r} cannot read {output_path}: "
                        f"{error}\n{format_completed(completed)}"
                    ) from error
                match = ASSEMBLY_INDEX_PATTERN.search(output_text)
                require(
                    match is not None and int(match.group(1)) == case.expected_index,
                    f"execution policy {name!r} did not preserve assembly index "
                    f"{case.expected_index}\n{format_completed(completed)}",
                )
                intermediate_path = working_directory / "inputIntermediateMAs"
                require(
                    intermediate_path.is_file() == expect_intermediate,
                    f"execution policy {name!r} produced an unexpected "
                    "intermediate-output state",
                )

    # Hydrogen removal leaves 24 processed bonds, below the retired fixed
    # 32-bond gate, while the prepared root/DAG frontier is large enough for
    # automatic parallel execution.
    high_work_case = SolverCase(
        name="automatic-high-work-sucrose",
        source=TEST_DIRECTORY / "sucrose.mol",
        expected_index=8,
        edges=24,
        active_mask_words=1,
    )
    with tempfile.TemporaryDirectory(
        prefix="assemblycpp-policy-automatic-high-work-"
    ) as temporary:
        working_directory = Path(temporary)
        input_name = "input.mol"
        shutil.copy2(high_work_case.source, working_directory / input_name)

        # The default is serial even for high estimated work. A malformed
        # parallel-only lease setting proves that no parallel setup was
        # attempted when --parallel is omitted.
        default_environment = environment.copy()
        default_environment["ASSEMBLYCPP_BRANCH_LEASE_SIZE"] = "invalid"
        try:
            default_completed = run_command(
                [
                    str(openmp),
                    input_name,
                    "--pathway=0",
                    "--threads=2",
                ],
                working_directory,
                default_environment,
                timeout,
            )
        except subprocess.TimeoutExpired as error:
            raise TestFailureError(
                f"{high_work_case.name}: default serial policy exceeded {timeout:g}s"
            ) from error
        require(
            default_completed.returncode == 0,
            f"{high_work_case.name}: default policy attempted parallel "
            f"execution\n{format_completed(default_completed)}",
        )
        require(
            fallback_prefix not in default_completed.stderr,
            f"{high_work_case.name}: default-off policy reported a fallback\n"
            f"{format_completed(default_completed)}",
        )

        try:
            completed = run_command(
                [
                    str(openmp),
                    input_name,
                    "--pathway=0",
                    "--parallel=auto",
                    "--threads=2",
                ],
                working_directory,
                environment,
                timeout,
            )
        except subprocess.TimeoutExpired as error:
            raise TestFailureError(
                f"{high_work_case.name}: automatic policy exceeded {timeout:g}s"
            ) from error
        require(
            completed.returncode == 0,
            f"{high_work_case.name}: automatic policy failed\n"
            f"{format_completed(completed)}",
        )
        require(
            fallback_prefix not in completed.stderr,
            f"{high_work_case.name}: high prepared work unexpectedly selected "
            f"serial execution\n{format_completed(completed)}",
        )
        output_path = working_directory / "inputOut"
        try:
            output_text = output_path.read_text(encoding="utf-8")
        except OSError as error:
            raise TestFailureError(
                f"{high_work_case.name}: cannot read {output_path}: {error}\n"
                f"{format_completed(completed)}"
            ) from error
        match = ASSEMBLY_INDEX_PATTERN.search(output_text)
        require(
            match is not None and int(match.group(1)) == high_work_case.expected_index,
            f"{high_work_case.name}: automatic parallel search did not preserve "
            f"assembly index {high_work_case.expected_index}\n"
            f"{format_completed(completed)}",
        )

    print(
        "PASS execution policy: default-off, auto/on/off, "
        "estimate/runtime/intermediate fallbacks, hard errors, and "
        "sub-32-bond high-work selection"
    )
    return len(scenarios) + 1


def run_invalid_lease_configuration_suite(
    openmp: Path,
    case: SolverCase,
    timeout: float,
) -> int:
    invalid_values = ("0", "", "nonnumeric")
    for value in invalid_values:
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "2",
                "OMP_THREAD_LIMIT": "2",
                "OMP_DYNAMIC": "FALSE",
                "ASSEMBLYCPP_BRANCH_LEASE_SIZE": value,
            }
        )
        with tempfile.TemporaryDirectory(
            prefix="assemblycpp-invalid-lease-"
        ) as temporary:
            working_directory = Path(temporary)
            input_name = "input.mol" if case.source.suffix == ".mol" else "input"
            shutil.copy2(case.source, working_directory / input_name)
            try:
                completed = run_command(
                    [
                        str(openmp),
                        input_name,
                        "--pathway=0",
                        "--parallel=on",
                        "--threads=2",
                    ],
                    working_directory,
                    environment,
                    timeout,
                )
            except subprocess.TimeoutExpired as error:
                raise TestFailureError(
                    f"invalid branch lease size {value!r} exceeded {timeout:g}s"
                ) from error
        display_value = repr(value)
        require(
            completed.returncode != 0,
            f"invalid branch lease size {display_value} was accepted",
        )
        require(
            "ASSEMBLYCPP_BRANCH_LEASE_SIZE" in completed.stderr,
            f"invalid branch lease size {display_value} did not name its "
            "environment variable in stderr",
        )
    print("PASS configuration: rejected lease sizes 0, empty, and nonnumeric")
    return len(invalid_values)


def run_distributed_parity_suite(
    serial: Path,
    target: ParallelTarget,
    cases: Sequence[SolverCase],
    timeout: float,
) -> int:
    runs = 0
    for case in cases:
        serial_index, _ = run_solver(serial, case, timeout)
        parallel_index, _ = run_solver(
            target.executable,
            case,
            timeout,
            topology=target.topology,
            mpiexec=target.mpiexec,
            numproc_flag=target.numproc_flag,
            branch_lease_size=target.branch_lease_size,
        )
        runs += 2
        require(
            serial_index == case.expected_index,
            f"{case.name}: serial index {serial_index}, expected {case.expected_index}",
        )
        require(
            parallel_index == serial_index,
            f"{case.name}: {target.label} index {parallel_index} does not match "
            f"serial index {serial_index}",
        )
        print(f"PASS parity {case.name}: {target.label} index {case.expected_index}")
    return runs


def run_telemetry_suite(
    serial: Path,
    telemetry_openmp: Path,
    cases: Sequence[SolverCase],
    repetitions: int,
    timeout: float,
) -> int:
    selected = [case for case in cases if case.name in TELEMETRY_CASE_NAMES]
    runs = 0
    for case in selected:
        serial_index, _ = run_solver(serial, case, timeout)
        runs += 1
        require(
            serial_index == case.expected_index,
            f"{case.name}: serial index {serial_index}, expected {case.expected_index}",
        )
        for workers in (2, 4):
            topology = ParallelTopology("openmp", (workers,))
            for lease_size in (1, DEFAULT_BRANCH_LEASE_SIZE):
                indices: list[int] = []
                for _ in range(repetitions):
                    index, document = run_solver(
                        telemetry_openmp,
                        case,
                        timeout,
                        topology=topology,
                        branch_lease_size=lease_size,
                        telemetry=True,
                    )
                    runs += 1
                    require(
                        document is not None,
                        f"{case.name}: telemetry document is absent",
                    )
                    require(
                        index == serial_index,
                        f"{case.name}: {workers}-worker/lease-{lease_size} "
                        f"telemetry index {index} does not match serial index "
                        f"{serial_index}",
                    )
                    validate_parallel_telemetry(
                        document,
                        case,
                        topology,
                        expected_lease_size=lease_size,
                    )
                    indices.append(index)
                require(
                    len(set(indices)) == 1,
                    f"{case.name}: repeated {workers}-worker/lease-{lease_size} "
                    f"telemetry runs returned {indices}",
                )
        print(
            f"PASS telemetry {case.name}: 2/4-worker lease-1/4 reductions "
            "and index parity"
        )
    return runs


def run_sparse_adaptive_telemetry_suite(
    serial: Path,
    telemetry_openmp: Path,
    timeout: float,
) -> int:
    """Force the adaptive depth-two path on one small, sparse molecule."""
    case = SPARSE_ADAPTIVE_TELEMETRY_CASE
    serial_index, _ = run_solver(serial, case, timeout)
    require(
        serial_index == case.expected_index,
        f"{case.name}: serial index {serial_index}, expected {case.expected_index}",
    )

    topology = ParallelTopology("openmp", (4,))
    index, document = run_solver(
        telemetry_openmp,
        case,
        timeout,
        topology=topology,
        telemetry=True,
    )
    require(document is not None, f"{case.name}: telemetry document is absent")
    require(
        index == serial_index,
        f"{case.name}: adaptive telemetry index {index} does not match serial "
        f"index {serial_index}",
    )
    validate_parallel_telemetry(
        document,
        case,
        topology,
        require_adaptive_splitting=True,
        require_depth_two_tasks=True,
    )
    print(
        f"PASS telemetry {case.name}: forced-parallel depth-two work and index parity"
    )
    return 2


def run_late_refill_adaptive_telemetry_suite(
    serial: Path,
    telemetry_openmp: Path,
    timeout: float,
) -> int:
    """Require shallow and telemetry-gated deep refill on an irregular tail."""
    case = LATE_REFILL_ADAPTIVE_TELEMETRY_CASE
    serial_index, _ = run_solver(serial, case, timeout)
    require(
        serial_index == case.expected_index,
        f"{case.name}: serial index {serial_index}, expected {case.expected_index}",
    )

    topology = ParallelTopology("openmp", (4,))
    index, document = run_solver(
        telemetry_openmp,
        case,
        timeout,
        topology=topology,
        telemetry=True,
    )
    require(document is not None, f"{case.name}: telemetry document is absent")
    require(
        index == serial_index,
        f"{case.name}: adaptive telemetry index {index} does not match serial "
        f"index {serial_index}",
    )
    validate_parallel_telemetry(
        document,
        case,
        topology,
        require_adaptive_splitting=True,
        require_depth_two_tasks=True,
    )

    parallel = require_mapping(document.get("parallel"), f"{case.name}: parallel")
    aggregate = require_mapping(
        parallel.get("aggregate"), f"{case.name}: parallel.aggregate"
    )
    root_candidates = require_nonnegative_integer(
        aggregate.get("branch_candidates"),
        f"{case.name}: parallel.aggregate.branch_candidates",
    )
    large_frontier_threshold = (
        4
        * topology.worker_count
        * ADAPTIVE_SPLITTING_POLICY["target_queued_tasks_per_worker"]
    )
    require(
        root_candidates >= large_frontier_threshold,
        f"{case.name}: root frontier {root_candidates} must start at or above "
        f"the {large_frontier_threshold}-task large-frontier threshold",
    )
    proactive_tail_refills = require_nonnegative_integer(
        aggregate.get("proactive_tail_refills"),
        f"{case.name}: parallel.aggregate.proactive_tail_refills",
    )
    require(
        proactive_tail_refills > 0,
        f"{case.name}: large root frontier did not trigger a proactive tail refill",
    )
    deeper_tasks = require_nonnegative_integer(
        aggregate.get("deeper_tasks_spawned"),
        f"{case.name}: parallel.aggregate.deeper_tasks_spawned",
    )
    task_steals = require_nonnegative_integer(
        aggregate.get("task_steals"),
        f"{case.name}: parallel.aggregate.task_steals",
    )
    scheduler_idle_waits = require_nonnegative_integer(
        aggregate.get("scheduler_idle_waits"),
        f"{case.name}: parallel.aggregate.scheduler_idle_waits",
    )
    deep_refills = require_nonnegative_integer(
        aggregate.get("deep_refill_activations"),
        f"{case.name}: parallel.aggregate.deep_refill_activations",
    )
    maximum_task_depth = require_nonnegative_integer(
        aggregate.get("maximum_task_depth_executed"),
        f"{case.name}: parallel.aggregate.maximum_task_depth_executed",
    )
    require(
        deeper_tasks > 0
        and task_steals > 0
        and scheduler_idle_waits > 0
        and deep_refills > 0
        and maximum_task_depth > 2,
        f"{case.name}: irregular tail did not prove starvation and expose "
        "stealable deeper work",
    )
    low_watermark = (
        topology.worker_count
        * ADAPTIVE_SPLITTING_POLICY["minimum_queued_tasks_per_worker"]
    )
    print(
        f"PASS telemetry {case.name}: late shallow/deep refill from "
        f"{root_candidates} roots across the {low_watermark}-task low "
        f"watermark and index parity"
    )
    return 2


def run_distributed_telemetry_suite(
    serial: Path,
    target: ParallelTarget,
    cases: Sequence[SolverCase],
    timeout: float,
) -> int:
    selected = [case for case in cases if case.name in TELEMETRY_CASE_NAMES]
    runs = 0
    for case in selected:
        serial_index, _ = run_solver(serial, case, timeout)
        index, document = run_solver(
            target.executable,
            case,
            timeout,
            topology=target.topology,
            mpiexec=target.mpiexec,
            numproc_flag=target.numproc_flag,
            branch_lease_size=target.branch_lease_size,
            telemetry=True,
        )
        runs += 2
        require(
            serial_index == case.expected_index,
            f"{case.name}: serial index {serial_index}, expected {case.expected_index}",
        )
        require(
            index == serial_index,
            f"{case.name}: {target.label} telemetry index {index} does not "
            f"match serial index {serial_index}",
        )
        require(document is not None, f"{case.name}: telemetry document is absent")
        validate_parallel_telemetry(
            document,
            case,
            target.topology,
            expected_lease_size=target.branch_lease_size,
        )
        print(f"PASS telemetry {case.name}: {target.label} reduction and index parity")
    return runs


def parse_arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--serial", type=Path, required=True, help="serial AssemblyCpp target"
    )
    parser.add_argument("--openmp", type=Path, help="OpenMP AssemblyCpp target")
    parser.add_argument(
        "--openmp-telemetry",
        type=Path,
        help="combined OpenMP and telemetry AssemblyCpp target",
    )
    parser.add_argument("--mpi", type=Path, help="two-rank MPI AssemblyCpp target")
    parser.add_argument(
        "--mpi-telemetry",
        type=Path,
        help="combined MPI and telemetry AssemblyCpp target",
    )
    parser.add_argument(
        "--hybrid",
        type=Path,
        help="two-rank, two-thread hybrid AssemblyCpp target",
    )
    parser.add_argument(
        "--hybrid-telemetry",
        type=Path,
        help="combined hybrid and telemetry AssemblyCpp target",
    )
    parser.add_argument(
        "--mpiexec",
        type=Path,
        help="MPI launcher required by MPI and hybrid targets",
    )
    parser.add_argument(
        "--mpiexec-numproc-flag",
        "--numproc-flag",
        dest="mpiexec_numproc_flag",
        default="-n",
        help="launcher flag placed immediately before the rank count (default: -n)",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--repetitions", type=positive_int, default=3)
    parser.add_argument("--telemetry-repetitions", type=positive_int, default=2)
    parser.add_argument("--timeout", type=positive_float, default=30.0)
    parsed = parser.parse_args(arguments)
    target_names = (
        "openmp",
        "openmp_telemetry",
        "mpi",
        "mpi_telemetry",
        "hybrid",
        "hybrid_telemetry",
    )
    if all(getattr(parsed, name) is None for name in target_names):
        parser.error("at least one parallel solver target is required")
    distributed_names = ("mpi", "mpi_telemetry", "hybrid", "hybrid_telemetry")
    distributed_requested = any(
        getattr(parsed, name) is not None for name in distributed_names
    )
    if distributed_requested and parsed.mpiexec is None:
        parser.error("--mpiexec is required with MPI or hybrid targets")
    if not parsed.mpiexec_numproc_flag:
        parser.error("--mpiexec-numproc-flag must not be empty")
    parsed.serial = resolve_path(parsed.serial)
    for name in target_names:
        value = getattr(parsed, name)
        if value is not None:
            setattr(parsed, name, resolve_path(value))
    if parsed.mpiexec is not None:
        parsed.mpiexec = resolve_command_path(parsed.mpiexec)
    return parsed


def main(arguments: Sequence[str] | None = None) -> int:
    try:
        options = parse_arguments(arguments)
        status = executable_status(
            (
                ("serial", options.serial),
                ("openmp", options.openmp),
                ("openmp-telemetry", options.openmp_telemetry),
                ("mpi", options.mpi),
                ("mpi-telemetry", options.mpi_telemetry),
                ("hybrid", options.hybrid),
                ("hybrid-telemetry", options.hybrid_telemetry),
                ("mpiexec", options.mpiexec),
            )
        )
        if status is not None:
            return status
        cases = load_cases(options.manifest)
        runs = 0
        if options.openmp is not None:
            openmp_target = ParallelTarget(
                "openmp-2",
                options.openmp,
                ParallelTopology("openmp", (2,)),
            )
            runs += run_parity_suite(
                options.serial,
                options.openmp,
                cases,
                options.repetitions,
                options.timeout,
            )
            runs += run_pathway_parity_suite(openmp_target, options.timeout)
            runs += run_openmp_execution_policy_suite(
                options.openmp,
                options.timeout,
            )
        if options.openmp_telemetry is not None:
            runs += run_telemetry_suite(
                options.serial,
                options.openmp_telemetry,
                cases,
                options.telemetry_repetitions,
                options.timeout,
            )
            runs += run_sparse_adaptive_telemetry_suite(
                options.serial,
                options.openmp_telemetry,
                options.timeout,
            )
            runs += run_late_refill_adaptive_telemetry_suite(
                options.serial,
                options.openmp_telemetry,
                options.timeout,
            )
        if options.mpi is not None:
            mpi_target = ParallelTarget(
                "mpi-2",
                options.mpi,
                MPI_TOPOLOGY,
                options.mpiexec,
                options.mpiexec_numproc_flag,
            )
            runs += run_distributed_parity_suite(
                options.serial,
                mpi_target,
                cases,
                options.timeout,
            )
            runs += run_pathway_parity_suite(mpi_target, options.timeout)
        if options.mpi_telemetry is not None:
            runs += run_distributed_telemetry_suite(
                options.serial,
                ParallelTarget(
                    "mpi-2",
                    options.mpi_telemetry,
                    MPI_TOPOLOGY,
                    options.mpiexec,
                    options.mpiexec_numproc_flag,
                ),
                cases,
                options.timeout,
            )
        if options.hybrid is not None:
            hybrid_target = ParallelTarget(
                "hybrid-2x2",
                options.hybrid,
                HYBRID_TOPOLOGY,
                options.mpiexec,
                options.mpiexec_numproc_flag,
            )
            runs += run_distributed_parity_suite(
                options.serial,
                hybrid_target,
                cases,
                options.timeout,
            )
            runs += run_pathway_parity_suite(hybrid_target, options.timeout)
        if options.hybrid_telemetry is not None:
            runs += run_distributed_telemetry_suite(
                options.serial,
                ParallelTarget(
                    "hybrid-2x2",
                    options.hybrid_telemetry,
                    HYBRID_TOPOLOGY,
                    options.mpiexec,
                    options.mpiexec_numproc_flag,
                ),
                cases,
                options.timeout,
            )
    except TestFailureError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    print(f"PASS: {runs} bounded parallel solver run(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
