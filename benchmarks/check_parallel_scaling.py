"""Report wall-time speedup and efficiency from paired parallel benchmarks."""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

PAIRED_COMPARISON_ORDER = (
    "baseline/candidate on odd rounds, candidate/baseline on even rounds"
)
ENVIRONMENT_KEY_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
MPI_RANK_FLAGS = frozenset(("-n", "-np", "--n", "--np", "--ntasks"))
MPI_COMPACT_RANK_PATTERN = re.compile(r"^-(?:n|np)([0-9]+)$")
MPI_EQUALS_RANK_PATTERN = re.compile(r"^--(?:n|np|ntasks)=([0-9]+)$")


class ScalingError(RuntimeError):
    """Raised when parallel scaling evidence is malformed or incompatible."""


@dataclass(frozen=True)
class TopologySpec:
    label: str
    workers: int
    path: Path


@dataclass(frozen=True)
class CorpusIdentity:
    manifest_sha256: str
    inputs: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class ExecutionIdentity:
    launcher: tuple[str, ...]
    environment: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class CaseScaling:
    name: str
    speedup: float


@dataclass(frozen=True)
class ScalingResult:
    spec: TopologySpec
    suite: str
    corpus: CorpusIdentity
    baseline_sha256: str
    baseline_execution: ExecutionIdentity
    candidate_sha256: str
    cases: tuple[CaseScaling, ...]
    suite_speedup: float


@dataclass(frozen=True)
class Regression:
    label: str
    workers: int
    case_name: str
    speedup: float


def parse_topology_spec(value: str) -> TopologySpec:
    parts = value.split(":", 2)
    if len(parts) != 3:
        raise ScalingError(
            f"invalid topology {value!r}: expected LABEL:WORKERS:PATH"
        )

    label, workers_text, path_text = (part.strip() for part in parts)
    if not label:
        raise ScalingError(f"invalid topology {value!r}: label must not be empty")
    if not workers_text.isdigit():
        raise ScalingError(
            f"invalid topology {value!r}: workers must be a positive integer"
        )
    workers = int(workers_text)
    if workers < 1:
        raise ScalingError(
            f"invalid topology {value!r}: workers must be a positive integer"
        )
    if not path_text:
        raise ScalingError(f"invalid topology {value!r}: path must not be empty")
    return TopologySpec(label, workers, Path(path_text).expanduser())


def string_at(document: object, keys: Sequence[str], context: str) -> str:
    value = document
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise ScalingError(f"missing {context}")
        value = value[key]
    if not isinstance(value, str) or not value.strip():
        raise ScalingError(f"invalid {context}: expected a non-empty string")
    return value.strip()


def positive_number_at(
    document: object,
    keys: Sequence[str],
    context: str,
) -> float:
    value = document
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise ScalingError(f"missing {context}")
        value = value[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScalingError(f"invalid {context}: expected a positive number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ScalingError(
            f"invalid {context}: expected a finite positive number"
        ) from error
    if not math.isfinite(result) or result <= 0:
        raise ScalingError(
            f"invalid {context}: expected a finite positive number"
        )
    return result


def execution_identity(
    document: dict[str, object],
    role: str,
    path: Path,
) -> ExecutionIdentity:
    execution = document.get("execution")
    if execution is None:
        # Early schema-v2 reports always launched both roles directly.
        return ExecutionIdentity((), ())
    if not isinstance(execution, dict):
        raise ScalingError(f"invalid execution configurations in {path}")
    config = execution.get(role)
    if not isinstance(config, dict):
        raise ScalingError(f"missing {role} execution configuration in {path}")

    launcher = config.get("launcher")
    if not isinstance(launcher, list) or any(
        not isinstance(argument, str) or not argument for argument in launcher
    ):
        raise ScalingError(f"invalid {role} launcher configuration in {path}")
    environment = config.get("environment")
    if not isinstance(environment, dict):
        raise ScalingError(f"invalid {role} environment configuration in {path}")

    normalized_environment = []
    for key, value in environment.items():
        if (
            not isinstance(key, str)
            or ENVIRONMENT_KEY_PATTERN.fullmatch(key) is None
            or not isinstance(value, str)
            or "\x00" in value
        ):
            raise ScalingError(f"invalid {role} environment configuration in {path}")
        normalized_environment.append((key, value))
    return ExecutionIdentity(
        tuple(launcher),
        tuple(sorted(normalized_environment)),
    )


def positive_worker_value(value: str, context: str) -> int:
    if not value.isdigit() or int(value) < 1:
        raise ScalingError(f"invalid {context}: expected a positive integer")
    return int(value)


def mpi_rank_count(execution: ExecutionIdentity, context: str) -> int:
    rank_count: int | None = None
    index = 0
    while index < len(execution.launcher):
        argument = execution.launcher[index]
        parsed: int | None = None
        if argument in MPI_RANK_FLAGS:
            if index + 1 >= len(execution.launcher):
                raise ScalingError(f"missing MPI rank count in {context}")
            parsed = positive_worker_value(
                execution.launcher[index + 1],
                f"MPI rank count in {context}",
            )
            index += 1
        else:
            match = MPI_COMPACT_RANK_PATTERN.fullmatch(argument)
            if match is None:
                match = MPI_EQUALS_RANK_PATTERN.fullmatch(argument)
            if match is not None:
                parsed = positive_worker_value(
                    match.group(1),
                    f"MPI rank count in {context}",
                )
        if parsed is not None:
            if rank_count is not None:
                raise ScalingError(f"multiple MPI rank counts in {context}")
            rank_count = parsed
        index += 1
    return 1 if rank_count is None else rank_count


def execution_worker_count(execution: ExecutionIdentity, context: str) -> int:
    environment = dict(execution.environment)
    thread_text = environment.get("OMP_NUM_THREADS")
    threads = (
        1
        if thread_text is None
        else positive_worker_value(thread_text, f"OMP_NUM_THREADS in {context}")
    )
    thread_limit_text = environment.get("OMP_THREAD_LIMIT")
    if thread_limit_text is not None:
        thread_limit = positive_worker_value(
            thread_limit_text,
            f"OMP_THREAD_LIMIT in {context}",
        )
        threads = min(threads, thread_limit)
    return mpi_rank_count(execution, context) * threads


def load_result(path: Path) -> dict[str, object]:
    try:
        with path.open(encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ScalingError(
            f"could not read benchmark result {path}: {error}"
        ) from error
    if not isinstance(document, dict):
        raise ScalingError(f"invalid benchmark result {path}: expected an object")
    if document.get("schema_version") != 2:
        raise ScalingError(
            f"invalid benchmark result {path}: expected schema_version 2"
        )
    return document


def corpus_identity(document: dict[str, object], path: Path) -> CorpusIdentity:
    corpus = document.get("corpus")
    if not isinstance(corpus, dict):
        raise ScalingError(f"missing corpus fingerprints in {path}")
    manifest_sha256 = string_at(
        corpus,
        ("manifest", "sha256"),
        f"manifest SHA-256 in {path}",
    )
    inputs = corpus.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ScalingError(f"missing corpus input fingerprints in {path}")

    fingerprints: dict[str, str] = {}
    for index, entry in enumerate(inputs):
        if not isinstance(entry, dict):
            raise ScalingError(f"invalid corpus input fingerprint {index} in {path}")
        name = string_at(entry, ("name",), f"corpus input {index} name in {path}")
        sha256 = string_at(
            entry,
            ("sha256",),
            f"corpus input {name!r} SHA-256 in {path}",
        )
        if name in fingerprints:
            raise ScalingError(f"duplicate corpus input fingerprint {name!r} in {path}")
        fingerprints[name] = sha256
    return CorpusIdentity(manifest_sha256, tuple(sorted(fingerprints.items())))


def parse_wall_samples(
    case: dict[str, object],
    role: str,
    runs: int,
    expected_assembly_index: int,
    case_name: str,
    path: Path,
) -> tuple[float, ...]:
    container = case.get(role)
    if not isinstance(container, dict):
        raise ScalingError(f"missing {role} measurements for {case_name!r} in {path}")
    samples = container.get("measurements")
    if not isinstance(samples, list):
        raise ScalingError(f"missing {role} measurements for {case_name!r} in {path}")
    if len(samples) != runs:
        raise ScalingError(
            f"wrong {role} sample count for {case_name!r} in {path}: "
            f"expected {runs}, got {len(samples)}"
        )

    wall_samples = []
    for expected_round, sample in enumerate(samples, start=1):
        if not isinstance(sample, dict):
            raise ScalingError(
                f"invalid {role} sample {expected_round} for {case_name!r} in {path}"
            )
        if sample.get("round") != expected_round or isinstance(
            sample.get("round"), bool
        ):
            raise ScalingError(
                f"invalid {role} round order for {case_name!r} in {path}"
            )
        wall_seconds = positive_number_at(
            sample,
            ("wall_seconds",),
            f"{role} wall sample for {case_name!r} in {path}",
        )
        assembly_index = sample.get("assembly_index")
        if (
            isinstance(assembly_index, bool)
            or not isinstance(assembly_index, int)
            or assembly_index != expected_assembly_index
        ):
            raise ScalingError(
                f"wrong {role} assembly index for {case_name!r} in {path}"
            )
        wall_samples.append(wall_seconds)
    return tuple(wall_samples)


def paired_median(baseline: Sequence[float], candidate: Sequence[float]) -> float:
    return statistics.median(
        baseline_value / candidate_value
        for baseline_value, candidate_value in zip(baseline, candidate, strict=True)
    )


def require_recorded_median(
    document: object,
    keys: Sequence[str],
    expected: float,
    context: str,
) -> None:
    recorded = positive_number_at(document, keys, context)
    if not math.isclose(recorded, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise ScalingError(f"inconsistent {context}: raw samples give {expected:.12g}")


def evaluate_report(spec: TopologySpec) -> ScalingResult:
    document = load_result(spec.path)
    suite = string_at(document, ("suite",), f"suite in {spec.path}")
    baseline_sha256 = string_at(
        document,
        ("executables", "baseline", "sha256"),
        f"baseline SHA-256 in {spec.path}",
    )
    candidate_sha256 = string_at(
        document,
        ("executables", "candidate", "sha256"),
        f"candidate SHA-256 in {spec.path}",
    )
    baseline_execution = execution_identity(document, "baseline", spec.path)
    candidate_execution = execution_identity(document, "candidate", spec.path)
    if (
        baseline_sha256 == candidate_sha256
        and baseline_execution == candidate_execution
    ):
        raise ScalingError(
            "candidate and baseline binary and execution configurations must "
            f"differ in {spec.path}"
        )

    baseline_workers = execution_worker_count(
        baseline_execution,
        f"baseline execution configuration in {spec.path}",
    )
    if baseline_workers != 1:
        raise ScalingError(
            f"baseline in {spec.path} is not serial: found {baseline_workers} workers"
        )
    candidate_workers = execution_worker_count(
        candidate_execution,
        f"candidate execution configuration in {spec.path}",
    )
    if candidate_workers != spec.workers:
        raise ScalingError(
            f"worker count mismatch for {spec.label!r} in {spec.path}: topology "
            f"declares {spec.workers}, execution configuration provides "
            f"{candidate_workers}"
        )

    comparison_order = string_at(
        document,
        ("schedule", "comparison_order"),
        f"comparison order in {spec.path}",
    )
    if comparison_order != PAIRED_COMPARISON_ORDER:
        raise ScalingError(f"unsupported paired comparison order in {spec.path}")

    runs = document.get("runs")
    if isinstance(runs, bool) or not isinstance(runs, int) or runs < 1:
        raise ScalingError(
            f"invalid run count in {spec.path}: expected a positive integer"
        )

    identity = corpus_identity(document, spec.path)
    cases = document.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ScalingError(f"missing benchmark cases in {spec.path}")

    case_scaling = []
    wall_samples: list[tuple[tuple[float, ...], tuple[float, ...]]] = []
    seen_names: set[str] = set()
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ScalingError(f"invalid benchmark case {index} in {spec.path}")
        name = string_at(case, ("name",), f"benchmark case {index} name in {spec.path}")
        if name in seen_names:
            raise ScalingError(f"duplicate benchmark case {name!r} in {spec.path}")
        seen_names.add(name)
        expected_assembly_index = case.get("expected_assembly_index")
        if isinstance(expected_assembly_index, bool) or not isinstance(
            expected_assembly_index, int
        ):
            raise ScalingError(
                f"invalid expected assembly index for {name!r} in {spec.path}"
            )

        candidate = parse_wall_samples(
            case,
            "candidate",
            runs,
            expected_assembly_index,
            name,
            spec.path,
        )
        baseline = parse_wall_samples(
            case,
            "baseline",
            runs,
            expected_assembly_index,
            name,
            spec.path,
        )
        speedup = paired_median(baseline, candidate)
        require_recorded_median(
            case,
            ("comparison", "paired_wall_speedup", "median"),
            speedup,
            f"paired wall median for {name!r} in {spec.path}",
        )
        case_scaling.append(CaseScaling(name, speedup))
        wall_samples.append((baseline, candidate))

    corpus_names = {name for name, _ in identity.inputs}
    if seen_names != corpus_names:
        raise ScalingError(
            f"benchmark case coverage does not match corpus in {spec.path}"
        )

    round_ratios = []
    for round_index in range(runs):
        baseline_total = sum(samples[0][round_index] for samples in wall_samples)
        candidate_total = sum(samples[1][round_index] for samples in wall_samples)
        if not math.isfinite(baseline_total) or not math.isfinite(candidate_total):
            raise ScalingError(f"non-finite suite round total in {spec.path}")
        round_ratios.append(baseline_total / candidate_total)
    suite_speedup = statistics.median(round_ratios)
    require_recorded_median(
        document,
        ("comparison", "paired_round_wall_speedup", "median"),
        suite_speedup,
        f"suite paired round-total wall median in {spec.path}",
    )
    return ScalingResult(
        spec,
        suite,
        identity,
        baseline_sha256,
        baseline_execution,
        candidate_sha256,
        tuple(case_scaling),
        suite_speedup,
    )


def evaluate_specs(specs: Sequence[TopologySpec]) -> list[ScalingResult]:
    if not specs:
        raise ScalingError("at least one topology report is required")

    seen_topologies: set[tuple[str, int]] = set()
    seen_paths: set[Path] = set()
    for spec in specs:
        key = (spec.label.casefold(), spec.workers)
        if key in seen_topologies:
            raise ScalingError(
                f"duplicate topology worker count: {spec.label}:{spec.workers}"
            )
        seen_topologies.add(key)
        resolved_path = spec.path.resolve()
        if resolved_path in seen_paths:
            raise ScalingError(f"benchmark report used more than once: {spec.path}")
        seen_paths.add(resolved_path)

    results = [evaluate_report(spec) for spec in specs]
    reference = results[0]
    candidates_by_label: dict[str, tuple[str, Path]] = {}
    for result in results[1:]:
        if result.suite != reference.suite:
            raise ScalingError(
                f"incompatible suite in {result.spec.path}: expected "
                f"{reference.suite!r}, got {result.suite!r}"
            )
        if result.corpus != reference.corpus:
            raise ScalingError(
                f"incompatible corpus fingerprints in {result.spec.path}"
            )
        if result.baseline_sha256 != reference.baseline_sha256:
            raise ScalingError(
                f"incompatible baseline fingerprint in {result.spec.path}"
            )
        if result.baseline_execution != reference.baseline_execution:
            raise ScalingError(
                f"incompatible baseline execution configuration in {result.spec.path}"
            )
    for result in results:
        label = result.spec.label.casefold()
        previous = candidates_by_label.get(label)
        if previous is not None and result.candidate_sha256 != previous[0]:
            raise ScalingError(
                f"incompatible candidate fingerprint for topology "
                f"{result.spec.label!r} in {previous[1]} and {result.spec.path}"
            )
        candidates_by_label[label] = (
            result.candidate_sha256,
            result.spec.path,
        )
    return sorted(
        results,
        key=lambda result: (result.spec.label.casefold(), result.spec.workers),
    )


def regressions(results: Sequence[ScalingResult]) -> list[Regression]:
    return [
        Regression(result.spec.label, result.spec.workers, case.name, case.speedup)
        for result in results
        for case in result.cases
        if case.speedup <= 1.0
    ]


def print_report(results: Sequence[ScalingResult]) -> list[Regression]:
    print(f"Parallel scaling for suite {results[0].suite}")
    print(f"Baseline SHA-256: {results[0].baseline_sha256}")
    print()
    print("Topology        Workers  Subject                     Speedup  Efficiency")
    print("------------------------------------------------------------------------")
    for result in results:
        label = result.spec.label
        workers = result.spec.workers
        print(
            f"{label:<15} {workers:>7}  {'suite round-total':<26} "
            f"{result.suite_speedup:>7.3f}x "
            f"{100 * result.suite_speedup / workers:>9.1f}%"
        )
        for case in result.cases:
            print(
                f"{label:<15} {workers:>7}  {f'case {case.name}':<26} "
                f"{case.speedup:>7.3f}x {100 * case.speedup / workers:>9.1f}%"
            )

    failures = regressions(results)
    print()
    if failures:
        print("Regressions (paired wall speedup <= 1.0x):")
        for failure in failures:
            print(
                f"  {failure.label}/{failure.workers} case {failure.case_name}: "
                f"{failure.speedup:.3f}x "
                f"({100 * failure.speedup / failure.workers:.1f}% efficiency)"
            )
    else:
        print("No case wall regressions detected.")
    return failures


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Report paired wall-time speedup and parallel efficiency for one "
            "suite across labelled worker topologies."
        )
    )
    parser.add_argument(
        "topologies",
        nargs="+",
        metavar="LABEL:WORKERS:PATH",
        help="topology label, positive worker count, and paired schema-v2 report",
    )
    parser.add_argument(
        "--require-all-faster",
        "--require-every-case-faster",
        action="store_true",
        help="fail when any case's paired wall speedup is at or below 1.0",
    )
    return parser


def main(arguments: Sequence[str] | None = None) -> int:
    parsed = create_argument_parser().parse_args(arguments)
    try:
        specs = [parse_topology_spec(value) for value in parsed.topologies]
        results = evaluate_specs(specs)
    except ScalingError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    failures = print_report(results)
    if parsed.require_all_faster:
        if failures:
            print(f"FAIL: {len(failures)} case wall regression(s).")
            return 1
        print("PASS: every case paired wall speedup exceeds 1.0x.")
    elif failures:
        print(
            "REPORT ONLY: regressions do not affect exit status; use "
            "--require-all-faster to enable the wall-time gate."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
