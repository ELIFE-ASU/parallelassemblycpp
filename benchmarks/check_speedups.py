"""Check paired benchmark reports before promoting an optimization."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

if __package__:
    from . import benchmark
else:
    import benchmark


@dataclass(frozen=True)
class GateFailure:
    suite: str
    subject: str
    speedup: float


class GateError(RuntimeError):
    """Raised when benchmark reports cannot be evaluated safely."""


MINIMUM_RUNS = {
    "quick": 100,
    "full": 100,
    "profile": 6,
    "scaling": 30,
}
PAIRED_COMPARISON_ORDER = (
    "baseline/candidate on odd rounds, candidate/baseline on even rounds"
)
MAX_CLOCK_TICKS = (1 << 64) - 1
ExecutionIdentity = tuple[tuple[str, ...], tuple[tuple[str, str], ...]]


def number_at(document: object, keys: Sequence[str], context: str) -> float:
    value = document
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise GateError(f"missing {context}")
        value = value[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise GateError(f"invalid {context}: expected a number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise GateError(f"invalid {context}: expected a finite number") from error
    if not math.isfinite(result):
        raise GateError(f"invalid {context}: expected a finite number")
    return result


def string_at(document: object, keys: Sequence[str], context: str) -> str:
    value = document
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise GateError(f"missing {context}")
        value = value[key]
    if not isinstance(value, str) or not value:
        raise GateError(f"invalid {context}: expected a non-empty string")
    return value


def execution_identity(
    document: dict[str, object], role: str, path: Path
) -> ExecutionIdentity:
    """Return one canonical launcher/environment identity from a report."""
    execution = document.get("execution")
    if execution is None:
        # Schema-v2 reports written before execution configurations existed used
        # the ordinary direct-launch configuration for both roles.
        return (), ()
    if not isinstance(execution, dict):
        raise GateError(f"invalid execution configurations in {path}")
    config = execution.get(role)
    if not isinstance(config, dict):
        raise GateError(f"missing {role} execution configuration in {path}")

    launcher = config.get("launcher")
    if not isinstance(launcher, list) or any(
        not isinstance(value, str) or not value for value in launcher
    ):
        raise GateError(f"invalid {role} launcher configuration in {path}")
    environment = config.get("environment")
    if not isinstance(environment, dict):
        raise GateError(f"invalid {role} environment configuration in {path}")
    normalized_environment = []
    for key, value in environment.items():
        if (
            not isinstance(key, str)
            or benchmark.ENVIRONMENT_KEY_PATTERN.fullmatch(key) is None
            or not isinstance(value, str)
            or "\x00" in value
        ):
            raise GateError(f"invalid {role} environment configuration in {path}")
        normalized_environment.append((key, value))
    return tuple(launcher), tuple(sorted(normalized_environment))


def load_result(path: Path) -> dict[str, object]:
    try:
        with path.expanduser().open(encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise GateError(f"could not read benchmark report {path}: {error}") from error
    if not isinstance(document, dict):
        raise GateError(f"invalid benchmark report {path}: expected a JSON object")
    if document.get("schema_version") != 2:
        raise GateError(
            f"invalid benchmark report {path}: expected schema_version 2"
        )
    return document


def validate_corpus_identity(
    document: dict[str, object],
    path: Path,
    expected: dict[str, object],
) -> None:
    corpus = document.get("corpus")
    if not isinstance(corpus, dict):
        raise GateError(f"missing corpus fingerprints in {path}")

    recorded_manifest_sha256 = string_at(
        corpus,
        ("manifest", "sha256"),
        f"manifest SHA-256 in {path}",
    )
    expected_manifest = expected.get("manifest")
    if not isinstance(expected_manifest, dict):
        raise GateError("internal error: expected benchmark manifest metadata")
    if recorded_manifest_sha256 != expected_manifest.get("sha256"):
        raise GateError(f"stale benchmark manifest fingerprint in {path}")

    recorded_inputs = corpus.get("inputs")
    if not isinstance(recorded_inputs, list):
        raise GateError(f"missing corpus input fingerprints in {path}")
    recorded_by_name: dict[str, str] = {}
    for index, entry in enumerate(recorded_inputs):
        if not isinstance(entry, dict):
            raise GateError(f"invalid corpus input fingerprint {index} in {path}")
        name = string_at(entry, ("name",), f"corpus input {index} name in {path}")
        sha256 = string_at(
            entry,
            ("sha256",),
            f"corpus input {name!r} SHA-256 in {path}",
        )
        if name in recorded_by_name:
            raise GateError(f"duplicate corpus input fingerprint {name!r} in {path}")
        recorded_by_name[name] = sha256

    expected_inputs = expected.get("inputs")
    if not isinstance(expected_inputs, list):
        raise GateError("internal error: expected benchmark input metadata")
    expected_by_name = {
        str(entry["name"]): str(entry["sha256"])
        for entry in expected_inputs
        if isinstance(entry, dict) and "name" in entry and "sha256" in entry
    }
    if recorded_by_name.keys() != expected_by_name.keys():
        raise GateError(f"wrong corpus input fingerprint coverage in {path}")
    changed = sorted(
        name
        for name, sha256 in recorded_by_name.items()
        if sha256 != expected_by_name[name]
    )
    if changed:
        raise GateError(
            f"stale benchmark input fingerprint(s) in {path}: {', '.join(changed)}"
        )


def parse_measurements(
    case_result: dict[str, object],
    role: str,
    runs: int,
    expected_assembly_index: int,
    case_name: str,
    path: Path,
) -> tuple[benchmark.Measurement, ...]:
    container = case_result.get(role)
    if not isinstance(container, dict):
        raise GateError(f"missing {role} measurements for {case_name!r} in {path}")
    samples = container.get("measurements")
    if not isinstance(samples, list):
        raise GateError(f"missing {role} samples for {case_name!r} in {path}")
    if len(samples) != runs:
        raise GateError(
            f"wrong {role} sample count for {case_name!r} in {path}: "
            f"expected exactly {runs}, got {len(samples)}"
        )

    measurements = []
    for expected_round, sample in enumerate(samples, start=1):
        if not isinstance(sample, dict):
            raise GateError(
                f"invalid {role} sample {expected_round} for {case_name!r} in {path}"
            )
        round_number = sample.get("round")
        if type(round_number) is not int or round_number != expected_round:
            raise GateError(
                f"invalid {role} round order for {case_name!r} in {path}"
            )
        wall_value = sample.get("wall_seconds")
        try:
            wall_seconds = float(wall_value)
        except (OverflowError, TypeError, ValueError):
            wall_seconds = math.nan
        if (
            isinstance(wall_value, bool)
            or not isinstance(wall_value, (int, float))
            or not math.isfinite(wall_seconds)
            or wall_seconds <= 0
        ):
            raise GateError(
                f"invalid {role} wall sample for {case_name!r} in {path}"
            )
        clock_ticks = sample.get("clock_ticks")
        if (
            type(clock_ticks) is not int
            or clock_ticks <= 0
            or clock_ticks > MAX_CLOCK_TICKS
        ):
            raise GateError(
                f"invalid {role} clock sample for {case_name!r} in {path}"
            )
        assembly_index = sample.get("assembly_index")
        if type(assembly_index) is not int or assembly_index != expected_assembly_index:
            raise GateError(
                f"wrong {role} assembly index for {case_name!r} in {path}"
            )
        measurements.append(
            benchmark.Measurement(wall_seconds, clock_ticks, assembly_index)
        )
    return tuple(measurements)


def require_recorded_median(
    document: object,
    keys: Sequence[str],
    expected: float,
    context: str,
) -> None:
    recorded = number_at(document, keys, context)
    if not math.isclose(recorded, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise GateError(f"inconsistent {context}: raw samples give {expected:.12g}")


def evaluate_results(
    paths: Sequence[Path],
    threshold: float,
) -> list[GateFailure]:
    manifest_path, manifest_cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
    cases_by_suite = {
        suite: [case for case in manifest_cases if suite in case.suites]
        for suite in benchmark.KNOWN_SUITES
    }
    expected_by_suite = {
        suite: {case.name for case in cases}
        for suite, cases in cases_by_suite.items()
    }
    corpus_by_suite = {
        suite: benchmark.benchmark_corpus_metadata(manifest_path, cases)
        for suite, cases in cases_by_suite.items()
    }
    documents: dict[str, tuple[Path, dict[str, object], int]] = {}
    candidate_hash: str | None = None
    baseline_hash: str | None = None
    candidate_execution: ExecutionIdentity | None = None
    baseline_execution: ExecutionIdentity | None = None

    for path in paths:
        document = load_result(path)
        suite = string_at(document, ("suite",), f"suite in {path}")
        if suite not in expected_by_suite:
            raise GateError(f"unknown suite {suite!r} in {path}")
        if suite in documents:
            raise GateError(
                f"duplicate suite {suite!r} in {documents[suite][0]} and {path}"
            )

        runs = document.get("runs")
        if isinstance(runs, bool) or not isinstance(runs, int):
            raise GateError(f"invalid run count in {path}: expected an integer")
        if runs < MINIMUM_RUNS[suite]:
            raise GateError(
                f"insufficient run count in {path}: {suite} requires at least "
                f"{MINIMUM_RUNS[suite]}, got {runs}"
            )
        if runs % 2:
            raise GateError(f"unbalanced odd run count in {path}: got {runs}")
        comparison_order = string_at(
            document,
            ("schedule", "comparison_order"),
            f"paired comparison order in {path}",
        )
        if comparison_order != PAIRED_COMPARISON_ORDER:
            raise GateError(f"invalid paired comparison order in {path}")
        validate_corpus_identity(document, path, corpus_by_suite[suite])

        current_candidate = string_at(
            document,
            ("executables", "candidate", "sha256"),
            f"candidate SHA-256 in {path}",
        )
        current_baseline = string_at(
            document,
            ("executables", "baseline", "sha256"),
            f"baseline SHA-256 in {path}",
        )
        current_candidate_execution = execution_identity(document, "candidate", path)
        current_baseline_execution = execution_identity(document, "baseline", path)
        if candidate_hash is None:
            candidate_hash = current_candidate
            baseline_hash = current_baseline
            candidate_execution = current_candidate_execution
            baseline_execution = current_baseline_execution
        elif (
            current_candidate != candidate_hash
            or current_baseline != baseline_hash
            or current_candidate_execution != candidate_execution
            or current_baseline_execution != baseline_execution
        ):
            raise GateError(
                "all reports must use the same candidate and baseline binary "
                "and execution configurations"
            )
        documents[suite] = (path, document, runs)

    missing_suites = [suite for suite in benchmark.KNOWN_SUITES if suite not in documents]
    if missing_suites:
        raise GateError(f"missing benchmark suites: {', '.join(missing_suites)}")
    if candidate_hash == baseline_hash and candidate_execution == baseline_execution:
        raise GateError(
            "candidate and baseline binary/execution identities must differ"
        )

    failures: list[GateFailure] = []
    for suite in benchmark.KNOWN_SUITES:
        path, document, runs = documents[suite]
        case_documents = document.get("cases")
        if not isinstance(case_documents, list):
            raise GateError(f"invalid cases in {path}: expected a list")
        results_by_name: dict[
            str,
            tuple[tuple[benchmark.Measurement, ...], tuple[benchmark.Measurement, ...]],
        ] = {}
        case_results_by_name: dict[str, dict[str, object]] = {}
        for index, case_result in enumerate(case_documents):
            if not isinstance(case_result, dict):
                raise GateError(f"invalid case {index} in {path}: expected an object")
            name = string_at(case_result, ("name",), f"case {index} name in {path}")
            if name in case_results_by_name:
                raise GateError(f"duplicate case {name!r} in {path}")
            case_results_by_name[name] = case_result

        actual_names = set(case_results_by_name)
        expected_names = expected_by_suite[suite]
        if actual_names != expected_names:
            missing = sorted(expected_names - actual_names)
            extra = sorted(actual_names - expected_names)
            details = []
            if missing:
                details.append(f"missing {', '.join(missing)}")
            if extra:
                details.append(f"unexpected {', '.join(extra)}")
            raise GateError(f"wrong {suite} case coverage in {path}: {'; '.join(details)}")

        for current_case in cases_by_suite[suite]:
            name = current_case.name
            case_result = case_results_by_name[name]
            reported_expected = case_result.get("expected_assembly_index")
            if (
                type(reported_expected) is not int
                or reported_expected != current_case.expected_assembly_index
            ):
                raise GateError(f"stale expected assembly index for {name!r} in {path}")
            assert current_case.expected_assembly_index is not None
            candidate = parse_measurements(
                case_result,
                "candidate",
                runs,
                current_case.expected_assembly_index,
                name,
                path,
            )
            baseline = parse_measurements(
                case_result,
                "baseline",
                runs,
                current_case.expected_assembly_index,
                name,
                path,
            )
            results_by_name[name] = (candidate, baseline)
            clock_summary = benchmark.paired_speedup_summary(
                candidate, baseline, "clock_ticks"
            )
            assert clock_summary is not None
            require_recorded_median(
                case_result,
                ("comparison", "paired_clock_speedup", "median"),
                clock_summary.median,
                f"paired clock median for {name!r} in {path}",
            )
            if clock_summary.median <= threshold:
                failures.append(
                    GateFailure(suite, f"case {name} clock", clock_summary.median)
                )

        round_wall_ratios = []
        round_clock_ratios = []
        for round_index in range(runs):
            candidate_wall = sum(
                results_by_name[case.name][0][round_index].wall_seconds
                for case in cases_by_suite[suite]
            )
            baseline_wall = sum(
                results_by_name[case.name][1][round_index].wall_seconds
                for case in cases_by_suite[suite]
            )
            candidate_clock = sum(
                results_by_name[case.name][0][round_index].clock_ticks
                for case in cases_by_suite[suite]
            )
            baseline_clock = sum(
                results_by_name[case.name][1][round_index].clock_ticks
                for case in cases_by_suite[suite]
            )
            round_wall_ratios.append(baseline_wall / candidate_wall)
            round_clock_ratios.append(baseline_clock / candidate_clock)

        for metric, label, ratios in (
            ("paired_round_wall_speedup", "suite round-total wall", round_wall_ratios),
            (
                "paired_round_clock_speedup",
                "suite round-total clock",
                round_clock_ratios,
            ),
        ):
            speedup = benchmark.summarize(ratios).median
            require_recorded_median(
                document,
                ("comparison", metric, "median"),
                speedup,
                f"{label} median in {path}",
            )
            if speedup <= threshold:
                failures.append(GateFailure(suite, label, speedup))
    return failures


def threshold_value(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("must be a finite non-negative number")
    return parsed


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fail unless every benchmark case and suite aggregate is faster "
            "in paired v2 reports."
        )
    )
    parser.add_argument(
        "results",
        type=Path,
        nargs="+",
        help="paired JSON reports for quick, full, profile, and scaling",
    )
    parser.add_argument(
        "--threshold",
        type=threshold_value,
        default=1.0,
        help="require every median to exceed this value (default: 1)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)
    try:
        failures = evaluate_results(arguments.results, arguments.threshold)
    except (GateError, benchmark.BenchmarkError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if failures:
        for failure in failures:
            print(
                f"FAIL {failure.suite}: {failure.subject} "
                f"speedup {failure.speedup:.6f} <= {arguments.threshold:.6f}"
            )
        return 1

    print(
        "PASS: all case clock medians and suite wall/clock medians exceed "
        f"{arguments.threshold:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
