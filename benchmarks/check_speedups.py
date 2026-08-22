"""Validate paired benchmark JSON files as an optimization promotion gate."""

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
    """Raised when benchmark results cannot be evaluated safely."""


MINIMUM_RUNS = {
    "quick": 100,
    "full": 100,
    "profile": 6,
    "scaling": 30,
}
PAIRED_COMPARISON_ORDER = (
    "baseline/candidate on odd rounds, candidate/baseline on even rounds"
)


def number_at(document: object, keys: Sequence[str], context: str) -> float:
    value = document
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise GateError(f"missing {context}")
        value = value[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise GateError(f"invalid {context}: expected a number")
    result = float(value)
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


def load_result(path: Path) -> dict[str, object]:
    try:
        with path.expanduser().open(encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise GateError(f"could not read benchmark result {path}: {error}") from error
    if not isinstance(document, dict):
        raise GateError(f"invalid benchmark result {path}: expected a JSON object")
    if document.get("schema_version") != 2:
        raise GateError(
            f"invalid benchmark result {path}: expected schema_version 2"
        )
    return document


def evaluate_results(
    paths: Sequence[Path],
    threshold: float,
) -> list[GateFailure]:
    _, manifest_cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
    expected_by_suite = {
        suite: {case.name for case in manifest_cases if suite in case.suites}
        for suite in benchmark.KNOWN_SUITES
    }
    documents: dict[str, tuple[Path, dict[str, object]]] = {}
    candidate_hash: str | None = None
    baseline_hash: str | None = None

    for path in paths:
        document = load_result(path)
        suite = string_at(document, ("suite",), f"suite in {path}")
        if suite not in expected_by_suite:
            raise GateError(f"unknown suite {suite!r} in {path}")
        if suite in documents:
            raise GateError(
                f"duplicate suite {suite!r} in {documents[suite][0]} and {path}"
            )
        documents[suite] = (path, document)

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
        if candidate_hash is None:
            candidate_hash = current_candidate
            baseline_hash = current_baseline
        elif current_candidate != candidate_hash or current_baseline != baseline_hash:
            raise GateError(
                "all suite results must use the same candidate and baseline binaries"
            )

    missing_suites = [suite for suite in benchmark.KNOWN_SUITES if suite not in documents]
    if missing_suites:
        raise GateError(f"missing benchmark suite result(s): {', '.join(missing_suites)}")
    if candidate_hash == baseline_hash:
        raise GateError("candidate and baseline binaries must be different")

    failures: list[GateFailure] = []
    for suite in benchmark.KNOWN_SUITES:
        path, document = documents[suite]
        cases = document.get("cases")
        if not isinstance(cases, list):
            raise GateError(f"invalid cases in {path}: expected a list")
        actual_names: set[str] = set()
        for index, case_result in enumerate(cases):
            if not isinstance(case_result, dict):
                raise GateError(f"invalid case {index} in {path}: expected an object")
            name = string_at(case_result, ("name",), f"case {index} name in {path}")
            if name in actual_names:
                raise GateError(f"duplicate case {name!r} in {path}")
            actual_names.add(name)
            speedup = number_at(
                case_result,
                ("comparison", "paired_clock_speedup", "median"),
                f"paired clock median for {name!r} in {path}",
            )
            if speedup <= threshold:
                failures.append(GateFailure(suite, f"case {name} clock", speedup))

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

        for metric, label in (
            ("paired_round_wall_speedup", "suite round-total wall"),
            ("paired_round_clock_speedup", "suite round-total clock"),
        ):
            speedup = number_at(
                document,
                ("comparison", metric, "median"),
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
            "Require every maintained benchmark case and suite aggregate to be "
            "faster in paired schema-v2 JSON results."
        )
    )
    parser.add_argument(
        "results",
        type=Path,
        nargs="+",
        help="one paired JSON result for each of quick, full, profile, and scaling",
    )
    parser.add_argument(
        "--threshold",
        type=threshold_value,
        default=1.0,
        help="speedup medians must be strictly greater than this value (default: 1)",
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
        "PASS: every case clock median and every suite round-total wall/clock "
        f"median exceeds {arguments.threshold:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
