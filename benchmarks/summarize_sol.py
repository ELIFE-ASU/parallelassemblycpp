"""Combine completed Sol benchmark reports into a CSV for quick analysis."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

if __package__:
    from . import benchmark, check_parallel_scaling
else:
    import benchmark
    import check_parallel_scaling


FIELDNAMES = (
    "run_kind",
    "suite",
    "case",
    "topology",
    "workers",
    "runs",
    "wall_median_seconds",
    "wall_mad_seconds",
    "wall_p95_seconds",
    "baseline_wall_median_seconds",
    "paired_speedup",
    "efficiency",
    "report",
)


def report_rows(path: Path, root: Path, *, parallel: bool) -> list[dict[str, object]]:
    """Validate raw measurements and summarize case and suite round totals."""
    checker = check_parallel_scaling
    document = checker.load_result(path)
    suite = checker.string_at(document, ("suite",), f"suite in {path}")
    runs = document.get("runs")
    if type(runs) is not int or runs < 1:
        raise checker.ScalingError(f"invalid run count in {path}")
    execution = checker.execution_identity(document, "candidate", path)
    workers = checker.execution_worker_count(execution, str(path))
    if parallel:
        match = re.fullmatch(r"omp-([0-9]+)\.json", path.name)
        if match is None or int(match[1]) < 2:
            raise checker.ScalingError(f"invalid OpenMP report filename: {path}")
        if checker.mpi_rank_count(execution, str(path)) != 1:
            raise checker.ScalingError(f"OpenMP report must use one process: {path}")
        checker.evaluate_report(checker.TopologySpec("omp", int(match[1]), path))
    elif (
        workers != 1
        or "--parallel=on" in execution.arguments
        or document.get("baseline_aggregate") is not None
        or document.get("comparison") is not None
    ):
        raise checker.ScalingError(
            f"suite report must contain serial measurements: {path}"
        )

    cases = document.get("cases")
    if not isinstance(cases, list) or not cases:
        raise checker.ScalingError(f"missing benchmark cases in {path}")
    candidate_totals = [0.0] * runs
    baseline_totals = [0.0] * runs if parallel else []
    names: set[str] = set()
    rows: list[dict[str, object]] = []

    def make_row(
        name: str, candidate: Sequence[float], baseline: Sequence[float]
    ) -> dict[str, object]:
        wall = benchmark.summarize(candidate)
        speedup = checker.paired_median(baseline, candidate) if baseline else None
        return {
            "run_kind": "parallel" if parallel else "suite",
            "suite": suite,
            "case": name,
            "topology": "omp" if parallel else "serial",
            "workers": workers,
            "runs": runs,
            "wall_median_seconds": wall.median,
            "wall_mad_seconds": wall.mad,
            "wall_p95_seconds": wall.p95,
            "baseline_wall_median_seconds": (
                benchmark.summarize(baseline).median if baseline else None
            ),
            "paired_speedup": speedup,
            "efficiency": None if speedup is None else speedup / workers,
            "report": path.relative_to(root).as_posix(),
        }

    for case in cases:
        if not isinstance(case, dict):
            raise checker.ScalingError(f"invalid benchmark case in {path}")
        name = checker.string_at(case, ("name",), f"case name in {path}")
        if name in names or name == "__suite__":
            raise checker.ScalingError(
                f"duplicate or reserved case name {name!r}: {path}"
            )
        names.add(name)
        expected = case.get("expected_assembly_index")
        if type(expected) is not int:
            raise checker.ScalingError(
                f"invalid expected assembly index: {path}/{name}"
            )
        candidate = checker.parse_wall_samples(
            case, "candidate", runs, expected, name, path
        )
        baseline = (
            checker.parse_wall_samples(case, "baseline", runs, expected, name, path)
            if parallel
            else ()
        )
        if not parallel and (
            case.get("baseline") is not None or case.get("comparison") is not None
        ):
            raise checker.ScalingError(f"unexpected paired measurements: {path}/{name}")
        rows.append(make_row(name, candidate, baseline))
        for index, value in enumerate(candidate):
            candidate_totals[index] += value
        for index, value in enumerate(baseline):
            baseline_totals[index] += value

    if names != {name for name, _ in checker.corpus_identity(document, path).inputs}:
        raise checker.ScalingError(f"case coverage does not match corpus in {path}")
    if not all(math.isfinite(value) for value in candidate_totals + baseline_totals):
        raise checker.ScalingError(f"non-finite suite round total in {path}")
    rows.append(make_row("__suite__", candidate_totals, baseline_totals))
    return rows


def write_summary(root: Path) -> tuple[Path, int]:
    """Replace the CSV only after every discovered report validates."""
    if not root.is_dir():
        raise check_parallel_scaling.ScalingError(f"output directory not found: {root}")
    rows = []
    for directory, pattern, parallel in (
        ("suites", "*.json", False),
        ("parallel", "omp-*.json", True),
    ):
        for path in sorted((root / directory).glob(pattern)):
            rows.extend(report_rows(path, root, parallel=parallel))

    destination = root / "summary.csv"
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=root, delete=False
        ) as stream:
            temporary = Path(stream.name)
            writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return destination, len(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="Sol job results directory")
    arguments = parser.parse_args(argv)
    try:
        path, rows = write_summary(arguments.output_dir.expanduser().resolve())
    except (check_parallel_scaling.ScalingError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(f"Analysis CSV: {path} ({rows} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
