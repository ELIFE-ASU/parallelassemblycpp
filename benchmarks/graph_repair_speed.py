"""Compare exact search time with graph-RePair bound calculation time.

These methods produce different guarantees: a fast bound does not prove the
minimum. Every sample starts a fresh process, and paired sample order alternates.
Algorithm timing excludes input parsing; process timing includes startup, parsing,
result serialization, teardown and subprocess overhead. Censored exact searches
are excluded from completed-case speedups, never substituted with the deadline.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "benchmarks" / "cases.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_cases(names: list[str]) -> list[dict]:
    with MANIFEST.open(encoding="utf-8", newline="") as stream:
        cases = list(csv.DictReader(stream, delimiter="\t"))
    unknown = set(names) - {case["name"] for case in cases}
    if unknown:
        raise ValueError(f"unknown cases: {', '.join(sorted(unknown))}")
    selected = [case for case in cases if not names or case["name"] in names]
    for case in selected:
        path = (MANIFEST.parent / case["input"]).resolve(strict=True)
        case["input"] = str(path.relative_to(ROOT))
        case["input_sha256"] = sha256(path)
        case["expected_assembly_index"] = int(case["expected_assembly_index"])
    return selected


def parse_events(output: str | bytes | None) -> dict:
    if isinstance(output, bytes):
        output = output.decode("utf-8")
    result = {}
    for line in (output or "").splitlines():
        event = json.loads(line)
        kind = event.pop("event")
        if kind not in {"started", "finished"} or kind in result:
            raise ValueError(f"unexpected probe event: {kind}")
        result[kind] = event
    return result


def run_sample(probe: Path, case: dict, method: str, timeout: float) -> dict:
    started = time.perf_counter()
    try:
        # The user-selected local benchmark binary is invoked without a shell.
        process = subprocess.run(  # noqa: S603
            [str(probe), method, str(ROOT / case["input"])],
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        elapsed = time.perf_counter() - started
        events = parse_events(process.stdout)
        if process.returncode != 0:
            raise RuntimeError(f"{case['name']} {method}: {process.stderr.strip()}")
        if set(events) != {"started", "finished"}:
            raise ValueError(f"{case['name']} {method}: incomplete probe output")
        sample = events["started"] | events["finished"]
        if sample["method"] != method:
            raise ValueError("probe returned an unexpected method")
        if not sample["succeeded"]:
            raise ValueError(f"{case['name']} {method}: calculation failed")
        for field in ("algorithm_seconds", "cpu_seconds"):
            if sample[field] < 0:
                raise ValueError(f"invalid {field}: {sample[field]}")
        sample["end_to_end_seconds"] = elapsed
        return sample
    except subprocess.TimeoutExpired as error:
        elapsed = time.perf_counter() - started
        events = parse_events(error.stdout)
        if "finished" in events:
            raise RuntimeError("probe timed out after returning a result") from error
        return {
            **events.get("started", {}),
            "method": method,
            "status": "timeout",
            "exact_completed": False,
            "algorithm_seconds": None,
            "cpu_seconds": None,
            "end_to_end_seconds": elapsed,
            "process_timeout_seconds": timeout,
            "calculation_started": "started" in events,
        }


def summarize_case(case: dict, samples: list[dict], runs: int) -> dict:
    """Validate deterministic answers and keep partial/censored search separate."""
    measured = [row for row in samples if row["phase"] == "measured"]
    completed_exact = [row for row in samples if row.get("exact_completed")]
    exact_values = {row["assembly_index"] for row in completed_exact}
    bound_values = {
        row["assembly_index"]
        for row in samples
        if row["method"] == "graph" and row["status"] == "upper_bound"
    }
    if len(exact_values) > 1 or len(bound_values) > 1:
        raise ValueError(f"{case['name']}: nondeterministic assembly indices")
    exact = next(iter(exact_values), None)
    bound = next(iter(bound_values), None)
    reference = case["expected_assembly_index"]
    if exact is not None and case["expectation"] == "reviewed" and exact != reference:
        raise ValueError(
            f"{case['name']}: exact result {exact} != reviewed {reference}"
        )
    for row in samples:
        if row["method"] == "graph" and row["status"] == "upper_bound":
            if not 0 <= row["assembly_index"] <= row["trivial_upper_bound"]:
                raise ValueError(f"{case['name']}: bound outside the trivial range")
            if exact is not None and row["assembly_index"] < exact:
                raise ValueError(f"{case['name']}: bound below completed exact result")

    stats = {}
    for method in ("exact", "graph"):
        rows = [row for row in measured if row["method"] == method]
        successful = [
            row for row in rows if row["status"] in {"completed", "upper_bound"}
        ]
        stats[method] = {
            "measured_samples": len(rows),
            "completed_samples": len(successful),
            "all_measured_completed": len(successful) == runs,
        }
        for field in ("algorithm_seconds", "cpu_seconds", "end_to_end_seconds"):
            values = [row[field] for row in successful]
            stats[method][f"median_{field}"] = (
                statistics.median(values) if values else None
            )

    comparisons = {
        f"paired_median_{field}_speedup": None
        for field in ("algorithm", "cpu", "end_to_end")
    }
    # A single timeout/limit in either warmup or measurement excludes the case
    # from completed comparisons even if some exact repetitions did finish.
    all_finished = all(row["status"] in {"completed", "upper_bound"} for row in samples)
    eligible = all_finished and all(
        stats[method]["all_measured_completed"] for method in stats
    )
    if eligible:
        pairs = {(row["round"], row["method"]): row for row in measured}
        for field in ("algorithm", "cpu", "end_to_end"):
            denominators = [
                pairs[(run, "graph")][f"{field}_seconds"] for run in range(runs)
            ]
            if all(value > 0 for value in denominators):
                ratios = [
                    pairs[(run, "exact")][f"{field}_seconds"] / denominators[run]
                    for run in range(runs)
                ]
                comparisons[f"paired_median_{field}_speedup"] = statistics.median(
                    ratios
                )
    deadlines = [
        row["process_timeout_seconds"]
        for row in samples
        if row["method"] == "exact" and row["status"] == "timeout"
    ]
    process_bound = None
    graph_process = stats["graph"]["median_end_to_end_seconds"]
    if deadlines and graph_process and stats["graph"]["all_measured_completed"]:
        process_bound = min(deadlines) / graph_process
    return {
        "exact_assembly_index": exact,
        "graph_upper_bound": bound,
        "gap_to_exact": None if exact is None or bound is None else bound - exact,
        "gap_to_manifest": None if bound is None else bound - reference,
        "exact_matches_manifest": None if exact is None else exact == reference,
        "completed_comparison": eligible,
        "methods": stats,
        "comparison": comparisons,
        "censored_process_speedup_lower_bound": process_bound,
    }


def benchmark_case(
    probe: Path, case: dict, case_number: int, runs: int, warmup: int, timeout: float
) -> dict:
    samples = []
    exact_active = True
    for phase, repetitions in (("warmup", warmup), ("measured", runs)):
        for round_number in range(repetitions):
            order = (
                ("exact", "graph")
                if (case_number + round_number) % 2 == 0
                else ("graph", "exact")
            )
            for position, method in enumerate(order):
                if method == "exact" and not exact_active:
                    continue
                sample = run_sample(probe, case, method, timeout)
                sample.update(phase=phase, round=round_number, position=position)
                samples.append(sample)
                if method == "exact" and sample["status"] != "completed":
                    exact_active = False
                if method == "graph" and sample["status"] != "upper_bound":
                    raise RuntimeError(
                        f"{case['name']}: graph bound {sample['status']}"
                    )
    return {**case, "samples": samples, **summarize_case(case, samples, runs)}


def aggregate(cases: list[dict]) -> dict:
    completed = [case for case in cases if case["completed_comparison"]]
    result = {
        "cases": len(cases),
        "completed_comparisons": len(completed),
        "excluded_cases": [
            case["name"] for case in cases if not case["completed_comparison"]
        ],
        "bounds_equal_to_completed_exact": sum(
            case["gap_to_exact"] == 0 for case in completed
        ),
        "maximum_gap_to_completed_exact": max(
            (case["gap_to_exact"] for case in completed), default=None
        ),
    }
    for field in ("algorithm", "cpu", "end_to_end"):
        exact = sum(
            case["methods"]["exact"][f"median_{field}_seconds"] for case in completed
        )
        graph = sum(
            case["methods"]["graph"][f"median_{field}_seconds"] for case in completed
        )
        ratios = [
            case["comparison"][f"paired_median_{field}_speedup"] for case in completed
        ]
        ratios = [value for value in ratios if value is not None]
        result[f"completed_sum_of_medians_{field}_speedup"] = (
            exact / graph if graph else None
        )
        result[f"completed_median_case_{field}_speedup"] = (
            statistics.median(ratios) if ratios else None
        )
        result[f"completed_exact_sum_of_medians_{field}_seconds"] = exact
        result[f"completed_graph_sum_of_medians_{field}_seconds"] = graph
    return result


def write_outputs(document: dict, json_path: Path, csv_path: Path | None) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    if csv_path is None:
        return
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"name": case["name"], "expectation": case["expectation"], **sample}
        for case in document["cases"]
        for sample in case["samples"]
    ]
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=10)
    parser.add_argument("--runs", type=int, default=6)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument(
        "--cpu", type=int, help="pin driver and all samples to this CPU"
    )
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--csv-output", type=Path)
    arguments = parser.parse_args()
    if arguments.runs < 1 or arguments.warmup < 0 or arguments.timeout <= 0:
        parser.error("runs and timeout must be positive; warmup must be nonnegative")
    try:
        probe = arguments.probe.resolve(strict=True)
        cases = load_cases(arguments.case)
        if arguments.cpu is not None:
            if not hasattr(os, "sched_setaffinity"):
                raise ValueError("CPU affinity is unsupported on this platform")
            if arguments.cpu not in os.sched_getaffinity(0):
                raise ValueError(f"CPU {arguments.cpu} is outside the allowed affinity")
            os.sched_setaffinity(0, {arguments.cpu})
        source_paths = [
            *sorted((ROOT / "src").glob("*")),
            Path(__file__),
            ROOT / "benchmarks" / "graphRepairSpeedProbe.cpp",
        ]
        sources = {
            str(path.relative_to(ROOT)): sha256(path)
            for path in source_paths
            if path.is_file()
        }
        git = subprocess.run(
            ["git", "rev-parse", "HEAD"],  # noqa: S607
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        document = {
            "schema_version": 1,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "runs": arguments.runs,
            "warmup": arguments.warmup,
            "timeout_seconds": arguments.timeout,
            "schedule": (
                "adjacent paired samples alternating AB/BA; "
                "fresh process for every sample"
            ),
            "scope": (
                "exact minimum versus heuristic upper bound; "
                "not exact-search acceleration"
            ),
            "pathway_output": False,
            "remove_hydrogens": True,
            "compensate_disjoint": False,
            "probe": {"path": str(probe), "sha256": sha256(probe)},
            "manifest": {
                "path": str(MANIFEST.relative_to(ROOT)),
                "sha256": sha256(MANIFEST),
            },
            "sources_sha256": sources,
            "git_revision": git.stdout.strip(),
            "machine": {
                "platform": platform.platform(),
                "python": sys.version,
                "cpu_affinity": sorted(os.sched_getaffinity(0))
                if hasattr(os, "sched_getaffinity")
                else None,
            },
            "cases": [],
        }
        for index, case in enumerate(cases):
            result = benchmark_case(
                probe, case, index, arguments.runs, arguments.warmup, arguments.timeout
            )
            document["cases"].append(result)
            document["summary"] = aggregate(document["cases"])
            write_outputs(document, arguments.json_output, arguments.csv_output)
            label = (
                "completed" if result["completed_comparison"] else "censored/incomplete"
            )
            print(
                f"{index + 1}/{len(cases)} {case['name']}: {label}; "
                f"bound={result['graph_upper_bound']}, "
                f"exact={result['exact_assembly_index']}",
                flush=True,
            )
        document["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write_outputs(document, arguments.json_output, arguments.csv_output)
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"benchmark error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
