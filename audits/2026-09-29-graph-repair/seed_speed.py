"""Paired exact-search timing before and after the constructive Re-Pair seed.

The seeded implementation was removed; the default candidate is its preserved
local experimental binary. Every sample is a fresh process on one CPU.
Finish the current pair after an incomplete result, then skip the case;
incomplete cases never enter speed ratios.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.graph_repair_speed import load_cases, run_sample, sha256  # noqa: E402

HERE = Path(__file__).resolve().parent
VARIANTS = ("baseline", "seeded")
FIELDS = ("algorithm", "cpu", "end_to_end")


def completed(sample: dict) -> bool:
    return sample["status"] == "completed" and sample.get("exact_completed", False)


def summarize(case: dict, samples: list[dict], runs: int, warmup: int) -> dict:
    values = {row["assembly_index"] for row in samples if completed(row)}
    if len(values) > 1:
        raise ValueError(f"{case['name']}: exact answers disagree")
    exact = next(iter(values), None)
    if (
        exact is not None
        and case["expectation"] == "reviewed"
        and exact != case["expected_assembly_index"]
    ):
        raise ValueError(f"{case['name']}: exact answer != reviewed reference")
    measured = [row for row in samples if row["phase"] == "measured"]
    methods = {}
    for variant in VARIANTS:
        rows = [row for row in measured if row["variant"] == variant]
        successful = [row for row in rows if completed(row)]
        methods[variant] = {
            "measured_samples": len(rows),
            "completed_samples": len(successful),
        }
        for field in FIELDS:
            durations = [row[f"{field}_seconds"] for row in successful]
            methods[variant][f"median_{field}_seconds"] = (
                statistics.median(durations) if durations else None
            )
    eligible = (
        len(samples) == 2 * (runs + warmup)
        and all(completed(row) for row in samples)
        and all(methods[v]["completed_samples"] == runs for v in VARIANTS)
    )
    pairs = {(row["round"], row["variant"]): row for row in measured}
    comparison = {}
    for field in FIELDS:
        positive = eligible and all(
            pairs[(run, "seeded")][f"{field}_seconds"] > 0 for run in range(runs)
        )
        ratios = (
            [
                pairs[(run, "baseline")][f"{field}_seconds"]
                / pairs[(run, "seeded")][f"{field}_seconds"]
                for run in range(runs)
            ]
            if positive
            else []
        )
        comparison[f"paired_median_{field}_speedup"] = (
            statistics.median(ratios) if ratios else None
        )
    return {
        "exact_assembly_index": exact,
        "completed_comparison": eligible,
        "methods": methods,
        "comparison": comparison,
    }


def aggregate(cases: list[dict]) -> dict:
    eligible = [case for case in cases if case["completed_comparison"]]
    result = {
        "cases": len(cases),
        "completed_comparisons": len(eligible),
        "excluded_cases": [c["name"] for c in cases if not c["completed_comparison"]],
    }
    for field in FIELDS:
        totals = {
            v: sum(c["methods"][v][f"median_{field}_seconds"] for c in eligible)
            for v in VARIANTS
        }
        ratios = [c["comparison"][f"paired_median_{field}_speedup"] for c in eligible]
        ratios = [ratio for ratio in ratios if ratio is not None]
        result[f"completed_median_case_{field}_speedup"] = (
            statistics.median(ratios) if ratios else None
        )
        result[f"completed_sum_of_medians_{field}_speedup"] = (
            totals["baseline"] / totals["seeded"] if totals["seeded"] else None
        )
        for variant in VARIANTS:
            key = f"completed_{variant}_sum_of_medians_{field}_seconds"
            result[key] = totals[variant]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    probe_name = "parallelassemblycpp_graph_repair_speed_probe"
    parser.add_argument(
        "--baseline",
        type=Path,
        default=ROOT / "build/graph-repair-before-seed" / probe_name,
    )
    parser.add_argument(
        "--seeded",
        type=Path,
        default=ROOT / "build/graph-repair-seeded-experiment" / probe_name,
    )
    parser.add_argument("--runs", type=int, default=6)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=10)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument(
        "--json-output", type=Path, default=HERE / "seed_speed_repeat.json"
    )
    args = parser.parse_args()
    if args.runs < 1 or args.warmup < 0 or args.timeout <= 0:
        parser.error("runs/timeout must be positive and warmup nonnegative")
    probes = {v: getattr(args, v).resolve(strict=True) for v in VARIANTS}
    cases = load_cases(args.case)
    if not hasattr(os, "sched_setaffinity") or args.cpu not in os.sched_getaffinity(0):
        parser.error("requested CPU affinity is unavailable")
    os.sched_setaffinity(0, {args.cpu})
    previous = json.loads((HERE / "speed.json").read_text(encoding="utf-8"))
    probe_metadata = {
        v: {"path": str(p), "sha256": sha256(p)} for v, p in probes.items()
    }
    baseline_matches = (
        probe_metadata["baseline"]["sha256"] == previous["probe"]["sha256"]
    )
    historical_seed = json.loads((HERE / "seed_speed.json").read_text(encoding="utf-8"))
    seeded_matches = (
        probe_metadata["seeded"]["sha256"]
        == historical_seed["probes"]["seeded"]["sha256"]
    )
    document = {
        "schema_version": 1,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "exact search before versus exact search seeded by graph Re-Pair",
        "runs": args.runs,
        "warmup": args.warmup,
        "timeout_seconds": args.timeout,
        "schedule": (
            "fresh processes; adjacent AB/BA pairs; "
            "stop case after first incomplete pair"
        ),
        "pathway_output": False,
        "remove_hydrogens": True,
        "compensate_disjoint": False,
        "probes": probe_metadata,
        "baseline_source_provenance": {
            "report": str((HERE / "speed.json").relative_to(ROOT)),
            "report_sha256": sha256(HERE / "speed.json"),
            "probe_matches_previous_report": baseline_matches,
            "sources_sha256": previous["sources_sha256"] if baseline_matches else None,
        },
        "seeded_source_provenance": {
            "report": str((HERE / "seed_speed.json").relative_to(ROOT)),
            "report_sha256": sha256(HERE / "seed_speed.json"),
            "probe_matches_previous_report": seeded_matches,
        },
        "seeded_sources_sha256": (
            historical_seed["seeded_sources_sha256"] if seeded_matches else None
        ),
        "manifest_sha256": sha256(ROOT / "benchmarks/cases.tsv"),
        "machine": {
            "platform": platform.platform(),
            "python": sys.version,
            "cpu_affinity": sorted(os.sched_getaffinity(0)),
        },
        "cases": [],
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    for case_number, case in enumerate(cases):
        samples = []
        active = True
        for phase, repetitions in (("warmup", args.warmup), ("measured", args.runs)):
            for round_number in range(repetitions):
                if not active:
                    break
                order = (
                    VARIANTS
                    if (case_number + round_number) % 2 == 0
                    else VARIANTS[::-1]
                )
                for position, variant in enumerate(order):
                    sample = run_sample(probes[variant], case, "exact", args.timeout)
                    sample.update(
                        variant=variant,
                        phase=phase,
                        round=round_number,
                        position=position,
                    )
                    samples.append(sample)
                    active = active and completed(sample)
        result = {
            **case,
            "samples": samples,
            **summarize(case, samples, args.runs, args.warmup),
        }
        document["cases"].append(result)
        document["summary"] = aggregate(document["cases"])
        document["updated_utc"] = datetime.now(timezone.utc).isoformat()
        args.json_output.write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8"
        )
        speedup = result["comparison"]["paired_median_algorithm_speedup"]
        label = f"{speedup:.3f}x calculation" if speedup else "censored/incomplete"
        print(
            f"{case_number + 1}/{len(cases)} {case['name']}: {label}; "
            f"exact={result['exact_assembly_index']}",
            flush=True,
        )
    document["finished_utc"] = datetime.now(timezone.utc).isoformat()
    args.json_output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
