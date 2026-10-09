"""Paired baseline/scalar/vector exact string-search timing, using only stdlib.

Each sample starts a fresh stringRepairSpeedProbe process. Algorithm and CPU
times come from the probe; process time also includes startup and serialization.
All six variant orders rotate across rounds. A timed-out variant is skipped for
the remaining rounds of that case, and no censored value becomes a speedup.
Completed assembly indices must agree across every variant and repetition.
Witness replay is covered by the string unit tests, not this timing probe.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

if __package__:
    from . import benchmark
else:
    import benchmark

VARIANTS = ("baseline", "scalar", "vector")
ORDERS = tuple(itertools.permutations(VARIANTS))
METRICS = ("algorithm_seconds", "cpu_seconds", "process_seconds")


def fixtures() -> list[dict]:
    cases = [
        ("unary-64", "a" * 64, False, 6),
        ("unary-256", "a" * 256, False, 8),
        ("vector-6-2", "aaabaaab", False, 4),
        ("vector-12-4", "aaabaaab" * 2, False, 5),
        ("vector-four-symbols", "abcd" * 2, False, 4),
        ("vector-five-symbols", "abcde" * 2, False, 5),
        ("vector-unicode", "ééé😀ééé😀", False, 4),
        ("search-improves-seed", "bbbabba", False, 4),
        ("reversed-search-improves-seed", "baabaa", True, 3),
        ("loose-composition", "ababcdcd", False, 5),
        ("loose-unbalanced", "aaaaabb", False, None),
        ("repeated-block-32", "abcddcba" * 4, False, None),
        ("homogeneous-blocks-48", "a" * 16 + "b" * 16 + "c" * 16, False, None),
    ]
    rng = random.Random(20261006)  # noqa: S311 - reproducible benchmark inputs.
    for alphabet, length in ((2, 16), (3, 24), (4, 32), (6, 48), (10, 80)):
        value = "".join(chr(97 + rng.randrange(alphabet)) for _ in range(length))
        cases.append((f"random-{alphabet}-{length}", value, False, None))
    return [
        {
            "name": name,
            "input": value,
            "length": len(value),
            "composition": sorted(Counter(value).values()),
            "input_sha256": hashlib.sha256(value.encode()).hexdigest(),
            "accept_reversed": reversed_,
            "expected_assembly_index": expected,
        }
        for name, value, reversed_, expected in cases
    ]


def decode_text(value: str | bytes | None) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else value or ""


def parse_events(output: str | bytes | None) -> dict:
    events = {}
    for line in decode_text(output).splitlines():
        event = json.loads(line)
        kind = event.get("event")
        if kind not in {"started", "finished"} or kind in events:
            raise ValueError(f"unexpected probe event: {kind}")
        events[kind] = event
    return events


def run_sample(probe: Path, case: dict, path: Path, timeout: float) -> dict:
    command = [str(probe), "exact", str(path)]
    if case["accept_reversed"]:
        command.append("--accept-reversed")
    started = time.perf_counter()
    try:
        process = subprocess.run(  # noqa: S603 - explicitly selected local probe.
            command, capture_output=True, text=True, check=False, timeout=timeout
        )
    except subprocess.TimeoutExpired as error:
        events = parse_events(error.stdout)
        if "finished" in events:
            raise RuntimeError("probe timed out after reporting completion") from error
        return {
            "status": "timeout",
            "events": events,
            "stderr": decode_text(error.stderr),
            "exact_completed": False,
            "timeout_seconds": timeout,
            "process_seconds": time.perf_counter() - started,
        }
    elapsed = time.perf_counter() - started
    if process.returncode != 0:
        raise RuntimeError(f"{case['name']}: probe failed: {process.stderr.strip()}")
    events = parse_events(process.stdout)
    if set(events) != {"started", "finished"}:
        raise ValueError("successful probe omitted its start or finish event")
    initial, final = events["started"], events["finished"]
    if (
        initial.get("method") != "exact"
        or initial.get("length") != case["length"]
        or initial.get("accept_reversed") != case["accept_reversed"]
        or final.get("method") != "exact"
        or final.get("exact_completed") is not True
        or final.get("upper_bound_only") is not False
    ):
        raise ValueError("probe returned mismatched input or incomplete exact result")
    index = final.get("assembly_index")
    if type(index) is not int or not -1 <= index <= case["length"] - 1:
        raise ValueError(f"invalid exact assembly index: {index}")
    for metric in METRICS[:2]:
        value = final.get(metric)
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid {metric}: {value}")
    return {
        **final,
        "status": "completed",
        "events": events,
        "stderr": process.stderr,
        "process_seconds": elapsed,
    }


def summarize(case: dict, samples: list[dict], runs: int) -> dict:
    indices = {row["assembly_index"] for row in samples if row["exact_completed"]}
    if len(indices) > 1:
        raise ValueError(f"{case['name']}: baseline/scalar/vector exact indices differ")
    index = next(iter(indices), None)
    expected = case["expected_assembly_index"]
    if expected is not None and index is not None and index != expected:
        raise ValueError(f"{case['name']}: exact index {index} != expected {expected}")
    measured = [row for row in samples if row["phase"] == "measured"]
    variants = {}
    for variant in VARIANTS:
        rows = [row for row in measured if row["variant"] == variant]
        complete = [row for row in rows if row["exact_completed"]]
        variants[variant] = {
            "completed_samples": len(complete),
            "censored": any(
                row["status"] == "timeout"
                for row in samples
                if row["variant"] == variant
            ),
            **{
                f"median_{metric}": statistics.median(row[metric] for row in complete)
                if complete
                else None
                for metric in METRICS
            },
        }
    paired = {(row["round"], row["variant"]): row for row in measured}
    comparisons = {}
    for numerator, denominator in (
        ("baseline", "scalar"),
        ("baseline", "vector"),
        ("scalar", "vector"),
    ):
        eligible = all(
            variants[variant]["completed_samples"] == runs
            and not variants[variant]["censored"]
            for variant in (numerator, denominator)
        )
        ratios = {}
        for metric in METRICS:
            values = (
                [
                    paired[(run, numerator)][metric]
                    / paired[(run, denominator)][metric]
                    for run in range(runs)
                ]
                if eligible
                and all(paired[(run, denominator)][metric] > 0 for run in range(runs))
                else []
            )
            ratios[metric] = statistics.median(values) if values else None
        comparisons[f"{numerator}_over_{denominator}"] = ratios
    return {
        "assembly_index": index,
        "completed_indices_agree": bool(indices),
        "all_variants_completed": all(
            row["completed_samples"] == runs and not row["censored"]
            for row in variants.values()
        ),
        "variants": variants,
        "paired_median_speedups": comparisons,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for variant in VARIANTS:
        parser.add_argument(f"--{variant}", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--runs", type=int, default=6)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=3.0)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--case", action="append", default=[])
    args = parser.parse_args()
    if (
        args.runs < 1
        or args.warmup < 0
        or not math.isfinite(args.timeout)
        or args.timeout <= 0
    ):
        parser.error(
            "runs must be positive; warmup nonnegative; timeout finite and positive"
        )
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    probes = {name: getattr(args, name).resolve(strict=True) for name in VARIANTS}
    cases = fixtures()
    unknown = set(args.case) - {case["name"] for case in cases}
    if unknown:
        parser.error(f"unknown cases: {', '.join(sorted(unknown))}")
    cases = [case for case in cases if not args.case or case["name"] in args.case]
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs = output.parent / f"{output.stem}-inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    report = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python": sys.version,
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "runs": args.runs,
        "warmup": args.warmup,
        "timeout_seconds": args.timeout,
        "witness_validation": (
            "Not exposed by timing probe; covered by string unit tests."
        ),
        "probes": {
            name: {
                "path": str(path),
                "sha256": benchmark.file_sha256(path),
            }
            for name, path in probes.items()
        },
        "cases": [],
    }
    for case_number, case in enumerate(cases):
        path = inputs / f"{case['name']}.txt"
        path.write_text(case["input"], encoding="utf-8")
        samples = []
        active = set(VARIANTS)
        for phase, repetitions in (("warmup", args.warmup), ("measured", args.runs)):
            for run in range(repetitions):
                order = ORDERS[(case_number + run) % len(ORDERS)]
                for position, variant in enumerate(order):
                    if variant not in active:
                        continue
                    sample = run_sample(probes[variant], case, path, args.timeout)
                    sample.update(
                        variant=variant, phase=phase, round=run, position=position
                    )
                    samples.append(sample)
                    if sample["status"] == "timeout":
                        active.remove(variant)
        summary = summarize(case, samples, args.runs)
        report["cases"].append({**case, **summary, "samples": samples})
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        timings = []
        for variant, data in summary["variants"].items():
            seconds = data["median_algorithm_seconds"]
            timings.append(
                f"{variant}={seconds:.6g}s"
                if seconds is not None
                else f"{variant}=timeout"
            )
        print(
            f"{case['name']}: AI={summary['assembly_index']} " + " ".join(timings),
            flush=True,
        )
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
