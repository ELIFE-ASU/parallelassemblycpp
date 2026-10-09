"""Compare the graph-pair construction bound with the initial solver incumbent.

The probe computes bounds only. Reviewed manifest values are regression
references; this runner does not rerun or establish exact search results.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

if __package__:
    from . import benchmark
else:
    import benchmark

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent


def load_cases(corpus: str) -> list[dict]:
    """Load and deduplicate the two maintained manifests by resolved input."""
    cases: dict[Path, dict] = {}
    if corpus in {"regression", "all"}:
        manifest = REPOSITORY_ROOT / "unitTests" / "regression_cases.tsv"
        with manifest.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                source = manifest.parent / row["molecule"]
                if not source.is_file():
                    source = source.with_suffix(".mol")
                source = source.resolve(strict=True)
                cases[source] = {
                    "name": row["molecule"],
                    "input": str(source.relative_to(REPOSITORY_ROOT)),
                    "expected_assembly_index": int(row["expected_assembly_index"]),
                    "expectation": "reviewed",
                    "suites": ["regression"],
                }
    if corpus in {"benchmarks", "all"}:
        manifest = REPOSITORY_ROOT / "benchmarks" / "cases.tsv"
        _, benchmark_cases = benchmark.load_manifest(manifest)
        for case in benchmark_cases:
            source = case.source
            expected = case.expected_assembly_index
            if source in cases:
                if cases[source]["expected_assembly_index"] != expected:
                    raise ValueError(f"conflicting expected indices for {source}")
                cases[source]["suites"].extend(case.suites)
                continue
            cases[source] = {
                "name": case.name,
                "input": str(source.relative_to(REPOSITORY_ROOT)),
                "expected_assembly_index": expected,
                "expectation": case.expectation,
                "suites": list(case.suites),
            }
    return list(cases.values())


def summarize(rows: list[dict]) -> dict:
    """Keep quality comparisons separate from algorithm timing."""
    if not rows:
        return {"cases": 0}
    return {
        "cases": len(rows),
        "strictly_improved": sum(row["improvement"] > 0 for row in rows),
        "unchanged": sum(row["improvement"] == 0 for row in rows),
        "worse": sum(row["improvement"] < 0 for row in rows),
        "equal_to_reference": sum(row["gap_to_reference"] == 0 for row in rows),
        "below_reference": sum(row["gap_to_reference"] < 0 for row in rows),
        "mean_initial_bound": statistics.mean(
            row["trivial_upper_bound"] for row in rows
        ),
        "mean_graph_repair_bound": statistics.mean(
            row["graph_repair_upper_bound"] for row in rows
        ),
        "mean_improvement": statistics.mean(row["improvement"] for row in rows),
        "mean_gap_to_reference": statistics.mean(
            row["gap_to_reference"] for row in rows
        ),
        "maximum_gap_to_reference": max(row["gap_to_reference"] for row in rows),
        "graph_better_than_trail": sum(
            row["graph_repair_upper_bound"] < row["trail_repair_upper_bound"]
            for row in rows
        ),
        "graph_equal_to_trail": sum(
            row["graph_repair_upper_bound"] == row["trail_repair_upper_bound"]
            for row in rows
        ),
        "graph_worse_than_trail": sum(
            row["graph_repair_upper_bound"] > row["trail_repair_upper_bound"]
            for row in rows
        ),
        "mean_trail_repair_bound": statistics.mean(
            row["trail_repair_upper_bound"] for row in rows
        ),
        "mean_trail_gap_to_reference": statistics.mean(
            row["trail_repair_upper_bound"] - row["expected_assembly_index"]
            for row in rows
        ),
        "median_algorithm_seconds": statistics.median(
            row["median_algorithm_seconds"] for row in rows
        ),
        "total_algorithm_seconds": sum(row["median_algorithm_seconds"] for row in rows),
    }


def trail_repair(atoms: list[str], edges: list[list[int]]) -> dict:
    """String RePair over one deterministic partition into simple paths.

    Paths cannot revisit a vertex: repeated labelled strings then always mean
    isomorphic graph fragments, including when separate occurrences share
    vertices. General walks would lose vertex identifications and are unsafe.
    The whole path orientation is canonicalized before compression; digrams
    themselves use ordinary directed matching without reversal matching.
    """
    adjacency: list[list[tuple[int, int]]] = [[] for _ in atoms]
    for edge_id, (first, second, _) in enumerate(edges):
        adjacency[first].append((edge_id, second))
        adjacency[second].append((edge_id, first))
    unused = set(range(len(edges)))
    paths: list[tuple[list[tuple[str, int, str]], list[int]]] = []
    while unused:
        first_edge = min(unused)
        first, second, bond = edges[first_edge]
        unused.remove(first_edge)
        vertices = [first, second]
        bonds = [bond]
        edge_ids = [first_edge]
        visited = {first, second}
        # Extend right, then left. An edge closing a cycle starts another path.
        for right in (True, False):
            while True:
                vertex = vertices[-1] if right else vertices[0]
                options = [
                    (edge_id, neighbor)
                    for edge_id, neighbor in adjacency[vertex]
                    if edge_id in unused and neighbor not in visited
                ]
                if not options:
                    break
                edge_id, neighbor = min(options)
                unused.remove(edge_id)
                visited.add(neighbor)
                if right:
                    vertices.append(neighbor)
                    bonds.append(edges[edge_id][2])
                    edge_ids.append(edge_id)
                else:
                    vertices.insert(0, neighbor)
                    bonds.insert(0, edges[edge_id][2])
                    edge_ids.insert(0, edge_id)
        forward = [
            (atoms[vertices[index]], bond, atoms[vertices[index + 1]])
            for index, bond in enumerate(bonds)
        ]
        reverse = [(end, bond, start) for start, bond, end in reversed(forward)]
        if reverse < forward:
            paths.append((reverse, list(reversed(edge_ids))))
        else:
            paths.append((forward, edge_ids))

    symbols = {
        token: index
        for index, token in enumerate(
            sorted({token for path, _ in paths for token in path})
        )
    }
    sequences = [
        [(symbols[token], [edge_id]) for token, edge_id in zip(path, ids, strict=True)]
        for path, ids in paths
    ]
    terminals = {}
    for sequence in sequences:
        for symbol, edge_ids in sequence:
            terminals.setdefault(symbol, {"id": symbol, "edges": edge_ids})
    next_symbol = len(symbols)
    rules = []
    while True:
        frequencies: Counter[tuple[int, int]] = Counter()
        for sequence in sequences:
            last_end: dict[tuple[int, int], int] = {}
            for position in range(len(sequence) - 1):
                pair = sequence[position][0], sequence[position + 1][0]
                if last_end.get(pair, -1) < position:
                    frequencies[pair] += 1
                    last_end[pair] = position + 1
        candidates = [pair for pair, count in frequencies.items() if count > 1]
        if not candidates:
            break
        selected = min(candidates, key=lambda pair: (-frequencies[pair], pair))
        for index, sequence in enumerate(sequences):
            replaced = []
            position = 0
            while position < len(sequence):
                if (
                    position + 1 < len(sequence)
                    and (sequence[position][0], sequence[position + 1][0]) == selected
                ):
                    left, right = sequence[position], sequence[position + 1]
                    joined = sorted(left[1] + right[1])
                    if not rules or rules[-1]["id"] != next_symbol:
                        rules.append(
                            {
                                "id": next_symbol,
                                "left": left[0],
                                "right": right[0],
                                "left_edges": left[1],
                                "right_edges": right[1],
                                "edges": joined,
                            }
                        )
                    replaced.append((next_symbol, joined))
                    position += 2
                else:
                    replaced.append(sequence[position])
                    position += 1
            sequences[index] = replaced
        next_symbol += 1
    residual = sum(len(sequence) for sequence in sequences)
    remaining_vertices = {vertex for edge in edges for vertex in edge[:2]}
    components = 0
    while remaining_vertices:
        pending = [remaining_vertices.pop()]
        components += 1
        while pending:
            vertex = pending.pop()
            for _, neighbor in adjacency[vertex]:
                if neighbor in remaining_vertices:
                    remaining_vertices.remove(neighbor)
                    pending.append(neighbor)
    upper_bound = len(rules) + max(residual - 1, 0)
    return {
        "trail_repair_upper_bound": upper_bound,
        "trail_rules": len(rules),
        "trail_residual_tokens": residual,
        "trail_count": len(paths),
        "certificate": {
            "schema": "graph-repair-assembly-v1",
            "upper_bound": upper_bound,
            "trivial_upper_bound": max(len(edges) - 1, 0),
            "rule_count": len(rules),
            "remaining_fragments": residual,
            "components": components,
            "compensate_disjoint": False,
            "atoms": atoms,
            "edges": edges,
            "terminals": list(terminals.values()),
            "rules": rules,
            "residual": [
                {"symbol": symbol, "edges": edge_ids}
                for sequence in sequences
                for symbol, edge_ids in sequence
            ],
        },
    }


def run_probe(
    executable: Path, cases: list[dict], repeats: int, timeout: float
) -> list[dict]:
    """Run the corpus once per process and validate every repetition's certificate.

    The probe's per-case timer excludes parsing, startup, and validation.
    """
    checker_spec = importlib.util.spec_from_file_location(
        "graph_repair_validation",
        REPOSITORY_ROOT / "unitTests" / "graphRepairTester.py",
    )
    if checker_spec is None or checker_spec.loader is None:
        raise ImportError("graph-repair certificate checker could not be loaded")
    checker = importlib.util.module_from_spec(checker_spec)
    checker_spec.loader.exec_module(checker)
    all_samples = []
    for _ in range(repeats):
        completed = subprocess.run(  # noqa: S603 -- explicit benchmark probe
            [str(executable), "--files"],
            input="".join(
                str(REPOSITORY_ROOT / case["input"]) + "\n" for case in cases
            ),
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        samples = [json.loads(line) for line in completed.stdout.splitlines()]
        if len(samples) != len(cases):
            raise ValueError("probe output count does not match selected cases")
        all_samples.append(samples)
    rows = []
    for index, case in enumerate(cases):
        sample = all_samples[0][index]
        if "error" in sample:
            raise ValueError(f"probe failed for {case['name']}: {sample['error']}")
        bound = sample["upper_bound"]
        trivial = sample["trivial_upper_bound"]
        timings = []
        for samples in all_samples:
            repeated = samples[index]
            if "error" in repeated:
                raise ValueError(
                    f"probe failed for {case['name']}: {repeated['error']}"
                )
            try:
                validated_bound = checker.validate_certificate(repeated)
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f"invalid construction for {case['name']}: {error}"
                ) from error
            if validated_bound != repeated["upper_bound"]:
                raise ValueError(f"construction count differs for {case['name']}")
            if repeated["upper_bound"] != bound:
                raise ValueError(f"nondeterministic bound for {case['name']}")
            seconds = repeated.get("elapsed_seconds")
            if (
                type(seconds) not in (int, float)
                or not math.isfinite(seconds)
                or seconds < 0
            ):
                raise ValueError(f"invalid algorithm timing for {case['name']}")
            timings.append(seconds)
        started_at = time.perf_counter()
        trail = trail_repair(sample["atoms"], sample["edges"])
        trail_seconds = time.perf_counter() - started_at
        if (
            checker.validate_certificate(trail.pop("certificate"))
            != trail["trail_repair_upper_bound"]
        ):
            raise ValueError(f"trail construction count differs for {case['name']}")
        if bound > trivial:
            raise ValueError(f"bound exceeds the initial bound for {case['name']}")
        if (
            case["expectation"] == "reviewed"
            and bound < case["expected_assembly_index"]
        ):
            raise ValueError(f"bound is below reviewed index for {case['name']}")
        if trail["trail_repair_upper_bound"] > trivial:
            raise ValueError(f"trail bound exceeds the initial bound: {case['name']}")
        if (
            case["expectation"] == "reviewed"
            and trail["trail_repair_upper_bound"] < case["expected_assembly_index"]
        ):
            raise ValueError(f"trail bound is below reviewed index: {case['name']}")
        rows.append(
            {
                **case,
                **trail,
                "bonds": len(sample["edges"]),
                "trivial_upper_bound": trivial,
                "graph_repair_upper_bound": bound,
                "graph_repair_rules": sample["rule_count"],
                "graph_repair_remaining_fragments": sample["remaining_fragments"],
                "construction_valid": True,
                "trail_construction_valid": True,
                "algorithm_seconds": timings,
                "trail_python_seconds": trail_seconds,
                "input_sha256": benchmark.file_sha256(REPOSITORY_ROOT / case["input"]),
                "improvement": trivial - bound,
                "gap_to_reference": bound - case["expected_assembly_index"],
                "median_algorithm_seconds": statistics.median(timings),
            }
        )
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument(
        "--corpus", choices=("regression", "benchmarks", "all"), default="all"
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="select a case name; repeatable",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument(
        "--json-output",
        type=Path,
        default=REPOSITORY_ROOT / "build" / "graph-repair-bounds.json",
    )
    parser.add_argument("--csv-output", type=Path)
    args = parser.parse_args(argv)
    if args.repeats < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error(
            "--repeats must be positive; --timeout must be finite and positive"
        )
    cases = load_cases(args.corpus)
    if args.case:
        selected = set(args.case)
        unknown = selected - {case["name"] for case in cases}
        if unknown:
            parser.error(f"unknown case names: {', '.join(sorted(unknown))}")
        cases = [case for case in cases if case["name"] in selected]
    executable = args.executable.resolve(strict=True)
    rows = run_probe(executable, cases, args.repeats, args.timeout)
    summaries = {
        "all": summarize(rows),
        "reviewed": summarize(
            [row for row in rows if row["expectation"] == "reviewed"]
        ),
        "provisional": summarize(
            [row for row in rows if row["expectation"] == "provisional"]
        ),
    }
    for suite in sorted({suite for row in rows for suite in row["suites"]}):
        summaries[suite] = summarize([row for row in rows if suite in row["suites"]])
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "executable": str(executable),
        "executable_sha256": benchmark.file_sha256(executable),
        "source_sha256": {
            str(path): benchmark.file_sha256(REPOSITORY_ROOT / path)
            for path in (
                Path("src/graphRepair.h"),
                Path("unitTests/graphRepairProbe.cpp"),
                Path("unitTests/graphRepairTester.py"),
                Path("benchmarks/graph_repair_benchmark.py"),
            )
        },
        "repeats": args.repeats,
        "baseline": (
            "Original hydrogen-filtered bond count minus one; "
            "the solver's initial incumbent."
        ),
        "trail_baseline": (
            "Deterministic edge-disjoint simple paths, lowest-index edge first, "
            "extended right then left. Canonical whole-path orientation. "
            "Ordinary most-frequent nonoverlapping directed digram RePair, "
            "without reversal matching. Dictionary joins plus residual tokens "
            "minus one. This does not reproduce an optimized trail partition."
        ),
        "reference_caveat": (
            "Reviewed regression indices and provisional benchmark guards, "
            "not newly proved optima."
        ),
        "timing_scope": (
            "Graph-pair algorithm only; parsing and certificate validation "
            "excluded. Python trail "
            "baseline time is separate and not directly comparable to C++."
        ),
        "summaries": summaries,
        "cases": rows,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.csv_output:
        args.csv_output.parent.mkdir(parents=True, exist_ok=True)
        columns = [
            "name",
            "input",
            "expectation",
            "bonds",
            "expected_assembly_index",
            "trivial_upper_bound",
            "trail_repair_upper_bound",
            "graph_repair_upper_bound",
            "improvement",
            "gap_to_reference",
            "median_algorithm_seconds",
            "construction_valid",
        ]
        with args.csv_output.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, columns, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps(summaries, indent=2))
    print(f"Report: {args.json_output}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (
        OSError,
        ValueError,
        benchmark.BenchmarkError,
        subprocess.SubprocessError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)
