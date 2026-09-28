"""Standalone, deterministic molecular B&B audit; no solver internals imported.

Enumerates every connected edge subset and every legal duplicate pair, with
independent brute-force labeled graph isomorphism. The reference search uses
physical-mask memoization only: no production bounds, DAG or canonical ordinal.
Also checks the proposed edge-type vector bound at every visited state.
"""

import argparse
import functools
import hashlib
import itertools
import json
import random
import selectors
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import TypedDict


class BoundExample(TypedDict):
    fragment_sizes: list[int]
    distinct_edge_types: int
    maximum_duplicate_size: int
    vector_bound: int
    scalar_bound: int
    exact_remaining_saving: int


class Metrics(TypedDict):
    oracle_states: int
    vector_bound_states_checked: int
    vector_strictly_tighter_than_scalar: int
    largest_vector_improvement: int
    vector_examples: list[BoundExample]


Canonical = tuple[tuple[tuple[tuple[str, int], int], ...], tuple[int, ...] | None]


def oracle(
    n: int,
    edges: list[tuple[int, int]],
    labels: list[str],
    bonds: list[int],
    metrics: Metrics,
) -> int:
    m = len(edges)
    incident = [0] * n
    edge_types = []
    for i, (a, b) in enumerate(edges):
        incident[a] |= 1 << i
        incident[b] |= 1 << i
        edge_types.append((bonds[i], *sorted((labels[a], labels[b]))))

    @functools.cache
    def components(mask: int) -> tuple[int, ...]:
        remaining, output = mask, []
        while remaining:
            todo, component = remaining & -remaining, 0
            while todo:
                bit = todo & -todo
                todo -= bit
                if bit & component:
                    continue
                component |= bit
                a, b = edges[bit.bit_length() - 1]
                todo |= (incident[a] | incident[b]) & mask & ~component
            remaining &= ~component
            if component.bit_count() >= 2:
                output.append(component)
        return tuple(output)

    @functools.cache
    def canonical(mask: int) -> Canonical:
        degrees, adjacency = {}, {}
        for i, (a, b) in enumerate(edges):
            if mask >> i & 1:
                degrees[a] = degrees.get(a, 0) + 1
                degrees[b] = degrees.get(b, 0) + 1
                adjacency[a, b] = adjacency[b, a] = bonds[i]
        groups = {}
        for vertex, degree in degrees.items():
            groups.setdefault((labels[vertex], degree), []).append(vertex)
        partitions = [groups[key] for key in sorted(groups)]
        prefix = tuple((key, len(groups[key])) for key in sorted(groups))
        best = None
        for permutations in itertools.product(
            *(itertools.permutations(group) for group in partitions)
        ):
            vertices = sum(permutations, ())
            code = tuple(
                adjacency.get((a, b), 0)
                for j, a in enumerate(vertices)
                for b in vertices[j + 1 :]
            )
            if best is None or code < best:
                best = code
        return prefix, best

    connected = []
    canonical_masks = {}
    for mask in range(1, 1 << m):
        if mask.bit_count() >= 2 and components(mask) == (mask,):
            connected.append(mask)
            canonical_masks[mask] = canonical(mask)

    @functools.cache
    def duplicates(parent: int) -> tuple[int, ...]:
        return tuple(mask for mask in connected if mask & parent == mask)

    @functools.cache
    def solve(state: tuple[int, ...]) -> int:
        groups = {}
        for j, parent in enumerate(state):
            for mask in duplicates(parent):
                groups.setdefault(canonical_masks[mask], []).append((j, mask))
        best, largest_duplicate = 0, 0
        for occurrences in groups.values():
            for position, (j, first) in enumerate(occurrences):
                for k, second in occurrences[position + 1 :]:
                    if first & second:
                        continue
                    size = first.bit_count()
                    largest_duplicate = max(largest_duplicate, size)
                    residual = [
                        fragment & ~(first | second) if t in (j, k) else fragment
                        for t, fragment in enumerate(state)
                    ]
                    child = [first]
                    for fragment in residual:
                        child.extend(components(fragment))
                    best = max(best, size - 1 + solve(tuple(sorted(child))))

        metrics["oracle_states"] += 1
        if largest_duplicate:
            sizes = [mask.bit_count() for mask in state]
            union = functools.reduce(int.__or__, state, 0)
            distinct = len(
                {kind for i, kind in enumerate(edge_types) if union >> i & 1}
            )
            surplus = sum(sizes) - distinct
            vector_bound = (
                surplus - (surplus + largest_duplicate - 1) // largest_duplicate
            )
            scalar_bound = sum(size // 2 for size in sizes) - 1
            for duplicate_size in range(3, largest_duplicate + 1):
                packing = sum(
                    size - (size + duplicate_size - 1) // duplicate_size
                    for size in sizes
                )
                scalar_bound = max(
                    scalar_bound, packing - (duplicate_size - 1).bit_length()
                )
            if vector_bound < best or scalar_bound < best:
                raise AssertionError(
                    {
                        "state": state,
                        "vector_bound": vector_bound,
                        "scalar_bound": scalar_bound,
                        "true_saving": best,
                        "edges": edges,
                        "labels": labels,
                        "bonds": bonds,
                    }
                )
            metrics["vector_bound_states_checked"] += 1
            if vector_bound < scalar_bound:
                metrics["vector_strictly_tighter_than_scalar"] += 1
                metrics["largest_vector_improvement"] = max(
                    metrics["largest_vector_improvement"], scalar_bound - vector_bound
                )
                if len(metrics["vector_examples"]) < 5:
                    metrics["vector_examples"].append(
                        {
                            "fragment_sizes": sizes,
                            "distinct_edge_types": distinct,
                            "maximum_duplicate_size": largest_duplicate,
                            "vector_bound": vector_bound,
                            "scalar_bound": scalar_bound,
                            "exact_remaining_saving": best,
                        }
                    )
        return best

    return m - 1 - solve(tuple(sorted(components((1 << m) - 1))))


@functools.cache
def addition_chain_length(n: int) -> int:
    """Exact IDDFS over all increasing addition chains, used only through 34."""

    def visit(chain: list[int], depth: int) -> bool:
        if chain[-1] == n:
            return True
        if depth == 0 or chain[-1] * (1 << depth) < n:
            return False
        candidates = {
            chain[i] + chain[j]
            for i in range(len(chain))
            for j in range(i + 1)
            if chain[-1] < chain[i] + chain[j] <= n
        }
        return any(
            visit([*chain, value], depth - 1)
            for value in sorted(candidates, reverse=True)
        )

    depth = (n - 1).bit_length()
    while not visit([1], depth):
        depth += 1
    return depth


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solver", type=Path, required=True, help="Compiled oracle_runner executable"
    )
    parser.add_argument("--output", type=Path, help="Write JSON validation evidence")
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("build/dev/libparallelassemblycpp.a"),
        help="Library linked into the supplied runner (recorded as provenance)",
    )
    parser.add_argument("--seed", type=int, default=120928)
    parser.add_argument("--random-cases", type=int, default=400)
    parser.add_argument("--large-cases", type=int, default=100)
    parser.add_argument("--case-timeout", type=float, default=30)
    arguments = parser.parse_args()
    metrics: Metrics = {
        "oracle_states": 0,
        "vector_bound_states_checked": 0,
        "vector_strictly_tighter_than_scalar": 0,
        "largest_vector_improvement": 0,
        "vector_examples": [],
    }
    counts = {
        "exhaustive_small_graphs": 0,
        "random_colored_graphs": 0,
        "factorized_27_edge_graphs": 0,
        "addition_chain_path_graphs": 0,
    }
    started = time.monotonic()
    # The executable is explicitly selected by the audit caller.
    process = subprocess.Popen(  # noqa: S603
        [str(arguments.solver.resolve())],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    ready = selectors.DefaultSelector()
    ready.register(process.stdout, selectors.EVENT_READ)

    def check(
        n: int,
        edges: list[tuple[int, int]],
        labels: list[str],
        bonds: list[int],
        expected: int,
        category: str,
    ) -> None:
        payload = "|".join(
            [
                "audit",
                str(n),
                " ".join(str(v + 1) for edge in edges for v in edge),
                " ".join(labels),
                " ".join(map(str, bonds)),
            ]
        )
        process.stdin.write(payload + "\n")
        process.stdin.flush()
        if not ready.select(arguments.case_timeout):
            raise TimeoutError(f"Solver response exceeded {arguments.case_timeout}s")
        response = process.stdout.readline().rstrip("\n")
        fields = response.split("\t", 4)
        if len(fields) != 5 or fields[1:4] != ["1", "0", "0"]:
            raise RuntimeError(f"Solver failed or stopped incompletely: {response!r}")
        if int(fields[0]) != expected:
            raise AssertionError(
                {
                    "expected": expected,
                    "actual": int(fields[0]),
                    "edges": edges,
                    "labels": labels,
                    "bonds": bonds,
                }
            )
        counts[category] += 1

    try:
        for n in range(3, 6):
            possible = list(itertools.combinations(range(n), 2))
            for bits in range(1 << len(possible)):
                edges = [edge for i, edge in enumerate(possible) if bits >> i & 1]
                if len(edges) < 2:
                    continue
                labels, bonds = ["C"] * n, [1] * len(edges)
                check(
                    n,
                    edges,
                    labels,
                    bonds,
                    oracle(n, edges, labels, bonds, metrics),
                    "exhaustive_small_graphs",
                )
        rng = random.Random(arguments.seed)  # noqa: S311 - reproducible test corpus
        for _ in range(arguments.random_cases):
            n = rng.randrange(4, 9)
            possible = list(itertools.combinations(range(n), 2))
            rng.shuffle(possible)
            edges = possible[: rng.randrange(3, min(len(possible), 11) + 1)]
            labels = [rng.choice(["C", "C", "N"]) for _ in range(n)]
            bonds = [rng.choice([1, 1, 2]) for _ in edges]
            check(
                n,
                edges,
                labels,
                bonds,
                oracle(n, edges, labels, bonds, metrics),
                "random_colored_graphs",
            )
        rng = random.Random(arguments.seed)  # noqa: S311 - reproducible test corpus
        for _ in range(arguments.large_cases):
            edges, labels, bonds, saving = [], [], [], 0
            for group in range(3):
                possible = list(itertools.combinations(range(6), 2))
                rng.shuffle(possible)
                part = possible[:9]
                kinds = [f"{group}_" + rng.choice(["C", "C", "N"]) for _ in range(6)]
                orders = [rng.choice([1, 1, 2]) for _ in part]
                saving += len(part) - 1 - oracle(6, part, kinds, orders, metrics)
                offset = len(labels)
                edges.extend((a + offset, b + offset) for a, b in part)
                labels += kinds
                bonds += orders
            # Distinct labels between gadgets prevent every cross-gadget reuse.
            check(
                len(labels),
                edges,
                labels,
                bonds,
                len(edges) - 1 - saving,
                "factorized_27_edge_graphs",
            )
        for edge_count in range(27, 35):
            edges = [(i, i + 1) for i in range(edge_count)]
            for capped in [False, True]:
                labels = ["C"] * (edge_count + 1)
                if capped:
                    labels[0], labels[-1] = "X", "Y"
                # The two unique terminal edges cannot participate in reuse.
                expected = (
                    addition_chain_length(edge_count - 2) + 2
                    if capped
                    else addition_chain_length(edge_count)
                )
                check(
                    edge_count + 1,
                    edges,
                    labels,
                    [1] * edge_count,
                    expected,
                    "addition_chain_path_graphs",
                )
        process.stdin.close()
        if process.wait(timeout=5) != 0:
            raise RuntimeError(f"Solver runner exited with status {process.returncode}")
        git_executable = shutil.which("git")
        if git_executable is None:
            raise RuntimeError("git is required to record baseline provenance")
        repository = Path(__file__).resolve().parents[2]
        revision = subprocess.run(  # noqa: S603 - discovered git, fixed arguments
            [git_executable, "rev-parse", "HEAD"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        compile_command = shlex.join(
            [
                "c++",
                "-std=c++20",
                "-O2",
                "-Isrc",
                "audits/2026-09-28-branch-and-bound/oracle_runner.cpp",
                str(arguments.library),
                "-o",
                str(arguments.solver),
            ]
        )
        evidence = {
            "baseline_git_revision": revision,
            "runner_sha256": hashlib.sha256(arguments.solver.read_bytes()).hexdigest(),
            "library_path": str(arguments.library.resolve()),
            "library_sha256": hashlib.sha256(
                arguments.library.read_bytes()
            ).hexdigest(),
            "oracle_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "compile_command": compile_command,
            "run_command": shlex.join([sys.executable, *sys.argv]),
            "bound_comparison": (
                "Vector bound compared with generic scalar maxDupBonds(M), "
                "not production targeted/class/pair/post-fragmentation bounds."
            ),
            "status": "passed",
            "seed": arguments.seed,
            "total_solver_cases": sum(counts.values()),
            "case_counts": counts,
            "bound_validation": metrics,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "solver_runner": str(arguments.solver.resolve()),
            "scope": (
                "Serial count-only molecular calculation; "
                "no production code imported by Python oracle."
            ),
            "limits": [
                "Small exhaustive corpus covers simple graphs on 3-5 vertices only.",
                "Random exact graphs have 4-8 vertices and at most 11 edges.",
                "27-edge tests factor into three label-disjoint 9-edge gadgets.",
                "27-34-edge paths use independently exhaustive addition-chain lengths.",
                (
                    "Bounds checked against exact savings at visited oracle states, "
                    "not every production internal pruning route."
                ),
                (
                    "Does not establish exhaustive correctness for arbitrary larger "
                    "graphs or concurrent execution."
                ),
            ],
        }
        if arguments.output:
            arguments.output.write_text(json.dumps(evidence, indent=2) + "\n")
        print(json.dumps(evidence, indent=2))
    finally:
        ready.close()
        if process.poll() is None:
            process.kill()
            process.wait()


if __name__ == "__main__":
    main()
