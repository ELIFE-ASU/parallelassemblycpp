"""Independent exact molecular and addition-chain oracles for graph tests.

Enumerates every connected edge subset and every legal duplicate pair, with
independent brute-force labeled graph isomorphism. The reference search uses
physical-mask memoization only: no production bounds, DAG or canonical ordinal.
Also checks the edge-type vector bound at every visited state.
"""

import functools
import itertools
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
    """Find the shortest increasing addition chain for a positive integer."""

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
