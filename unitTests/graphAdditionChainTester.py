"""Compare exact graph search with an independent exhaustive duplicate oracle.

The oracle uses physical edge masks and permutation-based graph isomorphism;
it imports no production bounds, canonicalization, or enumeration code.
"""

from __future__ import annotations

import argparse
import functools
import itertools
import json
import random
import subprocess
from collections import Counter
from pathlib import Path

from molecular_oracle import addition_chain_length, oracle

Graph = tuple[tuple[str, ...], tuple[tuple[int, int, int], ...]]


def components(graph: Graph) -> int:
    """Include isolated atoms, matching the exact solver's compensation policy."""
    atoms, edges = graph
    neighbours = [set() for _ in atoms]
    for first, second, _ in edges:
        neighbours[first].add(second)
        neighbours[second].add(first)
    remaining = set(range(len(atoms)))
    count = 0
    while remaining:
        pending = [remaining.pop()]
        count += 1
        while pending:
            found = neighbours[pending.pop()] & remaining
            remaining -= found
            pending.extend(found)
    return count


def remove_hydrogens(graph: Graph) -> Graph:
    atoms, edges = graph
    retained = [index for index, atom in enumerate(atoms) if atom != "H"]
    indices = {old: new for new, old in enumerate(retained)}
    return (
        tuple(atoms[index] for index in retained),
        tuple(
            (indices[first], indices[second], bond)
            for first, second, bond in edges
            if first in indices and second in indices
        ),
    )


def record(graph: Graph) -> str:
    atoms, edges = graph
    return "|".join(
        (
            "addition-chain-regression",
            str(len(atoms)),
            " ".join(str(vertex + 1) for edge in edges for vertex in edge[:2]),
            " ".join(atoms),
            " ".join(str(edge[2]) for edge in edges),
        )
    )


def labelled_path_index(bonds: tuple[int, ...]) -> int:
    """Independent bound-free search for paths, whose fragments are substrings.

    All vertices have the same label. Reversing an edge-label substring is
    therefore precisely graph isomorphism; no production graph code is used.
    """

    def canonical(part: tuple[int, ...]) -> tuple[int, ...]:
        return min(part, part[::-1])

    @functools.cache
    def savings(state: tuple[tuple[int, ...], ...]) -> int:
        groups = {}
        for index, part in enumerate(state):
            for start in range(len(part) - 1):
                for end in range(start + 2, len(part) + 1):
                    groups.setdefault(canonical(part[start:end]), []).append(
                        (index, start, end)
                    )
        best = 0
        for fragment, occurrences in groups.items():
            for first, second in itertools.combinations(occurrences, 2):
                left, start, end = first
                right, other_start, other_end = second
                if left == right and end > other_start:
                    continue
                child = [fragment]
                for index, part in enumerate(state):
                    if index == left == right:
                        remnants = (
                            part[:start],
                            part[end:other_start],
                            part[other_end:],
                        )
                    elif index == left:
                        remnants = (part[:start], part[end:])
                    elif index == right:
                        remnants = (part[:other_start], part[other_end:])
                    else:
                        remnants = (part,)
                    child.extend(canonical(part) for part in remnants if len(part) > 1)
                best = max(best, len(fragment) - 1 + savings(tuple(sorted(child))))
        return best

    return len(bonds) - 1 - savings((canonical(bonds),))


def run(probe: Path, *, full: bool = False) -> dict:
    metrics = {
        "oracle_states": 0,
        "vector_bound_states_checked": 0,
        "vector_strictly_tighter_than_scalar": 0,
        "largest_vector_improvement": 0,
        "vector_examples": [],
    }

    @functools.cache
    def exact(graph: Graph) -> int:
        atoms, edges = graph
        return oracle(
            len(atoms),
            [(first, second) for first, second, _ in edges],
            list(atoms),
            [bond for _, _, bond in edges],
            metrics,
        )

    cases: list[tuple[Graph, list[int], str]] = []

    def add(graph: Graph, category: str, expected: int | None = None) -> None:
        values = []
        for remove in (False, True):
            processed = remove_hydrogens(graph) if remove else graph
            value = exact(processed) if expected is None else expected
            values.extend((value, value - max(components(processed) - 1, 0)))
        cases.append((graph, values, category))

    # Exhaust all simple topologies, including cycles, branches, disjoint
    # components, zero bonds, and isolated atoms.
    for size in range(6 if full else 5):
        possible = list(itertools.combinations(range(size), 2))
        for bits in range(1 << len(possible)):
            graph = (
                ("C",) * size,
                tuple(
                    (*edge, 1)
                    for index, edge in enumerate(possible)
                    if bits >> index & 1
                ),
            )
            add(graph, "exhaustive_topologies")

    rng = random.Random(261006)  # noqa: S311 - reproducible test corpus
    for _ in range(240 if full else 48):
        size = rng.randrange(4, 9)
        possible = list(itertools.combinations(range(size), 2))
        rng.shuffle(possible)
        chosen = possible[: rng.randrange(3, min(len(possible), 10) + 1)]
        graph = (
            tuple(rng.choice(("C", "C", "N", "H", "Å")) for _ in range(size)),
            tuple(
                (first, second, rng.choice((1, 1, 2, 3))) for first, second in chosen
            ),
        )
        add(graph, "labelled_and_hydrogen_graphs")
        permutation = list(range(size))
        rng.shuffle(permutation)
        atoms = [""] * size
        for old, new in enumerate(permutation):
            atoms[new] = graph[0][old]
        edges = [(permutation[b], permutation[a], kind) for a, b, kind in graph[1]]
        rng.shuffle(edges)
        add((tuple(atoms), tuple(edges)), "permuted_graphs")

    # Labelled paths with a repeated motif and a unique bond exercise the
    # whole-input count before unique-type pruning and retained-class bounds.
    for kinds in ((1, 2), (1, 2, 3), (1, 1, 2), (1, 2, 3, 4)):
        bonds = (*kinds, *kinds, 7)
        add(
            (
                ("C",) * (len(bonds) + 1),
                tuple((index, index + 1, kind) for index, kind in enumerate(bonds)),
            ),
            "repeated_motifs_with_unique_bond",
        )

    # These larger cases have independent analytic answers. Homogeneous
    # paths realize every scalar chain, including lengths above ceil(log2 n).
    for size in (7, 15, 23, 31, 47, 63):
        add(
            (("C",) * (size + 1), tuple((edge, edge + 1, 1) for edge in range(size))),
            "scalar_chain_gaps",
            addition_chain_length(size),
        )
    for distinct in (5, 8):
        bonds = tuple(range(1, distinct + 1)) * 2
        # Each primitive type must enter the motif, and one reuse doubles it:
        # distinct - 1 joins for the motif plus one join for its second copy.
        add(
            (
                ("C",) * (len(bonds) + 1),
                tuple((index, index + 1, kind) for index, kind in enumerate(bonds)),
            ),
            "vector_composition_gaps",
            distinct,
        )

    for motif in ((1, 2, 3, 4), (1, 2, 1, 3, 4)):
        for bonds in (
            (*motif, 5, *motif, 5, *motif, 6),
            (*motif, 5, *reversed(motif), 6, *motif),
            (*motif, *motif, 5, *motif[:3], 6, *motif),
        ):
            add(
                (
                    ("C",) * (len(bonds) + 1),
                    tuple((index, index + 1, kind) for index, kind in enumerate(bonds)),
                ),
                "nested_labelled_paths",
                labelled_path_index(bonds),
            )

    completed = subprocess.run(  # noqa: S603 - explicitly selected test executable
        [str(probe.resolve())],
        input="".join(record(graph) + "\n" for graph, _, _ in cases),
        capture_output=True,
        text=True,
        timeout=240 if full else 120,
        check=False,
    )
    if completed.returncode:
        raise AssertionError(f"probe failed: {completed.stderr}")
    output = completed.stdout.splitlines()
    if len(output) != len(cases):
        raise AssertionError(
            f"expected {len(cases)} probe rows, received {len(output)}"
        )
    for line, (graph, expected, category) in zip(output, cases, strict=True):
        actual = json.loads(line)
        if actual != expected:
            raise AssertionError(
                f"{category}: exact results {actual} differ from oracle {expected}: "
                f"{graph!r}"
            )
    return {
        "graphs": len(cases),
        "exact_calculations": 4 * len(cases),
        "categories": dict(Counter(category for _, _, category in cases)),
        "oracle_states": metrics["oracle_states"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("--full", action="store_true")
    arguments = parser.parse_args()
    print(json.dumps(run(arguments.probe, full=arguments.full), sort_keys=True))
