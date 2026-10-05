"""Independently replay graph-RePair certificates and compare with an exact oracle.

The checker imports no production graph or canonicalization implementation.
Every binary rule and residual occurrence is checked with exact labeled graph
isomorphism; joins and the final physical edge partition are checked directly.
Run with --full for all simple graphs on 3-5 vertices and 400 random graphs.
"""

from __future__ import annotations

import argparse
import copy
import functools
import hashlib
import importlib.util
import itertools
import json
import random
import subprocess
import time
from collections import Counter
from pathlib import Path
from queue import Empty, Queue
from threading import Thread
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from types import ModuleType

Graph = tuple[list[str], list[tuple[int, int, int]]]


def require(condition: bool, explanation: str) -> None:
    if not condition:
        raise AssertionError(explanation)


def vertex_sets(edges: list[tuple[int, int, int]], mask: set[int]) -> set[int]:
    return {vertex for index in mask for vertex in edges[index][:2]}


def component_count(edges: list[tuple[int, int, int]], mask: set[int]) -> int:
    adjacency: dict[int, set[int]] = {}
    for index in mask:
        first, second, _ = edges[index]
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    remaining, result = set(adjacency), 0
    while remaining:
        pending = [remaining.pop()]
        result += 1
        while pending:
            neighbours = adjacency[pending.pop()] & remaining
            remaining -= neighbours
            pending.extend(neighbours)
    return result


def isomorphic(
    atoms: list[str],
    edges: list[tuple[int, int, int]],
    first: set[int],
    second: set[int],
) -> bool:
    """Exact backtracking with vertex invariants, independent of solver hashes."""
    if len(first) != len(second):
        return False
    if first == second:
        return True

    def adjacency(mask: set[int]) -> dict[int, dict[int, int]]:
        output: dict[int, dict[int, int]] = {}
        for index in mask:
            left, right, bond = edges[index]
            output.setdefault(left, {})[right] = bond
            output.setdefault(right, {})[left] = bond
        return output

    left_graph, right_graph = adjacency(first), adjacency(second)
    if len(left_graph) != len(right_graph):
        return False

    def signatures(graph: dict[int, dict[int, int]]) -> dict[int, tuple]:
        return {
            vertex: (
                atoms[vertex],
                tuple(
                    sorted(
                        (bond, atoms[other], len(graph[other]))
                        for other, bond in neighbours.items()
                    )
                ),
            )
            for vertex, neighbours in graph.items()
        }

    left_signatures, right_signatures = (
        signatures(left_graph),
        signatures(right_graph),
    )
    if Counter(left_signatures.values()) != Counter(right_signatures.values()):
        return False
    candidates = {
        left: [
            right
            for right in right_graph
            if left_signatures[left] == right_signatures[right]
        ]
        for left in left_graph
    }
    mapping: dict[int, int] = {}
    available = set(right_graph)

    def visit() -> bool:
        if len(mapping) == len(left_graph):
            return True
        left = min(
            (vertex for vertex in left_graph if vertex not in mapping),
            key=lambda vertex: (
                -sum(other in mapping for other in left_graph[vertex]),
                sum(other in available for other in candidates[vertex]),
                vertex,
            ),
        )
        for right in candidates[left]:
            if right not in available or any(
                left_graph[left].get(other) != right_graph[right].get(mapped)
                for other, mapped in mapping.items()
            ):
                continue
            mapping[left] = right
            available.remove(right)
            if visit():
                return True
            available.add(right)
            del mapping[left]
        return False

    return visit()


def validate_certificate(certificate: dict[str, Any]) -> int:
    """Replay a complete construction certificate, returning its proved bound."""
    require(certificate["schema"] == "graph-repair-assembly-v1", "wrong schema")
    atoms = certificate["atoms"]
    edges = [tuple(edge) for edge in certificate["edges"]]
    require(all(isinstance(atom, str) for atom in atoms), "invalid atom label")
    endpoints = set()
    for edge in edges:
        require(len(edge) == 3, "invalid edge record")
        first, second, bond = edge
        require(
            all(type(value) is int for value in edge)
            and 0 <= first < len(atoms)
            and 0 <= second < len(atoms)
            and first != second
            and bond > 0,
            "invalid edge",
        )
        key = tuple(sorted((first, second)))
        require(key not in endpoints, "parallel edge")
        endpoints.add(key)

    def edge_set(values: list[int]) -> set[int]:
        require(
            isinstance(values, list)
            and all(type(index) is int and 0 <= index < len(edges) for index in values),
            "invalid edge reference",
        )
        result = set(values)
        require(len(result) == len(values), "repeated edge in one fragment")
        require(bool(result), "empty fragment")
        return result

    symbols: dict[int, set[int]] = {}

    def add_symbol(symbol: int, mask: set[int]) -> None:
        require(type(symbol) is int and symbol >= 0, "invalid symbol ID")
        require(symbol not in symbols, "duplicate symbol definition")
        require(component_count(edges, mask) == 1, "disconnected symbol")
        symbols[symbol] = mask

    for terminal in certificate["terminals"]:
        mask = edge_set(terminal["edges"])
        require(len(mask) == 1, "nonprimitive terminal")
        add_symbol(terminal["id"], mask)

    for rule in certificate["rules"]:
        symbol = rule["id"]
        left, right, union = (
            edge_set(rule[key]) for key in ("left_edges", "right_edges", "edges")
        )
        require(not left & right, "join reuses a physical edge")
        require(left | right == union, "join does not produce declared graph")
        require(
            bool(vertex_sets(edges, left) & vertex_sets(edges, right)),
            "nonincident join",
        )
        for child, occurrence in ((rule["left"], left), (rule["right"], right)):
            require(child in symbols and child < symbol, "undefined or cyclic rule")
            require(
                isomorphic(atoms, edges, symbols[child], occurrence),
                "child occurrence does not match its labeled graph",
            )
        add_symbol(symbol, union)

    covered: set[int] = set()
    residual_vertices = []
    for occurrence in certificate["residual"]:
        mask = edge_set(occurrence["edges"])
        symbol = occurrence["symbol"]
        require(symbol in symbols, "unknown residual symbol")
        require(not covered & mask, "overlapping residual occurrences")
        require(
            isomorphic(atoms, edges, symbols[symbol], mask),
            "residual occurrence does not match its labeled graph",
        )
        covered |= mask
        residual_vertices.append(vertex_sets(edges, mask))
    require(covered == set(range(len(edges))), "residual does not cover target")
    components = component_count(edges, covered)
    # Replay the completion joins on physical residual occurrences. Any two
    # fragments sharing a vertex can be joined, including multiple attachments.
    remaining = list(residual_vertices)
    completion_joins = 0
    changed = True
    while changed:
        changed = False
        for first, second in itertools.combinations(range(len(remaining)), 2):
            if remaining[first] & remaining[second]:
                remaining[first] |= remaining.pop(second)
                completion_joins += 1
                changed = True
                break
    require(len(remaining) == components, "completion cannot construct components")
    rule_count, fragment_count = len(certificate["rules"]), len(residual_vertices)
    compensation = certificate["compensate_disjoint"]
    require(type(compensation) is bool, "invalid compensation flag")
    cost = rule_count + completion_joins
    if not compensation:
        cost += max(0, components - 1)
    trivial = max(0, len(edges) - (components if compensation else 1))
    require(certificate["components"] == components, "incorrect component count")
    require(certificate["rule_count"] == rule_count, "incorrect rule count")
    require(
        certificate["remaining_fragments"] == fragment_count, "incorrect fragment count"
    )
    require(certificate["trivial_upper_bound"] == trivial, "incorrect trivial bound")
    require(certificate["upper_bound"] == cost, "incorrect construction cost")
    require(cost <= trivial, "heuristic worsens trivial bound")
    return cost


def native_record(graph: Graph) -> str:
    atoms, edges = graph
    return "|".join(
        [
            "graph-repair-test",
            str(len(atoms)),
            " ".join(str(vertex + 1) for edge in edges for vertex in edge[:2]),
            " ".join(atoms),
            " ".join(str(edge[2]) for edge in edges),
        ]
    )


class Probe:
    """Exchange one-line UTF-8 records with a probe, bounding response waits."""

    def __init__(self, executable: Path, *, compensate: bool = False) -> None:
        command = [str(executable.resolve())]
        if compensate:
            command.append("--compensate-disjoint")
        self.process = subprocess.Popen(  # noqa: S603 - explicitly selected test tool
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            encoding="utf-8",
            bufsize=1,
        )
        # Windows selectors cannot wait on subprocess pipes. A reader thread
        # also keeps the timeout effective if a response has no final newline.
        self.responses: Queue[str | Exception | None] = Queue()
        self.reader = Thread(target=self._read_responses, daemon=True)
        self.reader.start()
        self.queries = 0

    def _read_responses(self) -> None:
        try:
            for line in self.process.stdout:
                self.responses.put(line)
        except (OSError, UnicodeError) as error:
            self.responses.put(error)
        finally:
            self.responses.put(None)

    def query(self, graph: Graph) -> dict[str, Any]:
        self.process.stdin.write(native_record(graph) + "\n")
        self.process.stdin.flush()
        try:
            line = self.responses.get(timeout=30)
        except Empty as error:
            raise TimeoutError("graph repair probe exceeded 30 seconds") from error
        if isinstance(line, Exception):
            raise line
        require(bool(line), "probe exited without a certificate")
        certificate = json.loads(line)
        require(certificate["atoms"] == graph[0], "probe changed atom labels")
        normalized = sorted((min(a, b), max(a, b), bond) for a, b, bond in graph[1])
        actual = sorted(
            (min(a, b), max(a, b), bond) for a, b, bond in certificate["edges"]
        )
        require(actual == normalized, "probe changed target edges")
        validate_certificate(certificate)
        self.queries += 1
        return certificate

    def close(self) -> None:
        self.process.stdin.close()
        if self.process.poll() is None:
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        self.reader.join(timeout=5)
        require(not self.reader.is_alive(), "probe stdout did not close")
        self.process.stdout.close()
        require(self.process.returncode == 0, "probe failed")


def load_oracle() -> ModuleType:
    path = Path(__file__).resolve().with_name("molecular_oracle.py")
    spec = importlib.util.spec_from_file_location("molecular_oracle", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_rejected(certificate: dict[str, Any]) -> None:
    try:
        validate_certificate(certificate)
    except (AssertionError, KeyError, TypeError, IndexError):
        return
    raise AssertionError("checker accepted an intentionally broken certificate")


def run(
    probe_path: Path, *, full: bool = False, random_cases: int | None = None
) -> dict:
    oracle_module = load_oracle()
    metrics = {
        "oracle_states": 0,
        "vector_bound_states_checked": 0,
        "vector_strictly_tighter_than_scalar": 0,
        "largest_vector_improvement": 0,
        "vector_examples": [],
    }
    results = Counter()
    probes = [Probe(probe_path), Probe(probe_path, compensate=True)]
    started = time.monotonic()
    largest_gap = 0

    @functools.cache
    def exact(graph: tuple) -> int:
        atoms, edges = graph
        if not edges:
            return 0
        return oracle_module.oracle(
            len(atoms),
            [(a, b) for a, b, _ in edges],
            list(atoms),
            [bond for _, _, bond in edges],
            metrics,
        )

    def check(graph: Graph, category: str, expected: int | None = None) -> dict:
        nonlocal largest_gap
        if expected is None:
            expected = exact((tuple(graph[0]), tuple(graph[1])))
        certificate = probes[0].query(graph)
        bound = certificate["upper_bound"]
        require(bound >= expected, f"bound {bound} below exact {expected}: {graph!r}")
        largest_gap = max(largest_gap, bound - expected)
        results[category] += 1
        results["strict_improvements"] += bound < certificate["trivial_upper_bound"]
        results["exact_bounds"] += bound == expected
        return certificate

    try:
        # All small topologies include overlap-only repeats, cycles, branching,
        # multi-vertex attachments, disconnected targets and isolated atoms.
        for size in range(3, 6 if full else 5):
            possible = list(itertools.combinations(range(size), 2))
            for bits in range(1 << len(possible)):
                edges = [
                    (*edge, 1)
                    for index, edge in enumerate(possible)
                    if bits >> index & 1
                ]
                if len(edges) >= 2:
                    check((["C"] * size, edges), "exhaustive_graphs")
        rng = random.Random(260929)  # noqa: S311 - reproducible test corpus
        count = random_cases if random_cases is not None else (400 if full else 60)
        for _ in range(count):
            size = rng.randrange(4, 9)
            possible = list(itertools.combinations(range(size), 2))
            rng.shuffle(possible)
            chosen = possible[: rng.randrange(3, min(len(possible), 11) + 1)]
            graph = (
                [rng.choice(["C", "C", "N", "Å"]) for _ in range(size)],
                [(a, b, rng.choice([1, 1, 2])) for a, b in chosen],
            )
            expected = exact((tuple(graph[0]), tuple(graph[1])))
            check(graph, "random_labelled_graphs", expected)
            # Renumbering can change a greedy packing; feasibility and the
            # exact lower comparator must survive it, not necessarily its score.
            permutation = list(range(size))
            rng.shuffle(permutation)
            atoms = [""] * size
            for old, new in enumerate(permutation):
                atoms[new] = graph[0][old]
            shuffled = [
                (permutation[b], permutation[a], bond) for a, b, bond in graph[1]
            ]
            rng.shuffle(shuffled)
            check((atoms, shuffled), "permuted_graphs", expected)

        for size in (0, 1, 2, 3, 7, 15, 27, 31, 32, 33, 63, 64, 65, 127, 128, 129):
            graph = (["C"] * (size + 1), [(edge, edge + 1, 1) for edge in range(size)])
            expected = oracle_module.addition_chain_length(size) if size else 0
            check(graph, "path_graphs", expected)
        shared = (["C"] * 9, [(0, vertex, 1) for vertex in range(1, 9)])
        check(shared, "shared_vertex_graphs", 3)
        triangle_pair = (
            ["C"] * 6,
            [(0, 1, 1), (1, 2, 1), (0, 2, 1), (3, 4, 1), (4, 5, 1), (3, 5, 1)],
        )
        check(triangle_pair, "disconnected_graphs")
        labelled_pair = (
            ["X", "N", "Å", "C"] * 2,
            [
                (offset + first, offset + second, bond)
                for offset in (0, 4)
                for first, second, bond in (
                    (0, 1, 32767),
                    (1, 2, 3),
                    (2, 3, 1),
                    (3, 0, 4),
                )
            ],
        )
        check(labelled_pair, "label_and_cycle_cases")
        check(
            (
                ["C"] * 5,
                [(0, 1, 1), (1, 2, 1), (2, 0, 1), (0, 3, 1), (3, 4, 1), (4, 0, 1)],
            ),
            "label_and_cycle_cases",
        )
        check(
            (
                ["H", "C", "H", "H", "C", "H"],
                [(0, 1, 1), (1, 2, 1), (3, 4, 1), (4, 5, 1)],
            ),
            "label_and_cycle_cases",
        )
        for graph in [
            ([], []),
            (["C", "N"], []),
            (["C"] * 3, [(0, 1, 1)]),
            triangle_pair,
        ]:
            normal = probes[0].query(graph)
            compensated = probes[1].query(graph)
            adjustment = max(0, normal["components"] - 1)
            require(
                compensated["upper_bound"] == normal["upper_bound"] - adjustment,
                "component compensation differs from certificate construction",
            )
            results["compensation_cases"] += 1

        valid = probes[0].query((["C"] * 9, [(edge, edge + 1, 1) for edge in range(8)]))
        require(bool(valid["rules"]), "repeated path produced no rules")
        mutations = []
        corrupt = copy.deepcopy(valid)
        corrupt["upper_bound"] -= 1
        mutations.append(corrupt)
        corrupt = copy.deepcopy(valid)
        corrupt["rules"][0]["left_edges"] = corrupt["rules"][0]["right_edges"]
        mutations.append(corrupt)
        corrupt = copy.deepcopy(valid)
        corrupt["rules"][0]["left"] = corrupt["rules"][0]["id"]
        mutations.append(corrupt)
        corrupt = copy.deepcopy(valid)
        corrupt["residual"].append(copy.deepcopy(corrupt["residual"][0]))
        mutations.append(corrupt)
        corrupt = copy.deepcopy(valid)
        corrupt["residual"].pop()
        mutations.append(corrupt)
        corrupt = copy.deepcopy(valid)
        corrupt["atoms"][0] = "different"
        mutations.append(corrupt)
        for mutation in mutations:
            check_rejected(mutation)
        results["invalid_certificates_rejected"] = len(mutations)
    finally:
        for probe in probes:
            probe.close()
    return {
        "status": "passed",
        "full_corpus": full,
        "probe": str(probe_path.resolve()),
        "probe_sha256": hashlib.sha256(probe_path.read_bytes()).hexdigest(),
        "source_sha256": {
            name: hashlib.sha256(
                (Path(__file__).resolve().parents[1] / name).read_bytes()
            ).hexdigest()
            for name in (
                "src/graphRepair.h",
                "unitTests/graphRepairProbe.cpp",
                "unitTests/graphRepairTester.py",
                "unitTests/molecular_oracle.py",
            )
        },
        "counts": dict(results),
        "certificates_replayed": sum(probe.queries for probe in probes),
        "largest_gap_to_exact": largest_gap,
        "oracle_states": metrics["oracle_states"],
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "scope": (
            "Independent construction replay and exact small-graph oracle comparison"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--random-cases", type=int)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    evidence = run(
        arguments.probe, full=arguments.full, random_cases=arguments.random_cases
    )
    serialized = json.dumps(evidence, indent=2) + "\n"
    if arguments.output:
        arguments.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
