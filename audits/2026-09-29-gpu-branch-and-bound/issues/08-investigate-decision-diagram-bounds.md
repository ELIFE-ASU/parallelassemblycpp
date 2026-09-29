# GPU-08: Investigate a sound decision-diagram relaxation for assembly bounds

**Priority:** P3 — exploratory research. **Type:** Formulation and correctness investigation.
**Dependencies:** None for formulation. Use [01: profiling](01-measure-gpu-feasibility.md) before integration; GPU implementation requires a useful, proved CPU relaxation.

## Context and research question

This solver minimizes molecular assembly steps by reusing equivalent, edge-disjoint fragments. [assemblyState.h](../../../src/assemblyState.h) represents remaining fragments and duplicate-bond savings; [improvedBnB.h](../../../src/improvedBnB.h) uses admissible bounds. [dagEnumeration.h](../../../src/dagEnumeration.h) records connected fragment occurrences and expansion relationships. **That fragment DAG is not a decision diagram over assembly decisions.** Shared terminology does not supply a valid relaxation.

Can a bounded-width decision diagram provide materially stronger assembly-index lower bounds, or useful feasible pathways, at acceptable cost? Tardivo, Michel, and van Hoeve separate CPU restricted diagrams from GPU relaxed diagrams. Applying that method first requires assembly-specific dynamic-programming states, transitions, objectives, and sound merging.

## Bounded tasks

1. Define a finite decision ordering and sufficient state for a stated subset of labelled graphs. Explain how fragment production/reuse, remaining structure, and costs affect future choices; show why equal states have equivalent feasible completions.
2. Build an exact CPU diagram for tiny instances. Demonstrate correspondence between its paths and legal assembly pathways, including objective accounting and terminal cases, using [pathwayGenerator.h](../../../src/pathwayGenerator.h) as a witness-format reference.
3. Propose one bounded-width merge operator. Prove that each merged state retains a superset of represented feasible completions, with costs optimistic for minimization. Specify any accompanying dominance rule and prove its safety separately.
4. Implement a small CPU prototype and width sweep. Compare bounds with existing bounds over enumerated tiny states; record construction cost, memory, bound improvement, and resulting search reduction. Restricted diagrams may return incumbents only through feasible, independently checked paths.
5. If evidence supports continuation, sketch layered GPU expansion, filtering, merging, and CPU/GPU queues; estimate memory and batching needs without production integration.

## Correctness requirements

Retain atom/bond labels, graph connectivity semantics, legal reuse, and occurrence overlap constraints, or explicitly relax them in the safe direction. An attractive numerical bound without the merge proof cannot prune exact search. Cross-check tiny optimum values and witnesses against exhaustive enumeration and the existing solver; equality on a sample supplements rather than replaces the proof.

## Deliverables and acceptance

Deliver a formulation note, merge/dominance proofs, a tiny reference prototype, reproducible comparisons, and a go/no-go recommendation. Success may be a counterexample disproving a proposed merge, excessive state growth, or bounds too weak/costly to help, documented clearly enough to prevent an unsound implementation.

## Primary reference

Tardivo, Michel, and van Hoeve, [GPU-Accelerated Relaxed Decision Diagrams for Branch-and-Bound Optimization (CP 2026)](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CP.2026.53). Its knapsack, independent-set, and Golomb-ruler results motivate the architecture; they supply neither an assembly relaxation nor an expected speedup here.
