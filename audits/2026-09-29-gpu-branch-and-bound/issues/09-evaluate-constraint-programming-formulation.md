# GPU-09: Evaluate an exact constraint-programming formulation of tiny assembly instances

**Priority:** P3 — alternative formulation research. **Type:** Exact-model feasibility experiment.
**Dependencies:** None for modelling. Consult [01: profiling](01-measure-gpu-feasibility.md) before performance comparisons or integration; no GPU port is required.

## Context and research question

The repository minimizes joining steps for a labelled molecular graph with fragment reuse. Search tracks fragments and duplicate savings in [assemblyState.h](../../../src/assemblyState.h), branches over matching occurrences in [improvedBnB.h](../../../src/improvedBnB.h), and captures pathways through [pathwayGenerator.h](../../../src/pathwayGenerator.h).

Can a finite integer/Boolean model express those semantics exactly for tiny graphs, and can a GPU constraint solver solve it usefully? Talbot's Turbo implements GPU propagation and backtracking using ternary constraint networks and on-demand “dive and solve” subproblem generation. It requires a constraint formulation, not arbitrary C++ recursion, graph objects, or caches.

## Bounded tasks

1. Choose and publish a small, bounded graph domain and step horizon. Include labelled paths, branching trees, cycles, repeated fragments, and overlapping occurrence traps. State hydrogen preprocessing, bond/atom labels, connectivity, and input assumptions explicitly, matching the repository.
2. Write one finite model, such as bounded fragment-production/reuse decisions. Define variables, legal joins, identity, availability, reuse, and objective. Prove the horizon covers an optimum; any fixed candidate catalogue must be complete for the stated domain.
3. Establish both directions: legal pathways map to satisfying assignments of equal cost, and satisfying assignments decode to legal pathways of equal cost. Prove symmetry breaking retains an optimal representative.
4. Solve tiny cases first with an independently inspectable CPU reference. Export a Turbo-compatible model only after semantic checks pass. Record model-generation time, variables/constraints before and after decomposition, memory, solve/proof times, and witness validation separately.
5. Compare exact optima and decoded witnesses against tiny exhaustive enumeration and the current solver. Distinguish finding a feasible result from proving optimality; preserve timeout and unsupported-instance outcomes.

## Correctness requirements

Prevent free creation of fragments, reuse before production, illegal overlap or joins, and label erasure. Do not substitute disjoint-fragment packing for complete assembly unless equivalence is proved. Decode witnesses to original graph identifiers and validate them independently of the CP constraints. Any declared limitation must be enforced at model input, not silently ignored.

## Deliverables and acceptance

Deliver the specification, correspondence argument, exporter/reference implementation, fixtures with pathways, reproducible commands, and recommendation. Useful negative outcomes include a semantic counterexample, impractical model growth, or slower complete solves including generation/transfers. Kernel speed or objective values alone cannot justify integration without feasibility and optimality evidence.

## Primary references

Pierre Talbot, [A GPU-based Constraint Programming Solver (AAAI 2026)](https://ptal.github.io/papers/aaai2026.pdf), with [Turbo's paper-version source](https://github.com/ptal/turbo/tree/aaai2026). Its block-cooperative propagation, shared incumbent, and lazy subproblem generation are architectural references, not evidence that molecular assembly already has an efficient CP encoding.
