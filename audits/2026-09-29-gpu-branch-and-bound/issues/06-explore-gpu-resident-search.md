# GPU-06: Evaluate a bounded GPU-resident assembly search prototype

**Priority/type:** P2 / exploratory architecture investigation.

**Required dependencies:** [01: feasibility measurements](01-measure-gpu-feasibility.md), [02: device-state ABI](02-design-device-state.md), and [07: canonical-state caching assessment](07-preserve-canonical-state-caching.md). **Optional dependency:** measured results from [05: hybrid batched search](05-integrate-hybrid-batched-search.md).

## Context and question

Exact molecular assembly search removes duplicate subgraph occurrences to minimize construction joins. Branches carry fragment masks, duplicated-bond savings, bounds, and canonical identity. Branch sizes vary, and cache hits change exploration. [searchContext.h](../../../src/searchContext.h) provides CPU work donation and transferable descriptors; their host containers and canonical-ID validity domains prevent direct device use.

Can GPU-resident depth-first exploration outperform CPU or hybrid search while preserving pruning? No speedup is assumed. Distinguish scheduling improvements from extra search caused by weaker canonicalization or caching.

## Bounded investigation

1. Select one measured workload family and explicitly supported mask/depth limits. Define compact DFS frames using the ABI from issue 02: indices into immutable data, fragment masks or reversible deltas, branch cursors, and objective metadata. Account for every worker's stack and scratch memory before launch.
2. Compare warp-per-subtree and block-per-subtree mappings. Measure occupancy, divergence, memory traffic, operations per node, and remaining idle workers near completion. Avoid a full breadth-first frontier expansion.
3. Compare bounded full-state work queues against an idle-worker list. In the latter, a busy explorer donates directly into a waiting worker's storage, with release/acquire publication and a measured minimum donation size. Limit this comparison to identical task sets.
4. Implement the canonicalization/cache strategy selected in issue 07, including any CPU handoff boundary. Measure its synchronization and transfer costs. Keep the existing bounds from [assemblyState.h](../../../src/assemblyState.h) and candidate semantics from [fragmentation.h](../../../src/fragmentation.h).

## Correctness constraints

Every accepted task must complete, transfer ownership, or spill losslessly to CPU search. Unsupported masks, depth overflow, and exhausted buffers require explicit fallback. Stale valid incumbents only weaken pruning; incumbent improvements must remain monotonic. Hash collisions cannot establish canonical equality, and worker-local canonical IDs cannot silently become global. Do not disable caching without measuring the resulting search expansion. Preserve complete duplicate routes, cancellation semantics, and sufficient provenance for valid pathway reconstruction.

## Acceptance and deliverables

Provide a runnable bounded prototype, memory budget, ownership/termination description, and CPU differential checks including forced overflow and donation. Report completed exact results, node counts, cache behavior, transfer costs, and end-to-end time including reconstruction separately where enabled. Deliver a go/no-go recommendation; unfavorable performance or incompatible state/caching costs are acceptable conclusions.

## Primary literature

- [Yamout et al., GPU vertex cover (2022)](https://arxiv.org/abs/2204.10402): local DFS stacks and global worklists.
- [Almasri et al., GPU maximal-clique enumeration (PACT 2023)](https://arxiv.org/html/2212.01473v3): demand-driven worker-list donation.
- [Cardone et al., GPU maximum clique (2024)](https://www.scitepress.org/publishedPapers/2024/128527/pdf/index.html): warp/block granularity.
- [Gmys et al., GPU branch-and-bound (2016)](https://doi.org/10.1016/j.parco.2016.01.008) and [PBB code](https://github.com/jangmys/pbb): compact permutation search; factoradic intervals do not directly represent assembly states.
