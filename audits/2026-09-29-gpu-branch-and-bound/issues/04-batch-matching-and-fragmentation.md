# GPU-04: Investigate GPU batches for duplicate matching and residual fragmentation

**Priority/type:** P1 / bounded performance investigation.

**Required dependencies:** [01: measure feasibility](01-measure-gpu-feasibility.md) and [02: device-state ABI](02-design-device-state.md). **Optional dependency:** [03: DAG enumeration](03-batch-dag-enumeration.md), if it supplies device-resident occurrence batches.

## Context and question

This solver minimizes molecular assembly joins by finding reusable, isomorphic, edge-disjoint subgraphs. Matching occurrences produces child assembly states; bounds and canonical-state caching then eliminate unnecessary exploration. [duplicateMatching.h](../../../src/duplicateMatching.h) currently streams valid pairs rather than storing every pair. [fragmentation.h](../../../src/fragmentation.h) removes occurrence masks, splits residual connected components, and postpones canonicalization until after inexpensive bounds.

Can GPU batches accelerate mask overlap filtering enough to improve complete searches, and, only if justified, also accelerate residual connected-component extraction? Candidate generation and graph canonicalization remain on the CPU in the first experiment.

## Bounded investigation

1. Capture representative occurrence classes from completed CPU searches, including same-parent and different-parent pairs, sparse classes, highly repetitive molecules, and single-word/wide masks. Record pair rejection rates, active mask words, and incumbent refresh frequency.
2. Prototype tiled overlap checks with bounded output buffers. Generate pair indices within tiles; do not materialize a quadratic pair array. Return compact surviving occurrence indices and a continuation cursor. Compare lane-per-pair and cooperative wide-mask processing.
3. Measure packing, transfers, launches, compaction, and CPU consumption together. Compare against CPU streaming with identical inputs and incumbents. Sweep a documented batch-size range.
4. Proceed to residual-component extraction only if overlap filtering leaves sufficient work. Compare batched mask subtraction and connectivity against [ufds.h](../../../src/ufds.h); preserve fragment edge counts, connectivity metadata, and retained-copy identity. Return residual masks for CPU canonicalization.

## Correctness constraints

Preserve every eligible duplicate-removal route and the existing admissible bounds, including [matchingBoundRefresh.h](../../../src/matchingBoundRefresh.h). Delayed incumbent reads may cause extra work, never invalid pruning. Buffer exhaustion must resume or fall back without dropping pairs. Match full surviving pair sets and child fragment multisets against CPU results; compare identities where pathway reconstruction requires them. Resolve unknown canonical IDs before cache lookup or recursive enumeration. Residual connected components are state fragments: their assembly costs must not be assumed independent because reuse can couple components.

## Acceptance and deliverables

Deliver a reproducible prototype or measured rejection, differential results, bounded-memory behavior, and phase/end-to-end timings using [benchmark conventions](../../../benchmarks/README.md). Report node-count changes separately from kernel throughput. A negative result is acceptable if it identifies transfer, small batches, connectivity, or canonicalization as the limiting cost and recommends a CPU threshold or stopping this route.

## Primary literature

- [VanCompernolle et al., GPU maximum clique with bitsets (2016)](https://www.cse.unr.edu/~fredh/papers/conf/158-mcsubog/paper.pdf): word-parallel set operations and small-instance limitations.
- [Quer et al., parallel maximum common subgraph (2020)](https://www.mdpi.com/2079-3197/8/2/48): CPU/GPU graph-search batching.
- [Amro et al., component-aware GPU vertex cover (2026)](https://arxiv.org/abs/2512.18334): connectivity and compact state; its additive decomposition is not an assembly proof.
