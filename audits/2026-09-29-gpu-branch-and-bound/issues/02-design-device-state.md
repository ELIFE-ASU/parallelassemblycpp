# GPU-02: Design a flat GPU representation for the immutable DAG and variable search states

**Priority:** P0. **Type:** representation design and bounded prototype.

**Dependencies:** Required: [01 — feasibility measurements](01-measure-gpu-feasibility.md) before fixing capacities or making performance claims. Optional: [07 — canonical-state caching](07-preserve-canonical-state-caching.md) informs shared cache integration; initial serialization work can proceed independently.

## Context and research question

Molecular assembly search holds physical edge masks for fragments, their edge counts and connectivity/canonical metadata, plus accumulated duplicate savings. Its immutable [runtime DAG](../../../src/dagEnumeration.h) already uses compressed sparse row arrays of nodes and transitions, carrying child canonical IDs and retained wide-mask words. However, [owning masks](../../../src/activeWordMask.h) use thread-local copy-on-write arenas; `std::vector` state objects and host pointers cannot simply become device data.

**Question:** Can a documented, bounded device representation preserve all search semantics while keeping the shared DAG resident and moving only compact state batches?

## Bounded investigation

1. Specify a versioned host/device ABI using explicit-width fields, array offsets and lengths. Cover molecule identity, mask width, DAG levels, nodes, transitions, canonical-ID namespace, fragment ranges, mask-word ranges, edge counts, connectivity knowledge, savings and lower bounds. Start from the existing [transfer descriptors](../../../src/searchContext.h), without copying their vector objects.
2. Pack and upload one immutable DAG per molecule. Compare array-of-structures and structure-of-arrays layouts using measured access patterns. Compute one-word child masks inline or store wide masks by validated offsets, matching current behaviour.
3. Specify the distinguished first fragment explicitly. [Enumeration](../../../src/improvedBnB.h) uses its canonical ID as an ordinal ceiling; transitions above that ceiling are excluded. [Canonical keys](../../../src/assemblyState.h) keep this ID first and sort only the remaining IDs. Preserve that role and the comparable ID ordering across batching, compaction and reconstruction. Resolve unknown IDs before an operation requiring canonical identity; never reorder or renumber them casually.
4. Implement CPU pack/unpack and a small device round-trip/transition-filter harness. Exercise one-, two- and wider-word masks, boundary bits, variable fragment counts, empty output batches, disconnected residuals and repeated DAG reuse. Preserve distinct physical occurrences even when their canonical IDs match.

## Correctness and acceptance

Validate offset arithmetic, integer conversion and capacities before access. Require zeroed unused mask bits and explicit allocation ownership/lifetimes. A full output queue must retry, spill or fall back to CPU without dropping work or treating overflow as pruning. Compare reconstructed states and allowed transition sets against CPU results, including different first-fragment ceilings.

Deliver the ABI specification, memory/transfer estimates, round-trip evidence and a recommendation. Rejecting a layout because footprint or packing cost is excessive is acceptable; full device recursion is outside scope.

**Primary reference:** [Gmys et al., GPU branch-and-bound using an Integer-Vector-Matrix representation (2016)](https://doi.org/10.1016/j.parco.2016.01.008). Its representation discipline is relevant; a permutation-specific IVM is not a direct molecular-state encoding.
