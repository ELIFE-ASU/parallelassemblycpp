# GPU-03: Prototype batched DAG expansion and duplicate-class construction

**Priority:** P1 — first kernel experiment. **Type:** Bounded prototype and benchmark.
**Dependencies:** [GPU-02](02-design-device-state.md) supplies a flat device
representation; [GPU-01](01-measure-gpu-feasibility.md) supplies representative
states and a current CPU baseline. Kernel design can start alongside profiling.

## Context and research question

The molecular assembly solver searches for reusable, edge-disjoint isomorphic
fragments. After initial enumeration, each search state reuses an immutable DAG
of connected fragments. DAG expansion generates occurrences, groups them by
canonical identity, and constructs eligible-edge masks for subsequent bounds.

The current DAG already stores transitions and retained masks in contiguous
arrays. Historical profiles found substantial cost in enumeration and class
construction, but those measurements precede subsequent CPU optimizations.
Can processing many states together make this stage faster after counting input
packing, output compaction, transfers and additional search work?

## Investigation

1. Capture representative CPU states, their first-fragment ordinal cutoff, and
   the complete enumeration outputs. Include small and large occurrence classes,
   repeated labels, cycles, split residual fragments and 64-bit mask boundaries.
2. Upload the immutable DAG once. Expand batches of parent occurrences using
   indexed transitions, edge-membership checks and canonical-ID cutoff checks.
   Start with one expansion level before reproducing the full enumeration loop.
3. Group output by state and canonical class; retain parent-fragment and DAG-node
   identities. Implement segmented mask unions and survivor compaction. Compare
   stable ordering with a reordered variant only after correctness is established.
4. Reproduce occurrence viability, target-mask accumulation, ordinal-overflow
   handling and final-level behavior. Measure grouping and allocation costs
   separately from raw transition throughput.
5. Sweep batch sizes and occurrence counts. Use bounded output buffers with
   lossless chunking or CPU fallback. Record transferred bytes, buffer peaks,
   batch latency and effective CPU-versus-GPU processing time.

## Correctness requirements

No eligible occurrence may disappear during grouping, compaction or overflow.
Preserve the distinguished first fragment and its cutoff semantics, including
the `overweight`/`last` transition between enumeration levels. Compare complete
occurrence multisets and target masks against CPU outputs; equality of the final
assembly index alone is insufficient. Output reordering requires separate
measurement of its effect on incumbent discovery and total search work.

## Deliverables and acceptance

- A reproducible replay benchmark and experimental implementation, with exact
  output comparisons against CPU enumeration on the captured states.
- Results separating kernel time, packing/transfers, memory and batch latency,
  with a measured crossover range or evidence that none was found.
- A recommendation for integration through GPU-05, further representation work,
  or closure. A faster isolated kernel does not establish solver speedup.

## Starting points and sources

- [`dagLevel`](../../../src/dagEnumeration.h), `dagGenerate`,
  `dagDuplicateGenerator` and `duplicateClassLevel::seal` in
  [duplicateMatching.h](../../../src/duplicateMatching.h), and
  `dagRecursiveEnumeration` in [improvedBnB.h](../../../src/improvedBnB.h).
- Historical profile report at `build/mask-results/README.md`, a local
  unpublished artifact not included in a fresh checkout. Recollect evidence
  using the [benchmark tooling](../../../benchmarks/README.md).
- Helbecque et al., [Portable PGAS-Based GPU-Accelerated Branch-And-Bound
  Algorithms at Scale (2025)](https://doi.org/10.1002/cpe.70321): batching
  architecture; applying it to assembly enumeration is this issue's hypothesis.
