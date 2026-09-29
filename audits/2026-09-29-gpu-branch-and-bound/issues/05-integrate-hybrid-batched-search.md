# GPU-05: Integrate CPU-controlled search with batched GPU enumeration

**Priority:** P2 — gated implementation experiment. **Type:** Performance prototype.
**Dependencies:** [02: device representation](02-design-device-state.md), [03: enumeration/grouping](03-batch-dag-enumeration.md), and fresh [01: profiling](01-measure-gpu-feasibility.md). [04: GPU filters](04-batch-matching-and-fragmentation.md) is optional, subject to its measurements.

## Context and research question

Molecular assembly search explores reuse of equivalent, edge-disjoint fragments, minimizing joining steps. It combines DAG enumeration, bounds, canonicalisation, transposition caching, and OpenMP workers. Entry points are [improvedBnB.h](../../../src/improvedBnB.h), [searchContext.h](../../../src/searchContext.h), [dagEnumeration.h](../../../src/dagEnumeration.h), and [assemblyTranspositionTable.h](../../../src/assemblyTranspositionTable.h).

Can batching independent states expose enough GPU parallelism to improve complete solves after transfer, scheduling, and additional-search costs? Keep frontier ownership, canonicalisation, cache decisions, accepted child creation, and pathway witnesses on CPU. Retain the immutable fragment DAG and masks on GPU; submit batches for occurrence expansion, canonical-class grouping, and any separately validated filters. This adapts hybrid B&B architecture without assuming the current cheap arithmetic bounds deserve offloading.

## Bounded tasks

1. Add an opt-in coordinator behind the CPU interface; preserve CPU fallback. Define request/result ownership and stable state identifiers.
2. Batch independent states, grouping compatible mask widths or work sizes. Adapt thresholds using profiling and kernel measurements; keep small or latency-sensitive work on CPU.
3. Overlap pinned-buffer transfers and GPU execution with CPU search using bounded asynchronous queues. Account for packing, grouping, copying, waiting, and device initialization separately.
4. Bound host/device memory. Specify queue backpressure, lossless CPU spill, oversized-request fallback, and cancellation while work is queued or in flight. Avoid making multi-GPU support part of this experiment.
5. Compare complete runs against current serial and OpenMP builds on fresh small, irregular, and search-heavy cases. Record hardware, revisions, wall time, nodes, incumbent timing, memory, preprocessing, and reconstruction costs.

## Correctness requirements

Every submitted state must finish exactly once or return to CPU ownership; buffer overflow must never truncate candidates. Preserve label equivalence, occurrence multiplicity, disjointness, canonical state keys, and valid winning witnesses. A stale incumbent may cause extra work only: pruning must use a valid bound and an incumbent backed by a feasible pathway. Interrupted runs must retain existing incomplete-result semantics. Exercise delayed results, cancellation, memory pressure, and frequent incumbent improvements against the CPU reference.

## Deliverables and acceptance

Deliver the opt-in prototype, differential checks, reproducible benchmark commands, and a recommendation. Accept promotion only with equal exact answers and valid pathways plus repeatable end-to-end benefit. A measured negative result identifying the break-even batch size or dominant overhead is a successful outcome; retain CPU execution by default when acceleration is unjustified.

## Primary references

- Helbecque et al., [Portable PGAS-Based GPU-Accelerated Branch-And-Bound Algorithms at Scale (2025)](https://onlinelibrary.wiley.com/doi/10.1002/cpe.70321), with [source code](https://github.com/Guillaume-Helbecque/GPU-accelerated-tree-search-Chapel).
- Tagliaferro et al., [Performance and Portability in Multi-GPU Branch-and-Bound (IPDPSW 2025)](https://orbilu.uni.lu/handle/10993/66041): CPU work pools and CUDA/HIP batched evaluation; its PFSP timings do not establish an assembly speedup.
