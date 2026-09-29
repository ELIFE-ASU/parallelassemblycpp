# GPU-01: Measure whether molecular assembly search has enough GPU-suitable work

**Priority:** P0. **Type:** measurement and feasibility investigation.

**Dependencies:** Required: none. Optional: [02 — device-state design](02-design-device-state.md) supplies concrete transfer sizes; measurement can start independently.

## Context and research question

This solver minimises molecular assembly index through a reverse search that removes one of two eligible isomorphic, edge-disjoint fragments and credits reuse savings. For a graph with `N` bonds and accumulated duplicate savings `D`, a feasible index is `N - 1 - D` before optional disconnected-input compensation. Search combines duplicate enumeration, matching, residual fragmentation, bounds, exact canonicalisation and transposition pruning. GPU throughput on one operation need not improve completion time if batching changes incumbent discovery or cache reuse.

Historical profiles in the local, potentially untracked `build/mask-results/README.md` reported 58% of pre-optimisation serial Paclitaxel cycles in enumeration/class processing. Those measurements describe older binaries; they neither establish today's bottleneck nor predict GPU speedup. The present baseline must be recorded independently, starting from revision `04129ee13bdaf165001afddecc22ea66b75a5285` or a documented successor.

**Question:** Which complete stages can supply sufficiently large batches to overcome packing, transfer, scheduling and extra search costs on available GPU hardware?

## Bounded investigation

1. Freeze source, executable hashes, compiler flags, CPU/GPU topology, placement, cache policy and pathway setting. Select Paclitaxel, long amino-acid cases, one-word molecules, mask-boundary cases and short negative controls from the [manifest](../../../benchmarks/cases.tsv).
2. Extend existing [telemetry](../../../src/searchTelemetry.h) only where necessary. Measure exclusive or explicitly nested costs for enumeration, class grouping, overlap tests, bounds, fragmentation, canonicalisation, canonical-ID registry operations and transposition lookup. Do not sum overlapping timers.
3. Record distributions of active fragments, mask words, class sizes, valid pairs, transitions, surviving candidates and available queue batches, including tails. Count expanded/pruned states separately from classes, pairs and blocks. Record steady-clock strict incumbent improvements and final proof time.
4. Establish current serial and tuned OpenMP end-to-end controls using the [benchmark protocol](../../../benchmarks/README.md). Keep instrumented runs separate from timing samples. Measure GPU setup, packing, copies, kernel execution, synchronisation and CPU survivor work when a prototype exists. Report optimisation and optional pathway reconstruction separately, with matched settings.

## Deliverables and acceptance

Deliver reproducible commands, raw measurements, distributions and a feasibility report stating which stage to prototype and why. An Amdahl-style ceiling and measured transfer/launch costs should bound expectations. Require exact assembly-index agreement and explain changes in total search work. A supported conclusion that batches are too small, CPU canonicalisation dominates, or no end-to-end gain is likely is successful completion; no CUDA implementation is required.

**Primary reference:** [Helbecque et al., *Portable PGAS-Based GPU-Accelerated Branch-And-Bound Algorithms at Scale* (2025)](https://doi.org/10.1002/cpe.70321) motivates evaluating complete heterogeneous execution, without transferring its speedups to molecular assembly.
