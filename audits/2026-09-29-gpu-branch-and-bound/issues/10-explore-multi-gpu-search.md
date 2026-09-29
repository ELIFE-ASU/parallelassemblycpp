# GPU-10: Conditionally investigate multi-GPU assembly search scheduling

**Priority/type:** P3 / conditional scalability investigation.

**Required dependencies:** a successful measured single-GPU result from [06: GPU-resident search](06-explore-gpu-resident-search.md) **or** [05: hybrid batched search](05-integrate-hybrid-batched-search.md), plus [02: device-state ABI](02-design-device-state.md). The selected route must include [07: caching assessment](07-preserve-canonical-state-caching.md). **Optional dependency:** additional hardware/platform measurements extending [01](01-measure-gpu-feasibility.md).

## Context and question

Exact molecular assembly search explores unequal duplicate-removal subtrees, sharing an assembly-index incumbent for pruning. [searchContext.h](../../../src/searchContext.h) supports CPU-local donation; [main.cpp](../../../src/main.cpp) implements MPI root refill and incumbent exchange. Transferable descriptors encode canonical-ID ownership and require adaptation across processes or devices.

Can a single-GPU benefit scale across two GPUs despite migration, canonicalization, and coordination costs? Start only after the dependency gate passes. Without suitable hardware and a positive single-GPU result, document deferral.

## Bounded investigation

1. Begin with one host and two GPUs. Compare independent root leases against hierarchical balancing: device-local donation first, then coarse inter-device task movement through the host. Reuse immutable molecule/DAG data where practical; document replication and memory costs.
2. Specify transferable branch descriptors with globally meaningful structural identity or receiver-side canonicalization. Separate shared immutable state from private DFS state. Tune only lease size, refill threshold, and incumbent exchange interval; record bytes moved, imbalance, tail idle time, and stolen-task usefulness.
3. Extend MPI progress only if two-GPU measurements justify it. Respect funneled MPI threading. Compare host-staged transfers with an optional peer/device-aware path.
4. Define distributed work accounting covering local stacks, device queues, host spill queues, leases, pending transfers, and acknowledgements. Exercise empty-frontier, delayed-message, cancellation, and memory-pressure cases before timing larger instances.

## Correctness constraints

Assembly-index incumbents are feasible upper bounds and must monotonically decrease; late messages may not overwrite a better value. Task handoff requires explicit ownership so no branch disappears during transfer or spill. Queue/stack exhaustion must retain work or use CPU fallback. Global termination requires all workers idle and no queued, leased, spilled, or in-flight work, not merely empty GPU queues. Invalid canonical-ID comparisons and hash-only duplicate suppression are forbidden. Cancellation remains distinguishable from a completed optimality proof. Preserve reconstruction provenance and validate any returned construction.

## Acceptance and deliverables

Deliver the dependency decision, protocol/state diagram, reproducible two-GPU experiment, correctness checks, and scaling report against the successful single-GPU baseline and existing CPU/MPI baseline. Report node counts, communication, memory, reconstruction, and end-to-end timing. A no-go result documenting negative scaling or insufficient work completes this investigation; cluster deployment is a separate decision.

## Primary literature

- [Gmys et al., hierarchical multi-GPU work stealing (2017)](https://doi.org/10.1002/cpe.4019): transferable scheduling hierarchy, with permutation-specific task intervals.
- [Gmys, exact flowshop search at scale (2022)](https://doi.org/10.1287/ijoc.2022.1193) and [PBB implementation](https://github.com/jangmys/pbb): distributed GPU branch-and-bound.
- [Almasri et al., GPU clique enumeration (2023)](https://arxiv.org/html/2212.01473v3): useful device-local balancing; its inter-GPU shared worker list did not repay its overhead.
