# GPU-07: Preserve exact canonical-state dominance across CPU and GPU work

**Priority:** P1. **Type:** correctness design and cache-placement experiment.

**Dependencies:** Required: [01 — feasibility measurements](01-measure-gpu-feasibility.md) before performance conclusions. Optional: [02 — device-state design](02-design-device-state.md) provides a transfer ABI; the correctness analysis and CPU cache ablations can start independently.

## Context and research question

Molecular assembly search can reach equivalent fragment states through different matching choices. The [transposition table](../../../src/assemblyTranspositionTable.h) stores the largest duplicated-bond score `D` reaching each canonical key; with the same residual choices, an equal or larger stored score dominates a revisit. The [key](../../../src/assemblyState.h) retains the first fragment's canonical ID and sorts only the remainder, because that first ID restricts subsequent enumeration. Canonical equivalence is therefore more than an unordered set of fragment hashes.

**Question:** Where should canonicalisation and exact state caching live in a heterogeneous solver, and how much useful pruning would be lost by GPU-local or bounded caches?

## Bounded investigation

1. Document current key construction, equality, score updates and the [L1-before-L2 lookup sequence](../../../src/searchContext.h). Keep exact CPU canonicalisation and host-owned authoritative lookup as the first experimental control. Compare it with batched host lookup, bounded admission and, only if justified, an exact device-cache prototype.
2. Specify canonical-ID validity. [Transferred tasks](../../../src/improvedBnB.h) can reuse immutable root-seed IDs, IDs from the active process registry, or the originating worker's private IDs. Unrelated private IDs and independently created process IDs are not interchangeable integers. Include molecule/namespace identity in transport contracts, and remap or recanonicalise when validity cannot be established.
3. Distinguish score dominance from completed-subtree knowledge. Entries are considered before recursive completion; GPU batches or donated descendants can remain outstanding. A stored score does not certify an exact residual optimum. Any future residual-bound/completion cache needs separate proof status, child aggregation and cancellation handling.
4. Run controlled cache-placement/admission ablations with matched bounds and scheduling. Measure retained bytes for local tables, shared tables, canonical registries and device queues separately, alongside hit rates, canonicalisation costs, expanded states, matching work and total completion time. Include long amino-acid and Paclitaxel cases; reduced storage alone is not success.

## Correctness and acceptance

Hash equality must never imply key equality: compare full key length and contents before pruning, including deliberate collisions. Test namespace mismatch, concurrent equal/better-score arrivals, cache exhaustion, reordered batches and cancellation. Skipped admission or forgotten entries may cause re-expansion; they must not imply dominance. Preserve pending work when pruning equivalent arrivals.

Deliver a correctness contract, adversarial tests, ablation data and a placement recommendation. Keeping all canonicalisation and authoritative caching on CPU, or finding no beneficial cache change, is acceptable.

**Primary reference:** [Coppé, Gillard and Schaus, *Decision Diagram-Based Branch-and-Bound with Caching for Dominance and Suboptimality Detection* (2024)](https://arxiv.org/abs/2211.13118). Its expansion-threshold caching motivates investigation; its proofs do not automatically establish molecular-state cache correctness.
