# GPU branch-and-bound: investigation issues

Ten self-contained issue drafts for exploring GPU acceleration of
`parallelassemblycpp`. Each file can be copied into a separate issue tracker
entry. This package records investigations, not a commitment to implement every
approach. No external issues have been published.

Literature was researched on 28 September 2026 and packaged on 29 September
2026. Repository entry points were checked against commit
`04129ee13bdaf165001afddecc22ea66b75a5285`. Historical performance observations
belong to the revisions in their original reports; they are not measurements
of this commit or predictions of GPU performance.

## Issue list

| ID | Self-contained issue | Main question | Sequence |
| --- | --- | --- | --- |
| GPU-01 | [Measure workload and GPU feasibility](issues/01-measure-gpu-feasibility.md) | Which work is substantial, batchable and worth accelerating today? | Start here |
| GPU-02 | [Design device state and a resident DAG](issues/02-design-device-state.md) | Can states and immutable fragment data be represented compactly and transferred exactly? | Foundation; design alongside GPU-01 |
| GPU-03 | [Batch DAG enumeration and class construction](issues/03-batch-dag-enumeration.md) | Does batched expansion, grouping and mask reduction repay its overhead? | First kernel experiment after GPU-02 |
| GPU-04 | [Batch matching filters and residual fragmentation](issues/04-batch-matching-and-fragmentation.md) | Can pair filtering and residual processing reduce CPU work economically? | Additional kernel experiment |
| GPU-05 | [Integrate CPU-controlled, GPU-batched search](issues/05-integrate-hybrid-batched-search.md) | Does batching improve complete solve time with CPU canonicalization and caching? | Integrate GPU-02/03; add GPU-04 if useful |
| GPU-06 | [Explore GPU-resident subtree search](issues/06-explore-gpu-resident-search.md) | Can cooperative DFS and dynamic donation outperform hybrid execution? | Larger backend experiment; use GPU-02/07 |
| GPU-07 | [Preserve canonical identity and state caching](issues/07-preserve-canonical-state-caching.md) | How can a GPU backend retain exact dominance and effective state reuse? | Design alongside GPU-06 |
| GPU-08 | [Investigate GPU decision-diagram bounds](issues/08-investigate-decision-diagram-bounds.md) | Can a sound assembly relaxation provide stronger bounds worth GPU evaluation? | Independent mathematical research |
| GPU-09 | [Evaluate an exact constraint-programming formulation](issues/09-evaluate-constraint-programming-formulation.md) | Can a compact assembly model use a general GPU constraint solver? | Independent alternative formulation |
| GPU-10 | [Explore multi-GPU search](issues/10-explore-multi-gpu-search.md) | Does a successful single-GPU design scale while preserving exact completion? | Conditional on GPU-05 or GPU-06 evidence |

The recommended first sequence is GPU-01, GPU-02, GPU-03, then GPU-05. GPU-04
can expand the useful offload region. GPU-06/07 form a separate, more substantial
backend direction. GPU-08 and GPU-09 explore alternative formulations and can
end with a mathematical feasibility result before any GPU implementation.

## Shared problem definition

A molecular input is a labelled graph. The solver minimizes the joining steps
needed to build it while reusing fragments already constructed. Its reverse
search removes one of two eligible, edge-disjoint isomorphic fragments of size
`k`, adding `k - 1` to accumulated savings `D`. For an original graph with `N`
edges, the search objective is `A = N - 1 - D` before any output compensation
for disconnected inputs. An upper bound `B` on further savings gives a lower
bound `N - 1 - D - B` on the attainable index.

Pruning must preserve that orientation. Independent bounds on remaining savings
combine with `min`, not addition. A feasible pathway establishes achievability;
proving minimum index also requires exhausting or soundly pruning all remaining
work. Assembly fragments can share reusable structures even across disconnected
components, so additive decomposition needs its own proof.

Existing CPU code already has task donation, work stealing, shared incumbents
and transposition caches. The research opportunity concerns GPU representations,
processing granularity and batching. The fragment-enumeration DAG is different
from the relaxed optimization decision diagrams considered in GPU-08.

## How to close an investigation

Every issue should produce a reproducible artifact and a proceed, revise or
stop recommendation. A demonstrated bottleneck, failed cost model or negative
benchmark is a valid result. Exploration need not meet production promotion
counts, but any promotion claim must follow the current
[benchmark protocol](../../benchmarks/README.md).

Record the code revision, GPU and CPU models, memory capacities, compiler and
driver versions, concurrency settings, cases, limits and repetitions. Separate
kernel/replay measurements from full solver measurements. Complete solve time
includes input preparation, enumeration, transfers, optimization and requested
pathway reconstruction. Report changed search work as well as time; a different
incumbent discovery order can change the size of the explored tree.

Use current serial/OpenMP results as controls, retain exactness and completion
status, and validate pathways when enabled. Use the
[independent molecular oracle](../2026-09-28-branch-and-bound/molecular_oracle.py)
within its supported small-instance scope. Do not claim that small-instance
agreement proves every production bound or scheduling rule correct.

## Literature map

These sources motivate the drafts; their performance results on other problems
do not establish a speedup for molecular assembly.

| Source | Reusable contribution | Related issues |
| --- | --- | --- |
| [Gmys et al., GPU B&B with Integer–Vector–Matrix, 2016](https://doi.org/10.1016/j.parco.2016.01.008) | Compact GPU-resident exploration state; permutation-specific representation | 02, 06 |
| [Gmys et al., hierarchical work stealing, 2017](https://doi.org/10.1002/cpe.4019) | Balancing within and between GPUs and CPU workers | 06, 10 |
| [Gmys, exact flowshop on GPU supercomputers, 2022](https://doi.org/10.1287/ijoc.2022.1193), [PBB code](https://github.com/jangmys/pbb) | Warp-based explorers and distributed GPU search | 02, 06, 10 |
| [Yamout et al., vertex cover, 2022](https://arxiv.org/abs/2204.10402) | Local DFS stacks, global worklist, compact graph state | 04, 06 |
| [Almasri et al., maximal clique enumeration, PACT 2023](https://arxiv.org/html/2212.01473v3) | Idle-worker lists and demand-driven subtree donation; this is enumeration rather than maximum-clique optimization | 06 |
| [Quer et al., maximum common subgraph, 2020](https://www.mdpi.com/2079-3197/8/2/48) | Compact iterative graph search and complementary CPU/GPU execution | 02, 05, 06 |
| [Helbecque et al., portable GPU B&B, 2025](https://doi.org/10.1002/cpe.70321), [code](https://github.com/Guillaume-Helbecque/GPU-accelerated-tree-search-Chapel) | CPU work pools and batched GPU evaluation | 01, 03, 05, 10 |
| [Talbot, GPU constraint programming, AAAI 2026](https://ojs.aaai.org/index.php/AAAI/article/download/38448/42410), [Turbo code](https://github.com/ptal/turbo/tree/aaai2026) | GPU propagation and on-demand subproblem generation through a CP formulation | 09 |
| [Tardivo, Michel and van Hoeve, relaxed decision diagrams, CP 2026](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CP.2026.53) | GPU-generated relaxed diagrams for bounds | 08 |
| [Amro et al., component-aware vertex cover, TPDS 2026](https://arxiv.org/abs/2512.18334) | Compact states and dynamically scheduled component-result joins; assembly independence is unproven | 04, 06 |

Each issue repeats the context, constraints and sources needed to explore it
without the original conversation. Dependencies identify supporting work rather
than silently assuming it has been completed.
