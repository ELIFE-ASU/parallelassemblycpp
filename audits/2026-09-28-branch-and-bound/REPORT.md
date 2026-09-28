**Branch-and-bound audit and literature assessment — 28 September 2026**

Audited algorithm revision: `815e5e6` on `main`. Scope: molecular enumeration,
matching, bounds, canonical-state dominance, parallel task transfer, and the
separate string search. Primary literature was checked through the audit date,
including 2026 publications and preprints. This audit adds evidence and
recommendations; it does not change the production solver.

The solver already implements a strong specialized exact search. No incorrect
assembly index or demonstrably unsound pruning rule was found in the inspected
code or the checks below. This is bounded evidence, not a formal proof or a
complete concurrency audit. The most promising next molecular experiment is a
cheap bond-type bound alongside the existing size bounds. For strings, a
witness-producing Re-Pair warm start is a distinct, practical opportunity.

| Priority | Action | Assessment |
| --- | --- | --- |
| First | Maintain an independent molecular oracle and test bounds against exact residual optima | Needed to validate further pruning changes |
| First molecular experiment | Add a bond-type/vector savings bound | Small adaptation with an independent admissibility argument |
| First string experiment | Produce a feasible grammar-based incumbent before exact search | Supported by recent string-assembly literature; retain a valid pathway |
| Before search-order changes | Add parallel incumbent histories and full prune/expansion counters | Current counters cannot separate search quality from execution cost |
| Later | Cache proved residual bounds or expansion thresholds | Potentially useful, but requires proof-state lifecycle changes |
| Conditional | Bound local cache retention; strengthen overlap-aware packing | Profile memory or loose bounds before implementation |
| Defer | GPU decision-diagram rewrite or learned branching | Substantial new machinery without assembly-specific evidence |

The algorithm minimizes `A = N - 1 - D`, where `N` is the original edge count
and `D` is accumulated duplicate savings. Removing one of two edge-disjoint
isomorphic fragments of size `k` increases `D` by `k - 1`. A bound `B` on
additional savings gives `LB = N - 1 - D - B`. Pruning when `LB >= incumbent`
is correct when only one optimum is required. A stronger bound on savings is
**smaller**; a stronger lower bound on the index is **larger**.

This orientation matters when combining new formulas. Independent bounds on
the same remaining savings combine with `min`, or their assembly-index lower
bounds combine with `max`. Their values must not simply be added.

The current implementation already has the following important properties:

- The first canonical fragment ID is preserved in the state key; only the
  remaining IDs are sorted. This preserves the descendant enumeration cutoff
  encoded by the first fragment. Sorting the entire key would require a new
  proof. See [state key](/home/louie/skunkworks/parallelassemblycpp/src/assemblyState.h:250)
  and [ordinal restriction](/home/louie/skunkworks/parallelassemblycpp/src/improvedBnB.h:162).
- The transposition tables compare complete keys after hashing. Larger
  `sumDupBonds` scores dominate smaller scores, and a strict improvement reopens
  the state. Hash collisions alone do not prune. See
  [local table](/home/louie/skunkworks/parallelassemblycpp/src/assemblyTranspositionTable.h:68).
- Class, pair and post-fragmentation bounds avoid expensive work before
  canonicalization. `maxDupBondsPrefix` retains a maximum across possible
  future duplicate sizes. See
  [prefix bound](/home/louie/skunkworks/parallelassemblycpp/src/assemblyState.h:116)
  and [child bound](/home/louie/skunkworks/parallelassemblycpp/src/improvedBnB.h:1543).
- Shared incumbents decrease monotonically; observing an older value causes
  extra work rather than unsafe pruning. Matching loops already refresh the
  value at class boundaries and every 64 bound evaluations.
- Donated tasks are queued before the producer inserts them in the
  transposition table. The receiver therefore does not discard a task because
  the producer prematurely recorded it as explored. Transfer storage also
  avoids moving masks owned by another thread's arena. See
  [donation](/home/louie/skunkworks/parallelassemblycpp/src/improvedBnB.h:1562).
- Adaptive root leases, local deques, stealing, deeper task donation and
  measured transfer amortization are already present. They are not missing
  features to introduce under a new name.

The main assurance gap is specific: the molecular bounds lack a permanent
independent exhaustive reference solver and direct admissibility tests.
Regression agreement and serial/parallel parity can preserve a shared pruning
error. String mode already has an independent reference search in
[stringAssemblyTester.cpp](/home/louie/skunkworks/parallelassemblycpp/unitTests/stringAssemblyTester.cpp:110).
The standalone molecular oracle retained beside this report is a starting
point. It checks the generic scalar and proposed vector bounds against exact
remaining savings. Extend it to cover every production pruning route,
including ordinal-restricted states, same-parent versus different-parent
matches, split residuals, and pruning exactly at the incumbent. Its small
graphs alone do not exercise the production pair-bound threshold of 27 edges.
The additional factorized and path cases below reach that threshold, but do
not replace direct admissibility checks for every targeted/pair route.

**Recommended molecular experiment.** The January 2026 JOSS publication of
Vimal et al.'s `assembly-theory` provides a useful independent implementation.
Its current source contains `VecSimple` and `VecSmallFrags` bounds using bond
types and small fragments. The reported Rust/v5 comparison uses different
thread counts, so it does not establish the speedup these bounds would produce
in this repository. See the [published paper](https://joss.theoj.org/papers/10.21105/joss.09318)
and [authors' bound implementation](https://raw.githubusercontent.com/DaymudeLab/assembly-theory/main/src/bounds.rs).
The latter is a moving `main` URL, inspected on the audit date; pin a revision
before implementing a port.

Start with the simpler bound. Let `S` be the sum of edges in the current state's
active fragments, `Z` the number of distinct edge types across those fragments,
and `M >= 2` a conservative upper bound on every future removable duplicate's
size. An edge type is `(bond label, unordered pair of endpoint labels)`.
Define:

```text
Q    = S - Z
Bvec = Q - ceil(Q / M)
LB   = N - 1 - D - min(Bexisting, Bvec)
```

The following is this audit's admissibility argument. Let `R` be the number
of edges discarded as duplicate copies over a continuation, and `t` its number
of removals. At least one representative edge of each initial type remains in
a retained copy or an irreducible remnant, so `R <= S - Z`. Since each removed
copy has at most `M` edges, `t >= ceil(R/M)`. Thus additional savings
`R - t <= Q - ceil(Q/M)`. The last expression is monotone in integer `Q`.
Singleton remnants discarded from active storage do not invalidate the
argument: they cannot subsequently be removed as reusable duplicates.

Implementation should precompute masks for edge types, count their intersection
with the current fragment union, and use integer ceiling division. Do not use
the original molecule's `Z` after fragmentation; that can overstate the types
that remain. Cache the metadata only after the fragmentation is valid. Use a
separate exact type ID scheme so collecting this bound cannot perturb canonical
ordinals. Skip the extra calculation when an existing bound already prunes.
Begin with state-level pruning, then consider class-level integration only
after validating its `M`. The insertion points are
[assembly state bounds](/home/louie/skunkworks/parallelassemblycpp/src/assemblyState.h:85)
and [post-fragmentation bound](/home/louie/skunkworks/parallelassemblycpp/src/improvedBnB.h:357).
The audit checked this formula against exact residual savings in **45,014
oracle state instances** with legal duplicates, with no underestimates. It
was strictly tighter than the generic scalar bound in **2,229 cases (4.95%)**,
by at most two saved joins. This is not a comparison against the stronger
production targeted/pair bounds, and does not measure additional production
prunes or wall time. Defer the more complicated small-fragment formula until
this simpler experiment is understood.

**Recommended string experiment.** Siebert, Chowdhury, Slocombe and Walker's
June 2026 preprint develops grammar-based assembly bounds. Bieniawski's 2026
preprint evaluates Re-Pair and introduces a tie-resolving branch-and-bound
variant that produces upper bounds. These support trying a bounded Re-Pair
pass as an incumbent generator. Neither makes a heuristic grammar an exact
assembly answer. See [Assembly Spaces](https://arxiv.org/abs/2606.15499)
and [Assembly Theory and the Smallest Grammar Problem](https://arxiv.org/abs/2608.19228).

This code already has an LZ-style residual bound at
[stringAssembly.h:705](/home/louie/skunkworks/parallelassemblycpp/src/stringAssembly.h:705);
adding an undifferentiated “compression bound” would miss that fact. The new
piece would be a feasible initial grammar and its assembly witness, replacing
the initial `n - 1` incumbent when better. Translate rule lengths to binary
join counts and replay the witness. Retain it when setting the incumbent:
the existing `>=` pruning can eliminate every equally good search branch, so
an external numeric bound alone cannot supply the required pathway. Ordinary
oriented concatenation remains feasible when reversal is allowed, although it
may give a weaker warm start. Respect Unicode scalar semantics and cancellation.
Bound tie exploration by work, and measure its cost separately.

The June preprint also suggests constructing molecular upper bounds from
trails. That is a later experiment for chain-rich molecules, with explicit
graph reconstruction. String compression of an arbitrary SMILES serialization
is not a justified graph bound, and formal hyperedge-grammar correspondence
does not by itself provide an efficient exact molecular solver.

**Useful but lower-priority bound work.** A certified table of shortest
addition-chain lengths can strengthen the `ceilLog2(k)` construction cost.
However, Seet et al.'s own algorithm paper explicitly discusses this option
and reports little timing benefit, potentially offset by cache cost. It is a
controlled benchmark candidate, not a newly discovered improvement. The
published version is December 2025; see
[JCIM](https://doi.org/10.1021/acs.jcim.5c01964) and the
[preprint, printed page 15](https://arxiv.org/pdf/2410.09100).
The December 2025 [Assembly Addition Chains preprint](https://arxiv.org/abs/2512.18030)
provides further formal context for size projections.

Use only certified exact lengths or proved lower bounds, never the length of
a merely found chain. Preserve the maximum over candidate sizes: shortest
chain length is not monotone. This audit's reasoning supports subtracting the
stronger construction cost once, including the separate extra operation for
a second distinct equal-sized fragment; production adoption still needs
route-by-route bound checks and timing.

A more substantial research experiment is to tighten the current eligible-edge
packing relaxation using occurrence conflicts. At a fixed size, occurrences
sharing a physical edge conflict. A valid clique cover of that conflict graph
upper-bounds how many occurrences can coexist; a greedy feasible packing gives
the opposite bound and is unsafe for this purpose. The mapping is this audit's
proposal, inspired by conflict-based optimization, including
[Dai and Chen, March 2026](https://doi.org/10.1007/s12532-026-00307-4).
Start on small, overlap-heavy classes near the pruning threshold. The cover
must include every occurrence permitted by the bound's class/ordinal route;
building a cover for one class and applying it to a union of classes is unsafe.
Prove the integration with retained copies before enabling pruning. Constructing
conflict graphs everywhere may cost more than it saves.

**Caching and search control.** Coppé, Gillard and Schaus's 2024 work caches
expansion thresholds that detect suboptimal revisits beyond ordinary score
dominance. That distinction is relevant here: the current table stores only
the best score reaching each key. See their
[paper, Sections 4.1–4.4](https://arxiv.org/html/2211.13118v3).

A possible adaptation stores an admissible bound `R(key)` on remaining
savings and tests `N - 1 - D - R(key) >= incumbent`. Exact completion values
are stronger still. However, current entries are inserted before recursion,
and donated children may be pending when a parent returns. A visited score is
not an exact residual value. Reusing completed-subtree information requires
separate proof status, correct aggregation of every child, and handling
cancellation and work still queued. Preserve the first-fragment cutoff in
the key. Start serial and do not replace the current safe in-progress
dominance mechanism with a completion assumption.

The existing telemetry records matching visits and cache outcomes, but lacks
full state expansion/prune counts and a parallel incumbent trajectory. Its
matching-refresh prune counters attribute only immediate extra savings caused
by a refresh; they are not totals for all pruning. Moreover,
`--write-intermediate-mas` forces serial execution at
[main.cpp:449](/home/louie/skunkworks/parallelassemblycpp/src/main.cpp:449).
Add optional per-worker counters and strict-improvement events using a steady
clock, keeping telemetry out of performance binaries. Separate time to a good
incumbent from time to prove optimality, and distinguish state, class, pair
and block counts.

There is already local evidence against assuming generic greedy diving helps.
Historical six-pair experiments on three long cases gave aggregate wall-time
speedups of approximately **0.999× at 8 workers** and **1.002× at 28 workers**.
These are historical, essentially neutral results, not current-branch
benchmarks. See [methods](/home/louie/skunkworks/parallelassemblycpp/build/incumbent-results/report-methods.md),
[8-worker results](/home/louie/skunkworks/parallelassemblycpp/build/incumbent-results/final-long-six-pairs-8.json)
and [28-worker results](/home/louie/skunkworks/parallelassemblycpp/build/incumbent-results/final-long-six-pairs-28.json).
The retained multi-start Paclitaxel experiment found a best greedy value of
25 versus the benchmark target of 23. Better branch guidance should therefore
target a demonstrated failure mode, with separate ablations for search order
and bound strength. A 2026 theoretical result also shows that accurate local
branching scores need not yield small overall trees; this is a caution about
evaluation, not evidence of a defect here. See
[Cheng and Basu](https://arxiv.org/abs/2601.23249).

Worker-local cache memory is another practical target. Every worker owns an
unbounded local table with monotonic key storage, and L1 insertion precedes
shared L2 lookup. The shared-cache byte budget does not bound total process
memory. See [local storage](/home/louie/skunkworks/parallelassemblycpp/src/assemblyTranspositionTable.h:130)
and [L1/L2 sequence](/home/louie/skunkworks/parallelassemblycpp/src/searchContext.h:280).
Measure actual retained bytes before trying capped admission or generations.
Forgetting a cache entry can preserve exactness by re-expanding it; rejecting
an insertion or matching only a fingerprint must never imply domination.

The July 2026 [GPU relaxed-decision-diagram paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CP.2026.53)
is relevant longer term, but its layered optimization diagrams are different
from this code's fragment-enumeration DAG. Adopting it requires an admissible
state-merging relaxation and suitable batched workloads. Its published gains
on other problems do not predict a GPU gain here.

For independent certification, the 2025
[Assembly in Directed Hypergraphs](https://arxiv.org/abs/2505.22826)
offers an ILP formulation worth considering as a second small-instance oracle.
Its assembly space and joins must match this program, and enumerating that
space can itself be expensive. The March 2026
[Certified Branch-and-Bound MaxSAT Solving](https://doi.org/10.1609/aaai.v40i17.38449)
shows the value of checkable pruning proofs in another domain; its proof format
is not directly reusable here. A returned pathway proves achievability, not
minimality. Independent residual-bound checks are the immediate practical step.

Validation performed during this audit:

- Rebuilt `build/dev`; all **22/22 CTest checks passed**, including focused
  regression, string oracle, transposition, fragmentation and telemetry checks.
  This is the focused developer suite, not the full `ci` molecular manifest.
- Separately rebuilt and ran the **OpenMP-enabled task-transfer unit: 1/1
  passed**, including concurrent root leases and cross-worker wide masks.
- Checked **1,590 solver inputs** against independent expected results. The
  unpruned molecular oracle covers **1,474 graphs**: all simple homogeneous
  graphs on 3–5 vertices with
  at least two edges, plus 400 deterministic random labelled graphs on 4–8
  vertices with at most 11 edges. Every result agreed. The oracle uses its own
  connected-subset enumeration and permutation-based isomorphism, and no
  production bounds, DAG or ordinal restriction. The remaining inputs are
  **100 factorized 27-edge graphs** with label-disjoint gadgets and **16
  homogeneous or uniquely capped paths of 27–34 edges**, checked with an
  independent exact addition-chain search. Across the small graph/gadget
  reference searches, **62,160 state instances** were visited; the generic
  scalar and proposed vector bound were admissible in all **45,014** states
  containing a legal duplicate.

No full MPI parity, sanitizer run, formal proof or new performance promotion
study was performed. The standalone oracle and its machine-readable results
are retained in this directory: [oracle](/home/louie/skunkworks/parallelassemblycpp/audits/2026-09-28-branch-and-bound/molecular_oracle.py),
[library runner](/home/louie/skunkworks/parallelassemblycpp/audits/2026-09-28-branch-and-bound/oracle_runner.cpp),
and [results](/home/louie/skunkworks/parallelassemblycpp/audits/2026-09-28-branch-and-bound/validation.json).
From the repository root, reproduce the independent check after building the
developer library:

```bash
g++ -std=c++20 -O2 -Isrc \
  audits/2026-09-28-branch-and-bound/oracle_runner.cpp \
  build/dev/libparallelassemblycpp.a -o /tmp/bnb-audit-runner
python audits/2026-09-28-branch-and-bound/molecular_oracle.py \
  --solver /tmp/bnb-audit-runner --output /tmp/bnb-audit-validation.json
```

The first implementation experiment should
compare the bond-type bound against a frozen baseline, record actual pruning
and bound cost, and preserve both indices and pathway validity. Run correctness
checks before isolated paired timings; a reduction in node count alone is not
a speedup. The repository's existing
[promotion protocol](/home/louie/skunkworks/parallelassemblycpp/benchmarks/README.md:338)
requires 100 rounds for quick/full, 6 for profile and 30 for scaling.
