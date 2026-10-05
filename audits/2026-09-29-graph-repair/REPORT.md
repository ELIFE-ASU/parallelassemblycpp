# Graph-pair assembly upper-bound benchmark

> Historical measurements from 29 September 2026. Corpus and validation counts
> describe the recorded run. The current CLI uses `--algorithm=re-pair`; the
> `--upper-bound=graph-repair` spelling below remains a legacy alias. See the
> [CLI reference](../../docs/cli.md) and [audit index](../README.md).

The GraphRePair-inspired bound improves the solver's initial bound on **1,054
of 1,067 reviewed molecular fixtures**, ties on 13, and never worsens it. It is
also tighter than the fixed simple-trail string RePair baseline on **872
fixtures**, ties on 161, and is worse on 34. It is a useful upper-bound heuristic,
but does not dominate string RePair on every molecule or establish exact indices.

## Results

All comparisons use the existing default convention: explicit hydrogen atoms
are removed, atom and bond labels are retained, and disconnected components
include the repository's extra component-joining cost. Smaller bounds are better.

| Corpus | Cases | Mean initial bound | Mean trail RePair | Mean graph-pair bound | Mean graph-pair gap to reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Reviewed regression fixtures | 1,067 | 15.643 | 11.738 | 9.976 | 0.933 |
| Additional provisional benchmarks | 22 | 67.227 | 29.727 | 20.045 | 2.000 |
| Deduplicated union | 1,089 | 16.685 | 12.101 | 10.179 | 0.954 |

For the reviewed corpus, the initial incumbent is the hydrogen-filtered bond
count minus one. The mean bound drops by 5.667 joins (36.2%). The graph-pair bound
equals the reviewed assembly index on 396 fixtures; its largest gap is six joins.
No graph-pair or trail bound falls below a reference value. The trail baseline's
mean gap on the reviewed corpus is 2.694 joins.

| Molecule | Initial bound | Trail RePair | Graph-pair | Reference MA | Reference status |
| --- | ---: | ---: | ---: | ---: | --- |
| Icosane | 18 | 6 | 6 | 6 | reviewed |
| Bisphenylmaleimide | 29 | 18 | 10 | 10 | reviewed |
| Sucrose | 23 | 16 | 10 | 8 | reviewed |
| Ketoconazole | 39 | 28 | 26 | 22 | reviewed |
| Ceftiofur | 36 | 28 | 29 | 26 | reviewed |
| Erythromycin | 52 | 32 | 27 | 20 | provisional |
| Paclitaxel | 67 | 36 | 30 | 23 | provisional |

Ceftiofur illustrates that a locally greedy graph choice can lose to a fixed
path decomposition. Of the 34 reviewed losses to trail RePair, 33 are by one
join and one is by two. Taking the minimum of both valid heuristics would lower
the reviewed mean to 9.943 and match 414 reference indices, although the production
option currently computes the graph-pair bound alone.

The median graph-pair calculation takes approximately **25 microseconds** on
this machine; the sum of the per-case median times is approximately **35 ms** for
the 1,089-case corpus. Each case has three samples. These measurements time the
C++ bound calculation, exclude parsing, certificate serialization and validation,
and use one process per complete corpus pass. They are local indicative timings,
not a performance promotion result. Python trail timings are retained separately;
they do not support a fair C++ versus Python speed comparison.

## Algorithms and validity

The starting point is Maneth and Peternek's
[Grammar-Based Graph Compression](https://arxiv.org/abs/1704.05254), considered for
the molecular upper-bound approach discussed in
[Assembly Spaces: Formal Definitions and Fast Methods for Approximating Assembly Indices](https://arxiv.org/abs/2606.15499).
The graph algorithm starts from individual labeled bonds and repeatedly groups
incident fragment pairs by exact isomorphism of their expanded graphs, using the
solver's existing exact graph canonicalization. It chooses
a class with the largest immediate join saving using a deterministic greedy set
of nonoverlapping occurrences. A new fragment costs one binary construction;
reuse of an existing fragment costs zero. Full expanded graph equivalence allows
cycles and joins sharing several vertices, and retains all attachment vertices.
This is an assembly adaptation inspired by GraphRePair, not an implementation or
runtime claim for the original graph compressor.

If there are `p` newly constructed fragments and `r` residual occurrences,
the default convention gives the bound `p + r - 1` for a nonempty graph. Residual
occurrences partition the physical bonds. Within each connected component they
can be joined along shared vertices; the default disconnected convention then
adds one join per extra component.

The comparison baseline partitions edges into simple trails deterministically:
start from the lowest-index unused edge, extend right and then left with the
lowest-index available edge, and stop an extension before revisiting a vertex.
Each whole trail is oriented to its lexicographically smaller labeled edge
sequence. Ordinary directed string RePair then repeatedly replaces the most
frequent nonoverlapping digram across those sequences. Ties are lexicographic.
Rules do not match reversed strings. Its bound is the number of binary dictionary
rules plus the number of remaining tokens minus one.

Cutting trails at repeated vertices matters: the same label sequence can describe
an open path and a cycle, so unrestricted walk strings do not alone establish
isomorphic graph fragments. This baseline's simple paths preserve the required
topology. It is a reproducible baseline rather than a reproduction of the paper's
trail optimization or a best-possible string RePair bound. Both heuristics can
depend on atom and edge ordering.

For all 1,089 inputs, **both algorithms' construction certificates** were checked
with the independent Python checker in `unitTests/graphRepairTester.py`. It
verifies labeled graph isomorphism, each binary rule, disjoint physical edge use,
coverage of the original graph, and the final residual joins. Thus 2,178
constructions were replayed. These checks establish valid upper bounds for these
inputs independently of whether their reference indices are optimal. The
separate exhaustive and randomized small-graph validation is recorded in
`validation.json` in this directory.

The small-graph audit completed **1,895 comparisons with an independent exact
oracle** and **1,904 certificate replays**. It covered all 1,074 simple graphs
with at least two edges on three through five vertices, 400 random labeled
graphs, 400 vertex permutations, 16 addition-chain path cases, shared-vertex and disconnected cases,
and three label/cycle cases. The audit also checked disconnected compensation
and rejection of six deliberately invalid certificates. Its largest gap to an
exact index was two joins.

The focused CMake suite passed all 24 tests, including the new graph-pair CLI,
library API, and construction checks alongside existing regressions. Separate
OpenMP, two-rank MPI, and hybrid smoke checks passed for the upper-bound mode.

## Reproduction and scope

Build configuration used here: CMake Release, GCC `/usr/bin/c++`, `-O3 -DNDEBUG`,
portable architecture (`PARALLELASSEMBLYCPP_X86_64_V3=OFF`). The benchmark stores
input, executable, and relevant source SHA-256 fingerprints in `bounds.json`.

```bash
cmake -S . -B build/graph-repair -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/graph-repair \
  --target parallelassemblycpp_graph_repair_probe
python3 benchmarks/graph_repair_benchmark.py \
  --executable build/graph-repair/parallelassemblycpp_graph_repair_probe \
  --corpus all --repeats 3 \
  --json-output audits/2026-09-29-graph-repair/bounds.json \
  --csv-output audits/2026-09-29-graph-repair/bounds.csv
```

`bounds.csv` contains one row per input, and `bounds.json` retains all samples,
fingerprints, and summaries. The runner reads `unitTests/regression_cases.tsv`
and `benchmarks/cases.tsv`, deduplicating by source path; the 15 reviewed full-suite
benchmarks are already included in the 1,067 regression cases. The other 22
benchmark cases remain explicitly provisional.

The new `--upper-bound=graph-repair` option is opt-in; default exact calculations
retain their existing search and behavior. This experiment assesses bound quality
and heuristic cost. It does not measure
exact-search acceleration, reproduce the linked paper's optimized trail method,
or replace the existing exact search. Corpus comparisons use the maintained
reference values; this run does not independently re-prove all 1,067 optima.
