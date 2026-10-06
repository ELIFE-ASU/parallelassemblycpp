# Addition-chain lower bounds

[Algorithm and provenance](algorithm.md) · [Development](development.md)

Addition chains relax the construction problem by forgetting topology or string
order. The exact searches retain their structural duplicate matching and
pathway construction. A bound may discard a search branch only when every
structural construction in that branch costs at least the reported bound.

## Scalar and vector relaxations

Let `l(n)` be the shortest scalar addition chain from `1` to `n`, counting one
operation for the sum of two previously available entries. Reuse is free. Any
construction of an `n`-edge fragment induces such a chain by replacing each
fragment with its edge count. Consequently its assembly index is at least
`l(n)`, which is at least `ceil(log2(n))`. For example, `l(15) = 5`, whereas
`ceil(log2(15)) = 4`.

For a labelled graph, associate each fragment with the counts of its primitive
bond types. A bond type comprises its bond label and its unordered pair of
endpoint atom labels. Binary joins add these vectors. Atom counts would not
work: graph joins identify vertices. For a string, use counts of Unicode scalar
values; reversing a string does not change this vector.

Let `L(v)` be the shortest vector addition chain from the free unit vectors to
the target count vector `v`. Every structural pathway projects to a vector
chain, so `L(v)` is a lower bound on structural assembly. Equal vectors do not
imply equal structures, and a vector chain does not provide a structural
pathway. Permuting vector coordinates does preserve the relaxed problem, so a
cache may sort the nonzero counts in its keys.

There is also a cheap support bound. With `d` distinct primitive types, at
least `d - 1` additions are necessary. Equality leaves exactly a tree with one
leaf of each type, so its vector is all ones. If any count exceeds one, at
least `d` additions are necessary. The implementation takes the maximum of
this bound, the total-count scalar bound, and each coordinate's scalar bound.

Only a certified lower bound on `l` or `L` may prune structural search. A chain
found by a heuristic is a constructive upper bound on the relaxed optimum; its
length is not automatically a structural lower bound. A search stopped after
exhausting depths below `d` may safely report `d`, even if it has not yet found
an optimal chain. An interrupted, partially searched depth does not justify
raising that certificate.

A caller may supply an incumbent as a sufficient proof threshold. Refinement
stops as soon as the certified bound reaches it; finding a relaxed construction
after that point cannot strengthen the structural optimality proof. Cached
certificates remain valid across thresholds, although a cache hit may return a
weaker bound than another search budget or threshold could establish.

The implementation uses an immutable table of exact scalar lengths through
256 and the logarithmic bound above that range. Vector bounds first take the
maximum of the total-count scalar bound, each individual coordinate's scalar
bound, and the support bound above. Targets with total count at most 16 can
additionally use bounded exhaustive search, with a default budget of 4,000 pair
evaluations per cache miss. Exhausting that budget leaves the last certified
bound intact; it does not claim to have solved that vector-chain problem. The
vector cache belongs to one search or worker and holds at most 256 targets.

These bounds are enabled by default. Configuring
`-DPARALLELASSEMBLYCPP_VECTOR_CHAIN_BOUNDS=OFF` retains the scalar improvements
while disabling vector refinement and composition penalties. This provides an
ablation and a fallback for workloads where vector analysis costs more than it
saves.

## Whole-input bound

Compute the graph count vector after the requested hydrogen removal and before
the exact search removes bond types occurring once. Those unique types still
contribute to the construction cost. The string counterpart counts the
original decoded symbols, before preprocessing removes nonreusable regions.

The bound is a floor on the complete assembly index. Combine it with other
lower bounds using `max`, and finish once a witnessed upper bound equals it.
Subtracting it from a residual count, or adding it to a branch score, would
charge some operations twice.

The graph search internally uses its uncompensated index. Optional disconnected
component compensation applies to the final index, after the internal bound
comparison. The same fixed adjustment must apply if a bound is ever displayed
in those output units. Empty-input conventions require separate handling; a
nonempty vector-chain floor must not change the legacy empty-graph result.

Graph root enumeration still honors its existing explicit enumeration cap.
Parallel graph workers skip calculation when the incumbent reaches the floor,
but continue draining root leases so the scheduler and MPI completion protocol
remain unchanged. String search can certify its complete Re-Pair seed before
exact enumeration begins.

## Stronger scalar penalties in graph savings bounds

The recursive graph search tracks an incumbent construction cost
`U = N - S - 1`, where `N` is the original bond count and `S` is the savings
already selected. A branch bound estimates the largest additional savings in
the remaining fragment forest.

Consider a completion whose largest reused fragment has `k` edges. Cut the
construction of each remaining target fragment at occurrences of reused
fragments; the resulting tokens have at most `k` edges. A target of size `n`
therefore requires at least `ceil(n / k)` tokens. The joins above these tokens
form a construction skeleton. The dictionary constructing the reused tokens
contains a `k`-edge object and needs at least `l(k)` operations. Those dictionary
operations are disjoint from the skeleton joins.

Compared with separately joining primitive edges, the total remaining savings
are therefore at most

```text
sum over target fragments [n - ceil(n / k)] - l(k).
```

Take the maximum over permitted largest duplicate sizes, including the case
with no further savings. This proves that each existing `ceil(log2(k))`
dictionary penalty can be strengthened to a certified scalar-chain bound.

The eligible-edge variant has the same proof. If only `e` edges may occur in
size-`k` tokens, at most `g = floor(e / k)` such tokens are possible; all other
tokens have size at most `k - 1`. Using the maximum possible number of large
tokens minimizes the required token count, giving the savings bound

```text
sum [n - g - ceil((n - g*k) / (k - 1))] - l(k).
```

This is the expression evaluated by `fixedSizeDupBondsForFragment`, followed
by the dictionary penalty. Eligibility masks change the packing argument, not
the required construction of the dictionary.

Pair filters evaluate the same generic expression before residual parent masks
are split into connected components. The function
`n - ceil(n / k)` is superadditive, so this unsplit calculation overestimates
possible savings. Replacing its size-dependent dictionary penalty does not
change that argument.

## Retained-fragment vector penalty

After selecting a duplicate pair, `postFragmentationCutoff` keeps one occurrence
as its first fragment `P`, of size `K`. This required fragment permits two
additional, narrowly scoped vector penalties:

1. If the next largest reused fragments belong to the same canonical class as
   `P`, their dictionary must construct `P`. Subtract a certified bound on
   `L(counts(P))` from that route's packing estimate.
2. If a different canonical class `Q` of the same size is also required, the
   dictionary must construct both `P` and `Q`. Neither can be an ancestor of the
   other because every join strictly increases edge count. A construction of
   `P` therefore excludes the final operation constructing `Q`, proving the
   penalty `L(counts(P)) + 1`. This remains valid when the two classes have
   identical count vectors.

These routes treat the retained first fragment as a dictionary target, even
when it has no further duplicate occurrence. The cutoff explicitly gives it
all `K` eligible edges. Other, smaller-size routes retain their scalar `l(k)`
penalties.

Applying the vector bound for `P` to the smaller-size routes would double-count
some skeleton operations. For example, a homogeneous path `P` with eight edges
can be constructed in three joins. With largest future duplicate size four,
its generic remaining savings are `8 - ceil(8/4) - l(4) = 4`. Incorrectly using
`l(8) = 3` as the dictionary penalty lowers those savings to three and could
prune an optimal construction.

Likewise, independent lower bounds for several remaining fragments must not be
summed: their constructions can share intermediate fragments. A future stronger
relaxation would have to represent all required targets in one shared vector
chain.

## Regression coverage

`graphAdditionChainTester.py` drives the public library through a batch probe
and compares exact results with the independent physical-mask graph oracle.
The corpus covers all tiny simple graph topologies, labelled and permuted
graphs, unique-type preprocessing, explicit hydrogen removal, isolated atoms,
disconnected compensation, and repeated input domains in one process.

Larger homogeneous paths have independently computed scalar-chain answers.
Repeated motifs containing distinct primitive types supply known vector-bound
gaps: a motif with `d` unique edge types repeated twice has index `d`, requiring
`d - 1` joins to combine the types and one to duplicate the motif. These cases
exercise improvements beyond the scalar logarithmic floor without relying on
the production chain solver for their expected answers.

Longer paths with nested, reversed, and overlapping labelled motifs are checked
by a separate bound-free substring search. Uniform atom labels make substring
reversal exactly the graph-isomorphism rule for these paths. This supplies
independent answers for deeper retained-fragment cases without the factorial
vertex-permutation cost of the general graph oracle.

The full graph corpus contains 1,598 inputs and performs 6,392 exact calculations
across hydrogen-removal and disconnected-compensation settings. Scalar table
entries through 256 are independently checked, and tiny vector targets are
compared with an independent breadth-first search. String tests check exact
indices, pathway replay, Unicode, cancellation, target status, and parallel
execution. The existing 1,067-case regression corpus and serial, OpenMP, MPI,
and hybrid checks also pass with the bounds enabled.

## Measured results

The following measurements compare the original commit
`67fd2460b9576ec75c0cfded7984a32d42d55217`, the scalar-only configuration, and the
default scalar-plus-vector configuration. Builds used GCC 15.2, Release mode,
and `x86-64-v3` on an Intel Core i7-14700KF. Each benchmark was pinned to CPU 0
and run without competing builds or tests, with one warmup and six measured
rounds. Graph runs alternate executable order; string runs rotate all six
three-executable orders. Speedups are medians of paired round ratios, so values
below one indicate a slowdown.

| Workload | Comparison | Metric | Speedup |
| --- | --- | --- | ---: |
| 15 full-suite graph fixtures | Original / default | Suite wall time | 1.014× |
| Same graph fixtures | Scalar / default | Suite wall time | 1.006× |
| Six homogeneous paths, 63–129 edges | Original / default | Suite wall time | 4.362× |
| Unary string, 64 symbols | Original / default | Algorithm time | 384.7× |
| Same unary string | Original / default | Process wall time | 8.58× |
| String composition `(6, 2)` | Scalar / default | Algorithm time | 1.91× |
| String composition `(12, 4)` | Scalar / default | Algorithm time | 4.22× |
| Five distinct symbols repeated twice | Scalar / default | Algorithm time | 2.11× |
| Random binary string, 16 symbols | Scalar / default | Algorithm time | 0.74× |

The path and unary improvements are largely due to the scalar bounds and early
optimality proof. The string composition comparisons isolate additional vector
benefits. General molecules show only a small aggregate improvement; individual
cases vary. Vector search adds overhead when composition is too weak to certify
optimality: the random binary case rises from about 53 to 74 microseconds. Such
small algorithm-time differences are often hidden by process startup, which is
why both timing scopes are retained in the reports. These results do not predict
a universal speedup or establish performance on the separate profiling corpus.

All completed string measurements agreed on exact indices. The original solver
timed out on the 256-symbol unary warmup at three seconds; that configuration
was not repeatedly retried and no speedup was assigned to the censored result.
Timing runs omit pathway output; regression tests independently validate
pathway witnesses. Parallel tests establish correctness, not parallel speedup.

Reproduction commands are in the [benchmark guide](../benchmarks/README.md#addition-chain-ablation).
The local experiment reports include raw samples and executable/input hashes:

- `build/addition-chains-final-full.json`
- `build/addition-chains-final-ablation.json`
- `build/addition-chains-final-paths.json`
- `build/addition-chains-strings-final.json`

These generated reports are build artifacts, not source-controlled fixtures.
