# Duplication consolidation for issue #146

This records the disposition of the repository-wide duplication audit and the
local implementation stages for [issue #146](https://github.com/ELIFE-ASU/parallelassemblycpp/issues/146).
The starting revision was `c478f7b`. The objective is fewer independently
maintained implementations while preserving scientific results, reference
oracles, search order, resource ownership and useful specialisations.

## Scope

The audit covered first-party C++ sources and headers, Python utilities and
benchmark drivers, shell/Slurm scripts, CMake/CI configuration, test programs,
test helpers and package-consumer code. Generated build trees, solver outputs,
binary assets and installed external libraries were excluded. Fixture data and
golden results were treated as evidence, not implementations to deduplicate.
No first-party source area was inaccessible. Runtime validation here is on
Linux with GCC; Windows/MSVC, Clang, multi-node MPI, Sol jobs and PGO training
still require their normal platform checks.

## Completed consolidations

The finding numbers group related occurrences from the audit. Narrow adapters
remain where they preserve an existing interface or policy.

| Finding | Canonical responsibility and changed callers | Preserved differences and validation |
| --- | --- | --- |
| 1. Virtual-child bound | `src/improvedBnB.h::pairSpecificGenericBound` takes explicit fragment metadata; its matching overload forwards to it. | Integer ceiling arithmetic, same-parent and different-parent cases, scalar lower bounds and the existing non-inlined kernel. `unitTests/consolidationKernelsTester.cpp` independently constructs child edge counts and checks the result. |
| 2. Process lifecycle and solver launch | `tools/process_utils.py::run_command` serves `benchmarks/benchmark.py` and `unitTests/parallelSolverTester.py`; string parity tests retain their existing runner import. `prepare_solver_launch` and `run_solver_command` serve timed and telemetry runs. | Existing runner signatures, direct-launch defaults, environment overlays, POSIX process-group cleanup, timeout output, interruption propagation and timing boundaries. Existing benchmark/process and parallel tests cover these paths. |
| 3. MOL rejection checks | `unitTests/molfileParserTester.cpp` reuses `expectRejected` and `withLastBond`. | All malformed bond/version cases, diagnostic checks and unchanged destination checks remain. Parser suite passes. |
| 4. Graph fixtures | `unitTests/graphTestFixtures.h::makeGraph` replaces the local builders in `treeCanonTester.cpp` and `cyclicCanonTester.cpp`. Cyclic collision fixtures reuse the prism/bipartite factories. | Permutations, atom/bond labels and edge insertion order remain explicit. The shared builder retains cyclic's always-on validation and strengthens the tree fixture checks; independent isomorphism checks remain separate. |
| 5. Telemetry test plumbing | `benchmarks/test_benchmark.py` uses one serialization/parse helper for telemetry fixtures. | Mutations, malformed schema cases and expected errors remain visible at each test site. No production validator is used to construct expected scientific answers. |
| 6. Hydrogen CLI matrix | `unitTests/unitTester.py::run_cli_checks` runs the same six scenarios over the independent MOL and native fixtures. | Literal input fixtures, option order, expected atom labels and edge counts, and cross-format equality remain. The full CLI suite exercises the matrix. |
| 7. Clock subtraction | `src/clockTicks.h` supplies `difference` and `budgetTicks` to `globalPrimitives.h`, `searchTelemetry.h`, `stringAssembly.h` and `stringApi.cpp`. | Unavailable-clock handling and unsigned rollover are shared. Budget saturation and raw API/telemetry conversion remain distinct. The public API still includes decoding in its timing; internal search timing does not. Compile-time boundary checks and string-library tests cover the change. |
| 8. Scheduler counters | `src/searchContext.h::resetSearchSchedulerStatistics` serves the two reset sites in `src/main.cpp::runParallelSearch`. | Only scheduler statistics are reset; topology, pointers and graph state remain explicit. Worker telemetry is captured before teardown reset. OpenMP/MPI telemetry tests exercise the lifecycle. |
| 9. Bonded components | `src/graphRepair.h::implementation::componentCount` reuses `molGraph::disjointFragments` and subtracts isolated vertices. | The graph API includes isolates; the Re-Pair bond model excludes them. Independent component tests and graph certificate/oracle suites cover empty, isolated and bonded cases. |
| 10. Touched-atom traversal | `src/ufds.h::ufdsSplit::forEachTouchedAtom` serves the small, two-word and wide split methods. | Storage-specific accumulation stays separate. Ascending traversal, sparse overflow-word sorting, component order, extra cycle edges, reset and generation wrap remain. Boundary tests cover 63/64/65, 127/128/129 and 513 edges with atoms beyond 512. |
| 11. Matching and pairability | `src/duplicateMatching.h::duplicateSet` shares reverse traversal and `visitOccurrencePairability`; initial DAG population and later DAG expansion supply their own callbacks. | Reverse matching order, interleaved scan/expansion, multi-fragment shortcut, cancellation polls, telemetry counts, partial results and frontier cleanup remain. Initial and later DAG expansion algorithms are still independent. |
| 12. Report structure and reuse | `benchmarks/report_reader.py` owns execution identity and corpus-fingerprint parsing. `check_parallel_scaling.py::evaluate_document` returns validated samples for `summarize_sol.py::report_rows`. | Historical identity return types, whitespace policies, empty-corpus policies and wall-only versus CPU-clock requirements remain. Summary reads and validates once but retains its original sequential floating-point accumulation. Adversarial `1e16, 1, 1` samples guard that distinction. |
| 13. Manifest/hash metadata | Graph benchmark readers reuse `benchmark.load_manifest` and `file_sha256`; PGO metadata reuses `benchmark_corpus_metadata` with a path adapter; the string timing driver reuses file hashing. | Existing report field shapes, repository-relative PGO paths, graph source reconciliation, conflict detection and suite ordering remain. Invalid benchmark manifests now consistently use the canonical reader's validation. Valid repository records compare identically. |
| 14. Scaling setup and plots | `paclitaxel_scaling.py::paired_paclitaxel_arguments` and `finish_scaling_plot` also serve `hybrid_scaling.py`. | MPI launcher placement, rank/thread layouts, OpenMP controls, measured curves, labels and topology-specific report names stay with their callers. Plot/setup tests cover both drivers. |
| 15. Exact string children | `src/stringAssembly.h::Search::exploreMatching` serves recursive serial matching and parallel root workers. | Parent saved-symbol counts versus the zero root count, worker-local canonical IDs/cache/pathway, pruning, ordering and stop polling remain. Exhaustive serial/OpenMP oracles and pathway replay cover the change. |
| 16. Collective success | `src/main.cpp::allRanksSucceeded` replaces the string helper and equivalent MPI success reductions for lease validation, graph preparation, scheduler preparation, thread setup, input loading and CLI parsing. | Same collective order, communicator, `MPI_INT`/`MPI_MIN`, and coordinating thread. Broadcasts, count/sum reductions, extrema and winner selection remain explicit. |
| 17. Molecular initialisation | `src/improvedBnB.h::prepareMolecularSearchInput` serves serial search and root production; `configureSearchMaskDomains` also serves worker binding. | Owning caches clear before mask destruction/reconfiguration. Context publication, seed registries and worker setup remain separate. Repeated 0/2/65/129/64/0-width tests check cache cleanup, unique bonds and isolates. |
| 18. Owner/view scans | `src/activeWordMask.h::activeWordMaskDetail::{findFirst,findNext}` share ascending scans through each representation's own accessor. | Zero width, end sentinels, null borrowed storage and owning-mask arenas remain unchanged. Owner-to-view delegation was rejected after a slowdown; the accepted form retains specialised accessors. Existing mask tests and pinned scan comparisons passed. |
| Small confirmed helpers | `stringAssemblyTester.cpp` shares exhaustive-string enumeration and a file-cleanup guard. `benchmark.cpu_description` reuses `cpu_topology.cpu_model` after its existing preferred query. `unitTester.py` removes one duplicate `--threads=2147483647` help case. | Corpus limits/order, temporary-file lifetime and CPU-description precedence remain. The thread boundary still runs in `run_flag_matrix_checks`. |

## Retained implementations and follow-up work

These items were evaluated after the high-confidence changes. They are not
silently counted as completed extractions.

| Candidate | Disposition and evidence | Requirement before further change |
| --- | --- | --- |
| `src/improvedBnB.h::{pairBoundFiltersMatching,pairBoundFiltersFragmentPair}` | Deferred. Occurrence and block paths use different refresh cadence, telemetry categories and deliberate hot-loop call boundaries. Their common arithmetic already uses the shared bound. | Workload-specific timing and counter equivalence for the complete callback adapter, including cancellation and refreshed incumbent cases. |
| Exact string canonical text versus `src/stringRepair.h::Substrings` | Retained. Exact-search canonical ordinals constrain descendant enumeration; Re-Pair ranks support a different algorithm. Replacing either risks changing pruning or deterministic choices. | A proof of ordering/identity compatibility and independent oracle, witness and reversal tests before any shared representation. |
| Remaining mask count/intersection/`any` implementations | Retained. Owning and borrowed storage have different null handling, small/wide fast paths and inlining requirements. Only the demonstrated scan kernel is shared. | Separate microbenchmarks for each storage mode and whole-search evidence, not just fewer source lines. |
| `benchmarks/benchmark.py` telemetry counter/summary validation | Deferred. Historical schema versions and conditional counter presence have distinct errors and acceptance rules. Test fixture plumbing was consolidated without collapsing production policy. | Explicit schema-version coverage map, malformed reports and overflow/absent-counter cases. |
| Promotion versus scaling sample validators | Retained. `check_speedups.py` requires integer CPU clocks, including values beyond `2**53`, and the current corpus; `check_parallel_scaling.py` needs wall samples and compatible paired reports. | Extract a smaller proven common layer only if it avoids flags and preserves validation order and diagnostics. |
| Small report-loading adapters | Retained where path expansion and error wording differ. `check_speedups.load_result` remains straightforward; a textual exception-rewriting adapter was rejected as harder to maintain. | No further action for line reduction alone. |
| Graph versus string NDJSON decoders | Retained. The graph decoder consumes event keys and has different missing-event errors; the string decoder retains the event in its record. | A common event contract adopted by both callers, rather than a compatibility layer with special cases. |
| Broad CLI-case removal in `unitTests/unitTester.py` | Only the exact repeated help-boundary case was removed. Disabled-output tests seed telemetry sentinels even on non-telemetry builds; output-failure cases require precise diagnostics and cover extra targets; mixed-line-ending tests assert timing records; malformed-option cases check more than the flag name. | Preserve these assertions in an explicit coverage map before removing any larger block. Relevant functions are `run_cli_checks`, `run_input_output_matrix_checks`, `run_flag_matrix_checks` and `run_string_record_checks`. |
| Statistical producers and acceptance gates | Retained. Independent recomputation detects bad reports; censoring, integer precision and sequential versus compensated sums differ. | Preserve independent validation and exact numeric policy. Do not have a gate trust producer summaries. |
| Scientific references | Retained: `ReferenceSearch`, string Re-Pair oracle/naive implementation, `unitTests/molecular_oracle.py`, addition-chain oracles, brute-force isomorphism and certificate/pathway replay. | These validate production algorithms; sharing their implementation would weaken scientific assurance. |
| Other lookalikes | Retained: cache strategies, initial versus later DAG expansion, tree versus cyclic canonicalisation, fixed-width MOL versus native parsing, production versus test telemetry validators, and specialised build configuration. Tiny CLI integer/suffix parsers were not replaced with a generic utility. | Consolidate only with a coherent shared contract and evidence that algorithm, grammar, lifetime and performance requirements agree. |

## Commit stages and verification

Each implementation stage is a separate local commit referencing `#146`:

1. `1cf07a3`: CPU clock arithmetic and rollover tests.
2. `ba634be`: Process lifecycle, benchmark launch setup and report identities.
3. `5a7c190`: Graph/parser/string/CLI test scaffolding.
4. `3b299aa`: Graph kernels and molecular input preparation, with independent kernel tests.
5. `5bd10b1`: Exact string child exploration.
6. `d13096f`: Owner/view mask scans, after focused performance evaluation.
7. `7afff51`: Benchmark manifests and fingerprints.
8. `c748764`: Paired scaling setup and plot finishing.
9. `981a099`: Telemetry fixture plumbing and the exact redundant CLI case.
10. `570d659`: Validated report reuse and policy-specific regression tests.
11. `ec6359c`: MPI success coordination and scheduler resets.

Production and tooling code has 277 fewer lines overall; tests have 160 more
lines and CMake has 14 more. The added independent boundary and policy tests
are intentional. These counts exclude this record and generated artifacts.

All 41 configured CTest tests passed with GCC 15.2, C++20, Release optimisation,
strict warnings, full regression coverage, OpenMP, MPI, hybrid and telemetry
enabled. Six MPI tests initially failed before solver execution because the
sandbox blocked sockets; all six passed when rerun with local socket access.
The 171 Python benchmark tests, repository text policy, Ruff lint/format and
`git diff --check` also passed. The final build reports no pending work.

`consolidationKernelsTester` and `activeWordMaskTester` passed AddressSanitizer,
UndefinedBehaviorSanitizer and leak detection with `-O1 -g
-fsanitize=address,undefined -fno-omit-frame-pointer`. LeakSanitizer required
execution outside the sandbox's ptrace restriction. All reference oracles and
golden pathway data remain independent and unchanged.

## Performance evidence

The untouched `c478f7b` source was archived into
`build/issue146/baseline-source` and compiled with the same portable GCC Release
settings as the candidate. Paired serial runs used `taskset -c 0`, alternating
AB/BA order and two warm-up rounds. The reports retain binary hashes and sample
measurements under the ignored `build/issue146/` directory.

| Corpus | Cases / measured rounds | Paired wall speedup | Paired CPU-clock speedup | Local report |
| --- | --- | --- | --- | --- |
| Full | 15 / 12 | 1.0165x | 1.0171x | `build/issue146/full.json` |
| Profile excluding Paclitaxel | 4 / 12 | 1.0253x | 1.0253x | `build/issue146/profile-bounded.json` |
| Mask boundaries plus four-component amino-acid input | 7 / 30 | 1.0006x | 0.9995x | `build/issue146/mask-boundaries.json` |

Speedup is baseline divided by candidate, using the median paired suite totals.
All completed measurements matched the expected assembly index. This is
regression evidence, not a complete performance promotion or a parallel
scaling claim. The boundary suite's essentially unchanged CPU median is within
its 0.0023x median absolute deviation. The complete profile attempt stopped because Paclitaxel exceeded
the 60-second warm-up timeout in the untouched baseline; its comparison remains
unverified (`build/issue146/profile.log`).

The owner/view scan experiment also compared matching checksums across mask
widths/densities on a pinned CPU. Direct owner-to-view delegation was rejected
after slowdowns. The retained word-accessor implementation passed the focused
mask tests; a longer 65-bit sparse check gave a 0.929 candidate/baseline runtime
ratio. Count, intersection and other specialised mask operations remain separate.
