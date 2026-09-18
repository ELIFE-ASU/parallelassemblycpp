# Search locality and incumbent experiment

The final candidate defaults to an adaptive worker-local mask cache. Each generation starts unordered; the 512th distinct mask promotes eligible one- or two-word keys into a flat table. The first allocation is 1,024 slots (24 KiB), growing at 70% load to at most 65,536 slots. Unordered overflow preserves every mask needed by DAG construction, so only the flat allocation is capped. Wider masks remain unordered. Clearing resets activation while retaining allocated flat capacity for reuse. Exact canonicalization is unchanged. Set `ASSEMBLY_CANONICAL_MASK_CACHE=unordered` for the previous backend.

The bounded greedy bootstrap and child-first sibling donation remain opt-in (`PARALLELASSEMBLYCPP_GREEDY_BOOTSTRAP=1` and `PARALLELASSEMBLYCPP_CHILD_FIRST=1`). The bootstrap is capped at 32 accepted steps, 65,536 root jobs and 262,144 parent-fragment checks. It publishes feasible incumbents without consuming root work or recording states as searched.

## Final adaptive-cache timings

All 37 unique manifest cases were checked in two alternating paired rounds at each worker budget (296 timed calculations). Every expected assembly index matched. Times below are medians; speedups are medians of paired baseline/candidate ratios, so they need not equal the ratio of medians. Values above 1 mean faster. Two pairs are exploratory evidence, especially for small differences; these are not the repository's longer promotion gate. This section uses only the completed `adaptive-full-8.json` and `adaptive-full-28.json` reports. Earlier fixed-flat experiments and the interrupted adaptive attempt are retained separately and excluded from these ratios. Two unrelated CPU-intensive `vasp_std` processes were active on the host, consuming over 1000% aggregate CPU near the start of the final pass; host snapshots are retained with the experiment. These contended exploratory timings cannot establish stable small speedups or regressions. The quiet rerun below concerns an older policy and does not remove this limitation from the final adaptive results.

| Workers | Whole-corpus paired speedup | MAD |
|---:|---:|---:|
| 8 | 1.056x | 0.130x |
| 28 | 0.950x | 0.013x |

Longer cases:

| Workers | Case | Baseline seconds | Candidate seconds | Paired speedup | MAD |
|---:|---|---:|---:|---:|---:|
| 8 | paclitaxel | 48.711 | 47.800 | 1.015x | 0.019x |
| 8 | amino-acid-scale-12c | 6.300 | 6.750 | 0.940x | 0.029x |
| 8 | amino-acid-scale-13c | 26.205 | 20.521 | 1.277x | 0.497x |
| 28 | paclitaxel | 15.489 | 19.092 | 0.830x | 0.114x |
| 28 | amino-acid-scale-12c | 2.225 | 2.594 | 0.869x | 0.109x |
| 28 | amino-acid-scale-13c | 12.729 | 10.473 | 1.182x | 0.192x |

Every case and raw sample is retained in [all-cases.csv](../build/search-locality-results/all-cases.csv), [adaptive 8-worker JSON](../build/search-locality-results/adaptive-full-8.json), and [adaptive 28-worker JSON](../build/search-locality-results/adaptive-full-28.json). All small cases remain in the CSV and raw reports; no timing gate or outlier exclusion was applied.

## Why activation became adaptive

The earlier fixed-flat v1 candidate activated its table immediately, beginning with 256 slots (6 KiB). Despite their filenames, `final-full-*` reports describe that older policy. The four-pair 8-worker full run was affected by an unrelated CPU-intensive computation. The four-pair 28-worker full run showed slower medians for the two long amino-acid cases, whose per-worker canonical-mask working sets were small. These observations motivated retaining the unordered path for small working sets and promoting only after 512 distinct masks. They do not identify the cause of every timing difference.

A completed quiet 8-worker rerun used the same fixed-flat v1 binary, four pairs and the three longer cases. It is supplemental evidence about v1, not the final adaptive policy. The earlier full runs total 592 timed calculations and the quiet rerun adds 24; all are separate from the 296-calculation final comparison.

| Earlier experiment | Workers | Case | Baseline seconds | Fixed-flat seconds | Paired speedup | MAD |
|---|---:|---|---:|---:|---:|---:|
| v1 full | 8 | paclitaxel | 37.202 | 36.751 | 1.020x | 0.005x |
| v1 full | 8 | amino-acid-scale-12c | 5.133 | 4.544 | 1.093x | 0.088x |
| v1 full | 8 | amino-acid-scale-13c | 19.026 | 18.885 | 1.002x | 0.013x |
| v1 full | 28 | paclitaxel | 15.288 | 14.871 | 1.023x | 0.004x |
| v1 full | 28 | amino-acid-scale-12c | 1.948 | 2.064 | 0.960x | 0.023x |
| v1 full | 28 | amino-acid-scale-13c | 7.829 | 8.067 | 0.973x | 0.021x |
| v1 quiet | 8 | paclitaxel | 28.019 | 27.723 | 1.008x | 0.004x |
| v1 quiet | 8 | amino-acid-scale-12c | 4.048 | 4.011 | 1.014x | 0.008x |
| v1 quiet | 8 | amino-acid-scale-13c | 15.909 | 15.715 | 1.011x | 0.002x |

The quiet v1 three-case aggregate paired speedup is 1.011x. Preserved reports: [v1 8-worker full](../build/search-locality-results/final-full-8.json), [v1 28-worker full](../build/search-locality-results/final-full-28.json), [v1 quiet rerun](../build/search-locality-results/quiet-long-8.json). The interrupted adaptive run is retained in `build/search-locality-results/adaptive-interrupted-contention/`; it is excluded as an incomplete experiment, not silently substituted or merged into the final comparison.

## Separate search-policy experiments

Both sides below use the earlier fixed-flat v1 executable, before adaptive activation was introduced. Each comparison has two alternating pairs, on Paclitaxel only; they are screening experiments, not evidence of final adaptive-cache performance or a stable small speedup. Host contention also limits interpretation of these screens.

| Workers | Experiment versus existing policy | Paired speedup | MAD |
|---:|---|---:|---:|
| 8 | Bounded greedy + child-first | 1.132x | 0.327x |
| 28 | Bounded greedy + child-first | 0.951x | 0.035x |
| 8 | Disable matching-batch refresh | 0.981x | 0.019x |
| 28 | Disable matching-batch refresh | 0.972x | 0.018x |

The existing relaxed refresh cadence remains enabled: class boundaries and every 64 matching-bound evaluations. No additional synchronization was added to inner matching iterations. `PARALLELASSEMBLYCPP_MATCHING_REFRESH=0` enables the comparison above.

The earlier `flat-pilot-*` reports also use an immediately active fixed-flat prototype and remain separate from final adaptive results.

## Untimed trajectories and work

Telemetry is excluded from timing aggregates. These are individual instrumented runs, whose instrumentation and scheduling can affect their elapsed times. An unrelated CPU-intensive job was active while the final adaptive traces were collected. Their work counters and retained-byte snapshots remain useful measurements of those runs, but incumbent discovery times are host-contended and must not be read as a controlled wall-time speedup.

Incumbent timestamps start before parallel root preparation/enumeration and use a steady clock; MPI uses barrier-relative rank-local epochs, not a globally synchronized clock. Indices are internal, before optional disconnected component compensation.

Expanded states entered duplicate enumeration, including the prepared root on each rank. Pruned states are materialized candidates rejected by a bound or transposition lookup. Skipped classes, occurrence pairs, and fragment-pair blocks have separate counters and are not added to state counts.

| Workers | Paclitaxel variant | Final incumbent seconds | Expanded states | Pruned states | Matching visits | Mask misses |
|---:|---|---:|---:|---:|---:|---:|
| 8 | Unordered control | 6.691 | 36,893,628 | 169,456,149 | 206,349,776 | 473,078 |
| 8 | Final adaptive (contended) | 3.784 | 37,289,217 | 170,799,700 | 208,088,916 | 474,917 |
| 28 | Unordered control | 2.871 | 39,771,056 | 179,041,876 | 218,812,931 | 825,711 |
| 28 | Final adaptive (contended) | 3.116 | 39,622,406 | 178,480,014 | 218,102,419 | 805,460 |

An earlier incumbent alone does not establish less total work. Differences in pruning, search order and local cache replication remain distinct possible contributors; the counters do not identify a unique cause.

![Paclitaxel incumbent trajectories](../build/search-locality-results/incumbents.png)

The older `traces-flat-*` and `traces-final-flat-*` files describe immediately active fixed-flat prototypes. `traces-greedy-child-*` uses unordered caches plus greedy bootstrap and child-first donation. Those historical trajectories are preserved but are not mixed into the final adaptive tables or plot.

### Earlier unordered search-policy traces

Both sides of this separate Paclitaxel experiment use unordered caches. These are earlier individual instrumented runs, not paired timings or measurements of the final adaptive policy.

| Workers | Earlier unordered policy | Final incumbent seconds | Expanded states | Matching visits |
|---:|---|---:|---:|---:|
| 8 | Existing search | 6.691 | 36,893,628 | 206,349,776 |
| 8 | Greedy + child-first | 3.872 | 37,196,484 | 207,647,688 |
| 28 | Existing search | 2.871 | 39,771,056 | 218,812,931 |
| 28 | Greedy + child-first | 4.984 | 41,813,520 | 227,269,092 |

At 8 workers, greedy plus child-first discovered the final incumbent earlier but expanded more states and performed more matching work. At 28 workers, it discovered the final incumbent later and also performed more work. These traces do not establish a consistent reduction in search work, so both switches remain opt-in. Their raw sources are `traces-control-8.json`, `traces-control-28.json`, `traces-greedy-child-8.json` and `traces-greedy-child-28.json`.

## Retained worker-local storage

Each worker is sampled before its caches are destroyed. The aggregate is a sum of those snapshots, not process RSS or a simultaneous memory peak. Assembly-state bytes count actual retained slot capacity and upstream monotonic-arena allocations. Unordered/container measurements are estimates excluding allocator headers and some arena slack. Shared cache and immutable seed storage are reported separately and are excluded from these local totals.

| Workers | Case | Control mask MiB | Adaptive mask MiB | Control state GiB | Adaptive state GiB |
|---:|---|---:|---:|---:|---:|
| 8 | paclitaxel | 25.997 | 17.008 | 2.544 | 2.625 |
| 8 | amino-acid-scale-12c | 0.061 | 0.077 | 0.382 | 0.382 |
| 8 | amino-acid-scale-13c | 0.066 | 0.083 | 1.529 | 1.529 |
| 28 | paclitaxel | 48.952 | 35.768 | 2.775 | 2.710 |
| 28 | amino-acid-scale-12c | 0.183 | 0.197 | 0.511 | 0.513 |
| 28 | amino-acid-scale-13c | 0.198 | 0.214 | 1.770 | 1.776 |

These retained-capacity measurements do not justify growing mask caches. Full worker records and trajectories are in `build/search-locality-results/traces-*.json`; placement sweep reports now also emit `search-profiles.json`.

## Reproduction and validation

Baseline: a frozen snapshot of the initial working tree, including its pre-existing uncommitted changes. Host: Intel Core i7-14700KF, 20 mixed P/E physical cores and 28 logical CPUs. Builds use GCC 15.2, Release/x86-64-v3, strict warnings, OpenMP, and neither LTO nor PGO. Placements use `OMP_PLACES=cores`, `OMP_PROC_BIND=close`, `OMP_DYNAMIC=FALSE`, and both `OMP_NUM_THREADS`/`OMP_THREAD_LIMIT` set to 8 or 28. These are worker budgets on heterogeneous cores; automatic search can choose fewer workers for small cases. Shared cache policy is `shared`, with unlimited entry admission.

Complete paired timing reports use sequential runs with no concurrent task build/test workload and zero separate warmups. No samples within those completed reports were discarded. The interrupted run and contended historical evidence are disclosed above. The final adaptive timings also ran under sustained unrelated host contention; process snapshots at the start and end are retained with the raw evidence. The runner rotates cases and alternates baseline/candidate order. Exact commands, executable fingerprints, environment and raw reports are retained under `build/search-locality-results/`. The runners are `run_comparison.py` and `collect_traces.py` in that directory.

Final-policy validation:

- 29 CTest checks passed: [ctest-adaptive.log](../build/search-locality-results/ctest-adaptive.log).
- 2 CTest checks passed: [ctest-adaptive-experiments.log](../build/search-locality-results/ctest-adaptive-experiments.log).

The main CTest configuration includes full serial regression, OpenMP/MPI/hybrid parity, telemetry invariants, library and canonicalization checks, the Python benchmark suite, and differential cache width/growth/overflow/promotion tests. The separate experimental log checks opt-in search-policy telemetry. Focused ASan/UBSan checks also passed on the final adaptive mask cache, including allocation-failure behavior. Earlier fixed-flat test logs remain supplemental and are not substituted for final adaptive-policy checks. See [benchmark instructions](README.md) for complete test and experiment commands.
