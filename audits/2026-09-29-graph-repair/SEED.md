# Historical molecular exact search seeded by the Re-Pair upper bound

**Historical molecular/graph experiment, subsequently removed.** Exact graph
calculations run without a graph Re-Pair prepass. The explicit
`--algorithm=re-pair` upper-bound-only mode remains available. Exact string
search is separate: it uses a Re-Pair seed and retains its witness, introduced
in commit `8d4585c`. See the [current CLI reference](../../docs/cli.md) and
[audit index](../README.md).

The implementation description, measurements, and validation below describe
the removed graph experiment. Its binaries and source patch were preserved
locally under `build/graph-repair-seeded-experiment/` before restoring the
unseeded graph solver; these artifacts are unpublished and are not included
in a repository checkout.

In the experimental version, `--algorithm=re-pair` returned the fast constructive upper bound alone. `--algorithm=full` (the default) computed that bound first, retained its construction, and used the incumbent to prune exact molecular search. The same initialization applied to serial, OpenMP, MPI, hybrid and library calls.

The initial incumbent is expressed in the exact solver’s original, uncompensated bond-count convention. Disconnected compensation is applied at output. The physical duplication witness is translated to the preprocessed edge universe after configuring each thread’s mask arena. Thus a seed-equal optimum or an early search limit still has a valid pathway in the existing exact JSON format. Seed preparation counts toward calculation time and the cooperative runtime budget.

The witness is obtained from the Re-Pair token forest by visiting larger fragments first, retaining one representative and removing the other equal fragments before descending into the representative’s children. It can omit unused grammar productions and therefore cannot worsen the original grammar bound. The standalone bound and certificate format are unchanged.

## Validation

- All 25 CTest checks passed, including the manifest and telemetry regressions.
- The new integration harness passed 56 CLI runs covering a seed already equal to optimum, improvement from seed, runtime/enumeration limits, pathway/no-pathway, disconnected compensation and empty/isolated graphs. Repeated library calls also passed without leaking options or graph state.
- Independent validation replayed 1,904 grammar certificates and 1,904 seed duplication paths, and checked 1,895 bounds against an exact small-graph oracle. Six corrupted duplication witnesses were rejected. See [seed_validation.json](seed_validation.json).
- OpenMP index/pathway and execution-policy checks passed. Thirty MPI/hybrid solver runs passed, including mask boundaries and improved-seed reconstruction. Additional equal-seed icosane pathway checks exercise the retained construction.
- All 13 reviewed pathway fixtures were independently replayed with unchanged exact indices. Four valid alternate optimum witnesses were updated (113, graphio_test, icosane and thioflavin); strict golden comparison remains.

## Exact-search timing

This integration did **not produce a general speedup** on the measured corpus. The median case had a calculation speed ratio of **0.968×** (about **3.3% longer**), including the Re-Pair prepass. Ratios are old exact time divided by seeded exact time; above 1 is faster.

Across the 34 completed cases, the sums of per-case median calculation times were **15.212 s before** and **15.243 s after** (about 0.2% more time). Median-case end-to-end speed ratio was **0.972×**. The new bound adds noticeable proportional overhead on the smallest problems and changes total exact-search cost little on the hardest completed ones. These are one-machine measurements, not a guarantee of speedup or slowdown on another workload.

Both versions completed the same 34 cases with identical exact results, all matching their manifest reference values. Both versions exceeded the 10-second process deadline on paclitaxel, amino-acid-scale-12c and amino-acid-scale-13c; all three are excluded from the completed-case timing summaries.

Times below are separate medians of six measurements per version. Speed ratios
are medians of the six paired old/seeded ratios; they need not equal the quotient
of the two displayed median times. The median-case summary gives each completed
case equal weight, whereas the sum of median times is dominated by the hardest
completed searches.

| Case | Old exact (ms) | Seeded exact (ms) | Paired calculation speed ratio | Exact index |
| --- | ---: | ---: | ---: | ---: |
| icosane | 0.178 | 0.281 | 0.739× | 6 |
| sr1001 | 0.447 | 0.539 | 0.823× | 22 |
| 5843 | 7.810 | 8.050 | 0.970× | 8 |
| bisphenylmaleimide | 11.012 | 11.668 | 0.950× | 10 |
| ketoconazole | 29.664 | 25.412 | 1.120× | 22 |
| THC | 4.054 | 4.259 | 0.931× | 12 |
| sucrose | 2.906 | 3.262 | 0.938× | 8 |
| graphio | 2.445 | 2.468 | 0.993× | 5 |
| Ceftiolene | 13.342 | 10.941 | 1.181× | 23 |
| Cefquinome | 7.280 | 7.198 | 0.999× | 25 |
| Cefpirome | 6.798 | 6.983 | 0.965× | 25 |
| dipyridamole | 7.892 | 8.022 | 0.991× | 15 |
| ceftiofur | 4.422 | 5.504 | 0.918× | 26 |
| dienogest | 4.717 | 5.033 | 0.974× | 11 |
| folic_acid | 3.047 | 3.192 | 0.956× | 20 |
| phosphatidylcholine | 459.487 | 454.307 | 1.009× | 22 |
| erythromycin | 1519.142 | 1530.669 | 0.998× | 20 |
| clarithromycin | 4857.277 | 4884.317 | 0.998× | 21 |
| amino-acid-scale-02c | 0.097 | 0.146 | 0.654× | 10 |
| amino-acid-scale-03c | 0.283 | 0.332 | 0.845× | 13 |
| amino-acid-scale-04c | 0.754 | 0.886 | 0.889× | 15 |
| amino-acid-scale-05c | 0.935 | 1.084 | 0.883× | 17 |
| amino-acid-scale-06c | 4.358 | 4.937 | 0.964× | 19 |
| amino-acid-scale-07c | 22.595 | 29.343 | 0.827× | 21 |
| amino-acid-scale-08c | 51.756 | 59.335 | 0.887× | 23 |
| amino-acid-scale-09c | 303.135 | 303.476 | 1.000× | 25 |
| amino-acid-scale-10c | 1465.162 | 1468.268 | 1.016× | 27 |
| amino-acid-scale-11c | 6244.474 | 6217.891 | 1.000× | 29 |
| mask-boundary-path-063b | 3.159 | 3.161 | 1.016× | 8 |
| mask-boundary-path-064b | 3.311 | 3.014 | 1.012× | 6 |
| mask-boundary-path-065b | 4.541 | 4.816 | 0.967× | 7 |
| mask-boundary-path-127b | 64.916 | 77.258 | 0.977× | 10 |
| mask-boundary-path-128b | 48.275 | 47.295 | 0.947× | 7 |
| mask-boundary-path-129b | 51.987 | 49.201 | 1.046× | 8 |

## Why the seed did not save much search work

Four separate instrumented diagnostic pairs (pathway output disabled) show that
the seed usually does not change the expensive search work in these examples:

| Molecule | Matching visits, unseeded → seeded | Exact canonicalisation calls, unseeded → seeded |
| --- | ---: | ---: |
| icosane | 107 → 103 | 121 → 121 |
| sucrose | 2,763 → 2,763 | 2,441 → 2,441 |
| ketoconazole | 35,789 → 35,789 | 21,855 → 21,855 |
| erythromycin | 2,217,419 → 2,217,419 | 469,692 → 469,692 |

The unseeded incumbent trajectories reach or beat the Re-Pair seed in the first
few improvements. For erythromycin, the seed is 27, but unseeded search progresses
through 39, 32, 29, 27, 26, 25, 24 and 23 almost immediately after enumeration.
It finds the eventual optimum 20 early and spends most remaining time proving
there is no better answer. A starting upper bound does not remove that proof.

The initial enumeration receives no incumbent and is unchanged by seeding.
Re-Pair also constructs and canonicalises fragments in a separate prepass and
derives its duplication witness. Its pair cache is local and the tree interner
is cleared before exact search, so those canonicalisation results are not reused.
The unchanged exact counters above therefore exclude additional Re-Pair work.

[seed_profile.json](seed_profile.json) contains raw telemetry, incumbent logs,
results and binary/input fingerprints. These diagnostic timings are instrumented
single runs and must not replace the paired benchmark timings. The very small
0.2% aggregate timing difference is not evidence of a meaningful slowdown by
itself; the counters explain why there is little opportunity for a speedup here.

## Method and reproduction

Same GCC 15.2.0 portable Release optimization (`-O3 -DNDEBUG`) and Intel Core i7-14700KF CPU 0. Six adjacent measured pairs follow one warmup pair, alternating old/seeded order. Each calculation uses a fresh process, removes explicit hydrogens, disables pathway output, uses the default disconnected convention and keeps the 50-million enumeration cap. The measured timer includes seed construction and exact setup/search, but excludes parsing and process startup. End-to-end process times are recorded separately. An incomplete sample excludes its entire case from paired aggregates.

```bash
python audits/2026-09-29-graph-repair/seed_speed.py \
  --baseline build/graph-repair-before-seed/parallelassemblycpp_graph_repair_speed_probe \
  --seeded build/graph-repair-seeded-experiment/parallelassemblycpp_graph_repair_speed_probe \
  --cpu 0 --runs 6 --warmup 1 --timeout 10
```

The old probe was preserved locally before integration; its executable hash matches [the previous timing data](speed.json), which records its source hashes. The candidate used the removed seeded implementation, preserved at the local path above. [seed_speed.json](seed_speed.json) retains all original samples, statuses, input and executable hashes, source provenance and summary calculations. At the end of the experiment, its original candidate path was replaced by the restored unseeded build. Reproduction requires the preserved historical binaries or reconstruction of the corresponding sources; a build of the current checkout does not recreate this candidate. The complete pair counts, exact results, ratios and fingerprints were independently checked at the time of the experiment.

The earlier **47× median speedup** compared obtaining a bound with solving exactly. It is not an exact-search acceleration claim. See [that separate experiment](SPEED.md).
