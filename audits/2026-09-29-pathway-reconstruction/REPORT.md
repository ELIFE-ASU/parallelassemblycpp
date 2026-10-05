# Deterministic pathway reconstruction

> Historical investigation from 29 September 2026. Implementation comparisons,
> validation counts, and timings below describe that change and its measured
> binaries. See the [audit index](../README.md) for scope and the
> [benchmark guide](../../benchmarks/README.md) for current commands.

## Measurement

The parallel solver proves an assembly index without recording a pathway,
then searches the prepared root frontier in serial order for the first witness
of that index. Previously, parallel telemetry ended before that second search
and disabled its counters. The benchmark runner also passed `--pathway=0`.

Telemetry now keeps `pathway_reconstruction` separate from optimization:

- `attempted` and `completed` describe the reconstruction attempt, including
  output failures.
- `elapsed_seconds` uses the steady clock and includes setup, witness search,
  output, and worker cleanup on the primary process.
- `cpu_seconds` uses the process CPU clock, so it can include activity from
  other process threads; it is not a thread CPU measurement.
- `counters` contains the same counter schema as a parallel worker, without
  contributing to the optimization aggregate or its legacy top-level counters.

The fields are zero when a second search is not attempted, including serial
searches that record their witness during optimization. They do not time
integrated serial witness generation. String assembly does not support this
telemetry. Parallel optimization timing and reconstruction timing are distinct
regions, but their sum excludes some reductions and final result output; use
the benchmark's process wall time for end-to-end comparisons.

The benchmark's opt-in `--pathways` enables pathways in baseline, candidate,
warm-up, and telemetry runs. The ordinary samples measure complete process
wall time; reconstruction counters and phase time come from a separate
instrumented run. Reports record `pathways_enabled`; comparison tools reject
mixed modes and treat historical reports without the field as disabled.

Example:

```bash
python benchmarks/benchmark.py \
  --input unitTests/Ceftiolene.mol --expected 23 \
  --executable build/parallel/ParallelAssemblyCppOMP \
  --telemetry-executable build/parallel/ParallelAssemblyCppOMPTelemetry \
  --candidate-parallel on --candidate-env OMP_NUM_THREADS=4 \
  --runs 9 --warmup 1 --pathways --telemetry \
  --json-output build/pathways-ceftiolene.json
```

## Changes to witness storage

The synchronous DFS records borrowed decision views into immutable masks.
It materializes an owning, replayable witness only when the incumbent improves.
Prepared roots borrow the immutable occurrence-word buffer, serial roots
borrow their active enumeration frame, and recursive decisions borrow the
immutable DAG. These sources survive the synchronous child search.

Snapshotting retains a common prefix of the previous witness and copies only
the changed suffix. Popping decisions truncates the known common prefix. This
avoids repeatedly copying ancestors when ordinary serial search improves at
successively greater depths. Owning snapshots survive producer-frame reuse and
unwinding; borrowed decisions never cross a thread or MPI boundary.

The known-target reconstruction normally retains one witness. Before this
change, it materialized two masks for every admitted state; now it materializes
only the winning steps. This removes wide-mask allocation/copy work without
changing branch order, pruning, matching equivalence, or the stopping rule.

## Replay and checkpoint investigation

An arbitrary parallel winner cannot replace the deterministic search:

1. Count-only search can traverse fragment-pair blocks in a different order
   from pathway search's legacy reverse-occurrence order.
2. The shared transposition table admits whichever worker arrives first. Its
   entries describe encountered states, not proof that all earlier serial
   branches have been exhausted.
3. Bounds prune equal optima once an incumbent is known. Recording the first
   parallel winner would therefore select a schedule-dependent witness.
4. Donated tasks carry a fragmentation state, not the complete decision
   ancestry needed to replay an original-edge pathway. Canonical IDs can also
   belong to a worker-local validity domain.

A complete replay design needs mask-free decision ancestry on donated tasks,
a serial-order key for each decision, and proof that all earlier branches
cannot contain the target. Checkpoints must preserve the appropriate canonical
ID domain or recanonicalize on import. Merely saving the winning masks or the
current transposition table does not provide the ordering proof.

A narrower follow-up experiment could search root job zero independently in
legacy order, with a private table and no donations or racing bound, and retain
its optimal witness. If it reaches the global optimum, its witness is a
candidate for deterministic replay after checking equivalence to known-target
traversal; otherwise the existing reconstruction remains necessary. Measure
checkpoint hit rate, duplicated work, lost parallelism, and total runtime
before enabling that policy. The current change keeps the deterministic search.

## Validation

- Pathway ownership tests exercise one-word and 129-edge masks, producer reuse,
  prefix extension, divergent suffix replacement, and empty witnesses.
  The witness tests also pass under AddressSanitizer and UndefinedBehaviorSanitizer.
- Repeated OpenMP, MPI, and hybrid pathway runs match the golden ketoconazole
  JSON. Telemetry tests cover success, disabled pathways, integrated serial
  generation, write failures, and isolation from worker counter reductions.
- Serial regression and telemetry regression pass.
- All 137 benchmark/tooling tests pass, including mode compatibility and
  malformed reconstruction telemetry; Ruff and repository text checks pass.

## Local timing evidence

Measured on an Intel Core i7-14700KF with GCC 15.2.0, C++20, `-O3 -DNDEBUG`,
and four OpenMP workers. The baseline is revision `2608a16` with the new
isolated timing/counters but the original owning decision stack and root-mask
copies; the candidate includes all storage changes described above. Both
binaries enable telemetry. There is one warm-up pair and nine measured pairs
per case, with execution order alternating each round. Compilation and other
task tests had finished before this comparison. Thread settings were:

```text
OMP_NUM_THREADS=4
OMP_DYNAMIC=FALSE
OMP_PLACES=cores
OMP_PROC_BIND=close
OMP_WAIT_POLICY=PASSIVE
```

Both binaries ran each input with
`--parallel=on --threads=4 --pathway=1 --telemetry=1`. The raw
[timing samples](timings.csv) include pathway SHA-256 digests. All 144 measured
executions and their warm-ups produced byte-identical pathways and identical
reconstruction counters for the same input across versions and repetitions.
The comparison script and binary metadata were retained locally in
`build/pathway-measurements/compare.py` and
`build/pathway-measurements/final-comparison.json`. These unpublished build
artifacts are not included in a repository checkout.

Times below are milliseconds, reported as median +/- median absolute deviation
(MAD); optimization is the candidate's median parallel elapsed time.

| Input | Reconstruction before | Reconstruction after | Optimization after |
| --- | ---: | ---: | ---: |
| ketoconazole | 0.564 +/- 0.035 | 0.573 +/- 0.051 | 6.766 |
| Ceftiolene | 3.219 +/- 0.131 | 3.192 +/- 0.105 | 3.845 |
| erythromycin | 10.399 +/- 0.553 | 10.348 +/- 0.520 | 401.401 |
| clarithromycin | 12.412 +/- 0.978 | 11.576 +/- 0.310 | 1207.602 |
| phosphatidylcholine | 1.173 +/- 0.057 | 1.126 +/- 0.023 | 123.881 |
| amino-acid scale 09c | 0.183 +/- 0.016 | 0.183 +/- 0.019 | 74.792 |
| homogeneous path 127b | 0.357 +/- 0.047 | 0.340 +/- 0.033 | 26.197 |
| homogeneous path 129b | 0.283 +/- 0.044 | 0.270 +/- 0.039 | 19.262 |

These samples establish the missing cost, not a reliable elapsed-time win.
The largest apparent reduction (about 7% for clarithromycin) remains within
the baseline's MAD; most differences are smaller. The reduction in mask
materializations is structural and covered by ownership tests, but broader
performance promotion is not justified from this desktop run. Ceftiolene's
reconstruction is comparable to its parallel optimization, while the longer
erythromycin/clarithromycin optimizations dominate their reconstruction.
A preliminary Paclitaxel run exceeded a 30-second wall timeout and is excluded;
no Paclitaxel reconstruction result or speedup is claimed.
