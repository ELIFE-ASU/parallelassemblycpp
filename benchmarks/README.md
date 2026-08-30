# Benchmarks

The benchmark tools run isolated AssemblyCpp calculations, verify their
assembly indices, and report wall time and program-reported `std::clock` ticks.

## Quick start

From the repository root:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite quick
```

Use `--list-cases` to inspect the selected cases and
`python benchmarks/benchmark.py --help` for all runner options.

## Corpus

`cases.tsv` is the maintained benchmark manifest.

| Column | Meaning |
| --- | --- |
| `name` | Unique command-line-safe case name. |
| `input` | Input path, relative to this directory unless absolute. |
| `expected_assembly_index` | Required result for every run. |
| `expectation` | `reviewed` regression value or `provisional` benchmark guard. |
| `suites` | Comma-separated suite names. |
| `workload` | Short label shown in reports. |

Malformed rows, duplicate names, missing files, and unknown suite or expectation
values are rejected.

| Suite | Scope |
| --- | --- |
| `quick` | Short, varied regression workloads for routine checks. |
| `full` | All quick cases plus larger reviewed inputs. |
| `profile` | Longer search-heavy inputs; most expectations are provisional. |
| `scaling` | Cumulative amino-acid and 64-bit mask-boundary series. |

The largest scaling case has 100 bonds and can take several minutes and over
1 GiB of memory. For a smoke run, use `--runs 1 --warmup 0`.

## Common runs

Run a suite and save the raw samples and summary:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite full \
  --json-output build/full.json
```

Run one custom input:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --input unitTests/sucrose.mol \
  --expected 8
```

Collect one additional, untimed telemetry run per case:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --telemetry-executable build/performance/AssemblyCppTelemetry \
  --suite scaling \
  --telemetry \
  --json-output build/scaling.json
```

Telemetry is excluded from timing aggregates. It records graph size, retained
masks, matching and canonicalisation activity, cache rates, phase clock ticks,
and phase memory. Parallel telemetry additionally records rank/thread topology,
distributed root-queue participation, dynamic branch leases, global root-branch
coverage, depth-two and deeper task transfers, local executions and steals,
scheduler idle waits, deep-refill activations, task-queue high-water marks,
maximum executed task depth, incumbent warm starts, steady-clock worker timing,
and all 31 raw search counters per worker plus their exact aggregate. The
parallel aggregate also reports shared-L2 hits, misses, collision-chain probe
steps, retained entry bytes, and contended shard-lock waits and wait time.
Parallel phase memory is disabled because `/proc` peak resets are process-wide.
A cache rate is `null` when no lookup occurred. Aggregate elapsed time is the
critical parallel-region wall time; worker elapsed and busy values are sums,
and each worker's busy time is its elapsed time minus measured scheduler-idle
wait time, not CPU utilization. The aggregate queue high-water mark and maximum
task depth are maxima across workers; the other scheduler event fields are
sums. Legacy VF2 counters remain in the schema but are zero with the exact
cyclic canonicaliser.

The aggregate and worker records expose `deeper_tasks_spawned`,
`deeper_tasks_executed`, `task_steal_attempts`, `task_steals`,
`local_task_executions`, `scheduler_idle_waits`,
`scheduler_idle_nanoseconds`, `deep_refill_activations`,
`task_queue_high_watermark`, and `maximum_task_depth_executed`. The reported
`busy_timing_method` is `elapsed_minus_scheduler_idle_time`.

## Paired comparisons

Keep the previous executable and pass it as the baseline:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/AssemblyCpp-before \
  --executable build/performance/AssemblyCpp \
  --suite quick \
  --json-output build/quick.json
```

Baseline and candidate runs are adjacent and alternate AB/BA order. Paired mode
defaults to six rounds. A speedup above `1.0` means the candidate is faster.
Use an even run count and compare on the same idle, pinned host.

`check_speedups.py` validates a complete four-suite promotion set:

```bash
python benchmarks/check_speedups.py \
  build/quick.json \
  build/full.json \
  build/profile.json \
  build/scaling.json
```

The gate requires every case clock median and every suite round-total wall and
clock median to exceed `1.0`. It recalculates medians from raw paired samples
and rejects stale corpus or executable fingerprints. Shared CI does not enforce
timing thresholds because host contention makes them unreliable.

Promotion reports require 100 rounds for `quick` and `full`, 6 for `profile`,
and 30 for `scaling`.

## LTO and PGO

Build the LTO candidate with:

```bash
cmake --preset performance-lto
cmake --build --preset performance-lto
```

PGO requires GCC. Generate, train, and consume the profile in order:

```bash
cmake --preset pgo-generate
cmake --build --preset pgo-train
cmake --preset performance-pgo
cmake --build --preset performance-pgo
```

`pgo-training.tsv` assigns a repetition count to every corpus case. Training
clears old GCC data, verifies each result, and records fingerprints for the
weights, manifest, and inputs. The profile-use build rejects incomplete or
stale training data. Re-run all four benchmark suites before promoting an LTO
or PGO build.

## Parallel scaling

The `parallel` preset builds serial, OpenMP, MPI, and hybrid executables plus a
telemetry-enabled sibling for each parallel topology:

```bash
cmake --preset parallel
cmake --build --preset parallel
```

Example launch commands:

```bash
OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./build/parallel/AssemblyCppOMP molecule.mol --pathway=0

mpirun --map-by slot --bind-to core -n 4 \
  ./build/parallel/AssemblyCppMPI molecule.mol --pathway=0

OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close \
  mpirun --map-by slot:PE=2 --bind-to core -n 2 \
  ./build/parallel/AssemblyCppHybrid molecule.mol --pathway=0

OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./build/parallel/AssemblyCppOMPTelemetry molecule.mol \
    --pathway=0 --telemetry=1

mpirun --map-by slot --bind-to core -n 4 \
  ./build/parallel/AssemblyCppMPITelemetry molecule.mol \
    --pathway=0 --telemetry=1
```

These placement flags use Open MPI syntax. Each process enumerates the root
once, then shares an immutable processed graph, canonical seed, runtime DAG,
and serialized root-job table. Workers own their post-seed canonical deltas,
fragmentation scratch, and search caches. MPI and hybrid ranks request disjoint
chunks from a rank-zero global queue, so work follows each rank's actual local
capacity while every root job is executed exactly once. Wide masks are rebuilt
inside the receiving worker from serialized words.

The adaptive MPI default uses one-root worker leases, with each rank-level
broker refill bundling one lease per local worker. This bounds tail imbalance
when a single root is much more expensive than its neighbours; faster ranks
simply issue more requests. Serial/OpenMP scheduling retains guided leases for
larger frontiers. Root jobs take priority over transferred work.
Within hybrid ranks, observed idle pressure can make a root search expose
immediate children as depth-two tasks. The idle trigger is at least half the
workers on ranks with fewer than eight local workers, and roughly one quarter
(with a minimum of two) on larger ranks.

Each rank has one lazily populated deque per local worker. Producers push to
their own deque and execute its newest task (owner LIFO); idle peers take the
oldest task from another worker (thief FIFO). After a worker fails to find
local or stealable work and the live idle count reaches the trigger, the
scheduler can request a refill before entering its 1 ms signal-poll wait.
Depths three and four are armed only through this observed-starvation path,
one level at a time while work at the preceding depth remains outstanding.
The estimated rank frontier must be below eight tasks per worker to request a
starvation refill. Donation stops around a target of sixteen tasks per worker,
with a hard rank-wide task-slot cap of thirty-two times the local worker count
and an absolute transferred-task depth cap of four. Within the process-local
scheduler, root/task outstanding counts, ready/slot/idle counts, and the
donation request occupy separate 64-byte-aligned storage; worker deque state is
aligned separately too. Its serial/OpenMP root-lease cursor is isolated there;
distributed searches use the broker's separate cursor.

The largest initial duplicate is evaluated first on each rank to publish a
valid first-step incumbent before concurrent searching begins. Improved bounds
are propagated periodically with a passive-target RMA minimum, allowing remote
progress to tighten local pruning before the final result reduction. Root work
uses the FUNNELED request broker, and the RMA heartbeat also propagates
cancellation so it stops issuing new chunks after an observed interrupt or
search failure. Task transfer is disabled in an MPI-only rank with one local
worker, but remains available within hybrid ranks.

Set the positive `ASSEMBLYCPP_BRANCH_LEASE_SIZE` environment variable to use a
fixed root lease size. Only MPI rank zero writes output.

The runner accepts separate launchers and environment variables for each role:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/parallel/AssemblyCpp \
  --baseline-launcher "taskset -c 0" \
  --executable build/parallel/AssemblyCppOMP \
  --candidate-launcher "taskset -c 0,2,4,6" \
  --candidate-env OMP_NUM_THREADS=4 \
  --suite profile --runs 6 \
  --json-output build/parallel-omp-4.json
```

Set `OMP_NUM_THREADS` explicitly for every OpenMP or hybrid report. The separate
telemetry executable inherits the candidate launcher and candidate environment,
so MPI rank counts and hybrid thread placement match the timed candidate. On
POSIX, a timeout or interrupt terminates the launcher's process group.

Compare topology reports with:

```bash
python benchmarks/check_parallel_scaling.py \
  omp:1:build/parallel-omp-1.json \
  omp:4:build/parallel-omp-4.json \
  mpi:4:build/parallel-mpi-4.json \
  hybrid:8:build/parallel-hybrid-8.json
```

Each value is `LABEL:WORKERS:PATH`. Reports must use the same suite, corpus, and
single-worker baseline configuration. The command reports paired wall-time
speedup and efficiency; add `--require-all-faster` to fail on a case at or below
`1.0`. Algorithm clock ticks are not compared across OpenMP and MPI layouts.

## Measurement notes

- Cases run serially in isolated temporary directories; launchers may add
  workers within a calculation.
- Case order rotates between rounds. Paired runs alternate AB and BA order.
- Wall time includes startup and parsing. `std::clock` starts after parsing and
  is platform-specific.
- Reports show median, median absolute deviation, p95, and raw samples. No
  outliers are removed automatically.
- JSON schema 2 stores schedules, platform data, executable and input
  fingerprints, summaries, raw samples, and optional telemetry.
- Scaling suites describe their selected workloads; they do not establish
  asymptotic complexity for arbitrary molecules.

Run the benchmark-tool tests with:

```bash
python -m unittest discover -s benchmarks -p 'test_*.py'
```

Build and run the solver-level parallel parity and telemetry checks with:

```bash
cmake --preset parallel-tests
cmake --build --preset parallel-tests
ctest --preset parallel-tests
```
