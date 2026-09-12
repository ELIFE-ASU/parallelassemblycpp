# Benchmarks

The benchmark tools run isolated ParallelAssemblyCpp calculations, verify their
assembly indices, and report wall time and program-reported `std::clock` ticks.
Plotting requires Matplotlib, included in `environment.yml`. For an existing
Conda environment, run `conda env update --file environment.yml`; for a separate
Python installation, run `python -m pip install matplotlib`.

## Quick start

From the repository root:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/ParallelAssemblyCpp \
  --suite quick
```

Use `--list-cases` to inspect the selected cases and
`python benchmarks/benchmark.py --help` for all runner options.

## ASU Sol batch job

After the [Sol environment setup](../README.md#quick-start) job succeeds,
submit from the repository root:

```bash
sbatch slurm/benchmark-sol.sbatch
```

The job builds optimized serial, OpenMP, and MPI/OpenMP hybrid executables in
its own result directory, runs the `full`, `profile`, and `scaling` suites, then
measures Paclitaxel OpenMP scaling at 2, 4, 8, 16, 32, 64, and 128 threads against
paired serial runs. It also measures hybrid layouts at the same total worker counts
from 4 through 128, trying each power-of-two MPI rank count with at least two
OpenMP threads per rank. The default includes 21 hybrid layouts; at 128 workers
these are `2x64`, `4x32`, `8x16`, `16x8`, `32x4`, and `64x2` (ranks × threads).
`full` includes every `quick` case. Workload scaling measures the amino-acid
and mask-boundary series; thread scaling holds Paclitaxel fixed while changing
the worker count and rank/thread split. Every calculation checks its expected
assembly index and disables pathway reconstruction. The default is six measured
rounds and one warm-up round per case or scaling layout, with a 600-second
timeout per calculation. These counts are for exploration, not the promotion
gate.

The allocation reserves one whole standard Sol CPU node exclusively: one task
with 128 physical cores, SMT excluded, and all node memory (`--mem=0`, normally
512 GiB). This matches
[Sol's standard compute hardware](https://docs.rc.asu.edu/supercomputer-hardware/).
The time limit is 24 hours in the `public` partition and QoS, using
[Sol's compute partitions](https://docs.rc.asu.edu/partitions-and-qos/).
Suite runs use one CPU; the thread sweep gets the full allocated CPU set via
`srun` and uses the drivers' CPU placement logic. All experiments stay on the
reserved node. Hybrid builds explicitly use `mpicxx` from the selected mamba
environment and launch with that environment's `mpirun`; no additional MPI
module is loaded.

Hybrid drivers run inside `srun --mpi=none`. Each candidate uses local Open MPI
with explicit node slots, `slot:PE=THREADS` mapping, core binding, and
oversubscription disabled. Only the MPI child environment has inherited Slurm
and MPI launcher identity variables removed, avoiding nested scheduler steps
and a dependency on Sol's site PMIx plugin. The parent retains job metadata and
all processes remain inside the allocated CPU/memory limits. Placement follows
[Open MPI's threaded-process mapping](https://docs.open-mpi.org/en/v5.0.7/man-openmpi/man1/mpirun.1.html).

Results default to `build/sol-benchmarks/<job-id>/`. Start with **`summary.csv`**:

| File | Contents |
| --- | --- |
| `summary.csv` | One row per case plus a suite-total row: MPI ranks, threads per rank, total workers, wall-time median, MAD, p95, paired speedup, and efficiency. |
| `status.txt` | `running`, `completed`, or a failure with its exit code. |
| `metadata.txt` | Job allocation, CPU hardware, Git revision/working-tree status, environment path, and run settings. |
| `suites/*.json` | Full timing samples, expected results, fingerprints, and summaries for each suite. |
| `suites/scaling.png`, `suites/scaling.pdf` | Wall time against workload size. |
| `parallel/omp-*.json` | Paired serial/OpenMP timing samples for each thread count. |
| `parallel/scaling.txt`, `parallel/scaling.png`, `parallel/scaling.pdf` | Validated thread-scaling summary and plots. |
| `parallel/cpu-topology.json` | CPU affinity, topology, and selected placement. |
| `hybrid/hybrid-RxT.json` | Paired serial/hybrid timing samples for each rank/thread layout. |
| `hybrid/scaling.txt`, `hybrid/scaling.png`, `hybrid/scaling.pdf` | Hybrid speedup and efficiency, grouped by MPI rank count. |
| `hybrid/cpu-topology.json` | Detected CPU allocation and MPI placement command for each layout. |
| `parallel-comparison.txt` | Validated comparison of OpenMP and hybrid results using the same serial reference. |
| `build.log`, `suites/*.txt`, `parallel.txt`, `hybrid.txt` | Build and runner logs. |

CSV speedup is the median of paired serial/parallel wall-time ratios;
efficiency is speedup divided by total workers (`mpi_ranks * threads_per_rank`),
expressed as a fraction. The `topology` column distinguishes `serial`, `omp`,
and `hybrid` rows; `mpi_ranks` is 1 for serial and OpenMP runs.
Values above 1 for speedup indicate improvement. Serial-only suite rows have
blank baseline, speedup, and efficiency fields. The `report` column points to
the raw JSON relative to the result directory. Rows with `case=__suite__`
summarize the per-round sum of all case wall times.

The CSV is refreshed after each suite and when the job exits, including on
ordinary failures. Completed scaling reports survive a later failed layout.
Check `status.txt` before treating a result set as complete; a forced kill may
leave it marked `running`. To regenerate the CSV from retained reports:

```bash
python benchmarks/summarize_sol.py build/sol-benchmarks/<job-id>
```

For a shorter first run on the reserved node (still including Paclitaxel):

```bash
sbatch --time=04:00:00 \
  slurm/benchmark-sol.sbatch --suites quick --threads 2,4 \
  --hybrid-layouts 2x2 --runs 2 --warmup 0
```

Customize storage, the installed environment, or the time limit:

```bash
sbatch --account=<your-account> --time=2-00:00:00 \
  slurm/benchmark-sol.sbatch \
  --env-prefix /data/your_group/envs/parallelassemblycpp \
  --output-dir /data/your_group/results/parallelassemblycpp-run-001 \
  --threads 2,4,8,16,32,64,128 \
  --hybrid-layouts 2x64,4x32,8x16,16x8,32x4,64x2 --runs 6
```

Use absolute paths for `--env-prefix`, `--output-dir`, and the optional
`--repo-dir` when submitting from outside the checkout. Existing output
directories are rejected. If `--threads` is omitted, the script selects powers
of two up to the allocated CPU count and adds that count if it is not a power
of two. Hybrid defaults derive from those selected totals; use `--hybrid-layouts`
for explicit rank/thread pairs or `--hybrid-layouts none` to run only the OpenMP
sweep. Each hybrid layout requires at least two ranks and two threads per rank,
and its total workers must fit the allocation. `--suites`, `--threads`,
`--hybrid-layouts`, `--runs`, `--warmup`, and `--timeout` follow the script path;
Slurm options such as `--account` and `--time` precede it.
Choosing fewer thread counts still reserves the whole node.
The job activates the existing mamba environment without updating it.

## Direct build shortcut

The runner's `--build` option is a direct `-O3 -DNDEBUG` x86-64-v3 shortcut.
Use the CMake presets for compiler-flag comparisons so target-specific GCC,
LTO, and PGO options are included.

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
  --executable build/performance/ParallelAssemblyCpp \
  --suite full \
  --json-output build/full.json
```

Run one custom input:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/ParallelAssemblyCpp \
  --input unitTests/sucrose.mol \
  --expected 8
```

Collect one additional, untimed telemetry run per case:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/ParallelAssemblyCpp \
  --telemetry-executable build/performance/ParallelAssemblyCppTelemetry \
  --suite scaling \
  --telemetry \
  --json-output build/scaling.json
```

Runs containing scaling cases automatically save PNG and PDF plots beside the
JSON report (`build/scaling.png` and `build/scaling.pdf` in this example), or at
those default paths when no JSON output is requested. This also applies when
selecting individual amino-acid cases with `--case`. Use `--plot-output PATH` to
choose another PNG path, with a matching `.pdf` sibling, or request plots for
another suite. The plot shows measured wall-time medians
with MAD error bars on a logarithmic axis, including the baseline for paired
runs. Amino-acid components and mask-boundary bond counts use separate panels.
Plots render without a display and are written only after successful validation.

Telemetry is excluded from timing aggregates. It records graph size, retained
masks, matching and canonicalisation activity, cache rates, phase clock ticks,
and phase memory. Parallel telemetry additionally records rank/thread topology,
distributed root-queue participation, dynamic branch leases, global root-branch
coverage, depth-two and deeper task transfers, local executions and steals,
scheduler idle waits, deep-refill activations, task-queue high-water marks,
maximum executed task depth, incumbent warm starts, steady-clock worker timing,
and all 37 raw search counters per worker plus their exact aggregate. The
parallel aggregate also reports shared-L2 hits, misses, collision-chain probe
steps, retained entry bytes, and contended shard-lock waits and wait time.
Parallel phase memory is disabled because `/proc` peak resets are process-wide.
A cache rate is `null` when no lookup occurred. Aggregate elapsed time is the
critical parallel-region wall time; worker elapsed and busy values are sums,
and each worker's busy time is its elapsed time minus measured scheduler-idle
wait time, not CPU utilization. The aggregate queue high-water mark and maximum
task depth and minimum useful task work are maxima across workers; the other
scheduler event fields are sums. Legacy VF2 counters remain in the schema but
are zero with the exact cyclic canonicaliser.

Matching scans refresh the shared best bound at class boundaries and
periodically before matching bounds. Their untimed telemetry records:

- `matching_bound_refresh_polls`: relaxed atomic loads at these scan refresh
  sites; recursion-entry refreshes are excluded.
- `matching_bound_refreshes`: polls that strictly improve the local best bound.
- `matching_bound_classes_pruned`: whole classes skipped by the immediately
  following class bound after a successful refresh.
- `matching_bound_pairs_pruned`: occurrence pairs rejected before fragmentation
  by the immediately following bound after a successful refresh.
- `matching_bound_blocks_pruned`: fragment-pair blocks rejected by the
  immediately following bound after a successful refresh.
- `matching_bound_candidates_pruned`: fragmented candidates rejected by the
  immediately following lower bound after a successful refresh, avoiding
  canonicalisation.

Each pruning counter requires that the previous local bound would have
admitted the work at that decision. These conservatively attribute immediate
avoided work to a refresh; later pruning from the improved bound is excluded.
Classes, pairs, blocks, and candidates have different units and must not be
summed as total work or interpreted as elapsed time saved. Current benchmarks
do not establish stale bounds during matching scans as a major bottleneck.
Historical schema-v1 reports without this counter group remain readable.

The aggregate and worker records expose `deeper_tasks_spawned`,
`deeper_tasks_executed`, `task_steal_attempts`, `task_steals`,
`local_task_executions`, `scheduler_idle_waits`,
`scheduler_idle_nanoseconds`, `deep_refill_activations`,
`task_queue_high_watermark`, and `maximum_task_depth_executed`. The reported
`busy_timing_method` is `elapsed_minus_scheduler_idle_time`.

Task-transfer measurements appear on each worker and in the aggregate:

- `task_serialization_nanoseconds`: steady-clock time acquiring reusable
  storage and serialising fragment metadata and mask words, before enqueueing.
- `task_execution_nanoseconds`: steady-clock time executing transferred tasks,
  including reconstruction, canonicalisation, and recursive search, before
  recycling their storage. Serialization of descendants can occur within this
  interval, so these two timings must not be added as disjoint costs.
- `tasks_immediately_pruned`: transferred tasks rejected by the initial
  incumbent bound or first transposition-table lookup before recursive search.
  Cancellation and exceptions are not counted as pruning.
- `task_buffers_created` and `task_buffers_reused`: donations that acquire new
  vector storage or reuse a returned descriptor's vector capacities. Each
  producer retains at most 32 returned buffers and 1 MiB of vector capacity;
  oversized and excess buffers are released. These counts describe buffer
  acquisitions, not allocator calls, and may include a subsequently cancelled
  donation.
- `tasks_rejected_as_too_small`: donations kept local by the measured minimum
  work threshold. They still execute within the current search.
- `task_minimum_work_units`: the largest threshold calculated by a worker,
  with the aggregate taking the maximum. Zero means no calibration window
  completed. A work unit is one possible edge pair within a fragment:
  `sum(edges * (edges - 1) / 2)`.

Every 16 completed transfers from a producer recalibrate its minimum work
threshold. The estimate combines serialization cost, execution time per work
unit for tasks that survive immediate pruning, and the fraction that survive.
A donation must offer estimated useful execution of at least eight times its
average serialization cost. A window containing only immediate prunes raises
the threshold above twice the largest observed task. Initial donations
calibrate without a minimum; every sixteenth undersized opportunity is still
transferred to refresh the estimate as the search changes. This calibration
runs in ordinary builds as well as telemetry builds. Historical schema-v1
reports without the transfer measurement group remain readable.

## Paired comparisons

Keep the previous executable and pass it as the baseline:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/ParallelAssemblyCpp-before \
  --executable build/performance/ParallelAssemblyCpp \
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
  ./build/parallel/ParallelAssemblyCppOMP molecule.mol \
    --pathway=0 --parallel=on --threads=4

mpirun --map-by slot --bind-to core -n 4 \
  ./build/parallel/ParallelAssemblyCppMPI molecule.mol \
    --pathway=0 --parallel=on --threads=1

OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close \
  mpirun --map-by slot:PE=2 --bind-to core -n 2 \
  ./build/parallel/ParallelAssemblyCppHybrid molecule.mol \
    --pathway=0 --parallel=on --threads=2

OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./build/parallel/ParallelAssemblyCppOMPTelemetry molecule.mol \
    --pathway=0 --parallel=on --threads=4 --telemetry=1

mpirun --map-by slot --bind-to core -n 4 \
  ./build/parallel/ParallelAssemblyCppMPITelemetry molecule.mol \
    --pathway=0 --parallel=on --threads=1 --telemetry=1
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

Set the positive `PARALLELASSEMBLYCPP_BRANCH_LEASE_SIZE` environment variable to use a
fixed root lease size. Only MPI rank zero writes output.

The runner accepts separate parallel policies, launchers, and environment
variables for each role:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/parallel/ParallelAssemblyCpp \
  --baseline-launcher "taskset -c 0" \
  --baseline-parallel off \
  --executable build/parallel/ParallelAssemblyCppOMP \
  --candidate-launcher "taskset -c 0,2,4,6" \
  --candidate-parallel on \
  --candidate-env OMP_NUM_THREADS=4 \
  --suite profile --runs 6 \
  --json-output build/parallel-omp-4.json
```

Set `OMP_NUM_THREADS` explicitly for every OpenMP or hybrid report. Use
`--candidate-parallel on` to require the candidate's parallel search path;
`on` fails instead of silently falling back when the requested topology cannot
be used. The selected parallel modes are recorded with the launcher and
environment in JSON reports. The separate telemetry executable inherits the
complete candidate execution configuration, so its mode, MPI rank count, and
hybrid thread placement match the timed candidate. On POSIX, a timeout or
interrupt terminates the launcher's process group.

### Paclitaxel thread sweep

`paclitaxel_scaling.py` detects the available CPUs and measures Paclitaxel at
every OpenMP thread count from 2 up to that limit. Each count gets its own
paired comparison with the same serial solver, with alternating AB/BA order
and assembly index validation against the manifest's provisional value of 23.
Only Paclitaxel runs; the other profile cases are excluded. Timings use
`--pathway=0`, excluding pathway reconstruction.

Use the parallel build above, or build just the serial and OpenMP executables:

```bash
cmake --preset performance -B build/paclitaxel \
  -DPARALLELASSEMBLYCPP_BUILD_OPENMP=ON
cmake --build build/paclitaxel --target ParallelAssemblyCpp ParallelAssemblyCppOMP
python benchmarks/paclitaxel_scaling.py \
  --build-dir build/paclitaxel \
  --output-dir build/paclitaxel-scaling
```

On Linux, detection intersects CPU topology with the driver's process affinity
and respects `SLURM_CPUS_PER_TASK` when set. It counts physical cores separately
from logical CPUs and orders one CPU from each physical core before adding
SMT siblings. For example, a host with 20 physical cores and 28 logical CPUs
runs 2, 3, 4, ..., 28 threads by default. `--physical-cores-only` stops at 20 in
that example. Use `--threads 2 4 8` to select a smaller subset. Requested counts
above the detected capacity are rejected.

Preview the detected CPU, selected counts, calculation count, and commands:

```bash
python benchmarks/paclitaxel_scaling.py --dry-run
python benchmarks/paclitaxel_scaling.py --physical-cores-only --dry-run
```

The default is six measured pairs plus one warm-up pair per thread count, or
`14 * (N - 1)` calculations for `N` available CPUs. A 28-CPU sweep makes 378
calculations. Paclitaxel can take several minutes per calculation, so a full
sweep can take many hours. `--timeout` sets the per-calculation limit in seconds
(default 600). For a smoke run:

```bash
python benchmarks/paclitaxel_scaling.py \
  --build-dir build/paclitaxel \
  --threads 2 --runs 1 --warmup 0 \
  --output-dir build/paclitaxel-smoke
```

Use an even `--runs` count for balanced timing evidence. Thread counts must be
distinct and at least two: the one-worker reference is the serial solver, and
the solver rejects forced parallel search with only one worker. The runner
sets `OMP_NUM_THREADS` and `OMP_THREAD_LIMIT` for each role, disables dynamic
teams, and uses `OMP_PROC_BIND=close`. When CPU IDs are available, each count
gets an explicit `OMP_PLACES` list with one logical CPU per place. This uses
physical cores before SMT siblings when physical topology is available,
including on CPUs with mixed core types. Otherwise it uses sorted CPU IDs.
When `taskset` is available, the serial reference is pinned to the first CPU
in that order. Otherwise it inherits the driver's affinity.

Run the driver itself inside a placement launcher or Slurm allocation so that
detection sees the intended CPUs, for example:

```bash
taskset -c 0,2,4,6 python benchmarks/paclitaxel_scaling.py --dry-run
```

The separate `--baseline-launcher` override remains available. A
`--candidate-launcher` requires explicit `--threads` and cannot be combined
with `--physical-cores-only`: the driver cannot inspect a launcher's target
allocation. In that mode, candidate CPU placement uses `OMP_PLACES=threads`
within the launcher's allocation; driver CPU IDs and capacity checks are not
applied to it. Keep serial placement fixed across the sweep.

When process affinity or physical topology is unavailable, the report states
the detection source and any unknown physical core count. CPU-count fallback
uses the process CPU count where supported, otherwise the system CPU count;
it does not invent CPU IDs. `--physical-cores-only` requires known topology.

`cpu-topology.json` records the detected topology, affinity, Slurm limit,
selected counts, baseline launcher, and CPU IDs for every planned run. Each
`omp-N.json` contains the usual schema-v2 paired samples and fingerprints.
After every count succeeds, the existing scaling checker validates the reports
and writes `scaling.txt` with CPU details, paired wall-time speedup and efficiency
(`speedup / threads`). It also saves `scaling.png` and `scaling.pdf` in the output directory,
plotting speedup and efficiency against thread count. Slowdowns are reported
without failing the benchmark; these measurements have no CI timing threshold.
Errors, incorrect assembly
indices, timeouts, and interrupts stop the sweep without writing a summary.
Existing reports are never overwritten; use a new output directory for a
repeat run. `--dry-run` prints the benchmark and checker commands without
building or running anything.

To collect a separate untimed telemetry calculation at each count, build
`ParallelAssemblyCppOMPTelemetry` in the same build directory and add `--telemetry`.
Telemetry uses the candidate's thread count and placement and is excluded from
the timing summary.

For a single four-thread comparison using the general runner:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/parallel/ParallelAssemblyCpp \
  --baseline-parallel off \
  --executable build/parallel/ParallelAssemblyCppOMP \
  --candidate-parallel on \
  --candidate-env OMP_NUM_THREADS=4 \
  --suite profile --case paclitaxel --runs 6 \
  --json-output build/parallel-paclitaxel-omp-4.json
```

### Comparing topology reports

Compare topology reports with:

```bash
python benchmarks/check_parallel_scaling.py \
  omp:2:build/parallel-omp-2.json \
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
