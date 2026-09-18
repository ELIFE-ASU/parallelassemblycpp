# Benchmarks

The [search locality experiment](search-locality-results.md) records the
8/28-worker comparison, incumbent trajectories, cache retention, and validation.

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
and all raw search counters per worker plus their exact aggregate. The
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

MPI refill measurements appear on each worker and in the aggregate:

- `mpi_refill_requests`: rank-level refill requests issued, including replies
  that report an exhausted global queue.
- `mpi_prefetched_refills`: refill requests issued while the rank still has
  local root work, so communication can overlap computation.
- `mpi_refill_replies`: refill exchanges completed, including empty replies
  and exchanges drained during termination. Completed searches have matching
  request and reply counts.
- `mpi_refill_wait_nanoseconds`: sampled steady-clock intervals with no
  leaseable local root slots while a refill is incomplete. These intervals
  can include computation on already leased roots and are not total worker
  idle time; `scheduler_idle_nanoseconds` measures worker idle waits.
- `mpi_progress_calls`: entries to the MPI controller's progress routine.
- `mpi_maximum_progress_gap_nanoseconds`: longest interval between the
  starts of controller progress calls, measured from search start through
  global completion. The interval includes computation, communication, and
  wait time; it does not measure only time outside MPI or prove that MPI made
  network progress during every call.
- `mpi_pending_refills_high_watermark`: largest number of simultaneously
  pending refill exchanges on a rank, bounded by one.

Each rank contributes its controller measurements through local worker zero
after global completion and request draining. Other local workers, and all
OpenMP workers, report zero for these MPI fields. The aggregate sums the
request, reply, prefetch, wait-time, and progress-call fields and takes the
maximum of the progress gap and pending-refill high-water mark across ranks.

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

## Shared transposition-cache experiments

`shared_cache_experiment.py` compares a retained baseline OpenMP executable with
four candidate configurations: shared lookup with unrestricted admission,
local-only lookup, shared lookup bounded to 256 MiB of entry/key storage, and
selective admission with the same bound. This separates storage-layout changes
from admission policy. Keep both executables and their instrumented siblings
from the respective builds:

```bash
python benchmarks/shared_cache_experiment.py \
  --baseline-executable build/shared-cache-baseline/ParallelAssemblyCppOMP \
  --baseline-telemetry-executable build/shared-cache-baseline/ParallelAssemblyCppOMPTelemetry \
  --executable build/shared-cache-candidate/ParallelAssemblyCppOMP \
  --telemetry-executable build/shared-cache-candidate/ParallelAssemblyCppOMPTelemetry \
  --threads 4 --runs 6 --warmup 1 \
  --output-dir build/shared-cache-experiment
```

The default runs each manifest case once per round, including cases belonging
to multiple suites. `--suite` or repeated `--case` options narrow the corpus.
Add `--dry-run` to inspect the selected cases, execution settings, and total
calculation count without creating files or requiring built executables.
`--parallel auto` is the default, preserving the solver's small-work fallback;
use `--parallel on` to require parallel search. All modes use the same thread
count and placement. Set `--launcher 'taskset -c 0-3'` to pin the experiment to
four CPUs available in the current allocation, or use the desired allocation
and override `OMP_PLACES` with `--env 'OMP_PLACES={0},{2},{4},{6}'`.

Timings reuse `benchmark.py`'s adjacent AB/BA pairs, rotating case order,
isolated temporary inputs, assembly-index checks, and executable/corpus
fingerprints. An even round count is required. Each variant has its own paired
schema-v2 JSON report. The separate `profiles.json` records one extra validated
run of every baseline and candidate case, with full telemetry and GNU time
peak RSS in KiB. The telemetry executable is optional; without it, the extra
run measures RSS using the ordinary executable. These profiles never enter
timing aggregates. RSS from an instrumented executable includes telemetry
overhead; for OpenMP it describes the whole process, including every worker.
GNU time does not report summed peak RSS across MPI ranks, so this driver is
intended for process-local OpenMP cache comparisons.

`summary.csv` combines paired wall/clock speedups, RSS, and every scalar shared
cache counter, including useful prunes, admissions/rejections, entry/key bytes,
arena and slot bytes, lock waits, and growth counts/times under the shard locks.
`profiles.json` has `complete: true` only after all requested comparisons finish.
Existing output directories are rejected. `--require-all-faster` exits with
failure when any case/variant wall or clock median ratio is at or below 1.0;
all raw samples and profiles remain available to inspect regressions. These
exploratory comparisons do not replace the four-suite promotion gate.

Use repeated `--variant NAME=POLICY[:BYTES]` options for a different selection:

```bash
--variant shared=shared --variant local=local \
--variant bounded=shared:268435456 --variant selective=selective:268435456
```

The driver sets `PARALLELASSEMBLYCPP_SHARED_CACHE_POLICY` and
`PARALLELASSEMBLYCPP_SHARED_CACHE_BYTES` explicitly; zero bytes means unlimited
entry/key admission. The cap excludes table slots, arena slack, the admission
filter, and worker L1 caches; compare those counters and process RSS as well.
The cap is divided equally across the 64 shards. Once a shard reaches its
budget, existing keys still support exact lookup and score improvement.
Selective admission uses a fixed-size, replaceable fingerprint filter: the
first recent sighting skips storage and a repeat permits admission. Fingerprint
collisions can admit extra keys, but pruning always compares the complete key.
Local mode keeps the canonical-ID registry and L1 behavior unchanged while
bypassing L2 lookups. The default solver policy remains unrestricted sharing.
`--env KEY=VALUE` applies to both roles, `--baseline-env KEY=VALUE` overrides only
the reference, and `--variant-env NAME:KEY=VALUE` overrides one candidate.
Evaluate runtime and pruning alongside hit rate: even an infrequent hit may
avoid a costly subtree.

### Index reservation and critical-section profiling

L2 keeps four-byte hash metadata contiguous and loads a key pointer only after
the metadata matches. The final comparison still checks the entire key and its
length. Expanded pointer arrays are allocated without a redundant zero fill:
hash metadata identifies occupied slots for lookup, growth, and destruction.
Telemetry reports `metadata_reject_count` and `key_comparison_count` separately
from `collision_chain_steps`, so probe counts need not be mistaken for pointer
loads.

`PARALLELASSEMBLYCPP_SHARED_CACHE_RESERVE_BYTES` optionally reserves additional
hash/pointer arrays before workers start. The default is zero. The budget is
split over 64 shards and capacities round down to powers of two. Reservations
no larger than the inline capacity do nothing; local-only policy also skips
reservation. The budget excludes the existing inline arrays and does not cap
later growth, entry/key admission, or arena slack. On a 64-bit build, 50,331,648
bytes reserves 65,536 slots per shard. Preallocation is included in wall time
and RSS; moving work out of a lock does not necessarily make the search faster.

Compare the same executable with and without reservation:

```bash
python benchmarks/shared_cache_experiment.py \
  --baseline-executable build/parallel/ParallelAssemblyCppOMP \
  --executable build/parallel/ParallelAssemblyCppOMP \
  --baseline-telemetry-executable build/parallel/ParallelAssemblyCppOMPTelemetry \
  --telemetry-executable build/parallel/ParallelAssemblyCppOMPTelemetry \
  --variant reserve48=shared \
  --variant-env reserve48:PARALLELASSEMBLYCPP_SHARED_CACHE_RESERVE_BYTES=50331648 \
  --threads 8 --runs 6 --warmup 1 \
  --output-dir build/reserve-experiment
```

The driver explicitly defaults both roles' reserve budgets to zero, protecting
comparisons from inherited settings. Use ordinary executables for timings and
separate telemetry executables for profiling. Lock-hold durations cover lookup,
admission, key allocation, and growth under the shard lock, including exception
exits. Growth includes index allocation, rehash, and old-index release; rehash
times only the relocation loop. Arena-refill durations measure worker arenas'
upstream allocation attempts, including failures, rather than every key
allocation.
Each duration has a total and maximum. These intervals overlap and totals sum
across workers; they do not measure CPU utilization or critical-path savings.
Instrumented timings include timer overhead and OS preemption.

Sharing immutable keys with L1 remains a separate ownership change. L1 currently
copies a missing key before consulting L2, so it would need a lookup/commit split
and a common immutable key representation, while keeping scores private. Any
borrowed L2 key must outlive every referring L1; rejected admissions still need
L1-owned storage. Index eviction alone can preserve search-lifetime arena keys,
but reclaiming their memory would require pins, reference counts, or retiring
whole generations only after all referring L1 caches are cleared.

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
simply issue more requests. A low watermark starts one asynchronous pending
refill while local root work remains. Replies carry the latest incumbent
bound, and completed work requests are drained before global termination.
All MPI progress stays on the initializing thread, preserving
`MPI_THREAD_FUNNELED`. Serial/OpenMP scheduling retains guided leases for
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
search failure. Ranks waiting for global completion continue servicing refill
requests and RMA, sleeping for 1 ms between progress calls while the completion
barrier remains pending. Task transfer is disabled in an MPI-only rank with
one local worker, but remains available within hybrid ranks.

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

With `--telemetry`, both the OpenMP sweep and `hybrid_scaling.py` also write
`search-profiles.json` and add search-work columns to `scaling.txt`. The hybrid
driver uses `ParallelAssemblyCppHybridTelemetry` and applies the same MPI launcher,
rank count, thread count, and binding as the timed candidate. Each placement's
profile records the full incumbent trajectory, expanded/pruned-state counts,
canonical-mask hits and misses, and retained bytes for worker-local caches.
The underlying paired report retains every worker's counters and trajectory.

For example, collect the requested 8-to-28-worker comparison on an allocation
with at least 28 available physical cores:

```bash
python benchmarks/paclitaxel_scaling.py \
  --build-dir build/parallel --threads 8 28 --physical-cores-only \
  --runs 6 --warmup 1 --telemetry --output-dir build/paclitaxel-work-profile
python benchmarks/hybrid_scaling.py \
  --build-dir build/parallel --cpus 28 --layouts 2x4 2x14 4x7 \
  --runs 6 --warmup 1 --telemetry --output-dir build/paclitaxel-hybrid-work-profile
```

`incumbent_trajectory` records strict improvements as
`{elapsed_nanoseconds, assembly_index, rank, global_worker_index}`. Its clock is
steady wall time since search initialization, including initial enumeration;
MPI rank starts are barrier aligned. The combined trajectory describes discovery
of feasible incumbents, so its last timestamp is when the final incumbent was
first found, rather than when optimality was proved or every MPI rank received
the new bound. `trajectory_index=internal_before_disjoint_compensation` identifies
internal search indices; for disconnected inputs these can differ from the
compensated assembly index in solver output. Compare timestamps with the
profile run duration to distinguish discovery time from proof work. Worker
trajectories may be empty when a worker publishes no improvement.

`states_expanded` counts materialized search states reaching expansion;
`states_pruned` counts materialized states rejected by the incumbent bound or
assembly-state cache, and `states_bound_pruned` identifies the bound subset.
`duplicate_classes_pruned`, `occurrence_pairs_pruned`, and
`fragment_pair_blocks_pruned` count avoided branches at different loop levels;
they must not be added to state counts. A rise in matching visits or canonical
misses alone cannot distinguish search order, pruning, and cache replication.

`local_caches` reports retained capacity at worker completion, with separate
canonical-mask, canonical-graph, canonical-tree, assembly-state, and residual
decomposition estimates. Totals sum across workers and MPI ranks. These are
capacity estimates excluding allocator overhead, shared storage, and scratch
space; they are neither process RSS nor a simultaneously sampled memory peak.
Compare them with the separate RSS profile from `shared_cache_experiment.py`
before increasing cache capacity. That experiment's CSV now includes state
counts, canonical-mask misses, final-incumbent time, and retained-byte columns
for both candidate and baseline profiles.

The local mask cache defaults to `ASSEMBLY_CANONICAL_MASK_CACHE=flat`;
`ASSEMBLY_CANONICAL_MASK_CACHE=unordered` selects the previous backend for
comparisons. The flat backend stores one- and two-word mask keys directly,
while each cache generation initially uses the unordered map. On the 512th
distinct mask, eligible domains promote all entries into a flat table whose
first allocation is 1,024 slots (24 KiB). Smaller working sets stay unordered;
wider masks always use the unordered path. Clearing the cache resets activation
while retaining any allocated flat capacity for reuse. The flat table grows
at 70% load up to 65,536 slots, with an unordered overflow table retaining
additional entries. Only the flat allocation is capped; total cache storage
retains the existing unbounded semantics because DAG construction needs every
admitted mask. Exact canonicalization is unchanged.
For a controlled same-binary comparison with raw timings and separate profiles:

```bash
python benchmarks/shared_cache_experiment.py \
  --baseline-executable build/parallel/ParallelAssemblyCppOMP \
  --executable build/parallel/ParallelAssemblyCppOMP \
  --baseline-telemetry-executable build/parallel/ParallelAssemblyCppOMPTelemetry \
  --telemetry-executable build/parallel/ParallelAssemblyCppOMPTelemetry \
  --threads 8 --parallel on --suite profile --case paclitaxel \
  --baseline-env ASSEMBLY_CANONICAL_MASK_CACHE=unordered \
  --variant flat=shared --variant-env flat:ASSEMBLY_CANONICAL_MASK_CACHE=flat \
  --runs 6 --warmup 1 --output-dir build/paclitaxel-flat-mask-comparison
```

The other experiment switches are `PARALLELASSEMBLYCPP_CHILD_FIRST=1`, which
keeps a promising child on the current worker while donating siblings, and
`PARALLELASSEMBLYCPP_GREEDY_BOOTSTRAP=1`, which tries a greedy initial solution
within fixed limits of 32 steps, 65,536 jobs, and 262,144 parent-fragment checks.
Both default to `0`. `PARALLELASSEMBLYCPP_MATCHING_REFRESH=0` disables the
existing batch-boundary incumbent refresh for an A/B comparison; its default
is `1`. Pass these through `--baseline-env` and `--variant-env` in the experiment
driver, changing one switch at a time. Use the incumbent trajectory to decide
whether bootstrap work is justified, then compare timings and total expanded
states with child-first donation. No switch adds a synchronization operation
to every matching iteration.

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

To cover every registered test, including the full serial regression manifest,
use an unfiltered CTest invocation. The `parallel-tests` test preset itself
selects only tests labelled `parallel`:

```bash
cmake --preset parallel-tests -DPARALLELASSEMBLYCPP_FULL_REGRESSION_TESTS=ON
cmake --build --preset parallel-tests
ctest --test-dir build/parallel-tests --output-on-failure
```

Omitting `--suite` selects every unique benchmark case once (37 cases across
`full`, `profile`, and `scaling`; `quick` is contained in `full`). A complete
paired benchmark check, with no timing gate on short cases, is:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/before/ParallelAssemblyCpp \
  --executable build/parallel/ParallelAssemblyCpp \
  --runs 6 --warmup 1 --timeout 600 \
  --json-output build/all-cases-comparison.json
```

Inspect individual longer cases and paired wall-time ratios as well as the suite
total. Use a fresh report path and an idle machine; telemetry runs are separate
from timings. Two rounds with no warm-up provide a fast correctness/performance
smoke check, while the repeated default above is more useful for judging gains.
