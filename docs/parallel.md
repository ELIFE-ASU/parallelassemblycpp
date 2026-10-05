# Parallel execution

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Run source-tree examples from the repository root unless stated otherwise.

The default `release` build is serial, and rejects `--parallel=on` with
`this executable was built without parallel support`. The `parallel` preset
adds the parallel executables and needs CMake 3.25 or newer, Ninja, a C++20
compiler, OpenMP, and MPI. It inherits the `performance` preset and therefore
targets **x86-64-v3**:

```bash
cmake --preset parallel
cmake --build --preset parallel
```

For a portable CPU target, override that optimization at configuration time:

```bash
cmake --preset parallel -DPARALLELASSEMBLYCPP_X86_64_V3=OFF
cmake --build --preset parallel
```

For OpenMP alone, without an MPI dependency, use a separate build directory:

```bash
cmake -S . -B build/openmp -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF -DPARALLELASSEMBLYCPP_BUILD_OPENMP=ON
cmake --build build/openmp
./build/openmp/ParallelAssemblyCppOMP unitTests/alanine.mol --parallel=auto
```

Every executable defaults to `--parallel=off`; choose `auto` or `on` to use
workers. These build options affect command-line targets; the installed library
remains serial. The [CLI guide](cli.md#options) describes automatic worker counts
and the restrictions imposed by finite runtime budgets and Re-Pair mode.

`build/parallel` then contains:

| Executable | Workers |
| --- | --- |
| `ParallelAssemblyCpp` | Serial search only. |
| `ParallelAssemblyCppOMP` | OpenMP threads in one process. |
| `ParallelAssemblyCppMPI` | One thread per MPI rank. |
| `ParallelAssemblyCppHybrid` | OpenMP threads inside each MPI rank. |

Each parallel executable has a `...Telemetry` sibling that also accepts
`--telemetry=1`.

For molecule and graph runs, telemetry reports the serial search for a
deterministic optimal pathway separately in `pathway_reconstruction`.
`elapsed_seconds` includes worker setup, witness search, cleanup, and pathway
file output; `cpu_seconds` uses the process CPU clock. Its search `counters`
are separate from the parallel optimization totals. `attempted` and
`completed` distinguish skipped reconstruction from a failed attempt.
Serial searches build their witness during optimization, so their separate
reconstruction timing is zero. Use the benchmark runner's `--pathways` option
to include pathways in end-to-end measurements (see
[benchmark instructions](../benchmarks/README.md)).

## OpenMP

Run the OpenMP executable directly and ask for a thread count:

```bash
OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./build/parallel/ParallelAssemblyCppOMP benchmarks/inputs/paclitaxel.mol \
    --parallel=on --threads=8
```

This writes `benchmarks/inputs/paclitaxelOut` and
`benchmarks/inputs/paclitaxelPathway`, using the same filenames as serial search.
Completed searches return the same index and pathway; timings differ.
An explicit `--threads` count applies to the
process as given, while `OMP_PLACES` and `OMP_PROC_BIND` pin the threads to
distinct cores so each worker keeps its caches local.

`--parallel=auto` chooses after preparing the root jobs and DAG, and explains a
serial choice on standard error:

```bash
./build/parallel/ParallelAssemblyCppOMP unitTests/alanine.mol --parallel=auto
```

```text
parallel: serial fallback: estimated work 0 (0 root jobs x 3 retained DAG nodes) is below 32768
```

Alanine offers too little search work to benefit from parallel coordination.
`--parallel=on` bypasses this workload threshold and runs the parallel search,
provided at least two workers are available and the selected options permit it.

## MPI

Launch the MPI executable with one rank per worker:

```bash
mpirun --map-by slot --bind-to core -n 8 \
  ./build/parallel/ParallelAssemblyCppMPI benchmarks/inputs/paclitaxel.mol \
    --parallel=on --threads=1
```

Every rank parses the same command line, compares its options with the others
during start-up, and the run stops with `MPI ranks received different
command-line options` if they disagree, so pass identical arguments to all of
them. Rank zero writes the output files and any diagnostic message. This
executable is built without OpenMP, so ranks are its only workers. Use
`--threads=1` or `auto`: a value above `1` makes `--parallel=on` fail and
`--parallel=auto` fall back to serial. With `--parallel=off`, the thread setting
is ignored. For molecular/native graph search, ranks claim disjoint chunks of
root work from a queue on rank zero, so a rank that finishes early asks for
more. String search uses fixed partitions of root branches across ranks.

## Hybrid MPI and OpenMP

The hybrid executable uses both: threads within a rank, ranks across sockets or
nodes. Give each rank the cores its threads need:

```bash
OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
  mpirun --map-by slot:PE=4 --bind-to core -n 4 \
  ./build/parallel/ParallelAssemblyCppHybrid benchmarks/inputs/paclitaxel.mol \
    --parallel=on --threads=4
```

That is sixteen workers as four ranks of four threads. With `--threads=auto`
the automatic budget is divided across the launched ranks, keeping at least one
thread per rank. These placement flags are Open MPI syntax; other launchers
spell them differently. Under a scheduler, launch with the `mpirun` from the
MPI installation the executable was built against: `slurm/benchmark-sol.sbatch`
runs the Conda environment's own `mpirun` inside a single-task
`srun --mpi=none` step.

## Notes

- Parallel search optimizes the index first, then deterministically
  reconstructs a winning pathway, so `--pathway=1` still works. Add
  `--pathway=0` to skip reconstruction when only the index is wanted.
- Finite `--runtime` budgets and `--write-intermediate-mas=1` require serial
  search. String mode supports parallel branch search but does not support
  telemetry or intermediate-index output.
- Speed-up depends on how much search a graph exposes, so measure rather than
  assume. [benchmarks/README.md](../benchmarks/README.md#parallel-scaling) covers
  paired serial/parallel comparisons, thread sweeps, and parallel telemetry.
