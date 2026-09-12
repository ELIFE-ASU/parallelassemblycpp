# parallelassemblycpp

parallelassemblycpp computes molecular assembly indices and recovers assembly
pathways. It provides a command-line tool, a reusable C++20 library, and optional
parallel search using OpenMP, MPI, or both together. The command-line tool also
computes string assembly indices for files containing one string per line.

This repository implements the algorithm described by Ian Seet, Keith Y.
Patarroyo, Gage Siebert, Sara I. Walker, and Leroy Cronin in [*Rapid Exploration
of Assembly Chemical Space of Molecular Graphs*](https://arxiv.org/abs/2410.09100).

## Quick start

Requirements: CMake 3.25 or newer, Ninja, and a C++20 compiler.

Install these tools directly, or create and activate the supplied Conda
environment:

```bash
conda env create --file environment.yml
conda activate parallelassemblycpp
```

Then configure, build, and run the release executable:

```bash
cmake --preset release
cmake --build --preset release
./build/release/ParallelAssemblyCpp unitTests/alanine.mol
```

The command writes `unitTests/alanineOut` and, by default,
`unitTests/alaninePathway`.

This runs ParallelAssemblyCpp directly from the build directory; installation is
optional.

<details>
<summary><strong>ASU Sol HPC setup</strong></summary>

From the repository root on Sol, submit the environment installation job:

```bash
sbatch slurm/install-sol.sbatch
```

The job loads Sol's `mamba/latest` module and creates
`$HOME/.conda/envs/parallelassemblycpp` from `environment.yml`. Resubmitting updates
that environment to satisfy the file. It then builds the `release` preset in
`build/sol-release` and checks `ParallelAssemblyCpp --help`. The environment
includes the compiler, CMake, Ninja, Open MPI, Python, Matplotlib, and Ruff; the
release executable uses serial search. See the development section for parallel
build presets.

The script requests one node, four CPUs, 16 GB RAM, and two hours in
`lightwork` with the `public` QoS, following ASU's guidance for
[environment creation and compilation](https://docs.rc.asu.edu/partitions-and-qos/#lightwork).
Output and errors go to `slurm-parallelassemblycpp-install-<job-id>.out` in the job's
working directory. An account can be selected with
`sbatch --account=<your-account> slurm/install-sol.sbatch`; use `myaccounts`
on Sol to list available accounts.

To choose another persistent environment location, pass an absolute prefix:

```bash
sbatch slurm/install-sol.sbatch /data/your_group/envs/parallelassemblycpp
```

An optional second argument supplies the absolute repository path when
submitting from another directory. Arguments are used because the script's
`--export=NONE` does not inherit custom variables from the submitting shell.
For example:

```bash
sbatch /path/to/parallelassemblycpp/slurm/install-sol.sbatch \
  /data/your_group/envs/parallelassemblycpp /path/to/parallelassemblycpp
```

After the setup job succeeds, activate the environment in subsequent jobs
using [ASU's supported activation syntax](https://docs.rc.asu.edu/mamba/):

```bash
module load mamba/latest
source activate "$HOME/.conda/envs/parallelassemblycpp"
./build/sol-release/ParallelAssemblyCpp unitTests/alanine.mol
```

Use your chosen prefix in `source activate` if you changed the default, and
run the executable from the repository root inside a compute allocation.
Finish jobs using the environment before resubmitting the installer, since
an update can replace their dependencies.

To run the benchmark suites and OpenMP plus MPI/OpenMP hybrid scaling sweeps on
a whole 128-core Sol node, submit `sbatch slurm/benchmark-sol.sbatch`. Results
include a combined `summary.csv`, raw JSON, and plots in
`build/sol-benchmarks/<job-id>/`. See the
[Sol benchmark job instructions](benchmarks/README.md#asu-sol-batch-job) for
resource overrides, thread counts, and output details.

</details>

<details>
<summary><strong>What the code does</strong></summary>

ParallelAssemblyCpp treats a molecule as a labelled graph: atoms are vertices and
bonds are edges. Its assembly index is the smallest number of joining steps
needed to build that graph when a fragment that has already been made can be
reused.

The program parses a V2000 MOL/SDF or native graph file, removes explicit
hydrogens by default, and enumerates connected fragments of the molecular
graph. It canonicalises those fragments so that structurally equivalent copies
can be recognised even when they use different atom or bond indices. A compact
directed acyclic graph (DAG) records the fragment relationships for reuse in
later search passes.

A branch-and-bound search then explores ways to reuse matching, disjoint
fragments. Lower bounds and a cache of previously visited canonical assembly
states eliminate branches that cannot improve the best result. The final
output contains the lowest assembly index found and, unless disabled, a JSON
description of a corresponding assembly pathway. If a runtime or enumeration
limit is reached, the reported index is the best result found so far and may
not be the proven minimum.

</details>

<details>
<summary><strong>Compared with the original AssemblyCpp v5</strong></summary>

This repository and the
[original AssemblyCpp v5](https://github.com/croningp/assemblycpp-v5) implement
the same molecular-assembly calculation. Both parse a labelled molecular graph,
enumerate connected fragments, identify structurally equivalent copies, search
possible fragment reuses with branch-and-bound, and recover a pathway for the
lowest assembly index found. This is therefore a re-engineering of the v5
implementation, not a different definition of the assembly index.

The comparison below is against the original repository's `main` branch at
[commit `f9209034`](https://github.com/croningp/assemblycpp-v5/commit/f9209034b0851d03282322bb6be697beaf030dda):

- **Graph representation and matching.** The original uses fixed 512-bit edge
  masks and Boost's VF2 implementation for cyclic graph isomorphism. This
  version uses compact, dynamically sized edge masks and in-project tree and
  cyclic canonicalisation, removing the vendored Boost dependency and the
  fixed 512-edge mask limit.
- **Search implementation.** This version adds compact DAG storage,
  frontier-driven enumeration, reusable canonical fragment identities,
  residual and transposition caches, tighter bounds, and allocation reuse.
  These changes reduce repeated graph and search work without changing the
  quantity being calculated.
- **Parallel execution.** The original solver is serial. Serial search remains
  the default here, while optional OpenMP, MPI, and hybrid builds can distribute
  independent search branches between workers and then deterministically
  reconstruct a winning pathway. This can reduce runtime for larger graphs,
  where the search exposes enough work to outweigh parallel coordination
  overhead; small graphs may see little or no speed-up.
- **Interfaces.** The original provides a C++17 `assembly` command with
  file-based results. This project provides a C++20 `ParallelAssemblyCpp` command
  plus an installable `ParallelAssemblyCpp::Library` with stream, file, and batch
  APIs that return results directly. Its command-line handling also validates
  options and reports interrupted or limited searches explicitly.
- **Project tooling.** This version expands the build and verification support
  with CMake presets, package installation and export, CI, focused and full
  regression suites, pathway and parallel-parity tests, telemetry, maintained
  benchmarks, and optional LTO and PGO builds.

When allowed to finish, both implementations target the same minimum assembly
index. A recovered pathway can differ when more than one optimal pathway
exists, and a runtime or enumeration limit can make this version return a
best-so-far result instead of a proven minimum.

</details>

<details>
<summary><strong>Installation</strong></summary>

Install ParallelAssemblyCpp when you want a standalone command, reusable library,
and CMake package outside the build directory. After the requirements above are
available, run these commands from the repository root:

```bash
cmake --preset release
cmake --build --preset release
cmake --install build/release --prefix build/install
```

The first two commands can be skipped after completing the quick start.
`build/install` is a user-writable installation prefix, so administrator access
is not required. Verify the installed command with:

```bash
./build/install/bin/ParallelAssemblyCpp --help
```

The installation contains:

- The command-line tool in `<prefix>/bin`.
- The public header in `<prefix>/include/parallelassemblycpp`.
- The static library and CMake package files in the platform's library
  directory, typically `<prefix>/lib`.

`--prefix` selects where the files are copied; it does not update `PATH`.
Replace `build/install` with another destination if needed, and add
`<prefix>/bin` to `PATH` to invoke `ParallelAssemblyCpp` from any directory.
Installing to a system location may require administrator privileges. On Windows,
the installed command is `build\install\bin\ParallelAssemblyCpp.exe`.

</details>

<details>
<summary><strong>Command line</strong></summary>

```text
ParallelAssemblyCpp INPUT [OPTIONS]
ParallelAssemblyCpp --help
```

By default, `INPUT` may be a V2000 MOL/SDF file or a ParallelAssemblyCpp native
graph file. With `--run-strings=1`, it is an exact text-file path instead. The
`.mol` and `.sdf` suffixes select MOL parsing case-insensitively; the `.mol`
suffix may be omitted when the file uses the lowercase `.mol` spelling. An
`.sdf` input reads its first V2000 structure. Native graph filenames must be
supplied in full. Options may appear before or after the input and use
`--name=value` syntax. Boolean values are `0` or `1`.

### Options

| Option | Default | Purpose |
| --- | --- | --- |
| `-h`, `--help` | — | Show command help. |
| `--runtime=<TICKS>` | Unlimited | Stop after the given `std::clock` budget. |
| `--enum-max=<COUNT>` | `50000000` | Limit retained connected masks in the initial DAG. |
| `--pathway=<0\|1>` | `1` | Write the recovered pathway. |
| `--run-strings=<0\|1>` | `0` | Treat `INPUT` as a file containing one string per line. |
| `--accept-palindromes=<0\|1>` | `0` | Identify a string fragment with its reversal. |
| `--parallel=<auto\|on\|off>` | `off` | Select parallel search automatically, require it, or disable it. |
| `--threads=<auto\|N>` | `auto` | Set the OpenMP thread count per process; `N` must be positive. |
| `--remove-hydrogens=<0\|1>` | `1` | Remove explicit hydrogens from MOL/SDF and native graph inputs. |
| `--verbose=<0\|1>` | `0` | Print the parsed input graph. |
| `--compensate-disjoint=<0\|1>` | `0` | Subtract one per processed component after the first. |
| `--memory-report=<0\|1>` | `0` | Write Linux peak virtual memory to `memUsage`. |
| `--telemetry=<0\|1>` | `0` | Write search telemetry. |
| `--write-intermediate-mas=<0\|1>` | `0` | Write each improved index and its clock tick. |

`--enum-max` includes one-edge masks. Runtime and enumeration limits return the
best index found so far, which may not be the proven minimum. Run
`ParallelAssemblyCpp --help` for full details and accepted legacy option names.
`--telemetry` is available only in telemetry-enabled executables.
`--threads=auto` treats the OpenMP runtime default (including `OMP_NUM_THREADS`)
as an upper limit. Once the prepared root jobs and DAG indicate enough work
for parallel search, it estimates one worker per 32,768 work units and rounds
the budget up to teams of 8, 16, 32, and so on. The eight-worker starting budget
leaves room for recursive work that this estimate can understate. The budget is
divided across launched MPI ranks, retaining at least one thread per rank.
Explicit thread counts apply to each process and are never reduced by this cap.

In a parallel-enabled executable, `--parallel=auto` prepares the root jobs and
DAG, then uses their estimated search work to choose parallel or serial
execution. A serial fallback reports its reason. `--parallel=on` forces parallel
search; automatic threads still use the workload cap, with at least two workers
in a single process. It fails if parallel execution cannot be honored, such as
when only one worker is available. `--parallel=off` always uses serial search.
Finite `--runtime` budgets and `--write-intermediate-mas=1` require serial
search, so `auto` reports a fallback and `on` reports an error. Pathway output
is supported with parallel optimization followed by bounded deterministic
reconstruction of a winning pathway.

### Outputs

For a MOL/SDF file, `INPUT` below excludes its recognised suffix.

| File | Contents |
| --- | --- |
| `INPUTOut` | Assembly index, search status, and `std::clock` ticks. |
| `INPUTPathway` | Recovered pathway JSON when `--pathway=1`. |
| `INPUTIntermediateMAs` | Improved indices when `--write-intermediate-mas=1`. |
| `INPUTTelemetry.json` | Search counters from a telemetry-enabled executable. |
| `memUsage` | Linux `VmPeak` value when `--memory-report=1`. |

### String assembly

Pass `--run-strings=1` to select the string algorithm. In this mode `INPUT` is
opened exactly as supplied and every line is processed as a separate string:

```bash
./build/release/ParallelAssemblyCpp strings.txt --run-strings=1
```

Results are written to `strings.txtOut`. With pathway output enabled, the
zero-based line number is included in each pathway name, such as
`strings.txt_0_Pathway`. `--accept-palindromes=1` treats a fragment and its
reversal as equivalent. String search is serial; `--parallel=on`, telemetry,
and intermediate-index output are unavailable in string mode. A finite
`--runtime` budget applies separately to each line.

The implementation is adapted from the standard-library string search in
[AssemblyCPP Public](https://gitlab.com/croningroup/public/assemblycpp-public)
at commit `2a87948`, authored by Stuart Marshall from work by Ian Seet and
Leroy Cronin. It ports the interval enumeration, Lempel-Ziv lower bound, and
recursive branch-and-bound search without importing the upstream Boost/VF2
graph implementation or vendored dependencies. The legacy spellings
`-runStrings`, `-acceptPalindromes`, and `-palindrome` remain accepted.

</details>

<details>
<summary><strong>Parallel execution</strong></summary>

The default `release` build is serial, and rejects `--parallel=on` with
`this executable was built without parallel support`. The `parallel` preset
adds the parallel executables and needs OpenMP and MPI alongside the
requirements above:

```bash
cmake --preset parallel
cmake --build --preset parallel
```

`build/parallel` then contains:

| Executable | Workers |
| --- | --- |
| `ParallelAssemblyCpp` | Serial search only. |
| `ParallelAssemblyCppOMP` | OpenMP threads in one process. |
| `ParallelAssemblyCppMPI` | One thread per MPI rank. |
| `ParallelAssemblyCppHybrid` | OpenMP threads inside each MPI rank. |

Each parallel executable has a `...Telemetry` sibling that also accepts
`--telemetry=1`.

### OpenMP

Run the OpenMP executable directly and ask for a thread count:

```bash
OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./build/parallel/ParallelAssemblyCppOMP benchmarks/inputs/paclitaxel.mol \
    --parallel=on --threads=8
```

Parallel search does not change the outputs: this writes
`benchmarks/inputs/paclitaxelOut` and `benchmarks/inputs/paclitaxelPathway`
exactly as a serial run would. An explicit `--threads` count applies to the
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

Alanine is far too small to repay coordination. `--parallel=on` turns that same
condition into an error, which is what a scaling script wants when a serial run
would be measured by mistake.

### MPI

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
executable is built without OpenMP, so ranks are its only workers and
`--threads` above `1` is rejected. Ranks claim disjoint chunks of root work
from a queue on rank zero, so a rank that finishes early asks for more instead
of waiting.

### Hybrid MPI and OpenMP

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

### Notes

- Parallel search optimizes the index first, then deterministically
  reconstructs a winning pathway, so `--pathway=1` still works. Add
  `--pathway=0` to skip reconstruction when only the index is wanted.
- Finite `--runtime` budgets and `--write-intermediate-mas=1` require serial
  search, and string mode (`--run-strings=1`) is serial throughout.
- Speed-up depends on how much search a graph exposes, so measure rather than
  assume. [benchmarks/README.md](benchmarks/README.md#parallel-scaling) covers
  paired serial/parallel comparisons, thread sweeps, and parallel telemetry.

</details>

<details>
<summary><strong>C++ library</strong></summary>

Installed packages export `ParallelAssemblyCpp::Library`. If ParallelAssemblyCpp
is installed to a non-system prefix, pass that prefix when configuring the
consuming project:

```bash
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH=/absolute/path/to/parallelassemblycpp/build/install
```

Then link the imported target in the consuming project's `CMakeLists.txt`:

```cmake
find_package(ParallelAssemblyCpp 0.1.0 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE ParallelAssemblyCpp::Library)
```

```cpp
#include <parallelassemblycpp.h>

#include <iostream>

int main()
{
    const auto result = parallelassemblycpp::calculate("molecule.mol");
    if (!result)
    {
        std::cerr << result.error << '\n';
        return 1;
    }
    std::cout << result.assemblyIndex << '\n';
}
```

`calculateMolfile` accepts a V2000 molfile stream, while `calculateGraph`
accepts a ParallelAssemblyCpp native graph stream. `calculateBatch` processes
several inputs sequentially without process startup between items. Library calls
do not create output files. Search state is process-global, so the API is
reusable but not thread-safe; use separate processes for concurrent work.

</details>

<details>
<summary><strong>Development</strong></summary>

Tests and benchmark tooling require Python 3.10 or newer. Parallel builds also
require OpenMP and MPI; PGO builds require GCC. The supplied Conda environment
includes the development dependencies.

### Tests

Run the focused developer suite:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the `ci` preset for the full regression suite and `parallel-tests` for
parallel parity and telemetry checks. See
[unitTests/README.md](unitTests/README.md) for targeted commands and fixture
details.

### Quality gates

All C++ targets compile with high-signal warnings treated as errors by default.
The `PARALLELASSEMBLYCPP_STRICT_WARNINGS` CMake option exists for toolchain diagnosis,
but changes should pass with it enabled. Check Python lint and formatting with:

```bash
ruff check .
ruff format --check .
python tools/check_repository.py
```

C++ variables, parameters, and data members use descriptive `lowerCamelCase`
names, while C++ macros use `UPPER_SNAKE_CASE`. Python follows PEP 8:
`snake_case` names, `UPPER_SNAKE_CASE` constants, single leading underscores
for private or intentionally unused names, and protocol-required double
underscores. Project CMake variables use the `PARALLELASSEMBLYCPP_UPPER_SNAKE_CASE`
prefix. C++ project-defined identifiers avoid leading underscores and
unexplained abbreviations. CI enforces the compiler and Python quality gates.

### Benchmarks

Build the optimized candidate, then run a maintained suite:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/ParallelAssemblyCpp \
  --suite quick
```

The `performance` preset targets x86-64-v3. See
[benchmarks/README.md](benchmarks/README.md) for the benchmark corpus, paired
comparisons, parallel builds, telemetry, LTO, PGO, and scaling guidance. The
[Paclitaxel thread sweep](benchmarks/README.md#paclitaxel-thread-sweep) automates
paired serial/OpenMP measurements across every available CPU count, with
automatic topology detection and an optional physical-core-only sweep.

### Packaging

These commands create binary and source archives for release distribution;
they are not required for a normal installation:

```bash
cpack --preset release
cmake --build --preset release --target package_source
```

Use CMake 4.3 or newer to normalize archive ownership. Create source archives
from a clean checkout because CPack includes the working tree.

</details>

## License

ParallelAssemblyCpp is licensed under
[Creative Commons Attribution-NonCommercial 4.0 International](License.md).
