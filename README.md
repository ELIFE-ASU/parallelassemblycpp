# parallelassemblycpp

Compute molecular and string assembly indices, recover assembly pathways, or
obtain fast constructive Re-Pair upper bounds. ParallelAssemblyCpp provides a
command-line tool and an installable C++20 library, with optional OpenMP, MPI,
and hybrid parallel search in the command-line executables.

**Exact search is the default.** A completed full search proves the minimum;
a runtime limit, enumeration limit, or interruption leaves a best-so-far result.
Re-Pair mode returns a heuristic upper bound and does not prove minimality.

[Command-line reference](docs/cli.md) · [Parallel execution](docs/parallel.md) ·
[C++ library](docs/library.md) · [Benchmarks](benchmarks/README.md) ·
[Tests](unitTests/README.md)

## Quick start

Requirements: **CMake 3.25 or newer, Ninja, and a C++20 compiler**. The serial
release build does not require Python, OpenMP, or MPI. CI covers GCC and Clang
on Linux and MSVC on Windows; the minimum-CMake job builds and installs with
CMake 3.25.0. The shell examples below use Bash. For native Windows builds,
see [installation](docs/installation.md#windows).

Clone the repository, or start in an existing checkout:

```bash
git clone https://github.com/ELIFE-ASU/parallelassemblycpp.git
cd parallelassemblycpp
```

Install the build tools directly, or use the supplied Conda environment, which
also includes the testing, benchmarking, MPI, and packaging dependencies:

```bash
conda env create --file environment.yml
conda activate parallelassemblycpp
```

Build and run:

```bash
cmake --preset release
cmake --build --preset release
./build/release/ParallelAssemblyCpp unitTests/alanine.mol
```

This writes `unitTests/alanineOut` and `unitTests/alaninePathway`. Add
`--pathway=0` when only the index is needed. Installation is optional.

ASU Sol users can submit `sbatch slurm/install-sol.sbatch` from the checkout;
see the [Sol setup and activation guide](docs/sol.md).

## Command line

```text
ParallelAssemblyCpp INPUT [OPTIONS]
ParallelAssemblyCpp [OPTIONS] -- INPUT
ParallelAssemblyCpp --help
```

| Calculation | Example |
| --- | --- |
| Exact molecular assembly | `ParallelAssemblyCpp molecule.mol` |
| Molecular Re-Pair upper bound | `ParallelAssemblyCpp molecule.mol --algorithm=re-pair` |
| Exact string assembly, one UTF-8 string per line | `ParallelAssemblyCpp strings.txt --run-strings=1` |
| String Re-Pair upper bound | `ParallelAssemblyCpp strings.txt --run-strings=1 --algorithm=re-pair` |
| Exact string assembly with reversal equivalence | `ParallelAssemblyCpp strings.txt --run-strings=1 --accept-palindromes=1` |

Use the executable's build path until it is installed on `PATH`. Options use
`--name=value`; booleans are `0` or `1`.

Molecular inputs may be V2000 MOL/SDF files or the
[five-line native graph format](docs/cli.md#native-graph-format). Only the first
SDF record is read. Explicit hydrogens are removed by default; charges,
isotopes, coordinates, and stereochemistry are not graph labels. String inputs
use Unicode code points without normalization.

Results go beside the input: `INPUTOut` contains the index and timing, and
`INPUTPathway` contains graph pathway JSON. MOL/SDF suffixes are removed from
the output stem. String pathways are named `INPUT_0_Pathway`, `INPUT_1_Pathway`,
and so on. Re-running the same input overwrites these outputs. A successful
exit alone does not prove minimality: check the result's status for limits,
interruption, or heuristic mode. See the [output and exit-code reference](docs/cli.md#outputs).

The [full CLI guide](docs/cli.md) covers options, limits, input validation,
telemetry, both Re-Pair certificate formats, and string behavior.

## Installation

After building the release preset:

```bash
cmake --install build/release --prefix build/install
./build/install/bin/ParallelAssemblyCpp --help
```

The installation includes the command, static library, public header, CMake
package, and user guides. Add `<prefix>/bin` to `PATH` to use the command outside
the checkout. See [installation details](docs/installation.md) for paths and
Windows instructions.

## Parallel execution

Serial search is the default in every executable. The `parallel` preset adds
`ParallelAssemblyCppOMP`, `ParallelAssemblyCppMPI`, and
`ParallelAssemblyCppHybrid`, plus telemetry variants. It requires OpenMP and
MPI and **targets x86-64-v3**. Use the portable override in the
[parallel build guide](docs/parallel.md) for other CPUs.

```bash
cmake --preset parallel
cmake --build --preset parallel
./build/parallel/ParallelAssemblyCppOMP benchmarks/inputs/paclitaxel.mol \
  --parallel=on --threads=4 --pathway=0
```

Full graph and string searches support parallel execution. Re-Pair-only runs
are serial. `--parallel=auto` may fall back to serial; `--parallel=on` reports
an error if the selected build or options cannot honor it. Speed-up depends on
the workload. See [OpenMP, MPI, and hybrid usage](docs/parallel.md) and the
[scaling benchmarks](benchmarks/README.md#parallel-scaling).

## C++ library

Installed packages export `ParallelAssemblyCpp::Library`:

```cmake
find_package(ParallelAssemblyCpp 0.1.0 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE ParallelAssemblyCpp::Library)
```

Include `<parallelassemblycpp.h>` to call `calculate`, `calculateMolfile`,
`calculateGraph`, `calculateBatch`, `calculateString`, or
`calculateStringBatch`. These return indices and status fields without creating
output files; the library does not return pathway JSON. Calls are serial and
use process-global state, so use separate processes for concurrent work.
See the [library guide](docs/library.md) for a complete consumer example,
options, and result semantics.

## Development

Python 3.10 or newer and Matplotlib are needed for the complete test suite;
both are included in the supplied Conda environment. Start with:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

| Guide | Contents |
| --- | --- |
| [Development and packaging](docs/development.md) | Presets, dependencies, quality checks, archives, and contribution guidance |
| [Tests](unitTests/README.md) | Focused/full regressions, parallel parity, and fixture maintenance |
| [Benchmarks](benchmarks/README.md) | Corpus, telemetry, paired comparisons, scaling, LTO, PGO, and Sol jobs |
| [ASU Sol](docs/sol.md) | Environment installation and job activation |
| [Algorithm and provenance](docs/algorithm.md) | Molecular search and comparison with AssemblyCpp v5 |
| [Dated audits](audits/README.md) | Historical measurements, validation evidence, and research proposals |

Report reproducible problems through
[GitHub issues](https://github.com/ELIFE-ASU/parallelassemblycpp/issues), including
the commit, build preset/compiler, command, input, and observed result.

## References

The molecular algorithm is described by Ian Seet, Keith Y. Patarroyo, Gage
Siebert, Sara I. Walker, and Leroy Cronin in
[*Rapid Exploration of Assembly Chemical Space of Molecular Graphs*](https://arxiv.org/abs/2410.09100).
The [provenance guide](docs/algorithm.md) describes this implementation's
relationship to AssemblyCpp v5; the [string guide](docs/cli.md#string-assembly)
records its upstream attribution. For reproducible research, record the
repository commit, build options, input handling, and whether results are
proven minima or upper bounds.

## License

ParallelAssemblyCpp is licensed under
[Creative Commons Attribution-NonCommercial 4.0 International](License.md).
