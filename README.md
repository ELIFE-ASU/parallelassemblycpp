# AssemblyCpp v5

AssemblyCpp computes molecular assembly indices and recovers assembly pathways.
It provides a command-line tool, a reusable C++20 library, and optional parallel
search using OpenMP, MPI, or both together. The command-line tool also computes
string assembly indices for files containing one string per line.

This repository implements the algorithm described by Ian Seet, Keith Y.
Patarroyo, Gage Siebert, Sara I. Walker, and Leroy Cronin in [*Rapid Exploration
of Assembly Chemical Space of Molecular Graphs*](https://arxiv.org/abs/2410.09100).

## Quick start

Requirements: CMake 3.25 or newer, Ninja, and a C++20 compiler.

Install these tools directly, or create and activate the supplied Conda
environment:

```bash
conda env create --file environment.yml
conda activate assemblycpp-v5
```

Then configure, build, and run the release executable:

```bash
cmake --preset release
cmake --build --preset release
./build/release/AssemblyCpp unitTests/alanine.mol
```

The command writes `unitTests/alanineOut` and, by default,
`unitTests/alaninePathway`.

This runs AssemblyCpp directly from the build directory; installation is
optional.

<details>
<summary><strong>What the code does</strong></summary>

AssemblyCpp treats a molecule as a labelled graph: atoms are vertices and
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
  file-based results. This project provides a C++20 `AssemblyCpp` command plus
  an installable `AssemblyCpp::Library` with stream, file, and batch APIs that
  return results directly. Its command-line handling also validates options
  and reports interrupted or limited searches explicitly.
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

Install AssemblyCpp when you want a standalone command, reusable library, and
CMake package outside the build directory. After the requirements above are
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
./build/install/bin/AssemblyCpp --help
```

The installation contains:

- The command-line tool in `<prefix>/bin`.
- The public header in `<prefix>/include/assemblycpp`.
- The static library and CMake package files in the platform's library
  directory, typically `<prefix>/lib`.

`--prefix` selects where the files are copied; it does not update `PATH`.
Replace `build/install` with another destination if needed, and add
`<prefix>/bin` to `PATH` to invoke `AssemblyCpp` from any directory. Installing
to a system location may require administrator privileges. On Windows, the
installed command is `build\install\bin\AssemblyCpp.exe`.

</details>

<details>
<summary><strong>Command line</strong></summary>

```text
AssemblyCpp INPUT [OPTIONS]
AssemblyCpp --help
```

By default, `INPUT` may be a V2000 MOL/SDF file or an AssemblyCpp native graph
file. With `--run-strings=1`, it is an exact text-file path instead. The
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
`AssemblyCpp --help` for full details and accepted legacy option names.
`--telemetry` is available only in telemetry-enabled executables.
`--threads=auto` uses the OpenMP runtime default and therefore honours settings
such as `OMP_NUM_THREADS`; an explicit thread count applies to each process.

In a parallel-enabled executable, `--parallel=auto` prepares the root jobs and
DAG, then uses their estimated search work to choose parallel or serial
execution. A serial fallback reports its reason. `--parallel=on` bypasses only
that work estimate: it fails if parallel execution cannot be honored, such as
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
./build/release/AssemblyCpp strings.txt --run-strings=1
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
<summary><strong>C++ library</strong></summary>

Installed packages export `AssemblyCpp::Library`. If AssemblyCpp is installed
to a non-system prefix, pass that prefix when configuring the consuming
project:

```bash
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH=/absolute/path/to/assemblycpp-v5/build/install
```

Then link the imported target in the consuming project's `CMakeLists.txt`:

```cmake
find_package(AssemblyCpp 5 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE AssemblyCpp::Library)
```

```cpp
#include <assemblycpp.h>

#include <iostream>

int main()
{
    const auto result = assemblycpp::calculate("molecule.mol");
    if (!result)
    {
        std::cerr << result.error << '\n';
        return 1;
    }
    std::cout << result.assemblyIndex << '\n';
}
```

`calculateMolfile` accepts a V2000 molfile stream, while `calculateGraph`
accepts an AssemblyCpp native graph stream. `calculateBatch` processes several
inputs sequentially without process startup between items. Library calls do not
create output files. Search state is process-global, so the API is reusable but
not thread-safe; use separate processes for concurrent work.

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
The `ASSEMBLYCPP_STRICT_WARNINGS` CMake option exists for toolchain diagnosis,
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
underscores. Project CMake variables use the `ASSEMBLYCPP_UPPER_SNAKE_CASE`
prefix. C++ project-defined identifiers avoid leading underscores and
unexplained abbreviations. CI enforces the compiler and Python quality gates.

### Benchmarks

Build the optimized candidate, then run a maintained suite:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite quick
```

The `performance` preset targets x86-64-v3. See
[benchmarks/README.md](benchmarks/README.md) for the benchmark corpus, paired
comparisons, parallel builds, telemetry, LTO, PGO, and scaling guidance.

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

AssemblyCpp is licensed under
[Creative Commons Attribution-NonCommercial 4.0 International](License.md).
