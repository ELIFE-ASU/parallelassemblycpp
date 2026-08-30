# AssemblyCpp v5

AssemblyCpp computes molecular assembly indices and recovers assembly pathways.
It provides a command-line tool and a reusable C++20 library.

This repository implements the algorithm described by Ian Seet, Keith Y.
Patarroyo, Gage Siebert, Sara I. Walker, and Leroy Cronin in [*Rapid Exploration
of Assembly Chemical Space of Molecular Graphs*](https://arxiv.org/abs/2410.09100).

## Quick start

Requirements: CMake 3.25 or newer, Ninja, and a C++20 compiler.

```bash
cmake --preset release
cmake --build --preset release
./build/release/AssemblyCpp unitTests/alanine.mol
```

The command writes `unitTests/alanineOut` and, by default,
`unitTests/alaninePathway`.

To use the supplied Conda environment:

```bash
conda env create --file environment.yml
conda activate assemblycpp-v5
```

## Command line

```text
AssemblyCpp INPUT [OPTIONS]
AssemblyCpp --help
```

`INPUT` may be a V2000 molfile or an AssemblyCpp native graph file. The `.mol`
suffix is optional; native graph filenames must be supplied in full. Options
may appear before or after the input and use `--name=value` syntax. Boolean
values are `0` or `1`.

### Options

| Option | Default | Purpose |
| --- | --- | --- |
| `-h`, `--help` | — | Show command help. |
| `--runtime=<TICKS>` | Unlimited | Stop after the given `std::clock` budget. |
| `--enum-max=<COUNT>` | `50000000` | Limit retained connected masks in the initial DAG. |
| `--pathway=<0\|1>` | `1` | Write the recovered pathway. |
| `--remove-hydrogens=<0\|1>` | `1` | Remove explicit hydrogens from molfiles. |
| `--verbose=<0\|1>` | `0` | Print the parsed input graph. |
| `--compensate-disjoint=<0\|1>` | `0` | Subtract one per processed component after the first. |
| `--memory-report=<0\|1>` | `0` | Write Linux peak virtual memory to `memUsage`. |
| `--telemetry=<0\|1>` | `0` | Write search telemetry. |
| `--write-intermediate-mas=<0\|1>` | `0` | Write each improved index and its clock tick. |

`--enum-max` includes one-edge masks. Runtime and enumeration limits return the
best index found so far, which may not be the proven minimum. Run
`AssemblyCpp --help` for full details and accepted legacy option names.
`--telemetry` is available only in telemetry-enabled executables.

### Outputs

For a molfile, `INPUT` below excludes the `.mol` suffix.

| File | Contents |
| --- | --- |
| `INPUTOut` | Assembly index, search status, and `std::clock` ticks. |
| `INPUTPathway` | Recovered pathway JSON when `--pathway=1`. |
| `INPUTIntermediateMAs` | Improved indices when `--write-intermediate-mas=1`. |
| `INPUTTelemetry.json` | Search counters from a telemetry-enabled executable. |
| `memUsage` | Linux `VmPeak` value when `--memory-report=1`. |

## C++ library

Installed packages export `AssemblyCpp::Library`:

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

`calculateMolfile` accepts a V2000 molfile stream. `calculateBatch` processes
several inputs sequentially without process startup between items. Library
calls do not create output files. Search state is process-global, so the API is
reusable but not thread-safe; use separate processes for concurrent work.

## Development

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

## Install and package

```bash
cmake --install build/release --prefix build/install
cpack --preset release
cmake --build build/release --target package_source
```

Use CMake 4.3 or newer to normalize archive ownership. Create source archives
from a clean checkout because CPack includes the working tree.

## License

AssemblyCpp is licensed under
[Creative Commons Attribution-NonCommercial 4.0 International](License.md).
