# AssemblyCpp v5

AssemblyCpp calculates molecular assembly indices and can recover assembly
pathways. It provides a command-line tool and a reusable C++20 library.

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
AssemblyCpp [OPTIONS] INPUT
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
| `--telemetry=<0\|1>` | `0` | Write search telemetry; available only in telemetry builds. |
| `--write-intermediate-mas=<0\|1>` | `0` | Write each improved index and its clock tick. |

`--enum-max` includes one-edge masks. Runtime and enumeration limits return the
best index found so far, which may not be the proven minimum. Run
`AssemblyCpp --help` for limit details and accepted legacy option names.

Use `--` before an input that starts with a dash:

```bash
./build/release/AssemblyCpp --pathway=0 -- -molecule.mol
```

### Outputs

For a molfile, `INPUT` below excludes the `.mol` suffix.

| File | Contents |
| --- | --- |
| `INPUTOut` | Assembly index, search status, and `std::clock` ticks. |
| `INPUTPathway` | Recovered pathway JSON when `--pathway=1`. |
| `INPUTIntermediateMAs` | Improved indices when `--write-intermediate-mas=1`. |
| `INPUTTelemetry.json` | Search counters and, for parallel runs, per-worker topology, work, and timing. |
| `memUsage` | Linux `VmPeak` value when `--memory-report=1`. |

An output error produces a non-zero exit status. After a handled Ctrl-C,
AssemblyCpp finalizes requested outputs when possible and exits with status 130
unless an output fails.

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

## Build presets

| Preset | Purpose |
| --- | --- |
| `release` | Portable release executable and library. |
| `dev` | Portable developer build with focused tests and telemetry. |
| `ci` | Portable release build with the full regression suite. |
| `performance` | x86-64-v3 build with a telemetry executable. |
| `performance-lto` | x86-64-v3 build with link-time optimization. |
| `parallel` | Experimental OpenMP, MPI, and hybrid executables. |
| `parallel-tests` | Portable parallel targets with correctness and telemetry tests. |
| `pgo-generate` / `performance-pgo` | GCC profile generation and profile-use builds. |

Performance presets require an x86-64-v3 processor. The parallel preset also
requires OpenMP and MPI. See [benchmarks/README.md](benchmarks/README.md) for
PGO, parallel-launch, telemetry, and comparison commands.

Train and build the PGO candidate in this order:

```bash
cmake --preset pgo-generate
cmake --build --preset pgo-train
cmake --preset performance-pgo
cmake --build --preset performance-pgo
```

Parallel search requires multiple workers, at least 32 bonds, pathway and
intermediate output disabled, and an unlimited runtime. Telemetry is available
from the `AssemblyCppOMPTelemetry`, `AssemblyCppMPITelemetry`, and
`AssemblyCppHybridTelemetry` targets. Other workloads run serially. Set
`ASSEMBLYCPP_PARALLEL_MIN_BONDS=0` only when testing smaller scaling inputs.
OpenMP workers dynamically lease root branches as they become available; hybrid
workers do the same within each MPI rank while ranks retain disjoint modulo
partitions. Root enumeration and DAG construction run once per process;
workers receive stable job indices and reconstruct masks in their own
thread-local arenas. `ASSEMBLYCPP_BRANCH_LEASE_SIZE` selects a positive lease
size and defaults to four branches, balancing lease overhead against tail
imbalance.

## Tests

Run the focused developer suite:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the `ci` preset instead to run every reviewed regression case. Run the
parallel parity and telemetry checks separately with:

```bash
cmake --preset parallel-tests
cmake --build --preset parallel-tests
ctest --preset parallel-tests
```

Test data and targeted commands are documented in
[unitTests/README.md](unitTests/README.md).

## Benchmarks

Build the optimized candidate, then run a maintained suite:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite quick
```

For a paired comparison:

```bash
python benchmarks/benchmark.py \
  --baseline-executable build/AssemblyCpp-before \
  --executable build/performance/AssemblyCpp \
  --suite quick \
  --json-output build/quick.json
```

See [benchmarks/README.md](benchmarks/README.md) for corpus, measurement, PGO,
and parallel-scaling details.

## Install and package

```bash
cmake --install build/release --prefix build/install
cpack --preset release
cmake --build build/release --target package_source
```

Packaging requires CMake 4.3 or newer. Create source archives from a clean
checkout because CPack includes the working tree. Packages include the README
and [CC BY-NC 4.0 license](License.md).
