# AssemblyCpp v5

https://arxiv.org/abs/2410.09100#

This repository contains the C++ implementation of v5 of the assembly algorithm,
with command-line and reusable C++ library interfaces.

School of Chemistry, The University of Glasgow, University Avenue, Glasgow G12 8QQ, United Kingdom

Authors Ian Seet, Leroy Cronin


## Build

AssemblyCpp uses CMake 3.25 or newer, Ninja, and a C++20 compiler. The release
preset builds a portable CPU baseline by default:

```bash
cmake --preset release
cmake --build --preset release
./build/release/AssemblyCpp --help
```

For an isolated development environment, install
[Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or another
Conda distribution, then run:

```bash
conda env create --file environment.yml
conda activate assemblycpp-v5
cmake --preset release
cmake --build --preset release
```

The environment provides CMake, Ninja, a C++ compiler, and Python. Recreate it
after changing `environment.yml` with:

```bash
conda env remove --name assemblycpp-v5
conda env create --file environment.yml
```

The performance preset preserves the former x86-64-v3 and POPCNT optimization
as an explicit opt-in. It also builds a separately instrumented executable:

```bash
cmake --preset performance
cmake --build --preset performance
./build/performance/AssemblyCpp --help
./build/performance/AssemblyCppTelemetry --help
```

Do not run the performance binaries on processors that lack the x86-64-v3
feature level. Portable and performance builds never mix in the same build
directory.

Search telemetry is compile-time opt-in so ordinary search loops contain no
counter branches. `--telemetry=1` is available only in
`AssemblyCppTelemetry`; the standard executable reports it as an unknown
option.

### Install and package

Install the portable executable, static library, public header, and documentation
to a chosen prefix:

```bash
cmake --install build/release --prefix build/install
./build/install/bin/AssemblyCpp --help
```

Create a platform-specific binary archive and a source archive with SHA-256
checksum sidecars:

```bash
cpack --preset release
cmake --build build/release --target package_source
```

Create source archives from a clean checkout because CPack packages the working
tree as well as tracked files.

Binary archive names state the operating system, architecture, and CPU
baseline. They are still specific to the platform and toolchain used to build
them; `portable` describes the CPU instruction baseline, not a universal
binary. CI builds its checked Linux archive on Ubuntu 22.04 with the system
compiler. Packages include the README and the CC BY-NC 4.0 license.

### C++ library and in-process batches

The CMake build exports the `AssemblyCpp::Library` target. An installed package
can be consumed with:

```cmake
find_package(AssemblyCpp 5 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE AssemblyCpp::Library)
```

`calculate` and `calculateMolfile` return assembly indices and status directly
without creating `INPUTOut` or pathway files. `calculateBatch` performs an
ordered sequence in one process, avoiding process startup for short molecules:

```cpp
#include <assemblycpp.h>

#include <iostream>
#include <string>
#include <vector>

int main()
{
    const std::vector<std::string> inputs = {
        "molecules/icosane.mol",
        "molecules/sucrose.mol",
    };
    for (const assemblycpp::CalculationResult& result :
         assemblycpp::calculateBatch(inputs))
    {
        if (!result)
        {
            std::cerr << result.input << ": " << result.error << '\n';
            continue;
        }
        std::cout << result.input << '\t' << result.assemblyIndex << '\n';
    }
}
```

The batch interface is sequential. Search storage is currently process-global,
so library calls are reusable but not thread-safe; use separate processes for
concurrent calculations. A runtime limit applies independently to each item.

## Running AssemblyCpp

```text
AssemblyCpp INPUT [OPTIONS]
AssemblyCpp [OPTIONS] INPUT
AssemblyCpp --help
```

`INPUT` can be either:

- A molfile path, with or without its `.mol` suffix. For example,
  `unitTests/alanine` and `unitTests/alanine.mol` read the same file.
- A file in AssemblyCpp's native graph format, such as
  `unitTests/graphio_test`. Pass the complete graph filename.

Options may appear before or after the input. Use `--` to stop option parsing
when an input name begins with a dash. Option values must follow the option in
the same argument using `=`. Boolean values must be exactly `0` (disabled) or
`1` (enabled).

### Options

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-h`, `--help` | — | — | Print usage, input, option, output, and compatibility details, then exit. |
| `--runtime=<TICKS>` | Non-negative integer | Unlimited | Cooperatively stop the search after this many elapsed `std::clock` ticks. `CLOCKS_PER_SEC` ticks represent one second according to the platform's C++ runtime; whether `clock()` measures processor or elapsed time is implementation-specific. Checks occur between search operations, so the budget can be exceeded by the duration of an in-progress operation. A value of `0` records the initial upper bound and then stops before enumeration. |
| `--enum-max=<COUNT>` | `1`–`INT_MAX` | `50000000` | Cap the number of unique connected subgraph masks retained while building the initial enumeration DAG, including one-edge masks. The search stops before retaining a mask beyond the limit. Larger limits can substantially increase memory use. |
| `--pathway=<0\|1>` | Boolean | `1` | Write the recovered assembly pathway to `INPUTPathway`. |
| `--remove-hydrogens=<0\|1>` | Boolean | `1` | Remove explicit hydrogen atoms from molfile inputs before calculation. This has no effect on native graph inputs. |
| `--verbose=<0\|1>` | Boolean | `0` | Print input summaries and parsed molecular graphs. |
| `--compensate-disjoint=<0\|1>` | Boolean | `0` | Subtract one from final and intermediate assembly indices for every component after the first in the processed graph. For molfiles, components are counted after optional hydrogen removal. Empty graphs receive no adjustment. |
| `--memory-report=<0\|1>` | Boolean | `0` | After the other calculation outputs are written successfully on Linux, write `/proc/self/status`'s peak virtual-memory value (`VmPeak`, in kB) to `memUsage` in the current working directory. The option has no output on other platforms. |
| `--telemetry=<0\|1>` | Boolean | `0` | In a telemetry-enabled build, write retained-mask and matching counts, canonicalisation activity, legacy VF2 counters, cache rates, and phase-specific memory to `INPUTTelemetry.json`. The current exact cyclic canonicaliser does not invoke VF2, so those compatibility counters remain zero. |
| `--write-intermediate-mas=<0\|1>` | Boolean | `0` | Write each newly found best assembly index and its elapsed `std::clock` tick to `INPUTIntermediateMAs`. |

The runtime and enumeration limits can stop an exhaustive search. `INPUTOut`
still begins with a numeric best assembly index, followed by a status line
naming the limit that was reached. The best value from a limited search is not
necessarily a proven minimum.

Canonical options use two leading dashes and kebab-case names. Both one and two
leading dashes are accepted for canonical and compatibility names, so existing
scripts can continue to use the former spellings:

| Canonical option | Compatibility names |
| --- | --- |
| `--runtime` | `-runTime` |
| `--enum-max` | `-enumMax` |
| `--pathway` | `-pathway` |
| `--remove-hydrogens` | `-removeHydrogens` |
| `--compensate-disjoint` | `-compensateDisjoint`, `-disjointCompensation` |
| `--memory-report` | `-memTest`, `-testMemory` |
| `--write-intermediate-mas` | `-writeIntermediateMAs` |

Unknown options, missing values, invalid values, duplicate value-taking options,
and missing inputs produce an error and a non-zero exit status.

### Outputs

For molfiles, `INPUT` below means the input path with `.mol` removed. For native
graph files it is the complete input path.

| File | When written | Contents |
| --- | --- | --- |
| `INPUTOut` | When the input is read and this file can be opened | Numeric assembly index, any limit or interruption status, and total elapsed `std::clock` ticks. |
| `INPUTPathway` | `--pathway=1` | Recovered pathway as JSON. |
| `INPUTIntermediateMAs` | `--write-intermediate-mas=1` | `std::clock` tick and improved assembly index pairs, using the selected disjoint-compensation setting. |
| `INPUTTelemetry.json` | `--telemetry=1` | Versioned search counters, cache statistics, processed mask width, phase clock ticks, and phase memory. |
| `./memUsage` | `--memory-report=1` on Linux, after other outputs succeed | Peak process virtual memory (`VmPeak`) reported by the kernel. |

Telemetry is opt-in so its accounting does not affect ordinary benchmark
timings. On Linux, each flat search phase resets and reads the process's
resident high-water mark through `/proc/self/clear_refs` and
`/proc/self/status`. Phase `peak_rss_kib` values are absolute resident sets and
must not be added together. If the kernel does not provide a resettable
high-water mark, phase peaks are `null`; start/end snapshots and all counters
remain available. Resetting `VmHWM` also resets `ru_maxrss`, so external
`time -v` peak-RSS measurements from a telemetry-enabled process are not valid.

If an enabled output cannot be opened or fully written, AssemblyCpp reports the
affected path and exits with a non-zero status.

On a Ctrl-C handled before final output finalization, the search unwinds
cooperatively and records an interrupted status. If all requested outputs can
be written, they are flushed (including the Linux memory report) and the process
exits with status `130`; an output failure still exits with status `1`.

### Examples

Calculate a molfile using the defaults:

```bash
./build/release/AssemblyCpp unitTests/alanine.mol
```

Skip pathway output and lower the enumeration budget:

```bash
./build/release/AssemblyCpp unitTests/alanine --pathway=0 --enum-max=1000000
```

Keep explicit hydrogen atoms and write intermediate best values:

```bash
./build/release/AssemblyCpp \
  --remove-hydrogens=0 \
  --write-intermediate-mas=1 \
  unitTests/alanine.mol
```

## Tests

The default developer preset builds the application, a telemetry sibling, and
the focused C++ tests. CTest checks both mask domains at every 64-bit boundary,
preserves the tree canonicaliser's constrained-stack regression, audits test
and benchmark tooling, checks the command-line interface, and runs a 20-case
regression sample:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the CI preset to run all reviewed regression cases:

```bash
cmake --preset ci
cmake --build --preset ci
ctest --preset ci
```

The CLI checks cover help content, validation errors,
compatibility names, input ordering, resource-limit
boundaries, hydrogen removal, disjoint compensation, output toggles and
failures, and Linux memory reporting. Each calculation runs in an isolated
temporary directory, so test artifacts do not modify `unitTests/`. Selected
cases also compare generated pathways with reviewed JSON golden files:

```bash
python unitTests/unitTester.py build/release/AssemblyCpp --jobs 4
```

For a quicker development check, limit the number of cases:

```bash
python unitTests/unitTester.py build/release/AssemblyCpp --limit 20
```

Audit the manifest for duplicate cases, missing fixtures, conflicting
expectations, and fixture coverage with:

```bash
python unitTests/unitTester.py --audit
```

Use `python unitTests/unitTester.py --help` to see all options, including custom
manifests, per-case timeouts, and verbose output.

GitHub Actions builds through CMake, runs the complete regression manifest,
smoke-tests a staged installation, and generates both package forms for every
push and pull request.

## Benchmark

The benchmark runner supports both the original single-input benchmark and
manifest-driven workload suites. It disables pathway generation, validates
configured assembly-index expectations on every calculation, and reports
wall-clock and algorithm clock-tick statistics. The no-argument workload is
five measured `ketoconazole` runs after one warm-up. Reuse the optimized CMake
build with:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp
```

The runner's `--build` option remains available as a direct-compiler shortcut.

Run the reviewed quick corpus for routine development, or the larger reviewed
corpus before merging an optimization:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp --suite quick
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp --suite full
```

Suite cases are run serially in rounds, with the starting case rotated each
round to reduce order and temperature bias. Summaries include the median,
median absolute deviation (MAD), and 95th percentile. List the available cases
and their expectation status with:

```bash
python benchmarks/benchmark.py --list-cases
python benchmarks/benchmark.py --suite profile --list-cases
```

The `profile` suite contains deliberately long inputs for profiling search hot
paths. Its non-ketoconazole expectations are marked provisional because those
fixtures are not part of the reviewed regression manifest. A one-pass profile
smoke run is:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite profile --runs 1 --warmup 0
```

The `scaling` suite measures a cumulative sequence of two through eight
disconnected amino-acid components (18 to 68 atoms, 16 to 60 bonds). Its table
places system size beside wall-time and algorithm clock-tick costs:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp --suite scaling
```

Use clock ticks to compare the smallest cases, where fixed process-startup and
input-parsing costs can dominate wall time. The sequence adds a different amino
acid at each step, so it characterizes this workload rather than an asymptotic
complexity law for arbitrary molecules.

For a paired before/after comparison, keep the old executable and pass it as
the baseline. Baseline and candidate calculations run adjacently, alternating
AB/BA order between rounds. Paired mode defaults to six measured rounds so both
executables occupy each position equally often. Reported speedups are paired
ratios; values greater than one mean the candidate is faster.

```bash
python benchmarks/benchmark.py \
  --suite quick \
  --baseline-executable build/AssemblyCpp-before \
  --executable build/performance/AssemblyCpp
```

Raw samples, executable SHA-256 fingerprints, platform details, and summaries
can be retained for later comparison:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --suite quick \
  --json-output benchmark-results.json
```

To reuse an existing build, increase the number of measured runs, or select a
single different input:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp --runs 10
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp \
  --input unitTests/sucrose.mol \
  --expected 8
```

The corpus format and measurement methodology are documented in
`benchmarks/README.md`. Use `python benchmarks/benchmark.py --help` for the
remaining options. Run comparisons on the same machine and under similar load;
the runner intentionally does not impose a pass/fail timing threshold.
