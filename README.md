# AssemblyCpp v5

https://arxiv.org/abs/2410.09100#

This repository contains the C++ implementation of v5 of the assembly algorithm, in a "script" style version.

School of Chemistry, The University of Glasgow, University Avenue, Glasgow G12 8QQ, United Kingdom

Authors Ian Seet, Leroy Cronin


## Build with Conda

Install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or another
Conda distribution, then run the following commands from the repository root:

```bash
conda env create --file environment.yml
conda activate assemblycpp-v5
mkdir -p build
c++ v5/main.cpp -std=c++23 -O3 -mpopcnt -march=x86-64-v3 \
  -I"${CONDA_PREFIX}/include" -o build/AssemblyCpp
```

The generated executable targets x86-64-v3 processors and requires POPCNT
support. The test and benchmark `--build` commands use the same target. For an
older x86-64 processor or a non-x86 platform, omit `-mpopcnt` and
`-march=x86-64-v3` to produce a portable `-O3` build instead.

The environment installs a C++ compiler, Boost (including the Boost Graph Library),
and Python for the test runner. Activating it places the Conda-provided compiler
on `PATH`. Recreate the environment after changing `environment.yml` with:

```bash
conda env remove --name assemblycpp-v5
conda env create --file environment.yml
```

To confirm that the executable was built successfully:

```bash
./build/AssemblyCpp --help
```

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
| `--compensate-disjoint=<0\|1>` | Boolean | `0` | Subtract one from final and intermediate assembly indices for every component after the first in the processed graph. For molfiles, components are counted after optional hydrogen removal. Empty graphs receive no adjustment. |
| `--memory-report=<0\|1>` | Boolean | `0` | After the other calculation outputs are written successfully on Linux, write `/proc/self/status`'s peak virtual-memory value (`VmPeak`, in kB) to `memUsage` in the current working directory. The option has no output on other platforms. |
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
| `./memUsage` | `--memory-report=1` on Linux, after other outputs succeed | Peak process virtual memory (`VmPeak`) reported by the kernel. |

If an enabled output cannot be opened or fully written, AssemblyCpp reports the
affected path and exits with a non-zero status.

On a Ctrl-C handled before final output finalization, the search unwinds
cooperatively and records an interrupted status. If all requested outputs can
be written, they are flushed (including the Linux memory report) and the process
exits with status `130`; an output failure still exits with status `1`.

### Examples

Calculate a molfile using the defaults:

```bash
./build/AssemblyCpp unitTests/alanine.mol
```

Skip pathway output and lower the enumeration budget:

```bash
./build/AssemblyCpp unitTests/alanine --pathway=0 --enum-max=1000000
```

Keep explicit hydrogen atoms and write intermediate best values:

```bash
./build/AssemblyCpp \
  --remove-hydrogens=0 \
  --write-intermediate-mas=1 \
  unitTests/alanine.mol
```

## Tests

The test runner compiles the executable, checks the command-line interface, and
runs the complete regression manifest in one command. The CLI checks cover help
content, validation errors, compatibility names, input ordering, resource-limit
boundaries, hydrogen removal, disjoint compensation, output toggles and
failures, and Linux memory reporting. Each calculation runs in an isolated
temporary directory, so test artifacts do not modify `unitTests/`. Selected
cases also compare generated pathways with reviewed JSON golden files:

```bash
python unitTests/unitTester.py --build --jobs 4
```

For a quicker development check, limit the number of cases:

```bash
python unitTests/unitTester.py --build --limit 20
```

Audit the manifest for duplicate cases, missing fixtures, conflicting
expectations, and fixture coverage with:

```bash
python unitTests/unitTester.py --audit
```

Use `python unitTests/unitTester.py --help` to see all options, including custom
manifests, per-case timeouts, and verbose output.

GitHub Actions audits the test data and runs the complete regression manifest
for every push and pull request.

## Benchmark

The benchmark runner supports both the original single-input benchmark and
manifest-driven workload suites. It disables pathway generation, validates
configured assembly-index expectations on every calculation, and reports
wall-clock and algorithm clock-tick statistics. The no-argument default remains
five measured `ketoconazole` runs after one warm-up:

```bash
python benchmarks/benchmark.py --build
```

Run the reviewed quick corpus for routine development, or the larger reviewed
corpus before merging an optimization:

```bash
python benchmarks/benchmark.py --suite quick
python benchmarks/benchmark.py --suite full
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
python benchmarks/benchmark.py --suite profile --runs 1 --warmup 0
```

For a paired before/after comparison, keep the old executable and pass it as
the baseline. Baseline and candidate calculations run adjacently, alternating
AB/BA order between rounds. Paired mode defaults to six measured rounds so both
executables occupy each position equally often. Reported speedups are paired
ratios; values greater than one mean the candidate is faster.

```bash
python benchmarks/benchmark.py \
  --suite quick \
  --baseline-executable build/AssemblyCpp-before \
  --executable build/AssemblyCpp
```

Raw samples, executable SHA-256 fingerprints, platform details, and summaries
can be retained for later comparison:

```bash
python benchmarks/benchmark.py \
  --suite quick \
  --json-output benchmark-results.json
```

To reuse an existing build, increase the number of measured runs, or select a
single different input:

```bash
python benchmarks/benchmark.py --runs 10
python benchmarks/benchmark.py \
  --input unitTests/sucrose.mol \
  --expected 8
```

The corpus format and measurement methodology are documented in
`benchmarks/README.md`. Use `python benchmarks/benchmark.py --help` for the
remaining options. Run comparisons on the same machine and under similar load;
the runner intentionally does not impose a pass/fail timing threshold.
