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
c++ v5/main.cpp -std=c++23 -O3 \
  -I"${CONDA_PREFIX}/include" -o build/AssemblyCpp
```

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
| `--runtime=<TICKS>` | Non-negative integer | Unlimited | Stop after this many raw `std::clock` ticks. `CLOCKS_PER_SEC` ticks represent one second according to the platform's C++ runtime; whether `clock()` measures processor or elapsed time is implementation-specific. |
| `--enum-max=<COUNT>` | `1`–`INT_MAX` | `50000000` | Limit the internal subgraph-enumeration state budget. Larger limits can substantially increase memory use. |
| `--pathway=<0\|1>` | Boolean | `1` | Write the recovered assembly pathway to `INPUTPathway`. |
| `--remove-hydrogens=<0\|1>` | Boolean | `1` | Remove explicit hydrogen atoms from molfile inputs before calculation. This has no effect on native graph inputs. |
| `--compensate-disjoint=<0\|1>` | Boolean | `0` | Subtract one from final and intermediate assembly indices for every disconnected input component after the first. |
| `--memory-report=<0\|1>` | Boolean | `0` | On Linux, write `/proc/self/status`'s peak virtual-memory value (`VmPeak`, in kB) to `memUsage` in the current working directory. The option has no output on other platforms. |
| `--write-intermediate-mas=<0\|1>` | Boolean | `0` | Write every improved assembly index and its `std::clock` tick to `INPUTIntermediateMAs`. |

The runtime and enumeration limits can stop an exhaustive search. When that
happens, the output may contain the best assembly index found so far rather
than a proven minimum.

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
| `INPUTOut` | Every successful calculation | Assembly index and total `std::clock` ticks. |
| `INPUTPathway` | `--pathway=1` | Recovered pathway as JSON. |
| `INPUTIntermediateMAs` | `--write-intermediate-mas=1` | `std::clock` tick and improved assembly index pairs, using the selected disjoint-compensation setting. |
| `./memUsage` | `--memory-report=1` on Linux | Peak process virtual memory (`VmPeak`) reported by the kernel. |

### Examples

Calculate a molfile using the defaults:

```bash
./build/AssemblyCpp unitTests/alanine.mol
```

Skip pathway generation and lower the enumeration budget:

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
content, validation errors, compatibility names, input ordering, output flags,
and Linux memory reporting. Each calculation runs in an isolated temporary
directory, so test artifacts do not modify `unitTests/`. Selected cases also
compare generated pathways with reviewed JSON golden files:

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
