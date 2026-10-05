# Development and packaging

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Run source-tree examples from the repository root unless stated otherwise.

Tests and benchmark tooling require Python 3.10 or newer. The full developer
suite includes plotting tests and therefore also needs Matplotlib. Parallel
builds require the selected OpenMP and/or MPI dependencies; PGO builds require
GCC. The supplied Conda environment includes these development dependencies.

Create and activate it as shown in the [quick start](../README.md#quick-start).
With your own Python environment, install Matplotlib there and activate it
before configuring. If CMake selected another interpreter, configure with
`-DPython3_EXECUTABLE=/absolute/path/to/python`. Ruff is needed for the separate
Python quality gates.

## Build presets

All presets use Ninja and write to `build/<preset>`. Configure before building;
`cmake --build --preset NAME` does not configure the directory for you.

| Preset | Configuration | Purpose |
| --- | --- | --- |
| `release` | Release, portable CPU target | Serial command and library; tests and telemetry disabled |
| `dev` | RelWithDebInfo, portable CPU target | Focused tests, 20 manifest cases, and telemetry |
| `ci` | Release, portable CPU target | Full manifest regression and telemetry |
| `performance` | Release, x86-64-v3 | Serial benchmarks and telemetry |
| `performance-lto` | Performance plus LTO | Link-time optimization experiments |
| `parallel` | Performance plus OpenMP, MPI, and hybrid | Parallel benchmarks and telemetry |
| `parallel-tests` | Release, portable CPU target | Parallel executables and tests; test preset selects the `parallel` label |
| `pgo-generate` | Performance plus LTO, GCC only | Profile instrumentation; `pgo-train` is the training build preset |
| `performance-pgo` | Performance plus LTO, GCC only | Use validated training profiles |

The `parallel` and performance presets require x86-64-v3 hardware. Override
`-DPARALLELASSEMBLYCPP_X86_64_V3=OFF` for a portable CPU target. PGO training and
profile-use instructions are in the [benchmark guide](../benchmarks/README.md#lto-and-pgo).
Use `cmake --list-presets=all` to list the available configure, build, test, and
package presets.

Enable individual optional executables with
`PARALLELASSEMBLYCPP_BUILD_OPENMP`, `PARALLELASSEMBLYCPP_BUILD_MPI`, or
`PARALLELASSEMBLYCPP_BUILD_HYBRID`. OpenMP alone does not require MPI; hybrid
requires both. See [parallel builds](parallel.md). In-source CMake builds are
rejected; use a separate build directory.

## Tests

Run the focused developer suite:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the `ci` preset for the full regression suite and `parallel-tests` for
parallel parity and telemetry checks. See
[unitTests/README.md](../unitTests/README.md) for targeted commands and fixture
details.

## Quality gates

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

## Benchmarks

Build the optimized candidate, then run a maintained suite:

```bash
cmake --preset performance
cmake --build --preset performance
python benchmarks/benchmark.py \
  --executable build/performance/ParallelAssemblyCpp \
  --suite quick
```

The `performance` preset targets x86-64-v3. See
[benchmarks/README.md](../benchmarks/README.md) for the benchmark corpus, paired
comparisons, parallel builds, telemetry, LTO, PGO, and scaling guidance. The
[Paclitaxel thread sweep](../benchmarks/README.md#paclitaxel-thread-sweep) automates
paired serial/OpenMP measurements across every available CPU count, with
automatic topology detection and an optional physical-core-only sweep.

## Packaging

These commands create binary and source archives for release distribution;
they are not required for a normal installation:

```bash
cmake --preset release
cmake --build --preset release
cpack --preset release
cmake --build --preset release --target package_source
```

Packaging only: CMake 4.3 or newer normalizes archive ownership, and the Conda
environment supplies it. Building the project still requires only CMake 3.25.
Create source archives from a clean checkout because CPack includes the working
tree.

The binary archive contains the serial executable, static library, public
header, CMake package, and guides. It does not bundle the compiler's runtime
libraries, so "portable" describes the CPU target rather than a statically
linked, universally runnable binary. CPack produces SHA-256 sidecars; binary
archives use TGZ on Unix and ZIP on Windows, and source archives use TGZ.

## Contributing and reporting problems

Use [GitHub issues](https://github.com/ELIFE-ASU/parallelassemblycpp/issues) for
bug reports and focused proposals. Include the commit, OS, compiler, CMake
preset/options, exact command, a reproducible input, and expected versus actual
results. For performance reports, include thread/rank placement, output mode,
time limits, and paired benchmark evidence.

For a code change, run the relevant tests and the quality gates above. Changes
to parsing, options, output formats, or the public API should update their
corresponding tests and [CLI](cli.md) or [library](library.md) documentation.
Performance claims should follow the benchmark guide's
[promotion criteria](../benchmarks/README.md#paired-comparisons).
