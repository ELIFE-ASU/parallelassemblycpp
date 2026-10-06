# Tests

Run commands below from the repository root. The CMake test builds require a
C++20 compiler, CMake 3.25 or newer, Ninja, and Python 3.10 or newer. The benchmark
tooling tests also require Matplotlib. Activate the supplied Conda environment,
or install Matplotlib into the Python environment selected by CMake. See the
[development guide](../docs/development.md) for setup and interpreter selection.

## Run the suites

The `dev` preset builds the portable solver, telemetry executable, library, and
unit tests. It runs CLI checks, tooling checks, and the first 20 manifest cases:

```bash
cmake --preset dev
cmake --build --preset dev --parallel 2
ctest --preset dev
```

The `ci` preset runs the complete regression manifest and the full exhaustive
and randomized graph Re-Pair corpus:

```bash
cmake --preset ci
cmake --build --preset ci --parallel 2
ctest --preset ci
```

Use CTest names and labels to select checks or list them without running them:

```bash
ctest --preset dev -N
ctest --preset dev -L unit
ctest --preset dev -R '^unit\.string-'
ctest --preset dev -L tooling
```

The portable `parallel-tests` build additionally requires OpenMP and MPI. Its
test preset selects only tests labelled `parallel`; run the `ci` suite as well
for complete serial regression coverage. See the
[parallel guide](../docs/parallel.md) for build and runtime requirements.

```bash
cmake --preset parallel-tests
cmake --build --preset parallel-tests --parallel 2
ctest --preset parallel-tests
```

The `release`, `performance`, and `parallel` presets have testing disabled.

## Select regression cases directly

[regression_cases.tsv](regression_cases.tsv) maps reviewed fixtures to their
expected assembly indices. [pathway_cases.tsv](pathway_cases.tsv) selects cases
whose pathway JSON must also match a file in
[expected_pathways/](expected_pathways/). Comparisons parse JSON, so formatting
and object-key order do not affect the result.

Pass the executable from your chosen build explicitly:

```bash
python3 unitTests/unitTester.py build/dev/ParallelAssemblyCpp --jobs 4
python3 unitTests/unitTester.py build/dev/ParallelAssemblyCpp --limit 20
python3 unitTests/unitTester.py build/dev/ParallelAssemblyCpp --pathways-only --verbose
```

`--jobs` controls independent test processes, not solver threads. `--limit N`
selects the first N cases after any `--pathways-only` filtering. `--timeout`
sets the regression timeout per case in seconds (default: 300). Every normal
run performs the CLI checks before the selected cases; limiting the manifest
does not skip those checks. CMake uses two regression processes and a
120-second timeout per case, plus a separate one-case telemetry regression.

The optional `--build` shortcut uses `$CXX` or `c++` to compile four standalone
tests (masks, tree canonicalization, cyclic canonicalization, and string
assembly), then an executable with telemetry. It requires a compiler accepting
GCC-style flags and an x86-64-v3 CPU. Its default output is
`build/ParallelAssemblyCpp`, which differs from the CMake preset paths. Use
CMake for portable builds and the complete test suite.

```bash
python3 unitTests/unitTester.py --build --pathways-only --verbose
```

## Coverage

| Area | Checks |
| --- | --- |
| Core graph search | Active-word masks and their lifetimes, transposition tables, tree and cyclic canonicalization, cancellation polling, matching-bound refresh, fragmentation metadata, task transfer, MOL/SDF and native graph parsing, and pathway generation. Telemetry variants exercise the cache and bound-refresh counters. |
| Addition-chain bounds | Exhaustive certification of the scalar table, independent vector-chain oracle, bounded-search exhaustion, cache isolation, and exact graph parity against independent graph and labelled-path oracles, including preprocessing and disconnected compensation. |
| CLI and manifest | Reviewed assembly indices and golden pathways; canonical and legacy flags; invalid values and duplicate aliases; input/output failures and unusual filenames; disabled-output preservation; limits; Linux memory reporting; string records, incompatible options, and malformed UTF-8 diagnostics. |
| Exact string search | Independent exhaustive short-string oracle, fragment-boundary pathway replay, interval merging, remnants, reversal, cancellation, target-index stopping, JSON escaping, and Unicode scalar/encoding boundaries. Independently replayed Re-Pair seeds cover nested, reversed, Unicode, and randomized constructions, optimal witness retention, and improvement of suboptimal seeds. |
| Re-Pair bounds | Independent graph and string certificate replay and small-instance exact oracles; deterministic choices, overlaps, Unicode and invalid encodings; CLI algorithm selectors, certificate formats, output failures, and heuristic status. Graph full-search checks verify that zero runtime retains the trivial initial bound without a graph Re-Pair prepass. |
| Public library | Exact and Re-Pair graph/string calls, file and stream inputs, batch recovery, budgets, and repeated calls without leaking calculation state or creating output files. |
| Tooling | Fixture audit, repository text policy, benchmark runner unit tests, and benchmark corpus validation. |

[parallelSolverTester.py](parallelSolverTester.py) repeats serial/OpenMP
calculations with 1, 2, and 4 workers. Cases cover the 63/64/65 and 127/128/129
edge-mask boundaries and a disconnected molecule. Telemetry checks cover
single-branch and chunked leases, adaptive task donation and tail refills,
worker aggregation, and complete branch coverage. MPI and hybrid variants also
exercise distributed work and refill reductions;
[mpiRefillTester.cpp](mpiRefillTester.cpp) checks pending replies, cancellation,
concurrent claims, and empty frontiers.

[parallelStringTester.py](parallelStringTester.py) checks exact string and
pathway parity across the enabled execution modes, including reversal and
Unicode. The string Re-Pair CLI tests also run through OpenMP, MPI, and hybrid
executables; this heuristic executes on the root process in distributed modes.

## Check the installed library

The package consumer exercises the exported `ParallelAssemblyCpp::Library`
target with both graph and string APIs. After building `release`, configure it
against a staged install (the example below uses a POSIX shell):

```bash
cmake --preset release
cmake --build --preset release --parallel 2
cmake --install build/release --prefix build/stage
cmake -S unitTests/packageConsumer -B build/package-consumer -G Ninja \
  -DCMAKE_PREFIX_PATH="$PWD/build/stage"
cmake --build build/package-consumer --parallel 2
ctest --test-dir build/package-consumer --output-on-failure --no-tests=error
```

[CI](../.github/workflows/test.yml) runs the serial and parallel suites, checks
Clang and MSVC unit builds, verifies the minimum CMake version, and tests staged
installs and binary/source packages. Its Python style checks use Ruff in
addition to the CTest tooling checks.

## Audit fixtures

Audit manifests, golden files, duplicate contents, and fixture coverage without
running solver calculations:

```bash
python3 unitTests/unitTester.py --audit
python3 unitTests/unitTester.py --audit --verbose
```

Fixture-only molfiles are reported by the audit and are not passing regression
cases. Add a reviewed expected index to the manifest before counting a fixture
as regression coverage. `--verbose` lists the fixture-only molecules.
