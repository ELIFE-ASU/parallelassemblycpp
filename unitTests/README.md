# Tests

`regression_cases.tsv` maps each reviewed fixture to its expected assembly
index. `pathway_cases.tsv` selects cases whose pathway JSON must also match a
file in `expected_pathways/`. Pathway comparisons ignore formatting.

Run the focused developer suite from the repository root:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the `ci` preset to run the complete regression manifest.

Parallel solver parity and per-worker telemetry invariants use a dedicated
preset so performance builds remain test-free:

```bash
cmake --preset parallel-tests
cmake --build --preset parallel-tests
ctest --preset parallel-tests
```

The parallel harness repeats serial/OpenMP calculations across 1, 2, and 4
workers, including the 63/64/65 and 127/128/129 edge-mask boundaries and a
disconnected molecule. It exercises both one-branch and chunked dynamic leases.
When MPI targets are available, the same harness also checks rank-partitioned
MPI and hybrid aggregation, exact lease reductions, and complete branch
coverage.

Run the Python harness directly when selecting cases or controlling parallelism:

```bash
python unitTests/unitTester.py build/release/ParallelAssemblyCpp --jobs 4
python unitTests/unitTester.py build/release/ParallelAssemblyCpp --limit 20
python unitTests/unitTester.py --build --pathways-only --verbose
```

Outside `--audit` mode, the harness runs command-line checks before the selected
regression cases. These cover validation, legacy names, limits, inputs, outputs,
and Linux memory reporting. They also run the five-case upstream string corpus, validate
line-ending and incompatible-option behavior, and check per-line pathway JSON.
Canonical and legacy flag spellings are exercised on actual inputs, with
invalid values, duplicate aliases, unusual filenames, output failures, and
disabled-output preservation checked separately. Graph-only flags in string
mode and enabled reversal matching in graph mode must produce diagnostics.
The focused `stringAssemblyTester` compares short binary and ternary inputs with
an independent exhaustive search and covers interval merging, remnants,
multi-step pathways, reversal, cancellation, target-index stopping, and JSON
escaping. Pathway checks replay fragment boundaries to reject impossible reuse.
Unicode strings also check scalar and encoding boundaries, embedded control
characters, malformed UTF-8 byte offsets, and parity with serial, OpenMP, MPI,
and hybrid execution. CLI checks preserve BOMs, combining marks, and Unicode
line separators as symbols and verify that malformed records stop processing
after preserving completed results.

The CMake suite also runs `graphRepairTester.py`, which independently replays
Re-Pair certificates and compares their bounds with an exact small-graph oracle.
The `ci` preset enables its full exhaustive and randomized corpus.
`graphRepairCliTester.py` checks the `--algorithm=full` and `--algorithm=re-pair`
selectors, the legacy `--upper-bound=graph-repair` selector, incompatible options,
and the distinct exact-search and heuristic certificate formats. Zero-runtime
checks ensure that full search retains its trivial initial bound without a
Re-Pair prepass.

`--build` compiles and runs four standalone C++ tests: masks, tree and cyclic
canonicalization, and string assembly. It then builds an x86-64-v3 executable
with telemetry as a test shortcut. Use CMake for the complete unit-test suite
and for portable builds on older x86-64 or non-x86 systems.

Audit manifests and fixture coverage without running calculations:

```bash
python unitTests/unitTester.py --audit
```

Fixture-only molfiles are reported by the audit and are not treated as passing
regression cases. Use `--verbose` to list them.
