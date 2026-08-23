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

Run the Python harness directly when selecting cases or controlling parallelism:

```bash
python unitTests/unitTester.py build/release/AssemblyCpp --jobs 4
python unitTests/unitTester.py build/release/AssemblyCpp --limit 20
python unitTests/unitTester.py --build --pathways-only --verbose
```

The harness always runs command-line checks before the selected regression
cases. These cover validation, legacy names, limits, inputs, outputs, and Linux
memory reporting. `--build` also compiles and runs the focused C++ tests, using
an x86-64-v3 executable with telemetry as a test shortcut. Use a CMake portable
build on older x86-64 or non-x86 systems.

Audit manifests and fixture coverage without running calculations:

```bash
python unitTests/unitTester.py --audit
```

Fixture-only molfiles are reported by the audit and are not treated as passing
regression cases. Use `--verbose` to list them.
