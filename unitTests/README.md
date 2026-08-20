# Test data

`regression_cases.tsv` is the canonical regression manifest. Each row names a
fixture in this directory and its reviewed assembly-index expectation. Keep the
fixture name unique and update an expectation only after confirming the change
is intentional.

`pathway_cases.tsv` selects regression cases that must also match a reviewed
JSON document in `expected_pathways/`. These comparisons parse JSON before
comparing it, so formatting changes do not invalidate a structurally identical
pathway. The selected cases cover numerical and named molfiles as well as the
plain graph-input format.

Before running molecule cases, `unitTester.py` also checks the executable's CLI:
help text, strict option validation, compatibility names, option/input ordering,
runtime and enumeration boundaries, hydrogen removal, disconnected-component
compensation, output toggles and write failures, `.mol` path handling, and Linux
memory reporting.

Run only the pathway golden cases with:

```bash
python unitTests/unitTester.py --build --pathways-only --verbose
```

`--build` targets x86-64-v3 with POPCNT and enables search telemetry so the CLI
and JSON schema checks can run. Aside from that diagnostic feature it matches
the repository's documented optimized build. A standard executable still runs
the ordinary CLI and regression checks; telemetry-only checks are detected and
skipped. Use a separately compiled portable executable on older x86-64 or
non-x86 systems.

The build path also runs focused active-word mask, iterative tree-
canonicalisation, and cyclic-canonicalisation tests. The cyclic tests compare
the exact coloured-2-core form with an independent exact labelled-multigraph
matcher across vertex permutations, atom and bond labels, pendant trees,
disconnected fallbacks, multiedges, and regular graphs that colour refinement
alone cannot distinguish. On Linux the tree tests use a 64 KiB stack and
include 32,766- and 32,767-node labelled paths, preventing recursive traversal
from returning unnoticed for either centroid shape.

Run the manifest and fixture integrity audit with:

```bash
python unitTests/unitTester.py --audit
```

Some `.mol` files are retained as a fixture library but do not yet have reviewed
expected values. The audit reports these as fixture-only molecules; they are not
silently treated as passing regression cases.
