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
output toggles, `.mol` path handling, and Linux memory reporting.

Run only the pathway golden cases with:

```bash
python unitTests/unitTester.py --build --pathways-only --verbose
```

Run the manifest and fixture integrity audit with:

```bash
python unitTests/unitTester.py --audit
```

Some `.mol` files are retained as a fixture library but do not yet have reviewed
expected values. The audit reports these as fixture-only molecules; they are not
silently treated as passing regression cases.
