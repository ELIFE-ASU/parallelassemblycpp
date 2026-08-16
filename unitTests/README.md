# Test data

`regression_cases.tsv` is the canonical regression manifest. Each row names a
fixture in this directory and its reviewed assembly-index expectation. Keep the
fixture name unique and update an expectation only after confirming the change
is intentional.

Run the manifest and fixture integrity audit with:

```bash
python unitTests/unitTester.py --audit
```

Some `.mol` files are retained as a fixture library but do not yet have reviewed
expected values. The audit reports these as fixture-only molecules; they are not
silently treated as passing regression cases.
