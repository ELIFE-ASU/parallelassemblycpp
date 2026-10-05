# Audits and experiments

These reports preserve evidence from specific dates, revisions, machines, and
test corpora. Their measurements and validation counts describe those runs.
Recommendations and references to the "current" implementation inside a dated
report refer to its original snapshot unless a later context note says otherwise.

Use the [project README](../README.md), [CLI reference](../docs/cli.md), and
[algorithm guide](../docs/algorithm.md) for current behavior. In particular:

- Exact molecular/graph search runs without a graph Re-Pair prepass. The
  explicit `--algorithm=re-pair` mode returns a constructive upper bound.
- Exact string search uses a Re-Pair construction as an initial incumbent and
  retains its pathway witness, added in commit `8d4585c`. This remains distinct
  from the removed molecular seeding experiment.
- The GPU documents are research issue drafts. They do not document an
  available GPU executable or promise implementation of the proposed work.

| Date | Report | Scope and status |
| --- | --- | --- |
| 28 September 2026 | [Branch-and-bound audit](2026-09-28-branch-and-bound/REPORT.md) | Audit of revision `815e5e6`, an independent small-graph oracle, and research recommendations. The recommendation to seed string search has since been implemented. |
| 29 September 2026 | [Graph Re-Pair bound quality](2026-09-29-graph-repair/REPORT.md) | Historical fixture comparison, independently checked certificates, and local heuristic timings. |
| 29 September 2026 | [Graph bound versus exact-search timing](2026-09-29-graph-repair/SPEED.md) | Cost of obtaining a bound compared with proving an exact result; these ratios do not measure an exact-solver speedup. |
| 29 September 2026 | [Molecular Re-Pair seeding experiment](2026-09-29-graph-repair/SEED.md) | Removed experiment that seeded exact graph search. Its negative timing result and validation data are retained. |
| 29 September 2026 | [Pathway reconstruction](2026-09-29-pathway-reconstruction/REPORT.md) | Investigation of reconstruction costs, witness storage, and local before/after timing evidence. |
| 29 September 2026 | [GPU investigation drafts](2026-09-29-gpu-branch-and-bound/README.md) | Ten proposed investigations, with repository entry points checked against revision `04129ee13bdaf165001afddecc22ea66b75a5285`. |

Tracked scripts, JSON, CSV, and figures beside the reports are the published
evidence. References to `build/` describe local, unpublished artifacts that are
not distributed with the repository; a fresh checkout cannot reproduce a removed
candidate merely by running those paths. Source links resolve to the current
checkout, so consult a report's recorded revision or source fingerprints when
reproducing historical results.

For new measurements, record fresh source and executable fingerprints, keep
the original audit data intact, and write results to a separate output directory.
Use the current [benchmark protocol](../benchmarks/README.md) and
[test instructions](../unitTests/README.md) when validating a new change.
