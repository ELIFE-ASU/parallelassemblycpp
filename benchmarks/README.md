# Benchmark corpus

`cases.tsv` defines the maintained AssemblyCpp performance corpus. Its columns
are:

| Column | Meaning |
| --- | --- |
| `name` | Unique CLI-safe case name. |
| `input` | Molfile or native graph path, relative to this directory unless absolute. |
| `expected_assembly_index` | Result required from every warm-up and measured calculation. |
| `expectation` | `reviewed` for regression-manifest values or `provisional` for profiling-only guards. |
| `suites` | Comma-separated membership in `quick`, `full`, `profile`, or `scaling`. |
| `workload` | Short description printed in the summary. |

The loader rejects malformed rows, duplicate names, unknown suites or
expectation statuses, missing inputs, and non-integer expected results.

## Suites

- `quick` covers an acyclic repeated chain, an automorphism-heavy graph, dense
  and symmetric cyclic graphs, and the existing ketoconazole baseline. One
  warm-up plus five rounds normally finishes in about a second.
- `full` contains the quick cases plus ten reviewed inputs covering native graph
  parsing and larger labelled or ring-rich molecules.
- `profile` contains a second-scale tree search and multi-second macrocycle
  searches. Except for ketoconazole, their expectations are provisional
  benchmark guards rather than reviewed regression expectations.
- `scaling` is a cumulative series of two through eight disconnected amino-acid
  components. It grows from 18 atoms and 16 bonds to 68 atoms and 60 bonds;
  every larger input contains the preceding graph and one additional component.
  The expectations are provisional benchmark guards.

Run the scaling series with:

```bash
python benchmarks/benchmark.py --suite scaling
```

The workload column puts atoms, bonds, and component counts next to the timing
statistics. Algorithm clock ticks are the clearest cost signal for the smaller
inputs because wall time also includes fixed process-startup and parsing costs.
The series measures this cumulative, changing-composition workload; it should
not be interpreted as an asymptotic complexity measurement for arbitrary
molecules.

## Measurement method

All calculations are serial and run in isolated temporary directories. Within
each round the first case is rotated to spread systematic order effects. In a
paired comparison, each baseline/candidate pair is adjacent; odd rounds use AB
order and even rounds use BA order. Paired mode defaults to six measured rounds
so each executable occupies each position equally often; explicitly requesting
an odd count emits a warning.

Wall time includes process startup and input parsing. The program-reported
`std::clock` values begin after parsing and remain in platform-specific ticks;
for the shortest cases, process overhead can dominate wall time.

The unpaired human-readable report gives median, MAD, p95, and raw clock-domain
medians. Paired reports calculate each baseline/candidate ratio before taking
the median, rather than dividing independent medians. Their primary corpus
result is the median paired round-total speedup, which naturally weights cases
by duration and is less sensitive to noise in very short processes. A
descriptive equal-weight geometric mean is also reported, with MAD; treat it
cautiously when a suite contains millisecond-scale cases.

JSON schema version 1 retains every sample, case metadata, aggregate summaries,
the deterministic schedule, platform information, and executable fingerprints.
Fingerprints are captured before measurement and verified afterward so a binary
replaced during a run cannot be misattributed. No outliers are removed
automatically, and no timing result is used as a CI pass/fail threshold.

Run the benchmark-runner tests with:

```bash
python -m unittest discover -s benchmarks -p 'test_*.py'
```
