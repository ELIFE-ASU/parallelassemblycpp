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
- `scaling` combines two complementary series. The cumulative amino-acid series
  grows from 18 atoms and 16 bonds to 113 atoms and 100 bonds; every larger
  input contains the preceding graph and one additional component. The five
  largest cases add isoleucine, leucine, lysine, methionine, and valine. A
  homogeneous-path boundary series covers 63, 64, 65, 127, 128, and 129 bonds.
  Those cases cross the residual-decomposition cache's scalar-to-wide key
  boundary and the first two 64-bit mask-word boundaries without introducing
  an exponential graph family. Wide caching remains eligible at later word
  boundaries because masks now grow to the processed graph's width rather than
  stopping at 512 bits. The expectations are provisional benchmark guards.

Build the optimized candidate and telemetry sibling from the repository root:

```bash
cmake --preset performance
cmake --build --preset performance
```

Run the scaling series with:

```bash
python benchmarks/benchmark.py \
  --executable build/performance/AssemblyCpp --suite scaling
```

The largest cases deliberately exercise expensive scaling behaviour. The
100-bond endpoint can take several minutes and more than 1 GiB of memory on a
performance build; use `--runs 1 --warmup 0` for a smoke run. The default
per-calculation timeout is 600 seconds and can be adjusted with `--timeout`.

Collect scaling telemetry with one additional untimed candidate run per case:

```bash
python benchmarks/benchmark.py --telemetry \
  --executable build/performance/AssemblyCpp \
  --telemetry-executable build/performance/AssemblyCppTelemetry \
  --suite scaling \
  --json-output scaling.json
```

The `performance` preset creates the ordinary timed candidate and its
instrumented sibling used above. The diagnostic run is performed after all
timed rounds and is excluded from wall-clock, algorithm-clock, and paired
speedup aggregates. In paired mode only the separate telemetry executable
receives the diagnostic run, so an older baseline can still be used for
timings. As an alternative, the runner's `--build` option creates both files
directly at the requested executable path.

## LTO and PGO candidates

`performance-lto` adds link-time optimization to the x86-64-v3 performance
configuration. `performance-pgo` adds a GCC profile trained by the maintained
corpus as well as LTO. Build them from the repository root with:

```bash
cmake --preset performance-lto
cmake --build --preset performance-lto

cmake --preset pgo-generate
cmake --build --preset pgo-train
cmake --preset performance-pgo
cmake --build --preset performance-pgo
```

The `pgo-train` target reads the strict, complete weight mapping in
`pgo-training.tsv`. Each training calculation runs serially in isolation and
must reproduce its expected assembly index. The weights deliberately emphasize
the short cases that can otherwise be lost among the long profile workloads.
Running the target again removes existing `.gcda` files first, so a profile
cannot silently mix code generations. A completion record is written atomically
only after the full schedule succeeds and fingerprints the weights, manifest,
and every input; the profile-use preset rejects partial, interrupted, or stale
training data.

Keep an unchanged `performance` executable as the baseline and compare every
suite before promoting either optimization preset. This Linux example pins all
work to one known performance core; choose an appropriate core for the host:

```bash
for candidate in performance-lto performance-pgo; do
  taskset -c 2 python benchmarks/benchmark.py \
    --baseline-executable build/performance/AssemblyCpp \
    --executable "build/${candidate}/AssemblyCpp" \
    --suite quick --runs 100 \
    --json-output "build/${candidate}-quick.json"
  taskset -c 2 python benchmarks/benchmark.py \
    --baseline-executable build/performance/AssemblyCpp \
    --executable "build/${candidate}/AssemblyCpp" \
    --suite full --runs 100 \
    --json-output "build/${candidate}-full.json"
  taskset -c 2 python benchmarks/benchmark.py \
    --baseline-executable build/performance/AssemblyCpp \
    --executable "build/${candidate}/AssemblyCpp" \
    --suite profile --runs 6 \
    --json-output "build/${candidate}-profile.json"
  taskset -c 2 python benchmarks/benchmark.py \
    --baseline-executable build/performance/AssemblyCpp \
    --executable "build/${candidate}/AssemblyCpp" \
    --suite scaling --runs 30 \
    --json-output "build/${candidate}-scaling.json"
  python benchmarks/check_speedups.py \
    "build/${candidate}-quick.json" \
    "build/${candidate}-full.json" \
    "build/${candidate}-profile.json" \
    "build/${candidate}-scaling.json"
done
```

Use an even run count so AB and BA ordering remains balanced. The gate requires
every case's paired algorithm-clock median and every suite's paired round-total
wall and clock medians to exceed 1.0. It recomputes those medians from the raw
paired samples and rejects reports whose manifest or input fingerprints no
longer match the maintained corpus. If a small case fails, repeat it with at
least 100 rounds. For a persistent PGO regression, increase the relevant
training weight, regenerate the profile, and repeat all suites; an LTO-only
failure instead requires flag tuning and the same full revalidation. Timing
remains a host-pinned manual release gate rather than a shared-CI assertion
because contention would make such a check unreliable. A failed gate blocks
promotion or default use, but does not prevent keeping the corresponding preset
available as an explicit experimental opt-in.

The workload column puts atoms, bonds, component counts, cache eligibility, and
active mask-word counts next to the relevant timing statistics. Algorithm clock
ticks are the clearest cost signal for the smaller inputs because wall time also
includes fixed process-startup and parsing costs. The two series measure a
cumulative changing-composition workload and isolated representation
boundaries; neither should be interpreted as an asymptotic complexity
measurement for arbitrary molecules.

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

JSON schema version 2 retains every timed sample, case metadata, aggregate
summaries, the deterministic schedule, platform information, executable,
manifest, and input fingerprints, and optional candidate telemetry. Telemetry
includes the raw
processed graph size, retained masks, matching visits, canonicalisation activity,
legacy VF2 counters, canonical/residual/assembly/pair-bound cache rates, and
absolute resident peak memory for setup, initial enumeration, DAG conversion,
assembly search, and output. A cache rate is `null` when no lookup occurred.

`retained_masks` counts unique initial-DAG masks accepted within the enumeration
limit, including one-edge roots. `matching_visits` counts valid materialised
duplicate pairs delivered to the search visitor. `canonicalisation_calls`
includes canonical-mask cache hits. `vf2_calls` and `vf2_matches` are retained for
schema compatibility and remain zero now that exact cyclic canonical codes replace
VF2. The residual lookup hit rate excludes first-occurrence,
small-residual, and disabled-cache bypasses; the request hit rate includes them.
Wide caches sample up to 6,144 residuals and return to direct decomposition when
no hit demonstrates reuse, so unique workloads do not keep paying lookup cost.
Fingerprints are captured before measurement and verified afterward so a binary
replaced during a run cannot be misattributed. No outliers are removed
automatically, and no timing result is used as a CI pass/fail threshold.

Run the benchmark-runner tests with:

```bash
python -m unittest discover -s benchmarks -p 'test_*.py'
```
