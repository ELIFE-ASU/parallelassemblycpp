# Task-transfer benchmark results

This archived report predates the rename to `parallelassemblycpp`. Executable
names, build options, and paths below retain the names used for these measurements.

**The all-benchmarks-faster requirement is not met by these measurements.**

The reports contain 5,140 timed executions across 37 distinct maintained cases. All completed measurements returned the expected assembly index. All four suites meet their documented promotion repetition counts.

## Correctness checks

All 23 CTest groups passed, including all 938 serial regression cases, the 13 golden pathway checks, 74 benchmark-tool tests, and OpenMP, MPI, and hybrid parity/telemetry tests. The focused transfer tests additionally passed with AddressSanitizer and UndefinedBehaviorSanitizer; leak detection was disabled because of the sandbox ptrace restriction. The final mixed-pruning calibration assertion passed after the full CTest run.

[CTest log](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/ctest-full.log)

## Method

Fresh GCC 15.2 Release baseline and candidate OpenMP executables use matching x86-64-v3 compiler flags. Both execute with eight threads pinned to CPUs 0,2,4,6,8,10,12,14, OMP_DYNAMIC=FALSE, OMP_PLACES=threads, and OMP_PROC_BIND=close. Baseline and candidate alternate AB/BA order and case order rotates between rounds. Telemetry runs are separate and excluded from timing. Wall time includes process startup and input parsing. Warm-up rounds per executable/case: quick=1, full=1, profile=1, scaling=1.

Speedup is baseline/candidate; values greater than 1 are quicker. Suite values are medians of paired round-total ratios. Case values are medians of paired case ratios. A ratio of the two displayed wall medians can differ from the median of the paired ratios. CPU-clock time is reported separately from elapsed time.

The serial executables have identical `.text` machine-code sections (SHA256 `8ffe87452064e62d71b96713984a784616006b5e76606756be86865788c7357a`) and identical section sizes. The transfer change does not change serial machine instructions. MPI and hybrid correctness checks do not establish performance on those topologies.

## Suite results

| Suite | Rounds | Wall speedup | CPU-clock speedup | Wall cases >1 | Clock cases >1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| quick | 100 | 0.9973× | 1.0007× | 2/5 | 2/5 |
| full | 100 | 1.0009× | 0.9980× | 8/15 | 6/15 |
| profile | 6 | 1.0024× | 1.0033× | 5/5 | 5/5 |
| scaling | 30 | 1.0040× | 1.0006× | 9/18 | 10/18 |

19/43 suite entries have wall speedup at or below 1; 20/43 have CPU-clock speedup at or below 1. Small differences may reflect timing variation; they do not justify an all-faster claim. The repository promotion gate requires every case CPU-clock median and every suite wall/CPU-clock median to exceed 1; case wall medians are also shown individually here.

[Exact gate output](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/promotion-v2/speedup-gate.log)

## Every measured case

| Suite | Case | Baseline wall median ms | Candidate wall median ms | Paired wall speedup | Paired CPU-clock speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| quick | icosane | 1.821 | 1.741 | 1.0190× | 0.9793× |
| quick | sr1001 | 1.834 | 1.866 | 0.9793× | 0.9760× |
| quick | 5843 | 7.470 | 7.437 | 1.0094× | 0.9991× |
| quick | bisphenylmaleimide | 6.440 | 6.512 | 0.9971× | 1.0008× |
| quick | ketoconazole | 6.178 | 6.233 | 0.9952× | 1.0028× |
| full | icosane | 1.726 | 1.707 | 0.9944× | 1.0110× |
| full | sr1001 | 1.830 | 1.829 | 1.0107× | 0.9609× |
| full | 5843 | 7.509 | 7.483 | 1.0016× | 0.9999× |
| full | bisphenylmaleimide | 6.531 | 6.473 | 1.0060× | 0.9998× |
| full | ketoconazole | 6.264 | 6.287 | 0.9974× | 1.0050× |
| full | THC | 2.864 | 2.880 | 1.0047× | 0.9958× |
| full | sucrose | 3.389 | 3.382 | 1.0102× | 0.9899× |
| full | graphio | 2.993 | 2.953 | 1.0173× | 1.0027× |
| full | Ceftiolene | 4.388 | 4.464 | 0.9846× | 0.9879× |
| full | Cefquinome | 3.533 | 3.544 | 0.9965× | 1.0019× |
| full | Cefpirome | 3.471 | 3.475 | 0.9943× | 1.0004× |
| full | dipyridamole | 4.694 | 4.709 | 1.0007× | 1.0018× |
| full | ceftiofur | 2.814 | 2.857 | 0.9845× | 0.9970× |
| full | dienogest | 4.009 | 4.021 | 1.0040× | 0.9988× |
| full | folic_acid | 2.709 | 2.746 | 0.9986× | 0.9894× |
| profile | ketoconazole | 7.009 | 6.767 | 1.0357× | 1.0086× |
| profile | phosphatidylcholine | 77.530 | 77.728 | 1.0044× | 1.0022× |
| profile | erythromycin | 255.882 | 256.319 | 1.0006× | 1.0018× |
| profile | clarithromycin | 777.670 | 772.484 | 1.0022× | 1.0023× |
| profile | paclitaxel | 42921.360 | 43010.314 | 1.0021× | 1.0030× |
| scaling | amino-acid-scale-02c | 2.277 | 2.272 | 0.9631× | 0.8723× |
| scaling | amino-acid-scale-03c | 2.090 | 2.297 | 0.9635× | 1.0056× |
| scaling | amino-acid-scale-04c | 2.524 | 2.419 | 1.0392× | 0.9950× |
| scaling | amino-acid-scale-05c | 2.557 | 2.571 | 1.0138× | 1.0163× |
| scaling | amino-acid-scale-06c | 3.137 | 3.159 | 0.9941× | 0.9837× |
| scaling | amino-acid-scale-07c | 5.372 | 5.422 | 0.9758× | 0.9966× |
| scaling | amino-acid-scale-08c | 11.780 | 11.606 | 1.0045× | 1.0055× |
| scaling | amino-acid-scale-09c | 52.800 | 52.627 | 0.9994× | 1.0001× |
| scaling | amino-acid-scale-10c | 222.082 | 216.901 | 1.0031× | 1.0037× |
| scaling | amino-acid-scale-11c | 960.456 | 966.785 | 0.9951× | 0.9980× |
| scaling | amino-acid-scale-12c | 6048.434 | 5946.716 | 1.0066× | 1.0067× |
| scaling | amino-acid-scale-13c | 23326.323 | 23280.502 | 1.0043× | 0.9989× |
| scaling | mask-boundary-path-063b | 4.344 | 4.401 | 0.9902× | 1.0054× |
| scaling | mask-boundary-path-064b | 4.238 | 3.925 | 1.0249× | 1.0052× |
| scaling | mask-boundary-path-065b | 4.290 | 4.428 | 0.9925× | 1.0001× |
| scaling | mask-boundary-path-127b | 24.394 | 24.034 | 1.0163× | 1.0075× |
| scaling | mask-boundary-path-128b | 18.266 | 18.017 | 0.9994× | 0.9945× |
| scaling | mask-boundary-path-129b | 18.653 | 18.649 | 1.0014× | 0.9912× |

## Secondary repeated checks

These two cases received 20 additional paired rounds after the exploratory run. They use the same binaries and placement but remain separate from the four promotion reports. Both showed slightly slower candidate medians; this does not support an improvement claim.

| Case | Rounds | Wall speedup | Wall speedup MAD | CPU-clock speedup |
| --- | ---: | ---: | ---: | ---: |
| clarithromycin | 20 | 0.9963× | 0.0040× | 0.9950× |
| amino-acid-scale-11c | 20 | 0.9962× | 0.0040× | 0.9960× |

[Secondary raw report](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/targeted-v2.json)

## Transfer measurements

Measurements below are separate untimed candidate executions. Task counts and thresholds vary with thread scheduling. The minimum-work column is the largest calibrated threshold observed during that run, not a fixed global constant. Work units sum edge-pair counts over fragments. Calibration samples 16 executions per producer and targets eight times the serialization cost in expected useful work, discounted by the immediate-pruning fraction. Buffer counts measure fresh or pooled acquisitions; a reused buffer can still need capacity growth, so they are not allocator-call counts.

| Case | Tasks | Serialization ns/task | Execution ns/task | Immediately pruned | Buffer reuse | Max minimum work units | Rejected small tasks |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| amino-acid-scale-13c | 5,553 | 139 | 5467147 | 39.2% | 49.8% | 524 | 2,351 |
| amino-acid-scale-12c | 2,535 | 128 | 2373951 | 36.3% | 53.0% | 434 | 1,226 |
| amino-acid-scale-11c | 2,257 | 134 | 406552 | 50.8% | 58.2% | 362 | 2,142 |
| amino-acid-scale-10c | 1,297 | 140 | 158574 | 54.0% | 62.5% | 463 | 317 |
| amino-acid-scale-09c | 1,017 | 124 | 85346 | 45.2% | 65.7% | 94 | 9 |
| amino-acid-scale-08c | 763 | 117 | 24602 | 55.0% | 58.8% | 182 | 85 |
| amino-acid-scale-07c | 247 | 193 | 20547 | 55.1% | 34.4% | 950 | 0 |
| amino-acid-scale-06c | 441 | 123 | 3564 | 69.6% | 66.9% | 188 | 76 |
| amino-acid-scale-05c | 99 | 283 | 4392 | 61.6% | 43.4% | 336 | 4 |
| amino-acid-scale-04c | 93 | 275 | 2784 | 55.9% | 33.3% | 143 | 14 |
| amino-acid-scale-03c | 2 | 1090 | 2582 | 0.0% | 0.0% | 0 | 0 |
| amino-acid-scale-02c | 3 | 401 | 368 | 66.7% | 0.0% | 0 | 0 |
| paclitaxel | 120,258 | 107 | 119894 | 41.3% | 70.5% | 1,606 | 7,057 |
| clarithromycin | 1,383 | 97 | 23845 | 38.3% | 66.3% | 1,496 | 30 |
| erythromycin | 539 | 103 | 14181 | 45.5% | 70.5% | 162 | 6 |
| phosphatidylcholine | 1,338 | 179 | 19177 | 46.9% | 62.9% | 669 | 1,233 |
| ketoconazole | 148 | 259 | 6082 | 50.7% | 32.4% | 620 | 24 |
| folic_acid | 137 | 193 | 10524 | 11.7% | 19.0% | 28 | 0 |
| dienogest | 17 | 387 | 2506 | 64.7% | 17.6% | 0 | 0 |
| ceftiofur | 154 | 275 | 9579 | 1.9% | 7.8% | 67 | 4 |
| dipyridamole | 3 | 651 | 378 | 100.0% | 66.7% | 0 | 0 |
| Cefpirome | 180 | 222 | 8938 | 1.1% | 24.4% | 72 | 0 |
| Cefquinome | 238 | 199 | 9661 | 5.0% | 13.4% | 85 | 10 |
| Ceftiolene | 301 | 252 | 4742 | 8.3% | 66.4% | 119 | 6 |
| THC | 53 | 229 | 3183 | 50.9% | 13.2% | 35 | 0 |
| sr1001 | 21 | 711 | 4541 | 23.8% | 33.3% | 0 | 0 |

## Canonicalisation comparison

Baseline and candidate telemetry below comes from separate untimed runs with the same eight-worker placement. Scheduling and pruning change the exact counts. These counters report canonicalisation activity; scheduling changes also affect them, so they do not isolate the effect of preserved IDs.

| Case | Baseline canonicalisation calls | Candidate calls | Baseline mask misses | Candidate mask misses | Fresh buffers | Reused buffers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| clarithromycin | 1,623,260 | 1,621,674 | 63,652 | 63,730 | 466 | 917 |
| amino-acid-scale-11c | 3,471,875 | 3,440,526 | 1,600 | 1,297 | 943 | 1,314 |
| amino-acid-scale-13c | 108,029,304 | 107,938,024 | 2,044 | 1,606 | 2,789 | 2,764 |

[Baseline telemetry](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-before/telemetry-omp8.json)

## Raw reports

Baseline source revision: `1ce92f95f37c6fb5c8ba848bc85049c327afd6f0`.

Reports include executable/input fingerprints, exact execution configuration, raw samples, medians, MAD, p95, and candidate telemetry.

- [quick](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/promotion-v2/quick.json)
- [full](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/promotion-v2/full.json)
- [profile](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/promotion-v2/profile.json)
- [scaling](/home/louie/skunkworks/assemblycpp-v5/build/task-transfer-validation/promotion-v2/scaling.json)

## Executables and reproduction

| Role | Executable | SHA256 |
| --- | --- | --- |
| candidate | `build/task-transfer-candidate/AssemblyCppOMP` | `4c2ae34c3271859218f0f6508a3bcd85974a15aebc5e2d843f84f7d3d8947c88` |
| baseline | `build/task-transfer-before/AssemblyCppOMP` | `290ab93f365c8bca8dc86768915d257bcccb9e115ec0bd50376add42e4ac9b79` |
| telemetry | `build/task-transfer-candidate/AssemblyCppOMPTelemetry` | `0eeaff9bf5a37a6cc86e88c697746400512f798ad66f86f492f913ecd0f462b7` |

Run from the repository root with the preserved baseline and candidate builds. The baseline sources were archived from the revision above; both builds use Release, GCC 15.2, and ASSEMBLYCPP_X86_64_V3=ON.

```bash
for taskBenchmarkSpec in quick:100 full:100 profile:6 scaling:30; do
  python3 benchmarks/benchmark.py \
    --baseline-executable build/task-transfer-before/AssemblyCppOMP \
    --executable build/task-transfer-candidate/AssemblyCppOMP \
    --telemetry-executable build/task-transfer-candidate/AssemblyCppOMPTelemetry \
    --baseline-parallel on --candidate-parallel on \
    --baseline-launcher "taskset -c 0,2,4,6,8,10,12,14" \
    --candidate-launcher "taskset -c 0,2,4,6,8,10,12,14" \
    --baseline-env OMP_NUM_THREADS=8 --candidate-env OMP_NUM_THREADS=8 \
    --baseline-env OMP_DYNAMIC=FALSE --candidate-env OMP_DYNAMIC=FALSE \
    --baseline-env OMP_PLACES=threads --candidate-env OMP_PLACES=threads \
    --baseline-env OMP_PROC_BIND=close --candidate-env OMP_PROC_BIND=close \
    --suite "${taskBenchmarkSpec%:*}" --runs "${taskBenchmarkSpec#*:}" \
    --warmup 1 --timeout 600 --telemetry \
    --json-output "build/task-transfer-validation/promotion-v2/${taskBenchmarkSpec%:*}.json"
done
python3 benchmarks/check_speedups.py \
  build/task-transfer-validation/promotion-v2/{quick,full,profile,scaling}.json
```
