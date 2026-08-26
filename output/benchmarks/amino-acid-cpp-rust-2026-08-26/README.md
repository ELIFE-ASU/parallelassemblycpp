# Amino-acid C++ / Rust benchmark

This report compares the shipped AssemblyCpp reference executable, the current
C++ branch, and the DaymudeLab Rust implementation on the maintained amino-acid
scaling inputs from 2 through 10 components.

## Main result

At 10 components (86 atoms, 76 bonds), all runs returned assembly index 27:

| Variant | Workers | Median wall time | MAD | Relative to shipped reference |
| --- | ---: | ---: | ---: | ---: |
| Shipped reference, serial | 1 | 5.729 s | 0.009 s | 1.00× |
| Current C++ serial | 1 | 2.392 s | 0.003 s | 2.39× faster |
| Rust serial | 1 | 39.950 s | 0.178 s | 6.97× slower |
| Current C++ OpenMP | 1 | 2.436 s | 0.009 s | 2.35× faster |
| Rust depth-one | 1 | 41.317 s | 0.081 s | 7.21× slower |
| Current C++ OpenMP | 2 | 1.963 s | 0.006 s | 2.92× faster |
| Rust depth-one | 2 | 17.952 s | 0.031 s | 3.13× slower |
| Current C++ OpenMP | 4 | 1.468 s | 0.003 s | 3.90× faster |
| Rust depth-one | 4 | 9.198 s | 0.035 s | 1.61× slower |
| Current C++ OpenMP | 8 | 1.154 s | 0.003 s | 4.97× faster |
| Rust depth-one | 8 | 5.227 s | 0.070 s | 1.10× faster |

The current C++ serial build is 16.70× faster than Rust serial at 10
components. At eight workers, current C++ OpenMP is 4.53× faster than Rust
depth-one. Rust's eight-worker mode nevertheless closes most of its serial gap
and narrowly beats the serial shipped reference.

Across the full 2–10 component range, using the sum of per-case medians:

- Current C++ serial is 2.37× faster than the shipped reference.
- Rust serial is 6.26× slower than the shipped reference and 14.84× slower
  than current C++ serial.
- Current C++ OpenMP at eight workers is 5.05× faster than the shipped
  reference.
- Rust depth-one at eight workers is 1.18× faster than the shipped reference.
- Current C++ OpenMP at eight workers is 4.29× faster than Rust depth-one at
  eight workers.

## Parallel scaling at 10 components

| Workers | Current C++ time | C++ speedup | C++ efficiency | Rust time | Rust speedup | Rust efficiency |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 2.436 s | 1.00× | 100.0% | 41.317 s | 1.00× | 100.0% |
| 2 | 1.963 s | 1.24× | 62.1% | 17.952 s | 2.30× | 115.1% |
| 4 | 1.468 s | 1.66× | 41.5% | 9.198 s | 4.49× | 112.3% |
| 8 | 1.154 s | 2.11× | 26.4% | 5.227 s | 7.90× | 98.8% |

Rust's apparent superlinear scaling at two and four workers is not pure CPU
efficiency: its parallel search shares a best-index bound, and scheduling can
change how much of the search tree is visited. The Rust project documents the
number of states searched as nondeterministic when parallelism is enabled.

## Method

- External CLI wall time from process spawn through exit, measured with
  `time.perf_counter_ns()`; input staging and result parsing are outside the
  timed region.
- One warm-up and five measured runs per case and variant.
- Variant order is cyclically rotated for each case and round.
- Every launch uses a fresh temporary working directory.
- Every sample requires a zero exit status and the exact expected assembly
  index.
- Serial variants are pinned to logical CPU 0. Parallel variants use one
  logical thread from each P-core: `0`, `0,2`, `0,2,4,6`, and
  `0,2,4,6,8,10,12,14` for 1, 2, 4, and 8 workers.
- Current OpenMP uses fixed worker counts, `OMP_DYNAMIC=FALSE`,
  `OMP_PROC_BIND=close`, and forces the parallel path on every case.
- Rust uses `--parallel none` for true serial and `--parallel depth-one` with
  `RAYON_NUM_THREADS` for worker scaling.
- Cases 2–5 are dominated by CLI startup and should not be used to judge solver
  scaling.
- The repeated matrix stops at 10 components because the Rust serial path is
  already about 40 seconds per launch there; balanced repeated 11-component
  serial runs would disproportionately extend the benchmark.

## Provenance and comparability

- Host: Intel Core i7-14700KF, Linux x86-64.
- Current C++ commit: `79696f6263e54aad0f1f10532de58f5885fafd19`.
- Current serial SHA-256: `ded5a7534a7c7693967d518837e1bcd2da1a607234c9523b49cfabde10f22be0`.
- Current OpenMP SHA-256: `7350a51060170ec7fb8e6829274730a37ccb613bfa2941355a47d365be51ee9e`.
- Shipped reference SHA-256: `571e70707783e637d65b4b4b35e96a6f0138e71c37ae09fec782752af54caedb`.
- Rust source: [DaymudeLab/assembly-theory](https://github.com/DaymudeLab/assembly-theory),
  commit `19370f85eed6784d4aa9a4e5be017d0c5126c2c7` on `main`, package version
  0.6.1.
- Rust binary SHA-256: `69194e7e6572699711445faca4c82f139e31c0fcec9f5ad205ef715fe95f2183`.
- Rust toolchain: rustc 1.98.0; build command `cargo build --locked --release`.
- Current C++ is a Release `x86-64-v3` build without LTO. Rust uses Cargo's
  generic release target with no `target-cpu` override. Consequently this is a
  comparison of the requested project builds, not a compiler-controlled
  language comparison.

Files:

- `benchmark.json`: raw samples, summaries, hashes, commands' runtime
  configuration, and source pins.
- `summary.csv`: one median/MAD/p95 row per case and variant.
- `run_benchmark.py`: reproducible neutral benchmark harness.
