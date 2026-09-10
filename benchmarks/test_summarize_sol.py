from __future__ import annotations

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from benchmarks import benchmark, summarize_sol


class SummarizeSolTests(unittest.TestCase):
    def write_report(
        self,
        root: Path,
        *,
        parallel: bool,
        mpi_ranks: int = 1,
        threads_per_rank: int = 4,
    ) -> Path:
        if mpi_ranks > 1:
            path = root / f"hybrid/hybrid-{mpi_ranks}x{threads_per_rank}.json"
        elif parallel:
            path = root / f"parallel/omp-{threads_per_rank}.json"
        else:
            path = root / "suites/full.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        suite = "profile" if parallel else "full"
        results = [
            benchmark.CaseResult(
                benchmark.BenchmarkCase(
                    name, root / f"{name}.mol", 23, "provisional", (suite,), "fixture"
                ),
                tuple(benchmark.Measurement(wall, 100, 23) for wall in candidate),
                (
                    tuple(benchmark.Measurement(wall, 100, 23) for wall in baseline)
                    if parallel
                    else ()
                ),
            )
            for name, candidate, baseline in (
                ("alpha", (1.0, 9.0), (9.0, 9.0)),
                ("beta,quoted", (2.0, 2.0), (4.0, 20.0)),
            )
        ]
        benchmark.write_json_report(
            path=path,
            candidate_metadata={"sha256": "candidate"},
            baseline_metadata={"sha256": "baseline"} if parallel else None,
            telemetry_metadata=None,
            corpus_metadata={
                "manifest": {"sha256": "manifest"},
                "inputs": [
                    {"name": result.case.name, "sha256": f"input-{result.case.name}"}
                    for result in results
                ],
            },
            manifest=root / "cases.tsv",
            suite=suite,
            runs=2,
            warmup=0,
            timeout=60.0,
            results=results,
            candidate_execution=benchmark.ExecutionConfig(
                launcher=("mpirun", "-n", str(mpi_ranks)) if mpi_ranks > 1 else (),
                environment=(
                    ("OMP_NUM_THREADS", str(threads_per_rank) if parallel else "1"),
                ),
                arguments=("--parallel=on" if parallel else "--parallel=off",),
            ),
            baseline_execution=(
                benchmark.ExecutionConfig(
                    environment=(("OMP_NUM_THREADS", "1"),),
                    arguments=("--parallel=off",),
                )
                if parallel
                else None
            ),
        )
        return path

    def read_summary(self, root: Path) -> list[dict[str, str]]:
        with (root / "summary.csv").open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(reader.fieldnames, list(summarize_sol.FIELDNAMES))
            return list(reader)

    def test_partial_and_combined_reports_use_paired_ratios_and_round_totals(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_report(root, parallel=False)
            self.assertEqual(summarize_sol.write_summary(root)[1], 3)
            serial = self.read_summary(root)
            self.assertEqual(serial[0]["run_kind"], "suite")
            self.assertEqual(serial[0]["topology"], "serial")
            self.assertEqual(serial[0]["workers"], "1")
            self.assertEqual(serial[0]["mpi_ranks"], "1")
            self.assertEqual(serial[0]["threads_per_rank"], "1")
            self.assertEqual(serial[0]["paired_speedup"], "")
            self.assertEqual(serial[0]["baseline_wall_median_seconds"], "")
            self.assertEqual(serial[1]["case"], "beta,quoted")
            self.assertEqual(serial[-1]["case"], "__suite__")
            self.assertEqual(float(serial[-1]["wall_median_seconds"]), 7.0)

            self.write_report(root, parallel=True)
            self.assertEqual(summarize_sol.write_summary(root)[1], 6)
            rows = self.read_summary(root)
            alpha, _, aggregate = rows[3:]
            self.assertEqual(alpha["run_kind"], "parallel")
            self.assertEqual(alpha["report"], "parallel/omp-4.json")
            self.assertEqual(alpha["topology"], "omp")
            self.assertEqual(alpha["workers"], "4")
            self.assertEqual(alpha["mpi_ranks"], "1")
            self.assertEqual(alpha["threads_per_rank"], "4")
            self.assertEqual(float(alpha["wall_median_seconds"]), 5.0)
            self.assertEqual(float(alpha["wall_mad_seconds"]), 4.0)
            self.assertAlmostEqual(float(alpha["wall_p95_seconds"]), 8.6)
            self.assertEqual(float(alpha["baseline_wall_median_seconds"]), 9.0)
            self.assertEqual(float(alpha["paired_speedup"]), 5.0)
            self.assertEqual(float(alpha["efficiency"]), 1.25)
            self.assertAlmostEqual(
                float(aggregate["paired_speedup"]), (13 / 3 + 29 / 11) / 2
            )

            self.write_report(root, parallel=True, mpi_ranks=2, threads_per_rank=64)
            self.assertEqual(summarize_sol.write_summary(root)[1], 9)
            combined = self.read_summary(root)
            self.assertEqual(combined[:6], rows)
            hybrid_alpha, _, hybrid_aggregate = combined[6:]
            self.assertEqual(hybrid_alpha["run_kind"], "parallel")
            self.assertEqual(hybrid_alpha["report"], "hybrid/hybrid-2x64.json")
            self.assertEqual(hybrid_alpha["topology"], "hybrid")
            self.assertEqual(hybrid_alpha["workers"], "128")
            self.assertEqual(hybrid_alpha["mpi_ranks"], "2")
            self.assertEqual(hybrid_alpha["threads_per_rank"], "64")
            self.assertEqual(float(hybrid_alpha["paired_speedup"]), 5.0)
            self.assertEqual(float(hybrid_alpha["efficiency"]), 5.0 / 128)
            self.assertEqual(
                hybrid_aggregate["paired_speedup"], aggregate["paired_speedup"]
            )
            self.assertAlmostEqual(
                float(hybrid_aggregate["efficiency"]),
                float(aggregate["paired_speedup"]) / 128,
            )
            self.assertEqual(list(root.glob("tmp*")), [])

    def test_parallel_only_and_empty_results_are_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.assertEqual(summarize_sol.write_summary(root)[1], 0)
            self.assertEqual(self.read_summary(root), [])
            self.write_report(root, parallel=True)
            self.assertEqual(summarize_sol.write_summary(root)[1], 3)

    def test_invalid_reports_preserve_previous_summary(self) -> None:
        for corruption in (
            "schema",
            "workers",
            "mpi",
            "samples",
            "assembly_index",
            "paired_median",
        ):
            with (
                self.subTest(corruption=corruption),
                tempfile.TemporaryDirectory() as tmp,
            ):
                root = Path(tmp)
                path = self.write_report(root, parallel=True)
                summarize_sol.write_summary(root)
                original = (root / "summary.csv").read_bytes()
                document = json.loads(path.read_text(encoding="utf-8"))
                if corruption == "schema":
                    document["schema_version"] = 1
                elif corruption == "workers":
                    document["execution"]["candidate"]["environment"][
                        "OMP_NUM_THREADS"
                    ] = "8"
                elif corruption == "mpi":
                    document["execution"]["candidate"]["launcher"] = [
                        "mpirun",
                        "-n",
                        "2",
                    ]
                elif corruption == "samples":
                    document["cases"][0]["candidate"]["measurements"].pop()
                elif corruption == "assembly_index":
                    document["cases"][0]["candidate"]["measurements"][0][
                        "assembly_index"
                    ] = 99
                else:
                    document["cases"][0]["comparison"]["paired_wall_speedup"][
                        "median"
                    ] = 1.8
                path.write_text(json.dumps(document), encoding="utf-8")
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    self.assertEqual(summarize_sol.main([str(root)]), 1)
                self.assertIn(str(path), stderr.getvalue())
                self.assertEqual((root / "summary.csv").read_bytes(), original)

    def test_serial_report_rejects_parallel_workers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.write_report(root, parallel=False)
            document = json.loads(path.read_text(encoding="utf-8"))
            document["execution"]["candidate"]["environment"]["OMP_NUM_THREADS"] = "2"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(
                summarize_sol.check_parallel_scaling.ScalingError, "serial measurements"
            ):
                summarize_sol.write_summary(root)

    def test_hybrid_uses_effective_thread_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = self.write_report(
                root, parallel=True, mpi_ranks=2, threads_per_rank=64
            )
            document = json.loads(path.read_text(encoding="utf-8"))
            document["execution"]["candidate"]["environment"].update(
                {"OMP_NUM_THREADS": "128", "OMP_THREAD_LIMIT": "64"}
            )
            path.write_text(json.dumps(document), encoding="utf-8")
            self.assertEqual(summarize_sol.write_summary(root)[1], 3)
            row = self.read_summary(root)[0]
            self.assertEqual(row["workers"], "128")
            self.assertEqual(row["mpi_ranks"], "2")
            self.assertEqual(row["threads_per_rank"], "64")

    def test_invalid_hybrid_reports_preserve_previous_summary(self) -> None:
        scenarios = (
            ("hybrid-1x64.json", 1, 64, 64, "invalid hybrid report filename"),
            ("hybrid-2x1.json", 2, 1, 1, "invalid hybrid report filename"),
            ("hybrid-invalid.json", 2, 64, 64, "invalid hybrid report filename"),
            ("hybrid-2x64.json", 1, 64, 64, "hybrid topology mismatch"),
            ("hybrid-2x64.json", 4, 32, 32, "hybrid topology mismatch"),
            ("hybrid-2x64.json", 2, 32, 32, "hybrid topology mismatch"),
            ("hybrid-2x64.json", 2, 64, 32, "hybrid topology mismatch"),
        )
        for filename, ranks, threads, limit, expected_error in scenarios:
            with (
                self.subTest(
                    filename=filename, ranks=ranks, threads=threads, limit=limit
                ),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                self.write_report(root, parallel=False)
                summarize_sol.write_summary(root)
                original = (root / "summary.csv").read_bytes()
                path = self.write_report(
                    root, parallel=True, mpi_ranks=2, threads_per_rank=64
                )
                document = json.loads(path.read_text(encoding="utf-8"))
                execution = document["execution"]["candidate"]
                execution["launcher"] = ["srun", f"--ntasks={ranks}"]
                execution["environment"].update(
                    {"OMP_NUM_THREADS": str(threads), "OMP_THREAD_LIMIT": str(limit)}
                )
                path = path.rename(path.with_name(filename))
                path.write_text(json.dumps(document), encoding="utf-8")
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    self.assertEqual(summarize_sol.main([str(root)]), 1)
                self.assertIn(expected_error, stderr.getvalue())
                self.assertIn(str(path), stderr.getvalue())
                self.assertEqual((root / "summary.csv").read_bytes(), original)
                self.assertEqual(list(root.glob("tmp*")), [])


if __name__ == "__main__":
    unittest.main()
