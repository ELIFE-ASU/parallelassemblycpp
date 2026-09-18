from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from benchmarks import benchmark
from benchmarks import shared_cache_experiment as experiment


class SharedCacheExperimentTests(unittest.TestCase):
    def arguments(self, *extra: str) -> argparse.Namespace:
        return experiment.create_argument_parser().parse_args(
            [
                "--baseline-executable",
                "before",
                "--executable",
                "after",
                "--output-dir",
                "unused",
                *extra,
            ]
        )

    def test_default_corpus_has_each_case_once(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = experiment.main(
                [
                    "--baseline-executable",
                    "before",
                    "--executable",
                    "after",
                    "--output-dir",
                    "unused",
                    "--dry-run",
                ]
            )
        self.assertEqual(status, 0)
        plan = json.loads(output.getvalue().split("\n", 2)[2])
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        self.assertEqual(plan["cases"], [case.name for case in cases])
        self.assertEqual(len(plan["cases"]), len(set(plan["cases"])))
        self.assertEqual(
            plan["execution"]["selective"]["environment"][experiment.BYTES_ENV],
            str(256 * 1024 * 1024),
        )

    def test_variant_parsing_and_environment_overrides(self) -> None:
        variant = experiment.variant_specification("bounded=shared:1024")
        self.assertEqual(variant, experiment.Variant("bounded", "shared", 1024))
        arguments = self.arguments(
            "--threads",
            "8",
            "--env",
            "OMP_PLACES={2},{4}",
            "--variant-env",
            "bounded:CUSTOM=value with spaces",
            "--variant-env",
            f"bounded:{experiment.RESERVE_BYTES_ENV}=1572864",
            "--baseline-env",
            "REFERENCE=1",
        )
        configs = experiment.execution_configs(arguments, [variant])
        candidate = dict(configs["bounded"].environment)
        self.assertEqual(candidate[experiment.POLICY_ENV], "shared")
        self.assertEqual(candidate[experiment.BYTES_ENV], "1024")
        self.assertEqual(candidate["CUSTOM"], "value with spaces")
        self.assertEqual(candidate[experiment.RESERVE_BYTES_ENV], "1572864")
        self.assertEqual(
            dict(configs["baseline"].environment)[experiment.RESERVE_BYTES_ENV], "0"
        )
        self.assertEqual(candidate["OMP_NUM_THREADS"], "8")
        self.assertEqual(candidate["OMP_PLACES"], "{2},{4}")
        self.assertNotIn("REFERENCE", candidate)
        self.assertEqual(dict(configs["baseline"].environment)["REFERENCE"], "1")
        self.assertEqual(configs["bounded"].arguments, ("--parallel=auto",))

    def test_invalid_variants_and_environment_targets(self) -> None:
        for value in (
            "baseline=shared",
            "profiles=shared",
            "../escape=local",
            "name=unknown",
            "name=shared:x",
            "name=shared:",
        ):
            with (
                self.subTest(value=value),
                self.assertRaises(argparse.ArgumentTypeError),
            ):
                experiment.variant_specification(value)
        variant = experiment.Variant("local", "local")
        with self.assertRaisesRegex(benchmark.BenchmarkError, "distinct"):
            experiment.execution_configs(self.arguments(), [variant, variant])
        with self.assertRaisesRegex(benchmark.BenchmarkError, "unknown variant"):
            experiment.execution_configs(
                self.arguments("--variant-env", "typo:KEY=value"), [variant]
            )
        with self.assertRaisesRegex(benchmark.BenchmarkError, "duplicate"):
            experiment.execution_configs(
                self.arguments("--env", "KEY=1", "--env", "KEY=2"), [variant]
            )

    def test_profiles_use_separate_launcher_and_keep_raw_telemetry(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "case.mol"
            source.write_text("fixture", encoding="utf-8")
            case = benchmark.BenchmarkCase("case", source, 7, "reviewed", (), "test")
            execution = benchmark.ExecutionConfig(
                launcher=("taskset", "-c", "2,4"), environment=(("POLICY", "local"),)
            )
            telemetry = {"raw": "retained"}

            def profile(
                _executable: Path,
                prepared: benchmark.PreparedCase,
                _timeout: float,
                *,
                execution: benchmark.ExecutionConfig,
            ) -> dict[str, str]:
                self.assertEqual(execution.launcher[:3], ("/usr/bin/time", "-f", "%M"))
                self.assertEqual(execution.launcher[5:], ("taskset", "-c", "2,4"))
                self.assertEqual(execution.environment, (("POLICY", "local"),))
                Path(execution.launcher[4]).write_text("8192\n", encoding="utf-8")
                self.assertEqual(prepared.case.expected_assembly_index, 7)
                return telemetry

            with (
                mock.patch.object(benchmark, "run_telemetry_once", side_effect=profile),
                mock.patch.object(benchmark, "run_once") as timed,
            ):
                profiles = experiment.collect_profiles(
                    Path("telemetry"),
                    [case],
                    execution,
                    30.0,
                    Path("/usr/bin/time"),
                    telemetry=True,
                )
            timed.assert_not_called()
            self.assertEqual(
                profiles["case"], {"peak_rss_kib": 8192, "telemetry": telemetry}
            )

    def test_summary_preserves_growth_and_uses_paired_ratios(self) -> None:
        case = benchmark.BenchmarkCase("case", Path("input"), 7, "reviewed", (), "test")
        result = benchmark.CaseResult(
            case,
            (benchmark.Measurement(1.0, 10, 7), benchmark.Measurement(4.0, 40, 7)),
            (benchmark.Measurement(2.0, 20, 7), benchmark.Measurement(6.0, 60, 7)),
        )
        profile = {
            "peak_rss_kib": 123,
            "telemetry": {
                "counters": {"states_expanded": 42, "states_pruned": 17},
                "local_caches": {"total_retained_bytes": 1024},
                "incumbent_trajectory": [
                    {"elapsed_nanoseconds": 0, "assembly_index": 9},
                    {"elapsed_nanoseconds": 12, "assembly_index": 7},
                ],
                "parallel": {
                    "aggregate": {
                        "shared_assembly_cache": {
                            "hits": 1,
                            "misses": 9,
                            "growth_count": 3,
                            "growth_nanoseconds": 80,
                            "max_growth_nanoseconds": 50,
                        }
                    }
                },
            },
        }
        rows = experiment.summary_rows(
            "shared", [result], {"case": profile}, {"case": profile}
        )
        self.assertEqual(rows[0]["paired_wall_speedup_median"], 1.75)
        self.assertEqual(rows[0]["paired_clock_speedup_median"], 1.75)
        self.assertEqual(rows[0]["candidate_cache_hit_rate"], 0.1)
        self.assertEqual(rows[0]["candidate_cache_growth_nanoseconds"], 80)
        self.assertEqual(rows[0]["candidate_states_expanded"], 42)
        self.assertEqual(rows[0]["candidate_states_pruned"], 17)
        self.assertEqual(rows[0]["candidate_total_retained_bytes"], 1024)
        self.assertEqual(rows[0]["candidate_final_incumbent_nanoseconds"], 12)
        self.assertEqual(rows[0]["baseline_peak_rss_kib"], 123)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.csv"
            experiment.write_summary(path, rows)
            with path.open(encoding="utf-8") as stream:
                saved = list(csv.DictReader(stream))
            self.assertEqual(saved[0]["candidate_cache_max_growth_nanoseconds"], "50")

    def test_odd_rounds_fail_before_running(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            status = experiment.main(
                [
                    "--baseline-executable",
                    "before",
                    "--executable",
                    "after",
                    "--output-dir",
                    "unused",
                    "--runs",
                    "3",
                    "--dry-run",
                ]
            )
        self.assertEqual(status, 1)
        self.assertIn("balanced AB/BA", stderr.getvalue())

    @unittest.skipUnless(Path("/usr/bin/time").exists(), "requires GNU time")
    def test_end_to_end_profiles_and_reports_validate_results(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "case.mol"
            source.write_text("fixture", encoding="utf-8")
            manifest = root / "cases.tsv"
            manifest.write_text(
                "\t".join(benchmark.MANIFEST_HEADER)
                + "\ncase\tcase.mol\t7\treviewed\tquick,full\tfixture\n",
                encoding="utf-8",
            )
            before = root / "before"
            before.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "Path('caseOut').write_text("
                "'has assembly index: 7\\ntime elapsed: 10\\n')\n",
                encoding="utf-8",
            )
            before.chmod(0o755)
            after = root / "after"
            shutil.copy2(before, after)
            with after.open("a", encoding="utf-8") as stream:
                stream.write("# different executable fingerprint\n")
            output = root / "result"
            arguments = [
                "--baseline-executable",
                str(before),
                "--executable",
                str(after),
                "--manifest",
                str(manifest),
                "--variant",
                "local=local",
                "--runs",
                "2",
                "--warmup",
                "0",
                "--output-dir",
                str(output),
                "--require-all-faster",
            ]
            with contextlib.redirect_stdout(io.StringIO()):
                # Equal clock ticks must fail the gate independently of noisy wall time.
                self.assertEqual(experiment.main(arguments), 1)
            report = json.loads((output / "local.json").read_text(encoding="utf-8"))
            self.assertEqual(len(report["cases"]), 1)
            self.assertEqual(len(report["cases"][0]["candidate"]["measurements"]), 2)
            profiles = json.loads(
                (output / "profiles.json").read_text(encoding="utf-8")
            )
            self.assertTrue(profiles["complete"])
            self.assertGreater(profiles["baseline"]["case"]["peak_rss_kib"], 0)
            self.assertGreater(profiles["variants"]["local"]["case"]["peak_rss_kib"], 0)
            self.assertIsNone(profiles["variants"]["local"]["case"]["telemetry"])
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()) as stderr,
            ):
                self.assertEqual(experiment.main(arguments), 1)
            self.assertIn("already exists", stderr.getvalue())
            broken = before.read_text(encoding="utf-8").replace("index: 7", "index: 8")
            before.write_text(broken, encoding="utf-8")
            arguments[arguments.index(str(output))] = str(root / "wrong-result")
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()) as stderr,
            ):
                self.assertEqual(experiment.main(arguments), 1)
            self.assertIn("expected assembly index 7, got 8", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
