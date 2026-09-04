from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from benchmarks import benchmark, check_speedups


class CheckSpeedupsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest, self.cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)

    def write_result(
        self,
        directory: Path,
        suite: str,
        *,
        case_speedups: dict[str, float] | None = None,
        wall_speedup: float = 1.02,
        clock_speedup: float = 1.03,
        candidate_hash: str = "candidate",
        baseline_hash: str = "baseline",
        candidate_execution: dict[str, object] | None = None,
        baseline_execution: dict[str, object] | None = None,
        runs: int | None = None,
        comparison_order: str = check_speedups.PAIRED_COMPARISON_ORDER,
    ) -> Path:
        overrides = case_speedups or {}
        effective_runs = check_speedups.MINIMUM_RUNS[suite] if runs is None else runs
        selected_cases = [case for case in self.cases if suite in case.suites]
        case_documents = []
        candidate_clock_ticks = 100_000
        baseline_clock_ticks: list[int] = []

        def measurements(
            wall_seconds: float,
            clock_ticks: int,
            assembly_index: int,
        ) -> list[dict[str, object]]:
            return [
                {
                    "round": round_number,
                    "wall_seconds": wall_seconds,
                    "clock_ticks": clock_ticks,
                    "assembly_index": assembly_index,
                }
                for round_number in range(1, effective_runs + 1)
            ]

        for case in selected_cases:
            assert case.expected_assembly_index is not None
            case_speedup = overrides.get(case.name, clock_speedup)
            case_baseline_ticks = round(candidate_clock_ticks * case_speedup)
            baseline_clock_ticks.append(case_baseline_ticks)

            actual_clock_speedup = case_baseline_ticks / candidate_clock_ticks
            case_documents.append(
                {
                    "name": case.name,
                    "expected_assembly_index": case.expected_assembly_index,
                    "candidate": {
                        "measurements": measurements(
                            1.0,
                            candidate_clock_ticks,
                            case.expected_assembly_index,
                        )
                    },
                    "baseline": {
                        "measurements": measurements(
                            wall_speedup,
                            case_baseline_ticks,
                            case.expected_assembly_index,
                        )
                    },
                    "comparison": {
                        "paired_wall_speedup": {"median": wall_speedup},
                        "paired_clock_speedup": {"median": actual_clock_speedup},
                    },
                }
            )

        aggregate_clock_speedup = sum(baseline_clock_ticks) / (
            len(selected_cases) * candidate_clock_ticks
        )
        document = {
            "schema_version": 2,
            "suite": suite,
            "runs": effective_runs,
            "schedule": {"comparison_order": comparison_order},
            "executables": {
                "candidate": {"sha256": candidate_hash},
                "baseline": {"sha256": baseline_hash},
            },
            "execution": {
                "candidate": candidate_execution or {"launcher": [], "environment": {}},
                "baseline": baseline_execution or {"launcher": [], "environment": {}},
            },
            "corpus": benchmark.benchmark_corpus_metadata(
                self.manifest, selected_cases
            ),
            "cases": case_documents,
            "comparison": {
                "paired_round_wall_speedup": {"median": wall_speedup},
                "paired_round_clock_speedup": {"median": aggregate_clock_speedup},
            },
        }
        path = directory / f"{suite}.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def write_all_results(self, directory: Path) -> list[Path]:
        return [self.write_result(directory, suite) for suite in benchmark.KNOWN_SUITES]

    def test_accepts_complete_faster_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            paths = self.write_all_results(Path(temp_directory))
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 0)
            self.assertEqual(
                stdout.getvalue(),
                "PASS: all case clock medians and suite wall/clock medians "
                "exceed 1.000000\n",
            )

    def test_rejects_case_clock_regression(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            paths[0] = self.write_result(
                directory,
                "quick",
                case_speedups={"icosane": 1.0},
            )
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 1)
            self.assertIn("case icosane clock", stdout.getvalue())

    def test_rejects_suite_round_total_regression(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            paths[3] = self.write_result(
                directory,
                "scaling",
                clock_speedup=0.99,
            )
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 1)
            self.assertIn("suite round-total clock", stdout.getvalue())

    def test_rejects_missing_suite_and_wrong_case_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths[:-1]])
            self.assertEqual(status, 2)
            self.assertIn("missing benchmark suite", stderr.getvalue())

            document = json.loads(paths[0].read_text(encoding="utf-8"))
            document["cases"].pop()
            paths[0].write_text(json.dumps(document), encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])
            self.assertEqual(status, 2)
            self.assertIn("wrong quick case coverage", stderr.getvalue())

    def test_rejects_mixed_binary_fingerprints(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            paths[1] = self.write_result(
                directory,
                "full",
                candidate_hash="different",
            )
            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 2)
            self.assertIn("same candidate and baseline", stderr.getvalue())

    def test_rejects_identical_candidate_and_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = [
                self.write_result(
                    directory,
                    suite,
                    candidate_hash="same",
                    baseline_hash="same",
                )
                for suite in benchmark.KNOWN_SUITES
            ]
            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 2)
            self.assertIn("must differ", stderr.getvalue())

    def test_accepts_same_binary_with_distinct_execution_configs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = [
                self.write_result(
                    directory,
                    suite,
                    candidate_hash="same-binary",
                    baseline_hash="same-binary",
                    candidate_execution={
                        "launcher": [],
                        "arguments": ["--mode=candidate"],
                        "environment": {},
                    },
                    baseline_execution={
                        "launcher": [],
                        "arguments": ["--mode=baseline"],
                        "environment": {},
                    },
                )
                for suite in benchmark.KNOWN_SUITES
            ]

            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 0)
            self.assertIn("PASS", stdout.getvalue())

    def test_legacy_missing_arguments_normalize_to_empty_identity(self) -> None:
        missing = {"execution": {"candidate": {"launcher": [], "environment": {}}}}
        explicit_empty = {
            "execution": {
                "candidate": {
                    "launcher": [],
                    "arguments": [],
                    "environment": {},
                }
            }
        }
        path = Path("report.json")

        self.assertEqual(
            check_speedups.execution_identity(missing, "candidate", path),
            check_speedups.execution_identity(explicit_empty, "candidate", path),
        )

    def test_rejects_mixed_execution_arguments_across_suites(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            paths[1] = self.write_result(
                directory,
                "full",
                candidate_execution={
                    "launcher": [],
                    "arguments": ["--candidate-mode=alternate"],
                    "environment": {},
                },
            )

            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

        self.assertEqual(status, 2)
        self.assertIn("execution configurations", stderr.getvalue())

    def test_rejects_malformed_execution_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            paths = self.write_all_results(Path(temp_directory))
            document = json.loads(paths[0].read_text(encoding="utf-8"))
            document["execution"]["candidate"]["arguments"] = [""]
            paths[0].write_text(json.dumps(document), encoding="utf-8")

            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

        self.assertEqual(status, 2)
        self.assertIn("invalid candidate arguments configuration", stderr.getvalue())

    def test_rejects_mixed_execution_configs_across_suites(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            paths = self.write_all_results(directory)
            paths[1] = self.write_result(
                directory,
                "full",
                candidate_execution={
                    "launcher": [],
                    "environment": {"OMP_NUM_THREADS": "4"},
                },
            )

            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 2)
            self.assertIn("execution configurations", stderr.getvalue())

    def test_rejects_undersampled_or_unbalanced_results(self) -> None:
        scenarios = (
            ({"runs": 2}, "requires at least"),
            ({"runs": 101}, "odd run count"),
            ({"comparison_order": "candidate first"}, "comparison order"),
        )
        for overrides, expected_error in scenarios:
            with (
                self.subTest(overrides=overrides),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                directory = Path(temp_directory)
                paths = self.write_all_results(directory)
                paths[0] = self.write_result(directory, "quick", **overrides)
                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    status = check_speedups.main([str(path) for path in paths])

                self.assertEqual(status, 2)
                self.assertIn(expected_error, stderr.getvalue())

    def test_rejects_stale_corpus_fingerprints(self) -> None:
        scenarios = (
            (("corpus", "manifest", "sha256"), "stale benchmark manifest"),
            (("corpus", "inputs", 0, "sha256"), "stale benchmark input"),
        )
        for keys, expected_error in scenarios:
            with (
                self.subTest(keys=keys),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                paths = self.write_all_results(Path(temp_directory))
                document = json.loads(paths[0].read_text(encoding="utf-8"))
                value = document
                for key in keys[:-1]:
                    value = value[key]
                value[keys[-1]] = "edited"
                paths[0].write_text(json.dumps(document), encoding="utf-8")

                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    status = check_speedups.main([str(path) for path in paths])

                self.assertEqual(status, 2)
                self.assertIn(expected_error, stderr.getvalue())

    def test_rejects_skeletal_or_edited_sample_evidence(self) -> None:
        scenarios = (
            ("sample count", "wrong candidate sample count"),
            ("round order", "round order"),
            ("stored summary", "inconsistent paired clock median"),
            ("raw sample", "inconsistent paired clock median"),
        )
        for scenario, expected_error in scenarios:
            with (
                self.subTest(scenario=scenario),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                paths = self.write_all_results(Path(temp_directory))
                document = json.loads(paths[0].read_text(encoding="utf-8"))
                first_case = document["cases"][0]
                if scenario == "sample count":
                    first_case["candidate"]["measurements"].pop()
                elif scenario == "round order":
                    first_case["candidate"]["measurements"][0]["round"] = 2
                elif scenario == "stored summary":
                    first_case["comparison"]["paired_clock_speedup"]["median"] = 9.0
                else:
                    for sample in first_case["baseline"]["measurements"]:
                        sample["clock_ticks"] += 1
                paths[0].write_text(json.dumps(document), encoding="utf-8")

                with contextlib.redirect_stderr(io.StringIO()) as stderr:
                    status = check_speedups.main([str(path) for path in paths])

                self.assertEqual(status, 2)
                self.assertIn(expected_error, stderr.getvalue())

    def test_rejects_noninteger_schema_version(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            paths = self.write_all_results(Path(temp_directory))
            document = json.loads(paths[0].read_text(encoding="utf-8"))
            document["schema_version"] = 2.0
            paths[0].write_text(json.dumps(document), encoding="utf-8")

            with contextlib.redirect_stderr(io.StringIO()) as stderr:
                status = check_speedups.main([str(path) for path in paths])

        self.assertEqual(status, 2)
        self.assertIn("schema_version 2", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
