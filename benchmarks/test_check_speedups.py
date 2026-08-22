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
        _, self.cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)

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
        runs: int | None = None,
        comparison_order: str = check_speedups.PAIRED_COMPARISON_ORDER,
    ) -> Path:
        overrides = case_speedups or {}
        effective_runs = (
            check_speedups.MINIMUM_RUNS[suite] if runs is None else runs
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
            "cases": [
                {
                    "name": case.name,
                    "comparison": {
                        "paired_clock_speedup": {
                            "median": overrides.get(case.name, 1.01)
                        }
                    },
                }
                for case in self.cases
                if suite in case.suites
            ],
            "comparison": {
                "paired_round_wall_speedup": {"median": wall_speedup},
                "paired_round_clock_speedup": {"median": clock_speedup},
            },
        }
        path = directory / f"{suite}.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def write_all_results(self, directory: Path) -> list[Path]:
        return [
            self.write_result(directory, suite)
            for suite in benchmark.KNOWN_SUITES
        ]

    def test_accepts_complete_faster_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            paths = self.write_all_results(Path(temp_directory))
            with contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = check_speedups.main([str(path) for path in paths])

            self.assertEqual(status, 0)
            self.assertIn("PASS", stdout.getvalue())

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
            self.assertIn("must be different", stderr.getvalue())

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


if __name__ == "__main__":
    unittest.main()
