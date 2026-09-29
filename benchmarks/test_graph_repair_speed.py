from __future__ import annotations

import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from benchmarks import graph_repair_speed as speed


class GraphRepairSpeedTests(unittest.TestCase):
    def setUp(self) -> None:
        self.case = {
            "name": "example",
            "input": "example.mol",
            "expected_assembly_index": 3,
            "expectation": "reviewed",
        }

    def sample(self, method: str, *, run: int = 0, seconds: float = 1) -> dict:
        return {
            "method": method,
            "status": "completed" if method == "exact" else "upper_bound",
            "exact_completed": method == "exact",
            "assembly_index": 3 if method == "exact" else 4,
            "trivial_upper_bound": 8,
            "algorithm_seconds": seconds,
            "cpu_seconds": seconds,
            "end_to_end_seconds": seconds + 0.01,
            "phase": "measured",
            "round": run,
        }

    def test_completed_uses_median_of_paired_ratios(self) -> None:
        samples = [
            self.sample("exact", seconds=2),
            self.sample("graph", seconds=1),
            self.sample("exact", run=1, seconds=100),
            self.sample("graph", run=1, seconds=10),
        ]
        result = speed.summarize_case(self.case, samples, 2)
        self.assertTrue(result["completed_comparison"])
        self.assertEqual(result["comparison"]["paired_median_algorithm_speedup"], 6)
        self.assertEqual(result["gap_to_exact"], 1)

    def test_one_measured_timeout_excludes_entire_case(self) -> None:
        timeout = {
            "method": "exact",
            "status": "timeout",
            "phase": "measured",
            "round": 1,
            "process_timeout_seconds": 10,
        }
        samples = [
            self.sample("exact"),
            self.sample("graph"),
            timeout,
            self.sample("graph", run=1),
        ]
        result = speed.summarize_case(self.case, samples, 2)
        self.assertFalse(result["completed_comparison"])
        self.assertIsNone(result["comparison"]["paired_median_algorithm_speedup"])
        self.assertAlmostEqual(
            result["censored_process_speedup_lower_bound"], 10 / 1.01
        )
        self.assertEqual(result["exact_assembly_index"], 3)

    def test_provisional_manifest_mismatch_is_flagged(self) -> None:
        self.case.update(expectation="provisional", expected_assembly_index=2)
        result = speed.summarize_case(
            self.case, [self.sample("exact"), self.sample("graph")], 1
        )
        self.assertFalse(result["exact_matches_manifest"])
        self.case["expectation"] = "reviewed"
        with self.assertRaisesRegex(ValueError, "reviewed"):
            speed.summarize_case(self.case, [self.sample("exact")], 1)

    def test_invalid_bound_is_rejected(self) -> None:
        for bound in (2, 9):
            graph = self.sample("graph")
            graph["assembly_index"] = bound
            with self.subTest(bound=bound), self.assertRaisesRegex(ValueError, "bound"):
                speed.summarize_case(self.case, [self.sample("exact"), graph], 1)

    def test_alternating_order_and_warmup(self) -> None:
        with patch.object(
            speed, "run_sample", side_effect=lambda _p, _c, m, _t: self.sample(m)
        ) as run:
            result = speed.benchmark_case(Path("probe"), self.case, 0, 2, 1, 10)
        self.assertEqual(
            [call.args[2] for call in run.call_args_list],
            ["exact", "graph", "exact", "graph", "graph", "exact"],
        )
        self.assertEqual(len(result["samples"]), 6)
        self.assertEqual(result["methods"]["exact"]["measured_samples"], 2)

    def test_warmup_timeout_skips_later_exact_samples(self) -> None:
        def sample(_probe: Path, _case: dict, method: str, _timeout: float) -> dict:
            if method == "exact":
                return {
                    "method": method,
                    "status": "timeout",
                    "process_timeout_seconds": 10,
                }
            return self.sample(method)

        with patch.object(speed, "run_sample", side_effect=sample) as run:
            result = speed.benchmark_case(Path("probe"), self.case, 0, 2, 1, 10)
        self.assertEqual(
            [call.args[2] for call in run.call_args_list],
            ["exact", "graph", "graph", "graph"],
        )
        self.assertFalse(result["completed_comparison"])
        self.assertIsNone(result["exact_assembly_index"])
        self.assertEqual(result["methods"]["graph"]["completed_samples"], 2)

    def test_enumeration_limit_is_not_completed(self) -> None:
        sample = self.sample("exact")
        sample.update(status="incomplete", exact_completed=False)
        result = speed.summarize_case(self.case, [sample, self.sample("graph")], 1)
        self.assertFalse(result["completed_comparison"])
        self.assertIsNone(result["exact_assembly_index"])

    def test_timeout_preserves_start_metadata(self) -> None:
        error = subprocess.TimeoutExpired(
            "probe", 10, output=b'{"event":"started","method":"exact","bonds":9}\n'
        )
        with patch.object(speed.subprocess, "run", side_effect=error):
            result = speed.run_sample(Path("probe"), self.case, "exact", 10)
        self.assertEqual(result["bonds"], 9)
        self.assertTrue(result["calculation_started"])
        self.assertIsNone(result["algorithm_seconds"])

    def test_manifest_cases_and_selection(self) -> None:
        self.assertEqual(len(speed.load_cases([])), 37)
        self.assertEqual(speed.load_cases(["icosane"])[0]["expected_assembly_index"], 6)
        with self.assertRaisesRegex(ValueError, "unknown cases"):
            speed.load_cases(["missing"])


if __name__ == "__main__":
    unittest.main()
