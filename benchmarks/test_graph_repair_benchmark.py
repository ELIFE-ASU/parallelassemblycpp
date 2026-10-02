from __future__ import annotations

import contextlib
import copy
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from benchmarks import graph_repair_benchmark as bounds


class GraphRepairBenchmarkTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.source = Path(directory.name) / "example"
        self.source.write_text("fixture", encoding="utf-8")
        self.case = {
            "name": "example",
            "input": str(self.source),
            "expected_assembly_index": 2,
            "expectation": "reviewed",
        }
        self.certificate = bounds.trail_repair(
            ["C", "C", "C", "C"], [[0, 1, 1], [1, 2, 1], [2, 3, 1]]
        )["certificate"]
        self.certificate["elapsed_seconds"] = 0.01

    def run_samples(self, samples: list[dict]) -> list[dict]:
        processes = [
            subprocess.CompletedProcess(
                [], 0, stdout=json.dumps(sample) + "\n", stderr=""
            )
            for sample in samples
        ]
        with patch.object(bounds.subprocess, "run", side_effect=processes):
            return bounds.run_probe(Path("probe"), [self.case], len(samples), 10)

    def test_all_repetitions_are_retained(self) -> None:
        second = copy.deepcopy(self.certificate)
        second["elapsed_seconds"] = 0.03
        rows = self.run_samples([self.certificate, second])
        self.assertEqual(rows[0]["algorithm_seconds"], [0.01, 0.03])
        self.assertEqual(rows[0]["median_algorithm_seconds"], 0.02)
        self.assertTrue(rows[0]["construction_valid"])
        self.assertTrue(rows[0]["trail_construction_valid"])

    def test_later_invalid_certificate_is_rejected(self) -> None:
        second = copy.deepcopy(self.certificate)
        second["remaining_fragments"] += 1
        with self.assertRaisesRegex(ValueError, "invalid construction for example"):
            self.run_samples([self.certificate, second])

    def test_later_probe_error_is_reported(self) -> None:
        with self.assertRaisesRegex(ValueError, "probe failed for example: failed"):
            self.run_samples([self.certificate, {"error": "failed"}])

    def test_invalid_algorithm_timings_are_rejected(self) -> None:
        for value in (float("nan"), float("inf"), -1, True, "1", None):
            second = copy.deepcopy(self.certificate)
            second["elapsed_seconds"] = value
            with (
                self.subTest(value=value),
                self.assertRaisesRegex(ValueError, "invalid algorithm timing"),
            ):
                self.run_samples([self.certificate, second])

    def test_nonfinite_timeout_is_rejected_before_running(self) -> None:
        for timeout in ("nan", "inf", "-inf"):
            with (
                self.subTest(timeout=timeout),
                contextlib.redirect_stderr(io.StringIO()) as stderr,
                self.assertRaises(SystemExit) as error,
            ):
                bounds.main(["--executable=probe", f"--timeout={timeout}"])
            self.assertEqual(error.exception.code, 2)
            self.assertIn("--timeout must be finite and positive", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
