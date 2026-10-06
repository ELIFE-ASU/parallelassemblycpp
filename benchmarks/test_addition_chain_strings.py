"""Correctness and censoring checks for addition-chain timing summaries."""

from __future__ import annotations

import unittest

from benchmarks.addition_chain_strings import fixtures, parse_events, summarize


class AdditionChainStringsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.case = next(case for case in fixtures() if case["name"] == "vector-6-2")
        self.samples = [
            {
                "variant": variant,
                "phase": "measured",
                "round": 0,
                "status": "completed",
                "exact_completed": True,
                "assembly_index": 4,
                "algorithm_seconds": seconds,
                "cpu_seconds": seconds,
                "process_seconds": seconds,
            }
            for variant, seconds in (
                ("baseline", 3.0),
                ("scalar", 2.0),
                ("vector", 1.0),
            )
        ]

    def test_parser_retains_event_fields_and_accepts_timeout_bytes(self) -> None:
        output = (
            b'{"event":"started","method":"exact","length":8}\n'
            b'{"event":"finished","assembly_index":4,"exact_completed":true}\n'
        )
        events = parse_events(output)
        self.assertEqual(events["started"]["length"], 8)
        self.assertEqual(events["finished"]["assembly_index"], 4)
        self.assertTrue(events["finished"]["exact_completed"])
        self.assertEqual(parse_events(None), {})

    def test_parser_rejects_duplicate_or_unknown_events(self) -> None:
        for output in (
            '{"event":"started"}\n{"event":"started"}\n',
            '{"event":"partial"}\n',
        ):
            with self.subTest(output=output), self.assertRaises(ValueError):
                parse_events(output)

    def test_paired_speedups_use_completed_rounds(self) -> None:
        result = summarize(self.case, self.samples, 1)
        self.assertTrue(result["all_variants_completed"])
        self.assertTrue(result["completed_indices_agree"])
        self.assertEqual(result["assembly_index"], 4)
        speedups = result["paired_median_speedups"]
        self.assertEqual(speedups["baseline_over_scalar"]["algorithm_seconds"], 1.5)
        self.assertEqual(speedups["baseline_over_vector"]["algorithm_seconds"], 3.0)
        self.assertEqual(speedups["scalar_over_vector"]["algorithm_seconds"], 2.0)

    def test_timeout_excludes_only_affected_pairs(self) -> None:
        self.samples[0].update(status="timeout", exact_completed=False)
        result = summarize(self.case, self.samples, 1)
        self.assertFalse(result["all_variants_completed"])
        speedups = result["paired_median_speedups"]
        self.assertIsNone(speedups["baseline_over_scalar"]["algorithm_seconds"])
        self.assertIsNone(speedups["baseline_over_vector"]["algorithm_seconds"])
        self.assertEqual(speedups["scalar_over_vector"]["algorithm_seconds"], 2.0)

    def test_warmup_timeout_cannot_be_hidden_by_completed_measurements(self) -> None:
        self.samples.append(
            {
                "variant": "vector",
                "phase": "warmup",
                "round": 0,
                "status": "timeout",
                "exact_completed": False,
            }
        )
        result = summarize(self.case, self.samples, 1)
        self.assertFalse(result["all_variants_completed"])
        self.assertIsNone(
            result["paired_median_speedups"]["baseline_over_vector"][
                "algorithm_seconds"
            ]
        )

    def test_missing_round_prevents_completed_comparison(self) -> None:
        result = summarize(self.case, self.samples, 2)
        self.assertFalse(result["all_variants_completed"])
        self.assertIsNone(
            result["paired_median_speedups"]["scalar_over_vector"]["algorithm_seconds"]
        )

    def test_cross_variant_index_disagreement_is_rejected(self) -> None:
        self.samples[1]["assembly_index"] = 5
        with self.assertRaisesRegex(ValueError, "exact indices differ"):
            summarize(self.case, self.samples, 1)

    def test_equal_but_wrong_indices_are_rejected_against_known_fixture(self) -> None:
        for sample in self.samples:
            sample["assembly_index"] = 5
        with self.assertRaisesRegex(ValueError, "expected 4"):
            summarize(self.case, self.samples, 1)


if __name__ == "__main__":
    unittest.main()
