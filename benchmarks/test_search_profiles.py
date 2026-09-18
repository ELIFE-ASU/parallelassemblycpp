from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

from benchmarks import search_profiles


class SearchProfileTests(unittest.TestCase):
    def fixture(self) -> dict[str, object]:
        def cache(size: int) -> dict[str, object]:
            return {
                "measurement": search_profiles.CACHE_MEASUREMENT,
                **dict.fromkeys(search_profiles.CACHE_COMPONENTS, size),
                "total_retained_bytes": size * len(search_profiles.CACHE_COMPONENTS),
            }

        def event(time: int, index: int, worker: int = 0) -> dict[str, int]:
            return {
                "elapsed_nanoseconds": time,
                "assembly_index": index,
                "rank": 0,
                "global_worker_index": worker,
            }

        counters = dict.fromkeys(search_profiles.STATE_COUNTERS, 2)
        trajectory = [event(0, 30), event(10, 25), event(20, 23, 1)]
        workers = [
            {
                "rank": 0,
                "global_worker_index": worker,
                "counters": dict.fromkeys(search_profiles.STATE_COUNTERS, 1),
                "local_caches": cache(100),
                "incumbent_trajectory": trajectory[:2]
                if worker == 0
                else trajectory[2:],
            }
            for worker in range(2)
        ]
        return {
            "counters": counters,
            "trajectory_clock": search_profiles.TRAJECTORY_CLOCK,
            "trajectory_index": "internal_before_disjoint_compensation",
            "trajectory_scope": (
                "feasible_incumbent_discovery; MPI rank starts barrier-aligned"
            ),
            "incumbent_trajectory": trajectory,
            "local_caches": cache(200),
            "parallel": {
                "aggregate": {"counters": dict(counters), "local_caches": cache(200)},
                "workers": workers,
            },
        }

    def test_legacy_and_complete_additive_profiles(self) -> None:
        search_profiles.validate_search_profile({"counters": {}})
        search_profiles.validate_search_profile(
            {
                "counters": {},
                "parallel": {
                    "aggregate": {"counters": {}},
                    "workers": [{"counters": {}}],
                },
            }
        )
        search_profiles.validate_search_profile(self.fixture())

    def test_state_counts_reject_partial_invalid_and_inconsistent_records(self) -> None:
        for location in ("top", "aggregate", "worker"):
            for value in (None, -1, True, 1.5):
                with self.subTest(location=location, value=value):
                    fixture = self.fixture()
                    record = (
                        fixture
                        if location == "top"
                        else (
                            fixture["parallel"]["aggregate"]
                            if location == "aggregate"
                            else fixture["parallel"]["workers"][0]
                        )
                    )
                    record["counters"]["states_expanded"] = value
                    with self.assertRaisesRegex(ValueError, "counter group"):
                        search_profiles.validate_search_profile(fixture)
        fixture = self.fixture()
        fixture["parallel"]["workers"][0]["counters"]["states_expanded"] += 1
        with self.assertRaisesRegex(ValueError, "does not match workers"):
            search_profiles.validate_search_profile(fixture)
        fixture = self.fixture()
        fixture["counters"]["states_bound_pruned"] += 1
        with self.assertRaisesRegex(ValueError, "exceed pruned"):
            search_profiles.validate_search_profile(fixture)

    def test_cache_retained_bytes_reject_invalid_components_and_sums(self) -> None:
        fixture = self.fixture()
        for name in (*search_profiles.CACHE_COMPONENTS, "total_retained_bytes"):
            for value in (None, -1, True, 1.0):
                with self.subTest(name=name, value=value):
                    malformed = copy.deepcopy(fixture)
                    malformed["local_caches"][name] = value
                    with self.assertRaisesRegex(ValueError, "retained bytes"):
                        search_profiles.validate_search_profile(malformed)
        malformed = copy.deepcopy(fixture)
        malformed["local_caches"]["total_retained_bytes"] += 1
        with self.assertRaisesRegex(ValueError, "total does not match"):
            search_profiles.validate_search_profile(malformed)
        malformed = copy.deepcopy(fixture)
        malformed["parallel"]["workers"][0]["local_caches"][
            "canonical_mask_retained_bytes"
        ] += 1
        malformed["parallel"]["workers"][0]["local_caches"]["total_retained_bytes"] += 1
        with self.assertRaisesRegex(ValueError, "does not match workers"):
            search_profiles.validate_search_profile(malformed)
        del fixture["parallel"]["workers"][0]["local_caches"]
        with self.assertRaisesRegex(ValueError, "retained bytes"):
            search_profiles.validate_search_profile(fixture)

    def test_trajectories_reject_invalid_order_clock_and_worker_identity(self) -> None:
        for name, value, message in (
            ("elapsed_nanoseconds", -1, "invalid incumbent"),
            ("elapsed_nanoseconds", 21, "time-ordered"),
            ("assembly_index", 30, "strict improvements"),
            ("assembly_index", True, "invalid incumbent"),
            ("global_worker_index", 2, "worker identity"),
            ("rank", 1, "worker identity"),
        ):
            fixture = self.fixture()
            fixture["incumbent_trajectory"][1][name] = value
            with (
                self.subTest(name=name, value=value),
                self.assertRaisesRegex(ValueError, message),
            ):
                search_profiles.validate_search_profile(fixture)
        fixture = self.fixture()
        fixture["trajectory_clock"] = "cpu_clock"
        with self.assertRaisesRegex(ValueError, "trajectory clock"):
            search_profiles.validate_search_profile(fixture)
        fixture = self.fixture()
        del fixture["parallel"]["workers"][1]["incumbent_trajectory"]
        with self.assertRaisesRegex(ValueError, "missing incumbent trajectory"):
            search_profiles.validate_search_profile(fixture)

    def test_placement_summary_preserves_trajectory_execution_and_separate_counts(
        self,
    ) -> None:
        telemetry = self.fixture()
        telemetry["counters"]["matching_visits"] = 25
        telemetry["caches"] = {
            "canonical_mask": {"hits": 99, "misses": 1, "hit_rate": 0.99}
        }
        execution = {"environment": {"OMP_NUM_THREADS": "2", "OMP_PLACES": "{0},{2}"}}
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "paired.json"
            report.write_text(
                json.dumps(
                    {
                        "execution": {"telemetry": execution},
                        "cases": [
                            {
                                "name": "paclitaxel",
                                "candidate": {"telemetry": telemetry},
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            summary = Path(directory) / "search-profiles.json"
            text = search_profiles.write_placement_profiles(
                summary, [("omp-2", report)]
            )
            profile = json.loads(summary.read_text())["placements"][0]
        self.assertEqual(profile["placement"], "omp-2")
        self.assertEqual(profile["execution"], execution)
        self.assertEqual(
            profile["incumbent_trajectory"], telemetry["incumbent_trajectory"]
        )
        self.assertEqual(profile["final_incumbent_nanoseconds"], 20)
        self.assertEqual(profile["final_incumbent_index"], 23)
        self.assertEqual(
            profile["trajectory_index"], "internal_before_disjoint_compensation"
        )
        self.assertEqual(profile["trajectory_scope"], telemetry["trajectory_scope"])
        self.assertEqual(profile["states_expanded"], 2)
        self.assertEqual(profile["states_pruned"], 2)
        self.assertEqual(profile["duplicate_classes_pruned"], 2)
        self.assertEqual(profile["total_retained_bytes"], 1000)
        self.assertEqual(profile["canonical_mask_misses"], 1)
        self.assertIn("omp-2 paclitaxel: states_expanded=2, states_pruned=2", text)


if __name__ == "__main__":
    unittest.main()
