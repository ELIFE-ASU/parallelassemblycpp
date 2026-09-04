from __future__ import annotations

import contextlib
import io
import json
import statistics
import tempfile
import unittest
from pathlib import Path

from benchmarks import check_parallel_scaling


class CheckParallelScalingTests(unittest.TestCase):
    def write_result(
        self,
        path: Path,
        *,
        suite: str = "quick",
        baseline_hash: str = "serial-baseline",
        candidate_hash: str = "parallel-candidate",
        manifest_hash: str = "manifest",
        walls: dict[str, tuple[list[float], list[float]]] | None = None,
        candidate_execution: dict[str, object] | None = None,
        baseline_execution: dict[str, object] | None = None,
    ) -> Path:
        case_walls = walls or {
            "alpha": ([1.0, 2.0], [2.0, 6.0]),
            "beta": ([2.0, 2.0], [4.0, 4.0]),
        }
        runs = len(next(iter(case_walls.values()))[0])
        cases = []
        round_baseline_totals = [0.0] * runs
        round_candidate_totals = [0.0] * runs

        for index, (name, (candidate, baseline)) in enumerate(case_walls.items()):
            self.assertEqual(len(candidate), runs)
            self.assertEqual(len(baseline), runs)
            assembly_index = index + 1

            def samples(
                values: list[float],
                result_index: int = assembly_index,
            ) -> list[dict[str, object]]:
                return [
                    {
                        "round": round_number,
                        "wall_seconds": value,
                        "assembly_index": result_index,
                    }
                    for round_number, value in enumerate(values, start=1)
                ]

            speedup = statistics.median(
                baseline_value / candidate_value
                for baseline_value, candidate_value in zip(
                    baseline, candidate, strict=True
                )
            )
            cases.append(
                {
                    "name": name,
                    "expected_assembly_index": assembly_index,
                    "candidate": {"measurements": samples(candidate)},
                    "baseline": {"measurements": samples(baseline)},
                    "comparison": {
                        "paired_wall_speedup": {"median": speedup},
                        # Deliberately nonsensical: clock data is not a parallel gate.
                        "paired_clock_speedup": {"median": 0.01},
                    },
                }
            )
            for round_index in range(runs):
                round_candidate_totals[round_index] += candidate[round_index]
                round_baseline_totals[round_index] += baseline[round_index]

        suite_speedup = statistics.median(
            baseline / candidate
            for baseline, candidate in zip(
                round_baseline_totals, round_candidate_totals, strict=True
            )
        )
        document = {
            "schema_version": 2,
            "suite": suite,
            "runs": runs,
            "schedule": {
                "comparison_order": check_parallel_scaling.PAIRED_COMPARISON_ORDER
            },
            "executables": {
                "baseline": {"sha256": baseline_hash},
                "candidate": {"sha256": candidate_hash},
            },
            "execution": {
                "candidate": (
                    {
                        "launcher": [],
                        "environment": {"OMP_NUM_THREADS": "4"},
                    }
                    if candidate_execution is None
                    else candidate_execution
                ),
                "baseline": (
                    {
                        "launcher": [],
                        "environment": {"OMP_NUM_THREADS": "1"},
                    }
                    if baseline_execution is None
                    else baseline_execution
                ),
            },
            "corpus": {
                "manifest": {"sha256": manifest_hash},
                "inputs": [
                    {"name": name, "sha256": f"input-{name}"} for name in case_walls
                ],
            },
            "cases": cases,
            "comparison": {
                "paired_round_wall_speedup": {"median": suite_speedup},
                "paired_round_clock_speedup": {"median": 0.01},
            },
        }
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def run_main(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = check_parallel_scaling.main(arguments)
        return status, stdout.getvalue(), stderr.getvalue()

    def test_reports_raw_speedups_and_efficiency_without_clock_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(Path(temp_directory) / "omp.json")
            status, stdout, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")
        self.assertIn("suite round-total", stdout)
        self.assertIn("2.250x", stdout)
        self.assertIn("56.2%", stdout)
        self.assertIn("case alpha", stdout)
        self.assertIn("2.500x", stdout)
        self.assertIn("62.5%", stdout)
        self.assertIn("No case wall-time regressions.", stdout)

    def test_accepts_recorded_forced_parallel_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "omp.json",
                candidate_execution={
                    "launcher": [],
                    "arguments": ["--parallel=on"],
                    "environment": {"OMP_NUM_THREADS": "4"},
                },
                baseline_execution={
                    "launcher": [],
                    "arguments": ["--parallel=off"],
                    "environment": {"OMP_NUM_THREADS": "1"},
                },
            )
            status, _, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")

    def test_legacy_missing_arguments_normalize_to_empty_identity(self) -> None:
        missing = {
            "execution": {
                "candidate": {
                    "launcher": [],
                    "environment": {"OMP_NUM_THREADS": "4"},
                }
            }
        }
        explicit_empty = {
            "execution": {
                "candidate": {
                    "launcher": [],
                    "arguments": [],
                    "environment": {"OMP_NUM_THREADS": "4"},
                }
            }
        }
        path = Path("report.json")

        missing_identity = check_parallel_scaling.execution_identity(
            missing, "candidate", path
        )
        explicit_identity = check_parallel_scaling.execution_identity(
            explicit_empty, "candidate", path
        )

        self.assertEqual(missing_identity, explicit_identity)
        self.assertEqual(missing_identity.arguments, ())

    def test_rejects_recorded_nonparallel_candidate_or_parallel_baseline(
        self,
    ) -> None:
        scenarios = (
            (["--parallel=auto"], ["--parallel=off"], "candidate arguments"),
            (
                ["--parallel=on", "--parallel=off"],
                ["--parallel=off"],
                "candidate arguments",
            ),
            (["--parallel=on"], ["--parallel=on"], "baseline arguments"),
        )
        for candidate_arguments, baseline_arguments, expected_error in scenarios:
            with (
                self.subTest(
                    candidate=candidate_arguments,
                    baseline=baseline_arguments,
                ),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                path = self.write_result(
                    Path(temp_directory) / "omp.json",
                    candidate_execution={
                        "launcher": [],
                        "arguments": candidate_arguments,
                        "environment": {"OMP_NUM_THREADS": "4"},
                    },
                    baseline_execution={
                        "launcher": [],
                        "arguments": baseline_arguments,
                        "environment": {"OMP_NUM_THREADS": "1"},
                    },
                )
                status, _, stderr = self.run_main([f"omp:4:{path}"])

            self.assertEqual(status, 2)
            self.assertIn(expected_error, stderr)

    def test_rejects_malformed_recorded_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "omp.json",
                candidate_execution={
                    "launcher": [],
                    "arguments": "--parallel=on",
                    "environment": {"OMP_NUM_THREADS": "4"},
                },
            )
            status, _, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 2)
        self.assertIn("invalid candidate arguments configuration", stderr)

    def test_report_only_lists_regressions_without_failing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "omp.json",
                walls={
                    "alpha": ([2.0, 2.0], [1.8, 1.8]),
                    "beta": ([1.0, 1.0], [2.0, 2.0]),
                },
            )
            status, stdout, _ = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 0)
        self.assertIn("omp/4 case alpha: 0.900x", stdout)
        self.assertIn(
            "REPORT ONLY: regressions do not change exit status. Use "
            "--require-all-faster to enforce the gate.",
            stdout,
        )

    def test_require_all_faster_fails_with_listed_regressions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "mpi.json",
                walls={"slow": ([1.0, 1.0], [1.0, 1.0])},
                candidate_execution={
                    "launcher": ["mpirun", "-n", "8"],
                    "environment": {},
                },
            )
            status, stdout, _ = self.run_main(["--require-all-faster", f"mpi:8:{path}"])

        self.assertEqual(status, 1)
        self.assertIn("mpi/8 case slow: 1.000x", stdout)
        self.assertIn("FAIL: 1 case wall regression.", stdout)

    def test_accepts_compatible_topologies_with_distinct_worker_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            omp_one = self.write_result(
                directory / "omp-1.json",
                candidate_execution={
                    "launcher": [],
                    "environment": {"OMP_NUM_THREADS": "1"},
                },
            )
            omp_four = self.write_result(directory / "omp-4.json")
            mpi_four = self.write_result(
                directory / "mpi-4.json",
                candidate_execution={
                    "launcher": ["mpirun", "-n", "4"],
                    "environment": {},
                },
            )
            status, stdout, stderr = self.run_main(
                [
                    f"omp:4:{omp_four}",
                    f"mpi:4:{mpi_four}",
                    f"omp:1:{omp_one}",
                ]
            )

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")
        omp_suite_rows = [
            line
            for line in stdout.splitlines()
            if line.startswith("omp") and "suite round-total" in line
        ]
        self.assertEqual([line.split()[1] for line in omp_suite_rows], ["1", "4"])

    def test_rejects_incompatible_suite_corpus_or_baseline(self) -> None:
        scenarios = (
            ({"suite": "full"}, "incompatible suite"),
            ({"manifest_hash": "different"}, "incompatible corpus"),
            ({"baseline_hash": "different"}, "incompatible baseline"),
            (
                {
                    "baseline_execution": {
                        "launcher": ["taskset", "-c", "2"],
                        "environment": {"OMP_NUM_THREADS": "1"},
                    }
                },
                "incompatible baseline execution",
            ),
        )
        for overrides, expected_error in scenarios:
            with (
                self.subTest(overrides=overrides),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                directory = Path(temp_directory)
                first = self.write_result(directory / "first.json")
                second = self.write_result(directory / "second.json", **overrides)
                status, _, stderr = self.run_main([f"omp:4:{first}", f"mpi:4:{second}"])

            self.assertEqual(status, 2)
            self.assertIn(expected_error, stderr)

    def test_rejects_different_candidate_builds_within_one_topology(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            first = self.write_result(
                directory / "omp-4.json",
                candidate_hash="omp-build-a",
            )
            second = self.write_result(
                directory / "omp-8.json",
                candidate_hash="omp-build-b",
                candidate_execution={
                    "launcher": [],
                    "environment": {"OMP_NUM_THREADS": "8"},
                },
            )
            status, _, stderr = self.run_main([f"omp:4:{first}", f"omp:8:{second}"])

        self.assertEqual(status, 2)
        self.assertIn("incompatible candidate fingerprint", stderr)

    def test_rejects_malformed_or_duplicate_worker_specs(self) -> None:
        malformed = (
            ("missing-fields", "expected LABEL:WORKERS:PATH"),
            ("omp:0:report.json", "workers must be a positive integer"),
            ("omp:many:report.json", "workers must be a positive integer"),
        )
        for value, expected_error in malformed:
            with self.subTest(value=value):
                status, _, stderr = self.run_main([value])
            self.assertEqual(status, 2)
            self.assertIn(expected_error, stderr)

        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            first = self.write_result(directory / "first.json")
            second = self.write_result(directory / "second.json")
            status, _, stderr = self.run_main([f"OMP:4:{first}", f"omp:4:{second}"])
        self.assertEqual(status, 2)
        self.assertIn("duplicate topology worker count", stderr)

    def test_rejects_edited_summaries_or_raw_samples(self) -> None:
        scenarios = (
            ("summary", "inconsistent paired wall median"),
            ("wall", "candidate wall sample"),
            ("round", "candidate round order"),
            ("round type", "candidate round order"),
            ("assembly", "wrong candidate assembly index"),
        )
        for scenario, expected_error in scenarios:
            with (
                self.subTest(scenario=scenario),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                path = self.write_result(Path(temp_directory) / "result.json")
                document = json.loads(path.read_text(encoding="utf-8"))
                first_case = document["cases"][0]
                if scenario == "summary":
                    first_case["comparison"]["paired_wall_speedup"]["median"] = 99
                elif scenario == "wall":
                    first_case["candidate"]["measurements"][0]["wall_seconds"] = 0
                elif scenario == "round":
                    first_case["candidate"]["measurements"][0]["round"] = 2
                elif scenario == "round type":
                    first_case["candidate"]["measurements"][0]["round"] = 1.0
                else:
                    first_case["candidate"]["measurements"][0]["assembly_index"] = 99
                path.write_text(json.dumps(document), encoding="utf-8")
                status, _, stderr = self.run_main([f"hybrid:4:{path}"])

            self.assertEqual(status, 2)
            self.assertIn(expected_error, stderr)

    def test_rejects_noninteger_schema_version(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(Path(temp_directory) / "result.json")
            document = json.loads(path.read_text(encoding="utf-8"))
            document["schema_version"] = 2.0
            path.write_text(json.dumps(document), encoding="utf-8")
            status, _, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 2)
        self.assertIn("schema_version 2", stderr)

    def test_derives_omp_mpi_and_hybrid_worker_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            omp = self.write_result(directory / "omp.json")
            mpi = self.write_result(
                directory / "mpi.json",
                candidate_execution={
                    "launcher": ["mpirun", "--np=4"],
                    "environment": {},
                },
            )
            hybrid = self.write_result(
                directory / "hybrid.json",
                candidate_execution={
                    "launcher": ["mpirun", "-n2"],
                    "environment": {"OMP_NUM_THREADS": "4"},
                },
            )
            status, stdout, stderr = self.run_main(
                [f"omp:4:{omp}", f"mpi:4:{mpi}", f"hybrid:8:{hybrid}"]
            )

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")
        self.assertIn("hybrid", stdout)

    def test_omp_thread_limit_caps_effective_workers(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "limited.json",
                candidate_execution={
                    "launcher": ["mpirun", "-n", "2"],
                    "environment": {
                        "OMP_NUM_THREADS": "8",
                        "OMP_THREAD_LIMIT": "4",
                    },
                },
            )
            status, _, stderr = self.run_main([f"hybrid:8:{path}"])

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")

        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "invalid-limit.json",
                candidate_execution={
                    "launcher": [],
                    "environment": {
                        "OMP_NUM_THREADS": "4",
                        "OMP_THREAD_LIMIT": "0",
                    },
                },
            )
            status, _, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 2)
        self.assertIn("OMP_THREAD_LIMIT", stderr)

    def test_rejects_worker_mismatch_or_nonserial_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(Path(temp_directory) / "mismatch.json")
            status, _, stderr = self.run_main([f"omp:8:{path}"])
        self.assertEqual(status, 2)
        self.assertIn("worker count mismatch", stderr)

        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "baseline.json",
                baseline_execution={
                    "launcher": [],
                    "environment": {"OMP_NUM_THREADS": "2"},
                },
            )
            status, _, stderr = self.run_main([f"omp:4:{path}"])
        self.assertEqual(status, 2)
        self.assertIn("baseline", stderr)
        self.assertIn("not serial", stderr)

    def test_allows_identical_binary_with_distinct_execution_configs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "same.json",
                baseline_hash="same",
                candidate_hash="same",
            )
            status, _, stderr = self.run_main([f"omp:4:{path}"])

        self.assertEqual(status, 0)
        self.assertEqual(stderr, "")

    def test_rejects_identical_binary_and_execution_configs(self) -> None:
        serial_execution = {
            "launcher": [],
            "environment": {"OMP_NUM_THREADS": "1"},
        }
        with tempfile.TemporaryDirectory() as temp_directory:
            path = self.write_result(
                Path(temp_directory) / "same.json",
                baseline_hash="same",
                candidate_hash="same",
                baseline_execution=serial_execution,
                candidate_execution=serial_execution,
            )
            status, _, stderr = self.run_main([f"omp:1:{path}"])

        self.assertEqual(status, 2)
        self.assertIn("binary and execution configurations must differ", stderr)


if __name__ == "__main__":
    unittest.main()
