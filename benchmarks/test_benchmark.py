from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import shlex
import signal
import tempfile
import textwrap
import time
import unittest
from itertools import pairwise
from pathlib import Path
from unittest import mock

from benchmarks import benchmark


def _find_component_root(parents: list[int], atom_index: int) -> int:
    while parents[atom_index] != atom_index:
        parents[atom_index] = parents[parents[atom_index]]
        atom_index = parents[atom_index]
    return atom_index


class BenchmarkTests(unittest.TestCase):
    def create_fixture(self, directory: Path, name: str = "input.mol") -> Path:
        fixture = directory / name
        fixture.write_text("fixture\n", encoding="utf-8")
        return fixture

    def create_manifest(
        self,
        directory: Path,
        rows: list[tuple[str, str, int, str, str, str]],
    ) -> Path:
        manifest = directory / "cases.tsv"
        contents = ["\t".join(benchmark.MANIFEST_HEADER)]
        contents.extend("\t".join(map(str, row)) for row in rows)
        manifest.write_text("\n".join(contents) + "\n", encoding="utf-8")
        return manifest

    def create_fake_executable(
        self,
        directory: Path,
        name: str,
        assembly_index: int,
        clock_ticks: int,
    ) -> Path:
        executable = directory / name
        executable.write_text(
            textwrap.dedent(
                f"""\
                #!/usr/bin/env python3
                import json
                import pathlib
                import sys

                input_name = sys.argv[-1]
                output_name = input_name.removesuffix(".mol") + "Out"
                pathlib.Path(output_name).write_text(
                    f"{{input_name}} has assembly index: {assembly_index}\\n"
                    "time elapsed: {clock_ticks}\\n",
                    encoding="utf-8",
                )
                if "--telemetry=1" in sys.argv:
                    phase = {{
                        "clock_ticks": 1,
                        "activations": 1,
                        "start_rss_kib": 10,
                        "peak_rss_kib": 12,
                        "end_rss_kib": 11,
                        "start_virtual_kib": 20,
                        "end_virtual_kib": 21,
                    }}
                    telemetry = {{
                        "schema_version": 1,
                        "processed_graph": {{
                            "atoms": 65,
                            "edges": 64,
                            "active_mask_words": 1,
                        }},
                        "counters": {{
                            "retained_mask_attempts": 3,
                            "retained_masks": 2,
                            "duplicate_mask_attempts": 1,
                            "rejected_masks": 0,
                            "matching_visits": 4,
                            "canonicalisation_calls": 5,
                            "vf2_calls": 1,
                            "vf2_matches": 1,
                        }},
                        "caches": {{
                            "canonical_mask": {{
                                "hits": 4, "misses": 1, "hit_rate": 0.8
                            }},
                            "canonical_class": {{
                                "insertions": 1, "reuses": 0, "reuse_rate": 0.0
                            }},
                            "residual_decomposition": {{
                                "eligible_for_processed_graph": True,
                                "requests": 5,
                                "eligible_requests": 5,
                                "small_molecule_bypasses": 0,
                                "wide_molecule_bypasses": 0,
                                "small_residual_bypasses": 1,
                                "first_occurrence_bypasses": 2,
                                "runtime_disabled_bypasses": 0,
                                "lookups": 2,
                                "hits": 1,
                                "misses": 1,
                                "admissions": 1,
                                "lookup_hit_rate": 0.5,
                                "request_hit_rate": 0.2,
                            }},
                            "assembly_state": {{
                                "lookups": 0,
                                "hits": 0,
                                "misses": 0,
                                "pruned_hits": 0,
                                "updated_hits": 0,
                                "hit_rate": None,
                            }},
                            "pair_bound": {{
                                "lookups": 0,
                                "hits": 0,
                                "misses": 0,
                                "hit_rate": None,
                            }},
                        }},
                        "memory": {{
                            "method": "linux_proc_vmhwm_reset",
                            "phase_peaks_are_absolute_not_additive": True,
                            "phase_peaks_complete": True,
                            "overall_peak_resident_kib": 12,
                            "process_peak_virtual_kib": 21,
                            "phases": {{
                                name: phase
                                for name in (
                                    "input_setup",
                                    "initial_enumeration",
                                    "dag_conversion",
                                    "assembly_search",
                                    "output",
                                )
                            }},
                        }},
                    }}
                    telemetry_name = (
                        input_name.removesuffix(".mol") + "Telemetry.json"
                    )
                    pathlib.Path(telemetry_name).write_text(
                        json.dumps(telemetry), encoding="utf-8"
                    )
                """
            ),
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return executable

    def add_valid_parallel_telemetry(
        self,
        telemetry: dict[str, object],
    ) -> None:
        counters = telemetry["counters"]
        caches = telemetry["caches"]
        assert isinstance(counters, dict)
        assert isinstance(caches, dict)
        canonical = caches["canonical_mask"]
        canonical_class = caches["canonical_class"]
        residual = caches["residual_decomposition"]
        assembly = caches["assembly_state"]
        pair_bound = caches["pair_bound"]
        assert all(
            isinstance(value, dict)
            for value in (canonical, canonical_class, residual, assembly, pair_bound)
        )
        aggregate_counters = {
            "retained_mask_attempts": counters["retained_mask_attempts"],
            "retained_masks": counters["retained_masks"],
            "duplicate_mask_attempts": counters["duplicate_mask_attempts"],
            "rejected_masks": counters["rejected_masks"],
            "matching_visits": counters["matching_visits"],
            "canonicalisation_calls": counters["canonicalisation_calls"],
            "canonicalisation_mask_cache_hits": canonical["hits"],
            "canonicalisation_mask_cache_misses": canonical["misses"],
            "canonical_class_insertions": canonical_class["insertions"],
            "canonical_class_reuses": canonical_class["reuses"],
            "vf2_calls": counters["vf2_calls"],
            "vf2_matches": counters["vf2_matches"],
            "residual_decomposition_requests": residual["requests"],
            "residual_cache_eligible_requests": residual["eligible_requests"],
            "residual_cache_small_molecule_bypasses": residual[
                "small_molecule_bypasses"
            ],
            "residual_cache_wide_molecule_bypasses": residual["wide_molecule_bypasses"],
            "residual_cache_small_residual_bypasses": residual[
                "small_residual_bypasses"
            ],
            "residual_cache_first_occurrence_bypasses": residual[
                "first_occurrence_bypasses"
            ],
            "residual_cache_runtime_disabled_bypasses": residual[
                "runtime_disabled_bypasses"
            ],
            "residual_cache_lookups": residual["lookups"],
            "residual_cache_hits": residual["hits"],
            "residual_cache_misses": residual["misses"],
            "residual_cache_admissions": residual["admissions"],
            "assembly_cache_lookups": assembly["lookups"],
            "assembly_cache_hits": assembly["hits"],
            "assembly_cache_misses": assembly["misses"],
            "assembly_cache_pruned_hits": assembly["pruned_hits"],
            "assembly_cache_updated_hits": assembly["updated_hits"],
            "pair_bound_cache_lookups": pair_bound["lookups"],
            "pair_bound_cache_hits": pair_bound["hits"],
            "pair_bound_cache_misses": pair_bound["misses"],
        }
        empty_counters = dict.fromkeys(benchmark.PARALLEL_TELEMETRY_COUNTERS, 0)
        graph = telemetry["processed_graph"]
        assert isinstance(graph, dict)
        worker_graph = {
            **graph,
            "residual_cache_eligible": residual["eligible_for_processed_graph"],
        }

        def phases(wall_nanoseconds: int) -> dict[str, dict[str, int]]:
            return {
                name: {"wall_nanoseconds": wall_nanoseconds, "activations": 1}
                for name in benchmark.SEARCH_TELEMETRY_PHASES
            }

        def scheduler_metrics(idle_waits: int) -> dict[str, int]:
            return {
                "depth_two_tasks_spawned": 0,
                "depth_two_tasks_executed": 0,
                "deeper_tasks_spawned": 0,
                "deeper_tasks_executed": 0,
                "task_steal_attempts": 0,
                "task_steals": 0,
                "local_task_executions": 0,
                "scheduler_idle_waits": idle_waits,
                "scheduler_idle_nanoseconds": 50,
                "deep_refill_activations": 0,
                "task_queue_high_watermark": 0,
                "maximum_task_depth_executed": 0,
                "proactive_tail_refills": 0,
                "warm_start_branches": 0,
            }

        workers = [
            {
                "rank": 0,
                "local_worker_index": 0,
                "global_worker_index": 0,
                "shard": {"index": 0, "count": 2},
                "branch_candidates": 2,
                "branch_assignments": 1,
                **scheduler_metrics(1),
                "elapsed_nanoseconds": 100,
                "busy_nanoseconds": 50,
                "processed_graph": worker_graph,
                "phases": phases(10),
                "counters": aggregate_counters,
            },
            {
                "rank": 0,
                "local_worker_index": 1,
                "global_worker_index": 1,
                "shard": {"index": 1, "count": 2},
                "branch_candidates": 2,
                "branch_assignments": 1,
                **scheduler_metrics(2),
                "elapsed_nanoseconds": 80,
                "busy_nanoseconds": 30,
                "processed_graph": worker_graph,
                "phases": phases(6),
                "counters": empty_counters,
            },
        ]
        telemetry["parallel"] = {
            "enabled": True,
            "mode": "openmp",
            "aggregation_scope": "process",
            "elapsed_timing_method": "parallel_region_steady_clock",
            "busy_timing_method": "elapsed_minus_scheduler_idle_time",
            "rank_count": 1,
            "local_threads": 2,
            "local_threads_per_rank": [2],
            "worker_count": 2,
            "shard_ownership": {
                "strategy": "root_branch_ordinal_modulo_worker_count",
                "complete": True,
                "shard_count": 2,
            },
            "branch_scan_complete": True,
            "aggregate": {
                "counters": aggregate_counters,
                "branch_candidates": 2,
                "branch_assignments": 2,
                "depth_two_tasks_spawned": 0,
                "depth_two_tasks_executed": 0,
                "deeper_tasks_spawned": 0,
                "deeper_tasks_executed": 0,
                "task_steal_attempts": 0,
                "task_steals": 0,
                "local_task_executions": 0,
                "scheduler_idle_waits": 3,
                "scheduler_idle_nanoseconds": 100,
                "deep_refill_activations": 0,
                "task_queue_high_watermark": 0,
                "maximum_task_depth_executed": 0,
                "proactive_tail_refills": 0,
                "warm_start_branches": 0,
                "elapsed_nanoseconds": 100,
                "worker_elapsed_nanoseconds": 180,
                "worker_busy_nanoseconds": 80,
                "shared_assembly_cache": {
                    "table_count": 0,
                    "hits": 0,
                    "misses": 0,
                    "collision_chain_steps": 0,
                    "allocated_bytes": 0,
                    "lock_acquisitions": 0,
                    "lock_waits": 0,
                    "lock_wait_nanoseconds": 0,
                },
            },
            "workers": workers,
        }

        memory = telemetry["memory"]
        assert isinstance(memory, dict)
        memory["method"] = "disabled_parallel"
        memory["phase_peaks_complete"] = False
        memory["overall_peak_resident_kib"] = None
        memory["process_peak_virtual_kib"] = None
        memory_phases = memory["phases"]
        assert isinstance(memory_phases, dict)
        for phase in memory_phases.values():
            assert isinstance(phase, dict)
            phase["wall_nanoseconds"] = 1
            for name in (
                "start_rss_kib",
                "peak_rss_kib",
                "end_rss_kib",
                "start_virtual_kib",
                "end_virtual_kib",
            ):
                phase[name] = None

    def add_valid_dynamic_parallel_telemetry(
        self,
        telemetry: dict[str, object],
        mode: str = "openmp",
    ) -> None:
        self.add_valid_parallel_telemetry(telemetry)
        parallel = telemetry["parallel"]
        assert isinstance(parallel, dict)
        del parallel["shard_ownership"]
        parallel["branch_scheduler"] = {
            "strategy": "distributed_global_root_queue",
            "lease_size": 2,
            "root_queue": {
                "participant_count": 1 if mode == "openmp" else 2,
            },
            "adaptive_splitting": dict(benchmark.ADAPTIVE_SPLITTING_POLICY),
            "complete": True,
        }
        aggregate = parallel["aggregate"]
        workers = parallel["workers"]
        assert isinstance(aggregate, dict)
        assert isinstance(workers, list)
        aggregate["branch_leases"] = 2

        if mode == "openmp":
            ranks = (0, 0)
            local_worker_indexes = (0, 1)
            rank_count = 1
        else:
            parallel["mode"] = mode
            parallel["aggregation_scope"] = "all_mpi_ranks"
            parallel["rank_count"] = 2
            parallel["local_threads"] = 1
            parallel["local_threads_per_rank"] = [1, 1]
            ranks = (0, 1)
            local_worker_indexes = (0, 0)
            rank_count = 2

        for worker, rank, local_worker_index in zip(
            workers,
            ranks,
            local_worker_indexes,
            strict=True,
        ):
            assert isinstance(worker, dict)
            del worker["shard"]
            worker["rank"] = rank
            worker["local_worker_index"] = local_worker_index
            worker["root_queue"] = {
                "participant_rank": rank,
                "participant_count": rank_count,
            }
            worker["branch_leases"] = 1

        first_worker = workers[0]
        second_worker = workers[1]
        assert isinstance(first_worker, dict)
        assert isinstance(second_worker, dict)
        first_worker.update(
            {
                "depth_two_tasks_spawned": 2,
                "depth_two_tasks_executed": 1,
                "deeper_tasks_spawned": 1,
                "deeper_tasks_executed": 0,
                "task_steal_attempts": 1,
                "task_steals": 0,
                "local_task_executions": 1,
                "deep_refill_activations": 1,
                "task_queue_high_watermark": 2,
                "maximum_task_depth_executed": 2,
                "proactive_tail_refills": 1,
                "warm_start_branches": 1,
            }
        )
        second_worker.update(
            {
                "depth_two_tasks_spawned": 0,
                "depth_two_tasks_executed": 1,
                "deeper_tasks_spawned": 0,
                "deeper_tasks_executed": 1,
                "task_steal_attempts": 3,
                "task_steals": 2,
                "local_task_executions": 0,
                "task_queue_high_watermark": 0,
                "maximum_task_depth_executed": 3,
                "warm_start_branches": 1 if mode != "openmp" else 0,
            }
        )
        aggregate.update(
            {
                "depth_two_tasks_spawned": 2,
                "depth_two_tasks_executed": 2,
                "deeper_tasks_spawned": 1,
                "deeper_tasks_executed": 1,
                "task_steal_attempts": 4,
                "task_steals": 2,
                "local_task_executions": 1,
                "deep_refill_activations": 1,
                "task_queue_high_watermark": 2,
                "maximum_task_depth_executed": 3,
                "proactive_tail_refills": 1,
                "warm_start_branches": rank_count,
            }
        )

    def create_forwarding_launcher(self, directory: Path) -> Path:
        launcher = directory / "forwarding launcher"
        launcher.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env python3
                import os
                import sys

                if len(sys.argv) < 5 or sys.argv[1] != "--require":
                    raise SystemExit(90)
                if os.environ.get(sys.argv[2]) != sys.argv[3]:
                    raise SystemExit(91)
                os.execv(sys.argv[4], sys.argv[4:])
                """
            ),
            encoding="utf-8",
        )
        launcher.chmod(0o755)
        return launcher

    def test_summary_statistics_support_one_sample_and_outliers(self) -> None:
        single = benchmark.summarize([4.0])
        self.assertEqual(single.median, 4.0)
        self.assertEqual(single.mad, 0.0)
        self.assertEqual(single.p95, 4.0)

        summary = benchmark.summarize([1.0, 2.0, 2.0, 3.0, 100.0])
        self.assertEqual(summary.median, 2.0)
        self.assertEqual(summary.mad, 1.0)
        self.assertGreater(summary.p95, summary.median)

    def test_geometric_mean_and_paired_speedup(self) -> None:
        self.assertAlmostEqual(benchmark.geometric_mean([2.0, 0.5]), 1.0)
        candidate = (
            benchmark.Measurement(2.0, 20, 7),
            benchmark.Measurement(4.0, 40, 7),
        )
        baseline = (
            benchmark.Measurement(4.0, 40, 7),
            benchmark.Measurement(2.0, 20, 7),
        )
        speedup = benchmark.paired_speedup_summary(candidate, baseline, "wall_seconds")
        self.assertIsNotNone(speedup)
        assert speedup is not None
        self.assertAlmostEqual(speedup.median, 1.25)

        zero_clock = (benchmark.Measurement(1.0, 0, 7),)
        self.assertIsNone(
            benchmark.paired_speedup_summary(zero_clock, baseline[:1], "clock_ticks")
        )

    def test_load_manifest_and_select_suite(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            self.create_fixture(directory, "one.mol")
            self.create_fixture(directory, "two")
            manifest = self.create_manifest(
                directory,
                [
                    (
                        "one",
                        "one.mol",
                        7,
                        "reviewed",
                        "quick,full",
                        "tree",
                    ),
                    (
                        "two",
                        "two",
                        8,
                        "provisional",
                        "profile",
                        "native graph",
                    ),
                ],
            )

            resolved, cases = benchmark.load_manifest(manifest)
            self.assertEqual(resolved, manifest.resolve())
            self.assertEqual([case.name for case in cases], ["one", "two"])
            self.assertEqual(
                [case.name for case in benchmark.select_cases(cases, "quick", [])],
                ["one"],
            )
            self.assertEqual(cases[1].expectation, "provisional")

    def test_manifest_rejects_invalid_header_and_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            self.create_fixture(directory)
            invalid = directory / "invalid.tsv"
            invalid.write_text("name\tinput\n", encoding="utf-8")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "invalid header"):
                benchmark.load_manifest(invalid)

            duplicate = self.create_manifest(
                directory,
                [
                    (
                        "same",
                        "input.mol",
                        7,
                        "reviewed",
                        "quick",
                        "one",
                    ),
                    (
                        "same",
                        "input.mol",
                        7,
                        "reviewed",
                        "full",
                        "two",
                    ),
                ],
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "duplicate"):
                benchmark.load_manifest(duplicate)

            empty = self.create_manifest(directory, [])
            with self.assertRaisesRegex(benchmark.BenchmarkError, "manifest is empty"):
                benchmark.load_manifest(empty)

            invalid_integer = self.create_manifest(
                directory,
                [("bad", "input.mol", "seven", "reviewed", "quick", "bad")],
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "assembly index"):
                benchmark.load_manifest(invalid_integer)

            unknown_suite = self.create_manifest(
                directory,
                [("bad", "input.mol", 7, "reviewed", "slow", "bad")],
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "unknown suites"):
                benchmark.load_manifest(unknown_suite)

            missing_fixture = self.create_manifest(
                directory,
                [("bad", "missing.mol", 7, "reviewed", "quick", "bad")],
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "input not found"):
                benchmark.load_manifest(missing_fixture)

    def test_case_rotation_is_deterministic(self) -> None:
        cases = ["a", "b", "c"]
        self.assertEqual(benchmark.rotate_cases(cases, 0), ["a", "b", "c"])
        self.assertEqual(benchmark.rotate_cases(cases, 1), ["b", "c", "a"])
        self.assertEqual(benchmark.rotate_cases(cases, 4), ["b", "c", "a"])

    def test_paired_schedule_rotates_cases_and_balances_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            first_source = self.create_fixture(directory, "first.mol")
            second_source = self.create_fixture(directory, "second.mol")
            cases = [
                benchmark.BenchmarkCase(
                    "first", first_source, 7, "reviewed", ("quick",), "first"
                ),
                benchmark.BenchmarkCase(
                    "second", second_source, 7, "reviewed", ("quick",), "second"
                ),
            ]
            calls: list[tuple[str, str]] = []

            def fake_run_once(
                executable: Path,
                prepared: benchmark.PreparedCase,
                timeout: float,
            ) -> benchmark.Measurement:
                del timeout
                calls.append((executable.name, prepared.case.name))
                return benchmark.Measurement(1.0, 100, 7)

            with (
                mock.patch.object(benchmark, "run_once", side_effect=fake_run_once),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                results = benchmark.run_benchmarks(
                    executable=Path("candidate"),
                    baseline_executable=Path("baseline"),
                    cases=cases,
                    runs=3,
                    warmup=0,
                    timeout=1.0,
                )

            self.assertEqual(
                calls,
                [
                    ("baseline", "first"),
                    ("candidate", "first"),
                    ("baseline", "second"),
                    ("candidate", "second"),
                    ("candidate", "second"),
                    ("baseline", "second"),
                    ("candidate", "first"),
                    ("baseline", "first"),
                    ("baseline", "first"),
                    ("candidate", "first"),
                    ("baseline", "second"),
                    ("candidate", "second"),
                ],
            )
            self.assertTrue(all(len(result.measurements) == 3 for result in results))
            self.assertTrue(
                all(len(result.baseline_measurements) == 3 for result in results)
            )

    def test_telemetry_run_reuses_candidate_launcher_and_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            source = self.create_fixture(directory)
            case = benchmark.BenchmarkCase(
                "sample", source, 7, "reviewed", ("quick",), "telemetry"
            )
            execution = benchmark.ExecutionConfig(
                launcher=("mpiexec", "-n", "2"),
                environment=(("OMP_NUM_THREADS", "3"),),
            )
            telemetry_document = {"schema_version": 1}

            with (
                mock.patch.object(
                    benchmark,
                    "run_once",
                    return_value=benchmark.Measurement(1.0, 100, 7),
                ),
                mock.patch.object(
                    benchmark,
                    "run_telemetry_once",
                    return_value=telemetry_document,
                ) as telemetry_run,
                contextlib.redirect_stdout(io.StringIO()),
            ):
                results = benchmark.run_benchmarks(
                    executable=Path("candidate"),
                    baseline_executable=None,
                    telemetry_executable=Path("candidate-telemetry"),
                    candidate_execution=execution,
                    cases=[case],
                    runs=1,
                    warmup=0,
                    timeout=1.0,
                )

            self.assertEqual(results[0].telemetry, telemetry_document)
            telemetry_run.assert_called_once()
            self.assertEqual(telemetry_run.call_args.kwargs.get("execution"), execution)

    def test_configured_telemetry_command_applies_launcher_and_environment(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            source = self.create_fixture(directory)
            working_directory = directory / "work"
            working_directory.mkdir()
            prepared = benchmark.PreparedCase(
                case=benchmark.BenchmarkCase(
                    "sample", source, 7, "reviewed", ("quick",), "telemetry"
                ),
                input_name=source.name,
                output_path=working_directory / "inputOut",
                telemetry_path=working_directory / "inputTelemetry.json",
                working_directory=working_directory,
            )
            execution = benchmark.ExecutionConfig(
                launcher=("mpiexec", "-n", "2"),
                environment=(("OMP_NUM_THREADS", "3"),),
            )
            telemetry_document = {"schema_version": 1}

            with (
                mock.patch.object(benchmark, "run_command") as run_command,
                mock.patch.object(benchmark, "parse_measurement"),
                mock.patch.object(
                    benchmark,
                    "parse_search_telemetry",
                    return_value=telemetry_document,
                ),
            ):
                result = benchmark.run_telemetry_once(
                    Path("candidate-telemetry"),
                    prepared,
                    1.0,
                    execution=execution,
                )

            self.assertEqual(result, telemetry_document)
            command, command_directory, timeout, environment = (
                run_command.call_args.args
            )
            self.assertEqual(
                command[:4],
                ["mpiexec", "-n", "2", "candidate-telemetry"],
            )
            self.assertEqual(command_directory, working_directory)
            self.assertEqual(timeout, 1.0)
            self.assertEqual(environment["OMP_NUM_THREADS"], "3")

    def test_role_execution_configs_launch_with_environment_and_reach_json(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            fixture = self.create_fixture(directory)
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=7, clock_ticks=200
            )
            telemetry = self.create_fake_executable(
                directory, "telemetry", assembly_index=7, clock_ticks=100
            )
            launcher = self.create_forwarding_launcher(directory)
            report_path = directory / "report.json"
            candidate_launcher = [
                str(launcher),
                "--require",
                "CANDIDATE_MODE",
                "candidate enabled",
            ]
            baseline_launcher = [
                str(launcher),
                "--require",
                "BASELINE_MODE",
                "baseline enabled",
            ]

            with contextlib.redirect_stdout(io.StringIO()):
                status = benchmark.main(
                    [
                        "--input",
                        str(fixture),
                        "--expected",
                        "7",
                        "--executable",
                        str(candidate),
                        "--baseline-executable",
                        str(baseline),
                        "--telemetry",
                        "--telemetry-executable",
                        str(telemetry),
                        "--candidate-launcher",
                        shlex.join(candidate_launcher),
                        "--baseline-launcher",
                        shlex.join(baseline_launcher),
                        "--candidate-env",
                        "CANDIDATE_MODE=candidate enabled",
                        "--baseline-env",
                        "BASELINE_MODE=baseline enabled",
                        "--runs",
                        "2",
                        "--warmup",
                        "0",
                        "--json-output",
                        str(report_path),
                    ]
                )

            self.assertEqual(status, 0)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(
                report["execution"]["candidate"],
                {
                    "launcher": candidate_launcher,
                    "environment": {"CANDIDATE_MODE": "candidate enabled"},
                },
            )
            self.assertEqual(
                report["execution"]["baseline"],
                {
                    "launcher": baseline_launcher,
                    "environment": {"BASELINE_MODE": "baseline enabled"},
                },
            )
            self.assertEqual(
                report["execution"]["telemetry"],
                report["execution"]["candidate"],
            )
            self.assertEqual(len(report["cases"][0]["candidate"]["measurements"]), 2)
            self.assertEqual(len(report["cases"][0]["baseline"]["measurements"]), 2)

    def test_same_executable_is_allowed_when_execution_configs_differ(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            fixture = self.create_fixture(directory)
            executable = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            report_path = directory / "report.json"

            with contextlib.redirect_stdout(io.StringIO()):
                status = benchmark.main(
                    [
                        "--input",
                        str(fixture),
                        "--expected",
                        "7",
                        "--executable",
                        str(executable),
                        "--baseline-executable",
                        str(executable),
                        "--candidate-env",
                        "OMP_NUM_THREADS=4",
                        "--baseline-env",
                        "OMP_NUM_THREADS=1",
                        "--runs",
                        "2",
                        "--warmup",
                        "0",
                        "--json-output",
                        str(report_path),
                    ]
                )

            self.assertEqual(status, 0)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(
                report["executables"]["candidate"]["sha256"],
                report["executables"]["baseline"]["sha256"],
            )
            self.assertEqual(
                report["execution"]["candidate"]["environment"],
                {"OMP_NUM_THREADS": "4"},
            )
            self.assertEqual(
                report["execution"]["baseline"]["environment"],
                {"OMP_NUM_THREADS": "1"},
            )

    def test_identical_binary_copies_require_distinct_execution_configs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=22, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=22, clock_ticks=100
            )
            self.assertNotEqual(candidate, baseline)
            self.assertEqual(candidate.read_bytes(), baseline.read_bytes())

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                status = benchmark.main(
                    [
                        "--executable",
                        str(candidate),
                        "--baseline-executable",
                        str(baseline),
                        "--runs",
                        "2",
                        "--warmup",
                        "0",
                    ]
                )

            self.assertEqual(status, 1)
            self.assertIn("same binary fingerprint", stderr.getvalue())

    def test_execution_config_rejects_invalid_or_ambiguous_settings(self) -> None:
        parser = benchmark.create_argument_parser()
        for arguments in (
            ["--candidate-env", "missing-separator"],
            ["--candidate-env", "1INVALID=value"],
            ["--candidate-launcher", "'unterminated"],
        ):
            with (
                self.subTest(arguments=arguments),
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit),
            ):
                parser.parse_args(arguments)

        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = self.create_fake_executable(
                directory, "candidate", assembly_index=22, clock_ticks=100
            )
            scenarios = (
                (
                    ["--candidate-launcher", str(directory / "missing-launcher")],
                    "candidate launcher not found",
                ),
                (
                    ["--baseline-env", "OMP_NUM_THREADS=1"],
                    "require --baseline-executable",
                ),
                (
                    [
                        "--candidate-env",
                        "OMP_NUM_THREADS=1",
                        "--candidate-env",
                        "OMP_NUM_THREADS=2",
                    ],
                    "duplicate candidate environment setting",
                ),
            )
            for arguments, expected_error in scenarios:
                with self.subTest(arguments=arguments):
                    stderr = io.StringIO()
                    with contextlib.redirect_stderr(stderr):
                        status = benchmark.main(
                            [
                                "--executable",
                                str(executable),
                                "--runs",
                                "1",
                                "--warmup",
                                "0",
                                *arguments,
                            ]
                        )
                    self.assertEqual(status, 1)
                    self.assertIn(expected_error, stderr.getvalue())

    @unittest.skipUnless(os.name == "posix", "POSIX process groups are required")
    def test_configured_timeout_terminates_launcher_process_group(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            fixture = self.create_fixture(directory)
            working_directory = directory / "working"
            working_directory.mkdir()
            marker = directory / "orphan-marker"
            launcher = directory / "timeout-launcher"
            launcher.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env python3
                    import subprocess
                    import sys
                    import time

                    child = (
                        "import pathlib,sys,time; time.sleep(0.4); "
                        "pathlib.Path(sys.argv[1]).write_text('orphan')"
                    )
                    subprocess.Popen([sys.executable, "-c", child, sys.argv[1]])
                    time.sleep(10)
                    """
                ),
                encoding="utf-8",
            )
            launcher.chmod(0o755)
            executable = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            case = benchmark.BenchmarkCase(
                "sample", fixture, 7, "reviewed", ("quick",), "timeout"
            )
            prepared = benchmark.PreparedCase(
                case=case,
                input_name=fixture.name,
                output_path=working_directory / "inputOut",
                telemetry_path=working_directory / "inputTelemetry.json",
                working_directory=working_directory,
            )

            with self.assertRaisesRegex(benchmark.BenchmarkError, "timed out"):
                benchmark.run_once(
                    executable,
                    prepared,
                    0.1,
                    execution=benchmark.ExecutionConfig(
                        launcher=(str(launcher), str(marker))
                    ),
                )
            time.sleep(0.6)
            self.assertFalse(marker.exists())

    @unittest.skipUnless(os.name == "posix", "POSIX process groups are required")
    def test_configured_interrupt_terminates_launcher_process_group(self) -> None:
        process = mock.Mock(pid=12345, returncode=-signal.SIGKILL)
        process.communicate.side_effect = [KeyboardInterrupt, ("", "")]

        with (
            mock.patch.object(benchmark.subprocess, "Popen", return_value=process),
            mock.patch.object(benchmark.os, "killpg") as kill_group,
            self.assertRaises(KeyboardInterrupt),
        ):
            benchmark.run_command(
                ["launcher", "candidate"],
                Path.cwd(),
                30.0,
                None,
            )

        kill_group.assert_called_once_with(12345, signal.SIGKILL)
        self.assertEqual(process.communicate.call_count, 2)

    def test_legacy_default_and_custom_input_resolution(self) -> None:
        parser = benchmark.create_argument_parser()
        manifest, cases = benchmark.resolve_requested_cases(parser.parse_args([]))
        self.assertIsNone(manifest)
        self.assertEqual(cases[0].source, benchmark.DEFAULT_INPUT.resolve())
        self.assertEqual(cases[0].expected_assembly_index, 22)

        with tempfile.TemporaryDirectory() as temp_directory:
            fixture = self.create_fixture(Path(temp_directory))
            manifest, cases = benchmark.resolve_requested_cases(
                parser.parse_args(["--input", str(fixture)])
            )
            self.assertIsNone(manifest)
            self.assertIsNone(cases[0].expected_assembly_index)
            self.assertEqual(cases[0].expectation, "unchecked")

        with self.assertRaisesRegex(benchmark.BenchmarkError, "cannot be combined"):
            benchmark.resolve_requested_cases(
                parser.parse_args(
                    ["--input", str(benchmark.DEFAULT_INPUT), "--suite", "quick"]
                )
            )
        with self.assertRaisesRegex(benchmark.BenchmarkError, "single-input mode"):
            benchmark.resolve_requested_cases(
                parser.parse_args(["--expected", "22", "--case", "ketoconazole"])
            )

    def test_legacy_single_case_summary_labels_remain(self) -> None:
        result = benchmark.CaseResult(
            case=benchmark.BenchmarkCase(
                "sample", Path("sample.mol"), 7, "supplied", (), "custom"
            ),
            measurements=(benchmark.Measurement(1.25, 125, 7),),
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            benchmark.print_summary([result])

        self.assertEqual(
            output.getvalue().splitlines()[:8],
            [
                "",
                "Summary",
                "  wall min:       1.250000 s",
                "  wall median:    1.250000 s",
                "  wall mean:      1.250000 s",
                "  clock min:      125 ticks",
                "  clock median:   125 ticks",
                "  clock mean:     125.0 ticks",
            ],
        )

    def test_unpaired_cli_still_defaults_to_five_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            executable = self.create_fake_executable(
                Path(temp_directory), "candidate", assembly_index=22, clock_ticks=100
            )
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = benchmark.main(
                    ["--executable", str(executable), "--warmup", "0"]
                )

            self.assertEqual(status, 0)
            self.assertIn("Runs: 5 measured, 0 warm-up", output.getvalue())

    def test_reviewed_corpus_values_match_regression_manifest(self) -> None:
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        regression_path = (
            benchmark.REPOSITORY_ROOT / "unitTests" / "regression_cases.tsv"
        )
        with regression_path.open(encoding="utf-8", newline="") as stream:
            rows = csv.DictReader(stream, delimiter="\t")
            expected = {
                row["molecule"]: int(row["expected_assembly_index"]) for row in rows
            }

        for case in cases:
            if case.expectation != "reviewed":
                continue
            fixture_name = case.source.name.removesuffix(".mol")
            self.assertIn(fixture_name, expected)
            self.assertEqual(
                case.expected_assembly_index,
                expected[fixture_name],
                case.name,
            )

        self.assertEqual(
            len(benchmark.select_cases(cases, "quick", [])),
            5,
        )
        self.assertEqual(
            len(benchmark.select_cases(cases, "full", [])),
            15,
        )
        self.assertEqual(
            len(benchmark.select_cases(cases, "profile", [])),
            4,
        )
        scaling = benchmark.select_cases(cases, "scaling", [])
        self.assertEqual(len(scaling), 18)
        self.assertTrue(all(case.expectation == "provisional" for case in scaling))
        self.assertEqual(
            [case.expected_assembly_index for case in scaling],
            [
                10,
                13,
                15,
                17,
                19,
                21,
                23,
                25,
                27,
                29,
                32,
                34,
                8,
                6,
                7,
                10,
                7,
                8,
            ],
        )
        self.assertEqual(
            [case.workload for case in scaling],
            [
                "18 atoms / 16 bonds / 2 comps",
                "27 atoms / 24 bonds / 3 comps",
                "36 atoms / 32 bonds / 4 comps",
                "43 atoms / 38 bonds / 5 comps",
                "53 atoms / 47 bonds / 6 comps",
                "63 atoms / 56 bonds / 7 comps",
                "68 atoms / 60 bonds / 8 comps",
                "77 atoms / 68 bonds / 9 comps",
                "86 atoms / 76 bonds / 10 comps",
                "96 atoms / 85 bonds / 11 comps",
                "105 atoms / 93 bonds / 12 comps",
                "113 atoms / 100 bonds / 13 comps",
                "64 atoms / 63 bonds / cache on / 1 word",
                "65 atoms / 64 bonds / cache on / 1 word",
                "66 atoms / 65 bonds / cache on / 2 words",
                "128 atoms / 127 bonds / cache on / 2 words",
                "129 atoms / 128 bonds / cache on / 2 words",
                "130 atoms / 129 bonds / cache on / 3 words",
            ],
        )

    def test_scaling_corpus_is_a_cumulative_size_series(self) -> None:
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        scaling = [
            case
            for case in benchmark.select_cases(cases, "scaling", [])
            if case.name.startswith("amino-acid-scale-")
        ]
        expected_sizes = [
            (18, 16, 2),
            (27, 24, 3),
            (36, 32, 4),
            (43, 38, 5),
            (53, 47, 6),
            (63, 56, 7),
            (68, 60, 8),
            (77, 68, 9),
            (86, 76, 10),
            (96, 85, 11),
            (105, 93, 12),
            (113, 100, 13),
        ]
        blocks: list[tuple[list[str], list[str]]] = []

        for case, (expected_atoms, expected_bonds, expected_components) in zip(
            scaling, expected_sizes, strict=True
        ):
            lines = case.source.read_text(encoding="utf-8-sig").splitlines()
            atom_count = int(lines[3][0:3])
            bond_count = int(lines[3][3:6])
            atom_lines = lines[4 : 4 + atom_count]
            bond_lines = lines[4 + atom_count : 4 + atom_count + bond_count]
            parents = list(range(atom_count))

            for line in bond_lines:
                first = _find_component_root(parents, int(line[0:3]) - 1)
                second = _find_component_root(parents, int(line[3:6]) - 1)
                parents[first] = second
            component_count = len(
                {
                    _find_component_root(parents, atom_index)
                    for atom_index in range(atom_count)
                }
            )

            self.assertEqual(
                (atom_count, bond_count, component_count),
                (expected_atoms, expected_bonds, expected_components),
                case.name,
            )
            self.assertEqual(bond_count, atom_count - component_count, case.name)
            self.assertNotIn("H", {line[31:34].strip() for line in atom_lines})
            self.assertEqual(lines[4 + atom_count + bond_count], "M  END")
            blocks.append((atom_lines, bond_lines))

        for previous, current in pairwise(blocks):
            previous_atoms, previous_bonds = previous
            current_atoms, current_bonds = current
            self.assertEqual(current_atoms[: len(previous_atoms)], previous_atoms)
            self.assertEqual(current_bonds[: len(previous_bonds)], previous_bonds)

    def test_scaling_corpus_covers_cache_and_word_boundaries(self) -> None:
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        boundary_cases = [
            case
            for case in benchmark.select_cases(cases, "scaling", [])
            if case.name.startswith("mask-boundary-path-")
        ]
        expected_bonds = [63, 64, 65, 127, 128, 129]

        self.assertEqual(
            [case.name for case in boundary_cases],
            [f"mask-boundary-path-{bonds:03d}b" for bonds in expected_bonds],
        )
        for case, bond_count in zip(boundary_cases, expected_bonds, strict=True):
            lines = case.source.read_text(encoding="utf-8").splitlines()
            atom_count = int(lines[1])
            edge_values = [int(value) for value in lines[2].split()]
            edges = list(zip(edge_values[::2], edge_values[1::2], strict=True))

            self.assertEqual(atom_count, bond_count + 1, case.name)
            self.assertEqual(len(edges), bond_count, case.name)
            self.assertEqual(
                edges,
                [(atom, atom + 1) for atom in range(1, atom_count)],
                case.name,
            )
            self.assertEqual(lines[3].split(), ["C"] * atom_count, case.name)
            self.assertEqual(lines[4].split(), ["1"] * bond_count, case.name)

        self.assertEqual(
            [(bonds + 63) // 64 for bonds in expected_bonds],
            [1, 1, 2, 2, 2, 3],
        )
        self.assertEqual(
            [bonds >= 31 for bonds in expected_bonds],
            [True, True, True, True, True, True],
        )

    def test_corpus_run_and_json_report(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            self.create_fixture(directory)
            manifest = self.create_manifest(
                directory,
                [
                    (
                        "sample",
                        "input.mol",
                        7,
                        "reviewed",
                        "quick",
                        "integration",
                    )
                ],
            )
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=7, clock_ticks=200
            )
            telemetry_executable = self.create_fake_executable(
                directory, "telemetry", assembly_index=7, clock_ticks=100
            )
            report_path = directory / "report.json"
            report_path.write_text("stale report\n", encoding="utf-8")

            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                status = benchmark.main(
                    [
                        "--manifest",
                        str(manifest),
                        "--suite",
                        "quick",
                        "--executable",
                        str(candidate),
                        "--baseline-executable",
                        str(baseline),
                        "--warmup",
                        "1",
                        "--json-output",
                        str(report_path),
                        "--telemetry",
                        "--telemetry-executable",
                        str(telemetry_executable),
                    ]
                )

            self.assertEqual(status, 0)
            self.assertIn(
                "Comparison (speedup > 1.0 = candidate faster)",
                output.getvalue(),
            )
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["schema_version"], 2)
            self.assertEqual(
                report["corpus"]["manifest"]["sha256"],
                benchmark.file_sha256(manifest),
            )
            self.assertEqual(
                report["corpus"]["inputs"][0]["sha256"],
                benchmark.file_sha256(directory / "input.mol"),
            )
            self.assertEqual(len(report["cases"]), 1)
            case = report["cases"][0]
            self.assertEqual(len(case["candidate"]["measurements"]), 6)
            self.assertEqual(len(case["baseline"]["measurements"]), 6)
            self.assertEqual(case["comparison"]["paired_clock_speedup"]["median"], 2.0)
            self.assertIn("sha256", report["executables"]["candidate"])
            self.assertIn("sha256", report["executables"]["telemetry"])
            self.assertTrue(report["platform"]["cpu_model"])
            self.assertEqual(report["runs"], 6)
            self.assertEqual(
                report["comparison"]["primary_metric"],
                "paired_round_wall_speedup",
            )
            self.assertTrue(report["telemetry"]["enabled"])
            telemetry = case["candidate"]["telemetry"]
            self.assertEqual(telemetry["counters"]["retained_masks"], 2)
            self.assertEqual(
                telemetry["caches"]["residual_decomposition"]["lookup_hit_rate"],
                0.5,
            )
            malformed = json.loads(json.dumps(telemetry))
            malformed["caches"]["canonical_mask"]["hit_rate"] = 0.999
            malformed_path = directory / "malformed-rate.json"
            malformed_path.write_text(json.dumps(malformed), encoding="utf-8")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "invalid cache rate"):
                benchmark.parse_search_telemetry(malformed_path)

            malformed = json.loads(json.dumps(telemetry))
            malformed["processed_graph"]["edges"] = 65
            malformed["processed_graph"]["active_mask_words"] = 2
            malformed["caches"]["residual_decomposition"][
                "eligible_for_processed_graph"
            ] = False
            malformed_path = directory / "malformed-cache-eligibility.json"
            malformed_path.write_text(json.dumps(malformed), encoding="utf-8")
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent residual cache eligibility",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            wider = json.loads(json.dumps(telemetry))
            wider["processed_graph"] = {
                "atoms": 514,
                "edges": 513,
                "active_mask_words": 9,
            }
            wider_path = directory / "dynamic-mask-width.json"
            wider_path.write_text(json.dumps(wider), encoding="utf-8")
            parsed = benchmark.parse_search_telemetry(wider_path)
            self.assertEqual(parsed["processed_graph"]["active_mask_words"], 9)

            self.assertNotIn("parallel", telemetry)
            parallel_telemetry = json.loads(json.dumps(telemetry))
            self.add_valid_parallel_telemetry(parallel_telemetry)
            parallel_path = directory / "parallel-telemetry.json"
            parallel_path.write_text(
                json.dumps(parallel_telemetry),
                encoding="utf-8",
            )
            parsed = benchmark.parse_search_telemetry(parallel_path)
            self.assertEqual(parsed["parallel"]["worker_count"], 2)
            self.assertEqual(
                parsed["parallel"]["aggregate"]["worker_busy_nanoseconds"],
                80,
            )
            summary_output = io.StringIO()
            with contextlib.redirect_stdout(summary_output):
                benchmark.print_telemetry_summary(
                    [
                        benchmark.CaseResult(
                            case=benchmark.BenchmarkCase(
                                "parallel-sample",
                                directory / "input.mol",
                                7,
                                "reviewed",
                                ("quick",),
                                "parallel telemetry",
                            ),
                            measurements=(benchmark.Measurement(1.0, 1, 7),),
                            telemetry=parsed,
                        )
                    ]
                )
            self.assertIn("Parallel worker telemetry", summary_output.getvalue())
            self.assertIn("openmp:2w", summary_output.getvalue())

            dynamic_telemetry = json.loads(json.dumps(telemetry))
            self.add_valid_dynamic_parallel_telemetry(dynamic_telemetry)
            dynamic_path = directory / "parallel-dynamic-openmp.json"
            dynamic_path.write_text(
                json.dumps(dynamic_telemetry),
                encoding="utf-8",
            )
            parsed = benchmark.parse_search_telemetry(dynamic_path)
            self.assertNotIn("shard_ownership", parsed["parallel"])
            self.assertEqual(
                parsed["parallel"]["branch_scheduler"]["lease_size"],
                2,
            )
            self.assertEqual(
                parsed["parallel"]["branch_scheduler"]["adaptive_splitting"][
                    "maximum_depth"
                ],
                4,
            )
            self.assertEqual(parsed["parallel"]["aggregate"]["branch_leases"], 2)
            self.assertEqual(
                parsed["parallel"]["aggregate"]["deeper_tasks_executed"],
                1,
            )
            self.assertEqual(
                parsed["parallel"]["aggregate"]["task_queue_high_watermark"],
                2,
            )
            self.assertEqual(
                parsed["parallel"]["aggregate"]["maximum_task_depth_executed"],
                3,
            )

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["branch_scheduler"]["adaptive_splitting"][
                "maximum_depth"
            ] = 2
            malformed_path = directory / "parallel-dynamic-task-depth-policy.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent branch scheduler",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][0]["scheduler_idle_waits"] = -1
            malformed_path = directory / "parallel-dynamic-negative-scheduler.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker measurement",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][0]["task_steals"] = 2
            malformed_path = directory / "parallel-dynamic-steal-attempts.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "worker task steals exceed attempts",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][0]["local_task_executions"] = 0
            malformed_path = directory / "parallel-dynamic-task-executions.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "worker task executions do not match transferred tasks",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["aggregate"]["scheduler_idle_waits"] += 1
            malformed_path = directory / "parallel-dynamic-scheduler-sum.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate scheduler field scheduler_idle_waits",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["aggregate"]["task_queue_high_watermark"] += 1
            malformed_path = directory / "parallel-dynamic-queue-maximum.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate task queue high-water mark",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][1][
                "maximum_task_depth_executed"
            ] = 5
            malformed_dynamic["parallel"]["aggregate"][
                "maximum_task_depth_executed"
            ] = 5
            malformed_path = directory / "parallel-dynamic-depth-bound.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker maximum task depth",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][1][
                "maximum_task_depth_executed"
            ] = 0
            malformed_path = directory / "parallel-dynamic-zero-task-depth.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker maximum task depth",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][0]["depth_two_tasks_spawned"] += 1
            malformed_dynamic["parallel"]["aggregate"]["depth_two_tasks_spawned"] += 1
            malformed_path = directory / "parallel-dynamic-depth-two-balance.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "spawned depth-two tasks were not each executed once",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["workers"][0]["busy_nanoseconds"] = 51
            malformed_parallel["parallel"]["aggregate"]["worker_busy_nanoseconds"] = 81
            malformed_path = directory / "parallel-scheduler-busy-time.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent worker busy time",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            mpi_dynamic_telemetry = json.loads(json.dumps(telemetry))
            self.add_valid_dynamic_parallel_telemetry(
                mpi_dynamic_telemetry,
                mode="mpi",
            )
            mpi_dynamic_path = directory / "parallel-dynamic-mpi.json"
            mpi_dynamic_path.write_text(
                json.dumps(mpi_dynamic_telemetry),
                encoding="utf-8",
            )
            parsed = benchmark.parse_search_telemetry(mpi_dynamic_path)
            self.assertEqual(parsed["parallel"]["rank_count"], 2)
            self.assertEqual(
                [
                    worker["root_queue"]["participant_rank"]
                    for worker in parsed["parallel"]["workers"]
                ],
                [0, 1],
            )

            malformed_shared_cache = json.loads(json.dumps(mpi_dynamic_telemetry))
            malformed_shared_cache["parallel"]["aggregate"]["shared_assembly_cache"][
                "lock_acquisitions"
            ] = 1
            malformed_path = directory / "parallel-shared-cache-lookups.json"
            malformed_path.write_text(
                json.dumps(malformed_shared_cache),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent shared assembly-cache lookups",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            legacy_mpi_dynamic = json.loads(json.dumps(mpi_dynamic_telemetry))
            del legacy_mpi_dynamic["parallel"]["aggregate"]["shared_assembly_cache"]
            legacy_scheduler = legacy_mpi_dynamic["parallel"]["branch_scheduler"]
            legacy_scheduler["strategy"] = (
                "dynamic_leases_with_static_mpi_rank_partition"
            )
            del legacy_scheduler["root_queue"]
            legacy_scheduler["rank_partition_count"] = 2
            for worker in legacy_mpi_dynamic["parallel"]["workers"]:
                rank = worker["rank"]
                del worker["root_queue"]
                worker["rank_partition"] = {"index": rank, "count": 2}
            legacy_mpi_path = directory / "parallel-dynamic-legacy-v1-mpi.json"
            legacy_mpi_path.write_text(
                json.dumps(legacy_mpi_dynamic),
                encoding="utf-8",
            )
            parsed_legacy = benchmark.parse_search_telemetry(legacy_mpi_path)
            self.assertEqual(parsed_legacy["schema_version"], 1)
            self.assertEqual(
                parsed_legacy["parallel"]["branch_scheduler"]["strategy"],
                "dynamic_leases_with_static_mpi_rank_partition",
            )
            self.assertEqual(
                [
                    worker["rank_partition"]["index"]
                    for worker in parsed_legacy["parallel"]["workers"]
                ],
                [0, 1],
            )

            malformed_legacy = json.loads(json.dumps(legacy_mpi_dynamic))
            malformed_legacy["parallel"]["workers"][0]["branch_assignments"] = 2
            malformed_legacy["parallel"]["workers"][1]["branch_assignments"] = 0
            malformed_legacy["parallel"]["workers"][1]["branch_leases"] = 0
            malformed_legacy["parallel"]["aggregate"]["branch_leases"] = 1
            malformed_legacy_path = (
                directory / "parallel-dynamic-legacy-rank-redistribution.json"
            )
            malformed_legacy_path.write_text(
                json.dumps(malformed_legacy),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "rank 0 branch assignments do not match partition",
            ):
                benchmark.parse_search_telemetry(malformed_legacy_path)

            redistributed_dynamic = json.loads(json.dumps(mpi_dynamic_telemetry))
            redistributed_dynamic["parallel"]["workers"][0]["branch_assignments"] = 2
            redistributed_dynamic["parallel"]["workers"][1]["branch_assignments"] = 0
            redistributed_dynamic["parallel"]["workers"][1]["branch_leases"] = 0
            redistributed_dynamic["parallel"]["aggregate"]["branch_leases"] = 1
            redistributed_path = directory / "parallel-dynamic-rank-redistribution.json"
            redistributed_path.write_text(
                json.dumps(redistributed_dynamic),
                encoding="utf-8",
            )
            redistributed = benchmark.parse_search_telemetry(redistributed_path)
            self.assertEqual(
                [
                    worker["branch_assignments"]
                    for worker in redistributed["parallel"]["workers"]
                ],
                [2, 0],
            )

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["branch_scheduler"]["lease_size"] = 0
            malformed_path = directory / "parallel-dynamic-lease-size.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent branch scheduler",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(mpi_dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][1]["root_queue"][
                "participant_rank"
            ] = 0
            malformed_path = directory / "parallel-dynamic-root-queue.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker root queue participant",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["aggregate"]["branch_leases"] = 3
            malformed_path = directory / "parallel-dynamic-aggregate-leases.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate branch leases do not match workers",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_dynamic = json.loads(json.dumps(dynamic_telemetry))
            malformed_dynamic["parallel"]["workers"][0]["branch_leases"] = 0
            malformed_path = directory / "parallel-dynamic-worker-leases.json"
            malformed_path.write_text(
                json.dumps(malformed_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent worker leases",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            idle_lease_dynamic = json.loads(json.dumps(dynamic_telemetry))
            idle_lease_dynamic["parallel"]["workers"][1]["branch_assignments"] = 0
            idle_lease_dynamic["parallel"]["aggregate"]["branch_assignments"] = 1
            idle_lease_dynamic["parallel"]["branch_scan_complete"] = False
            idle_lease_path = directory / "parallel-dynamic-idle-lease.json"
            idle_lease_path.write_text(
                json.dumps(idle_lease_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "inconsistent worker leases",
            ):
                benchmark.parse_search_telemetry(idle_lease_path)

            idle_lease_dynamic["parallel"]["workers"][1]["branch_leases"] = 0
            idle_lease_dynamic["parallel"]["aggregate"]["branch_leases"] = 1
            idle_lease_path.write_text(
                json.dumps(idle_lease_dynamic),
                encoding="utf-8",
            )
            parsed = benchmark.parse_search_telemetry(idle_lease_path)
            self.assertEqual(
                parsed["parallel"]["workers"][1]["branch_leases"],
                0,
            )

            incomplete_dynamic = json.loads(json.dumps(idle_lease_dynamic))
            incomplete_dynamic["parallel"]["branch_scan_complete"] = True
            incomplete_path = directory / "parallel-dynamic-incomplete-scan.json"
            incomplete_path.write_text(
                json.dumps(incomplete_dynamic),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "complete branch scan has incomplete assignments",
            ):
                benchmark.parse_search_telemetry(incomplete_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["enabled"] = False
            malformed_path = directory / "parallel-disabled.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "parallel telemetry is not enabled",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["workers"][1]["global_worker_index"] = 0
            malformed_path = directory / "parallel-duplicate-worker.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker identity",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["workers"][0]["busy_nanoseconds"] = 101
            malformed_path = directory / "parallel-worker-busy.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "worker busy time exceeds elapsed time",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["aggregate"]["counters"][
                "matching_visits"
            ] += 1
            malformed_path = directory / "parallel-counter-sum.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate counter matching_visits does not match workers",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["workers"][1]["branch_candidates"] = 3
            malformed_path = directory / "parallel-branch-disagreement.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate branch count does not match workers",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            del malformed_parallel["parallel"]["workers"][0]["phases"]["output"]
            malformed_path = directory / "parallel-worker-phases.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "invalid worker phases",
            ):
                benchmark.parse_search_telemetry(malformed_path)

            malformed_parallel = json.loads(json.dumps(parallel_telemetry))
            malformed_parallel["parallel"]["aggregate"]["counters"][
                "matching_visits"
            ] += 1
            malformed_parallel["parallel"]["workers"][0]["counters"][
                "matching_visits"
            ] += 1
            malformed_path = directory / "parallel-legacy-mismatch.json"
            malformed_path.write_text(
                json.dumps(malformed_parallel),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "aggregate counters do not match legacy telemetry",
            ):
                benchmark.parse_search_telemetry(malformed_path)

    def test_json_output_rejects_protected_paths_before_measurement(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            fixture = self.create_fixture(directory)
            manifest = self.create_manifest(
                directory,
                [
                    (
                        "sample",
                        fixture.name,
                        7,
                        "reviewed",
                        "quick",
                        "integration",
                    )
                ],
            )
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=7, clock_ticks=200
            )
            telemetry = self.create_fake_executable(
                directory, "telemetry", assembly_index=7, clock_ticks=100
            )
            baseline_alias = directory / "baseline-report.json"
            baseline_alias.hardlink_to(baseline)
            input_alias = directory / "input-report.json"
            input_alias.hardlink_to(fixture)

            scenarios = (
                (candidate, "candidate executable"),
                (baseline_alias, "baseline executable"),
                (telemetry, "telemetry executable"),
                (manifest, "benchmark manifest"),
                (input_alias, "benchmark input"),
            )
            protected_contents = {
                path: path.read_bytes()
                for path in (candidate, baseline, telemetry, manifest, fixture)
            }

            with mock.patch.object(benchmark, "run_benchmarks") as run_benchmarks:
                for report_path, expected_error in scenarios:
                    with self.subTest(report_path=report_path):
                        stderr = io.StringIO()
                        with contextlib.redirect_stderr(stderr):
                            status = benchmark.main(
                                [
                                    "--manifest",
                                    str(manifest),
                                    "--suite",
                                    "quick",
                                    "--executable",
                                    str(candidate),
                                    "--baseline-executable",
                                    str(baseline),
                                    "--telemetry",
                                    "--telemetry-executable",
                                    str(telemetry),
                                    "--runs",
                                    "2",
                                    "--warmup",
                                    "0",
                                    "--json-output",
                                    str(report_path),
                                ]
                            )

                        self.assertEqual(status, 1)
                        self.assertIn(expected_error, stderr.getvalue())

                run_benchmarks.assert_not_called()

            for path, contents in protected_contents.items():
                self.assertEqual(path.read_bytes(), contents)

    def test_search_telemetry_parser_rejects_malformed_data(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            path = Path(temp_directory) / "telemetry.json"
            path.write_text("not json", encoding="utf-8")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "could not parse"):
                benchmark.parse_search_telemetry(path)

            path.write_text(
                json.dumps({"schema_version": 1, "counters": {}}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "missing search"):
                benchmark.parse_search_telemetry(path)

    def test_unchecked_ab_run_rejects_index_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            fixture = self.create_fixture(directory)
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=7, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=8, clock_ticks=100
            )

            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                status = benchmark.main(
                    [
                        "--input",
                        str(fixture),
                        "--executable",
                        str(candidate),
                        "--baseline-executable",
                        str(baseline),
                        "--runs",
                        "1",
                        "--warmup",
                        "0",
                    ]
                )

            self.assertEqual(status, 1)
            self.assertIn("odd number of paired rounds", stderr.getvalue())
            self.assertIn("does not match baseline", stderr.getvalue())

    def test_build_rejects_baseline_alias_before_compiling(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = self.create_fake_executable(
                directory, "AssemblyCpp", assembly_index=22, clock_ticks=100
            )

            stderr = io.StringIO()
            with (
                mock.patch.object(benchmark, "build_executable") as build,
                contextlib.redirect_stderr(stderr),
            ):
                status = benchmark.main(
                    [
                        "--executable",
                        str(executable),
                        "--baseline-executable",
                        str(executable),
                        "--build",
                    ]
                )

            self.assertEqual(status, 1)
            build.assert_not_called()
            self.assertIn("same file", stderr.getvalue())

    def test_build_rejects_telemetry_hardlink_before_compiling(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            candidate = self.create_fake_executable(
                directory, "candidate", assembly_index=22, clock_ticks=100
            )
            baseline = self.create_fake_executable(
                directory, "baseline", assembly_index=22, clock_ticks=100
            )
            telemetry = directory / "telemetry"
            telemetry.hardlink_to(baseline)

            stderr = io.StringIO()
            with (
                mock.patch.object(benchmark, "build_executable") as build,
                contextlib.redirect_stderr(stderr),
            ):
                status = benchmark.main(
                    [
                        "--executable",
                        str(candidate),
                        "--baseline-executable",
                        str(baseline),
                        "--telemetry",
                        "--telemetry-executable",
                        str(telemetry),
                        "--build",
                    ]
                )

            self.assertEqual(status, 1)
            build.assert_not_called()
            self.assertIn("same file", stderr.getvalue())

    def test_executable_change_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            executable = self.create_fake_executable(
                Path(temp_directory), "candidate", assembly_index=7, clock_ticks=100
            )
            metadata = benchmark.executable_metadata(executable)
            executable.write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "changed"):
                benchmark.verify_executable_unchanged(executable, metadata, "candidate")


if __name__ == "__main__":
    unittest.main()
