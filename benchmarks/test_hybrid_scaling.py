from __future__ import annotations

import contextlib
import io
import json
import os
import shlex
import sys
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from typing import TYPE_CHECKING
from unittest import mock

from benchmarks import (
    benchmark,
    check_parallel_scaling,
    cpu_topology,
    hybrid_scaling,
    paclitaxel_scaling,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence


class HybridScalingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.topology = cpu_topology.CpuTopology(
            logical_cpus=16,
            physical_cores=8,
            cpu_order=(*range(8, 16), *range(24, 32)),
            affinity_cpus=(*range(8, 16), *range(24, 32)),
            model="Fixture CPU",
            slurm_cpus_per_task=None,
            source="fixture",
        )
        detection_patch = mock.patch.object(
            cpu_topology, "detect_cpu_topology", return_value=self.topology
        )
        self.detect = detection_patch.start()
        self.addCleanup(detection_patch.stop)
        real_which = hybrid_scaling.shutil.which

        def find_executable(command: str, path: str | None = None) -> str | None:
            if command in {"taskset", "mpirun"}:
                return sys.executable
            return real_which(command, path=path)

        tools_patch = mock.patch.object(
            hybrid_scaling.shutil, "which", side_effect=find_executable
        )
        tools_patch.start()
        self.addCleanup(tools_patch.stop)

    def create_build(self, directory: Path) -> Path:
        build = directory / "parallel build"
        build.mkdir()
        for name in ("ParallelAssemblyCpp", "ParallelAssemblyCppHybrid"):
            executable = paclitaxel_scaling.executable_path(build, name)
            executable.write_text(f"fixture for {name}\n", encoding="utf-8")
            executable.chmod(0o755)
        return build

    def run_main(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = hybrid_scaling.main(arguments)
        return status, stdout.getvalue(), stderr.getvalue()

    def benchmark_commands(self, output: str) -> list[list[str]]:
        commands = []
        for line in output.splitlines():
            if line.startswith("#"):
                continue
            command = shlex.split(line)
            for index, argument in enumerate(command):
                if Path(argument).name == "benchmark.py":
                    commands.append(command[index + 1 :])
                    break
        return commands

    def measurements(
        self,
        executable: Path,
        baseline_executable: Path | None,
        cases: Sequence[benchmark.BenchmarkCase],
        runs: int,
        warmup: int,
        timeout: float,
        telemetry_executable: Path | None = None,
        candidate_execution: benchmark.ExecutionConfig | None = None,
        baseline_execution: benchmark.ExecutionConfig | None = None,
    ) -> list[benchmark.CaseResult]:
        self.assertEqual(executable.stem, "ParallelAssemblyCppHybrid")
        self.assertEqual(baseline_executable.stem, "ParallelAssemblyCpp")
        self.assertEqual([case.name for case in cases], ["paclitaxel"])
        self.assertEqual((runs, warmup, timeout), (2, 0, 12.5))
        self.assertIsNone(telemetry_executable)
        self.assertEqual(candidate_execution.arguments, ("--parallel=on",))
        self.assertEqual(baseline_execution.arguments, ("--parallel=off",))
        threads = int(dict(candidate_execution.environment)["OMP_NUM_THREADS"])
        ranks = int(candidate_execution.launcher[-1])
        self.assertEqual(
            dict(candidate_execution.environment),
            {
                "OMP_NUM_THREADS": str(threads),
                "OMP_THREAD_LIMIT": str(threads),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            },
        )
        self.assertEqual(
            dict(baseline_execution.environment),
            {
                "OMP_NUM_THREADS": "1",
                "OMP_THREAD_LIMIT": "1",
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "{8}",
            },
        )
        self.assertEqual(baseline_execution.launcher, (sys.executable, "-c", "8"))
        workers = threads * ranks
        return [
            benchmark.CaseResult(
                cases[0],
                tuple(
                    benchmark.Measurement(value / workers, 800, 23)
                    for value in (16.0, 12.0)
                ),
                tuple(benchmark.Measurement(value, 800, 23) for value in (8.0, 12.0)),
            )
        ]

    def test_sweep_writes_valid_paired_results_and_separate_rank_curves(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build, output = self.create_build(root), root / "reports"
            with (
                mock.patch.object(
                    benchmark, "run_benchmarks", side_effect=self.measurements
                ) as run,
                mock.patch.object(
                    hybrid_scaling,
                    "save_scaling_plot",
                    wraps=hybrid_scaling.save_scaling_plot,
                ) as plot,
            ):
                status, stdout, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(output),
                        "--cpus",
                        "8",
                        "--layouts",
                        "2x2",
                        "2x4",
                        "4x2",
                        "--runs",
                        "2",
                        "--warmup",
                        "0",
                        "--timeout",
                        "12.5",
                    ]
                )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(run.call_count, 3)
            for layout in ("2x2", "2x4", "4x2"):
                document = json.loads((output / f"hybrid-{layout}.json").read_text())
                workers = hybrid_scaling.hybrid_layout(layout).workers
                self.assertEqual(
                    document["comparison"]["paired_round_wall_speedup"]["median"],
                    0.75 * workers,
                )
                result = check_parallel_scaling.evaluate_report(
                    check_parallel_scaling.TopologySpec(
                        f"hybrid-{layout}", workers, output / f"hybrid-{layout}.json"
                    )
                )
                self.assertEqual(result.suite_speedup, 0.75 * workers)
            saved = (output / "scaling.txt").read_text()
            self.assertIn("hybrid-2x4", saved)
            self.assertIn("hybrid-4x2", saved)
            self.assertIn("75.0%", saved)
            self.assertIn(saved.strip(), stdout)
            self.assertEqual(
                (output / "scaling.png").read_bytes()[:8], b"\x89PNG\r\n\x1a\n"
            )
            self.assertTrue((output / "scaling.pdf").read_bytes().startswith(b"%PDF-"))
            figure = plot.call_args.args[3]
            speedup_axis, efficiency_axis = figure.axes
            self.assertEqual(list(speedup_axis.lines[0].get_xdata()), [4, 8])
            self.assertEqual(list(speedup_axis.lines[0].get_ydata()), [3.0, 6.0])
            self.assertEqual(list(speedup_axis.lines[1].get_xdata()), [8])
            self.assertEqual(list(efficiency_axis.lines[0].get_ydata()), [75.0, 75.0])
            metadata = json.loads((output / "cpu-topology.json").read_text())
            self.assertEqual(
                metadata["topology"], json.loads(json.dumps(asdict(self.topology)))
            )
            self.assertEqual(metadata["allocated_physical_cores"], 8)
            self.assertEqual(
                [
                    (item["ranks"], item["threads"], item["workers"])
                    for item in metadata["runs"]
                ],
                [(2, 2, 4), (2, 4, 8), (4, 2, 8)],
            )

    def test_dry_run_needs_no_tools_or_build_and_prints_requested_layouts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "missing output"
            with (
                mock.patch.object(hybrid_scaling.shutil, "which", return_value=None),
                mock.patch.object(benchmark, "main") as run,
                mock.patch.dict(sys.modules, {"matplotlib.backends.backend_agg": None}),
            ):
                status, stdout, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(Path(temporary) / "missing build"),
                        "--output-dir",
                        str(output),
                        "--cpus",
                        "128",
                        "--layouts",
                        "2x64",
                        "64x2",
                        "--dry-run",
                    ]
                )
            self.assertEqual(status, 0, stderr)
            run.assert_not_called()
            self.assertFalse(output.exists())
            commands = self.benchmark_commands(stdout)
            self.assertEqual(len(commands), 2)
            for command, ranks, threads in zip(commands, (2, 64), (64, 2), strict=True):
                parsed = benchmark.create_argument_parser().parse_args(command)
                self.assertEqual(
                    (parsed.runs, parsed.warmup, parsed.timeout), (6, 1, 600.0)
                )
                self.assertEqual(parsed.candidate_launcher[-1], str(ranks))
                self.assertIn(f"slot:PE={threads}", parsed.candidate_launcher)
                self.assertIn("localhost:128", parsed.candidate_launcher)
                self.assertEqual(
                    parsed.json_output.name, f"hybrid-{ranks}x{threads}.json"
                )
            self.assertIn("hybrid-2x64:128:", stdout)
            self.assertIn("hybrid-64x2:128:", stdout)

    def test_launcher_cleans_child_scheduler_identity_and_binds_core_teams(
        self,
    ) -> None:
        identities = {
            "SLURM_JOB_ID": "42",
            "SLURM_NTASKS": "1",
            "PMI_RANK": "0",
            "PMIX_RANK": "0",
            "OMPI_COMM_WORLD_SIZE": "1",
        }
        with mock.patch.dict(
            os.environ, {**identities, "OMPI_MCA_btl": "self,vader"}, clear=True
        ):
            launcher = hybrid_scaling.mpi_launcher(
                Path("/env with spaces/bin/mpirun"), 128, hybrid_scaling.Layout(8, 16)
            )
            for name, value in identities.items():
                index = launcher.index(name)
                self.assertEqual(launcher[index - 1], "-u")
                self.assertEqual(os.environ[name], value)
            self.assertNotIn("OMPI_MCA_btl", launcher)
            self.assertIn("/env with spaces/bin/mpirun", launcher)
            self.assertEqual(
                launcher[-10:],
                (
                    "--host",
                    "localhost:128",
                    "--map-by",
                    "slot:PE=16",
                    "--bind-to",
                    "core",
                    "--nooversubscribe",
                    "--report-bindings",
                    "-n",
                    "8",
                ),
            )

    def test_invalid_layouts_and_cpu_budgets_are_parser_errors(self) -> None:
        for extra in (
            ["--layouts", "1x8"],
            ["--layouts", "8x1"],
            ["--layouts", "2X2"],
            ["--layouts", "two"],
            ["--layouts", "2x2", "2x2"],
            ["--layouts", "2x8"],
            ["--cpus", "0", "--layouts", "2x2"],
        ):
            with (
                self.subTest(extra=extra),
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit) as error,
            ):
                hybrid_scaling.main(["--cpus", "8", "--dry-run", *extra])
            self.assertEqual(error.exception.code, 2)

    def test_execution_requires_physical_capacity_and_affinity(self) -> None:
        for topology in (
            replace(self.topology, physical_cores=4),
            replace(self.topology, physical_cores=None),
            replace(self.topology, affinity_cpus=None),
            replace(self.topology, cpu_order=None),
        ):
            with (
                self.subTest(topology=topology),
                mock.patch.object(benchmark, "main") as run,
            ):
                self.detect.return_value = topology
                status, _, stderr = self.run_main(["--cpus", "8", "--layouts", "2x2"])
                self.assertEqual(status, 1, stderr)
                run.assert_not_called()

    def test_missing_launcher_or_solver_is_rejected_before_output(self) -> None:
        for missing in ("taskset", "mpirun", "ParallelAssemblyCppHybrid"):
            with (
                self.subTest(missing=missing),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                build, output = self.create_build(root), root / "reports"
                if missing == "ParallelAssemblyCppHybrid":
                    paclitaxel_scaling.executable_path(build, missing).unlink()
                real_which = hybrid_scaling.shutil.which

                def find_tool(
                    command: str,
                    path: str | None = None,
                    excluded: str = missing,
                    resolver: Callable[..., str | None] = real_which,
                ) -> str | None:
                    return None if command == excluded else resolver(command, path=path)

                with (
                    mock.patch.object(
                        hybrid_scaling.shutil, "which", side_effect=find_tool
                    ),
                    mock.patch.object(benchmark, "main") as run,
                ):
                    status, _, stderr = self.run_main(
                        [
                            "--build-dir",
                            str(build),
                            "--output-dir",
                            str(output),
                            "--cpus",
                            "8",
                            "--layouts",
                            "2x2",
                        ]
                    )
                self.assertEqual(status, 1, stderr)
                self.assertIn(missing, stderr)
                run.assert_not_called()
                self.assertFalse(output.exists())

    def test_invalid_report_cannot_produce_summary_or_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build, output = self.create_build(root), root / "reports"
            with (
                mock.patch.object(benchmark, "main", return_value=0),
                mock.patch.object(
                    check_parallel_scaling,
                    "evaluate_specs",
                    side_effect=check_parallel_scaling.ScalingError("invalid report"),
                ),
                mock.patch.object(hybrid_scaling, "save_scaling_plot") as plot,
            ):
                status, _, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(output),
                        "--cpus",
                        "8",
                        "--layouts",
                        "2x2",
                    ]
                )
            self.assertEqual(status, 1, stderr)
            self.assertIn("invalid report", stderr)
            plot.assert_not_called()
            self.assertEqual(
                {path.name for path in output.iterdir()}, {"cpu-topology.json"}
            )

    def test_failure_retains_completed_reports_and_skips_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build, output = self.create_build(root), root / "reports"
            calls = 0

            def run_benchmark(arguments: Sequence[str]) -> int:
                nonlocal calls
                calls += 1
                if calls == 1:
                    parsed = benchmark.create_argument_parser().parse_args(arguments)
                    parsed.json_output.write_text("completed evidence\n")
                    return 0
                return 17

            with (
                mock.patch.object(benchmark, "main", side_effect=run_benchmark),
                mock.patch.object(check_parallel_scaling, "evaluate_specs") as evaluate,
            ):
                status, _, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(output),
                        "--cpus",
                        "8",
                        "--layouts",
                        "2x2",
                        "2x4",
                        "4x2",
                    ]
                )
            self.assertEqual(status, 17)
            self.assertEqual(calls, 2)
            self.assertIn("completed JSON reports were retained", stderr)
            self.assertEqual(
                (output / "hybrid-2x2.json").read_text(), "completed evidence\n"
            )
            evaluate.assert_not_called()
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"hybrid-2x2.json", "cpu-topology.json"},
            )

    def test_existing_outputs_fail_before_measurement(self) -> None:
        for filename in (
            "hybrid-2x2.json",
            "scaling.txt",
            "scaling.png",
            "scaling.pdf",
            "cpu-topology.json",
        ):
            with (
                self.subTest(filename=filename),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                build, output = self.create_build(root), root / "reports"
                output.mkdir()
                existing = output / filename
                existing.write_text("existing evidence\n")
                with mock.patch.object(benchmark, "main") as run:
                    status, _, stderr = self.run_main(
                        [
                            "--build-dir",
                            str(build),
                            "--output-dir",
                            str(output),
                            "--cpus",
                            "8",
                            "--layouts",
                            "2x2",
                        ]
                    )
                self.assertEqual(status, 1, stderr)
                self.assertIn("refusing to overwrite", stderr)
                run.assert_not_called()
                self.assertEqual(existing.read_text(), "existing evidence\n")


if __name__ == "__main__":
    unittest.main()
