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
    paclitaxel_scaling,
)

if TYPE_CHECKING:
    from collections.abc import Sequence


class PaclitaxelScalingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.topology = cpu_topology.CpuTopology(
            logical_cpus=4,
            physical_cores=2,
            cpu_order=(0, 2, 1, 3),
            affinity_cpus=(0, 1, 2, 3),
            model="Fixture CPU",
            slurm_cpus_per_task=None,
            source="fixture",
        )
        detection_patch = mock.patch.object(
            cpu_topology, "detect_cpu_topology", return_value=self.topology
        )
        self.detect = detection_patch.start()
        self.addCleanup(detection_patch.stop)
        real_which = paclitaxel_scaling.shutil.which
        self.taskset: str | None = sys.executable

        def find_executable(command: str, path: str | None = None) -> str | None:
            return (
                self.taskset if command == "taskset" else real_which(command, path=path)
            )

        taskset_patch = mock.patch.object(
            paclitaxel_scaling.shutil, "which", side_effect=find_executable
        )
        self.find_taskset = taskset_patch.start()
        self.addCleanup(taskset_patch.stop)

    def create_build(self, directory: Path, *, telemetry: bool = False) -> Path:
        build = directory / "parallel build"
        build.mkdir()
        names = ["AssemblyCpp", "AssemblyCppOMP"]
        if telemetry:
            names.append("AssemblyCppOMPTelemetry")
        for name in names:
            executable = build / (name + (".exe" if os.name == "nt" else ""))
            executable.write_text(f"fixture for {name}\n", encoding="utf-8")
            executable.chmod(0o755)
        return build

    def run_main(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = paclitaxel_scaling.main(arguments)
        return status, stdout.getvalue(), stderr.getvalue()

    def benchmark_commands(self, output: str) -> list[list[str]]:
        commands = []
        for line in output.splitlines():
            if line.startswith("#"):
                continue
            command = shlex.split(line)
            scripts = [
                index
                for index, part in enumerate(command)
                if Path(part).name == "benchmark.py"
            ]
            if scripts:
                commands.append(command[scripts[0] + 1 :])
        return commands

    def thread_counts(self, output: str) -> list[int]:
        return [
            int(
                dict(
                    benchmark.create_argument_parser()
                    .parse_args(command)
                    .candidate_environment
                )["OMP_NUM_THREADS"]
            )
            for command in self.benchmark_commands(output)
        ]

    def test_sweep_records_topology_and_checks_real_paired_reports(self) -> None:
        launcher = shlex.quote(sys.executable)

        def measurements(
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
            self.assertEqual(executable.stem, "AssemblyCppOMP")
            self.assertIsNotNone(baseline_executable)
            self.assertEqual(baseline_executable.stem, "AssemblyCpp")
            self.assertIsNotNone(telemetry_executable)
            self.assertEqual(telemetry_executable.stem, "AssemblyCppOMPTelemetry")
            self.assertEqual([case.name for case in cases], ["paclitaxel"])
            self.assertEqual((runs, warmup, timeout), (2, 0, 12.5))
            self.assertIsNotNone(candidate_execution)
            self.assertIsNotNone(baseline_execution)
            self.assertEqual(candidate_execution.arguments, ("--parallel=on",))
            self.assertEqual(baseline_execution.arguments, ("--parallel=off",))
            self.assertEqual(candidate_execution.launcher, (sys.executable, "-B"))
            self.assertEqual(baseline_execution.launcher, (sys.executable, "-u"))
            candidate_environment = dict(candidate_execution.environment)
            workers = int(candidate_environment["OMP_NUM_THREADS"])
            shared_environment = {
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
            }
            self.assertEqual(
                candidate_environment,
                {
                    **shared_environment,
                    "OMP_PLACES": "threads",
                    "OMP_NUM_THREADS": str(workers),
                    "OMP_THREAD_LIMIT": str(workers),
                },
            )
            self.assertEqual(
                dict(baseline_execution.environment),
                {
                    **shared_environment,
                    "OMP_PLACES": "{0}",
                    "OMP_NUM_THREADS": "1",
                    "OMP_THREAD_LIMIT": "1",
                },
            )
            case = cases[0]
            self.assertEqual(case.expected_assembly_index, 23)
            return [
                benchmark.CaseResult(
                    case,
                    tuple(
                        benchmark.Measurement(8.0 / workers, 800, 23)
                        for _ in range(runs)
                    ),
                    tuple(benchmark.Measurement(8.0, 800, 23) for _ in range(runs)),
                )
            ]

        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = self.create_build(directory, telemetry=True)
            output = directory / "scaling reports"
            with (
                mock.patch.object(
                    benchmark, "run_benchmarks", side_effect=measurements
                ) as run,
                mock.patch.object(
                    check_parallel_scaling,
                    "evaluate_specs",
                    wraps=check_parallel_scaling.evaluate_specs,
                ) as evaluate,
                mock.patch.object(
                    check_parallel_scaling,
                    "print_report",
                    wraps=check_parallel_scaling.print_report,
                ) as report,
                mock.patch.object(
                    paclitaxel_scaling,
                    "save_scaling_plot",
                    wraps=paclitaxel_scaling.save_scaling_plot,
                ) as plot,
            ):
                status, stdout, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(output),
                        "--threads",
                        "2",
                        "4",
                        "--runs",
                        "2",
                        "--warmup",
                        "0",
                        "--timeout",
                        "12.5",
                        "--baseline-launcher",
                        f"{launcher} -u",
                        "--candidate-launcher",
                        f"{launcher} -B",
                        "--telemetry",
                    ]
                )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(stderr, "")
            self.assertEqual(run.call_count, 2)
            evaluate.assert_called_once()
            report.assert_called_once()
            plot.assert_called_once()
            saved_plot = output / "scaling.png"
            self.assertEqual(saved_plot.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertIn(f"Scaling plot: {saved_plot}", stdout)
            saved_pdf = output / "scaling.pdf"
            self.assertTrue(saved_pdf.read_bytes().startswith(b"%PDF-"))
            self.assertIn(f"Scaling plot: {saved_pdf}", stdout)
            self.assertEqual(plot.call_args.args[1], saved_plot)
            figure = plot.call_args.args[2]
            speedup_axis, efficiency_axis = figure.axes
            self.assertEqual(list(speedup_axis.lines[0].get_xdata()), [2, 4])
            self.assertEqual(list(speedup_axis.lines[0].get_ydata()), [2.0, 4.0])
            self.assertEqual(list(efficiency_axis.lines[0].get_xdata()), [2, 4])
            self.assertEqual(list(efficiency_axis.lines[0].get_ydata()), [100.0, 100.0])
            specs = evaluate.call_args.args[0]
            self.assertEqual(
                [(spec.label, spec.workers, spec.path) for spec in specs],
                [("omp", count, output / f"omp-{count}.json") for count in (2, 4)],
            )
            for workers in (2, 4):
                document = json.loads(
                    (output / f"omp-{workers}.json").read_text(encoding="utf-8")
                )
                self.assertEqual(document["suite"], "profile")
                self.assertEqual(
                    [case["name"] for case in document["cases"]], ["paclitaxel"]
                )
                self.assertEqual(
                    document["comparison"]["paired_round_wall_speedup"]["median"],
                    workers,
                )
            saved_report = (output / "scaling.txt").read_text(encoding="utf-8")
            self.assertIn("case paclitaxel", saved_report)
            self.assertIn("2.000x", saved_report)
            self.assertIn("4.000x", saved_report)
            self.assertIn("100.0%", saved_report)
            self.assertIn(saved_report.strip(), stdout)
            self.assertIn("Fixture CPU", saved_report.split("case paclitaxel")[0])
            topology_report = json.loads(
                (output / "cpu-topology.json").read_text(encoding="utf-8")
            )
            self.assertEqual(topology_report["schema_version"], 1)
            self.assertEqual(
                topology_report["topology"],
                json.loads(json.dumps(asdict(self.topology))),
            )
            self.assertEqual(topology_report["threads"], [2, 4])
            self.assertEqual(topology_report["selection"], "explicit")
            self.assertEqual(topology_report["cpu_kind"], "logical")
            self.assertEqual(topology_report["candidate_placement"], "launcher")
            self.assertEqual(
                topology_report["runs"],
                [
                    {
                        "threads": count,
                        "cpu_ids": None,
                        "report": str(output / f"omp-{count}.json"),
                    }
                    for count in (2, 4)
                ],
            )
            self.assertEqual(
                topology_report["baseline_launcher"], [sys.executable, "-u"]
            )

    def test_dry_run_prints_default_sweep_without_build_or_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = directory / "missing build"
            output = directory / "missing output"
            with (
                mock.patch.object(benchmark, "main") as run,
                mock.patch.dict(sys.modules, {"matplotlib.backends.backend_agg": None}),
            ):
                status, stdout, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(output),
                        "--dry-run",
                    ]
                )
            self.assertEqual(status, 0, stderr)
            run.assert_not_called()
            self.assertFalse(build.exists())
            self.assertFalse(output.exists())
            self.assertTrue(stdout.startswith("# "))
            self.assertIn("Fixture CPU", stdout)
            self.detect.assert_called_once()
            commands = self.benchmark_commands(stdout)
            self.assertEqual(len(commands), 3, stdout)
            for command, workers in zip(commands, (2, 3, 4), strict=True):
                arguments = benchmark.create_argument_parser().parse_args(command)
                self.assertEqual(arguments.suite, "profile")
                self.assertEqual(arguments.case, ["paclitaxel"])
                self.assertEqual((arguments.runs, arguments.warmup), (6, 1))
                self.assertEqual(arguments.timeout, 600)
                self.assertEqual(arguments.candidate_parallel, "on")
                self.assertEqual(arguments.baseline_parallel, "off")
                self.assertEqual(
                    dict(arguments.candidate_environment)["OMP_THREAD_LIMIT"],
                    str(workers),
                )
                self.assertEqual(
                    dict(arguments.candidate_environment)["OMP_PLACES"],
                    ",".join(f"{{{cpu}}}" for cpu in self.topology.cpu_order[:workers]),
                )
                self.assertEqual(
                    dict(arguments.baseline_environment)["OMP_PLACES"], "{0}"
                )
                self.assertEqual(
                    arguments.baseline_launcher, (sys.executable, "-c", "0")
                )
                self.assertEqual(arguments.json_output, output / f"omp-{workers}.json")
            self.assertIn("check_parallel_scaling.py", stdout)
            for workers in (2, 3, 4):
                self.assertIn(f"omp:{workers}:{output / f'omp-{workers}.json'}", stdout)

    def test_invalid_thread_counts_are_parser_errors(self) -> None:
        for threads in (["1"], ["0"], ["-2"], ["2", "2"], ["many"]):
            with (
                self.subTest(threads=threads),
                contextlib.redirect_stderr(io.StringIO()),
                mock.patch.object(benchmark, "main") as run,
                self.assertRaises(SystemExit) as raised,
            ):
                paclitaxel_scaling.main(["--dry-run", "--threads", *threads])
            self.assertEqual(raised.exception.code, 2)
            run.assert_not_called()

    def test_automatic_sweep_covers_all_available_hybrid_cpu_threads(self) -> None:
        cpu_order = (*range(8, 28), *range(40, 48))
        self.detect.return_value = replace(
            self.topology,
            logical_cpus=28,
            physical_cores=20,
            cpu_order=cpu_order,
            affinity_cpus=cpu_order,
            model="Hybrid fixture CPU",
        )
        status, stdout, stderr = self.run_main(["--dry-run"])
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.thread_counts(stdout), list(range(2, 29)))
        commands = self.benchmark_commands(stdout)
        for count in (20, 21, 28):
            arguments = benchmark.create_argument_parser().parse_args(
                commands[count - 2]
            )
            places = dict(arguments.candidate_environment)["OMP_PLACES"]
            self.assertEqual(
                places, ",".join(f"{{{cpu}}}" for cpu in cpu_order[:count])
            )
        first = benchmark.create_argument_parser().parse_args(commands[0])
        self.assertEqual(first.baseline_launcher, (sys.executable, "-c", "8"))
        self.assertEqual(dict(first.baseline_environment)["OMP_PLACES"], "{8}")

    def test_physical_core_sweep_excludes_smt_siblings(self) -> None:
        cpu_order = (*range(8, 28), *range(40, 48))
        self.detect.return_value = replace(
            self.topology,
            logical_cpus=28,
            physical_cores=20,
            cpu_order=cpu_order,
            affinity_cpus=cpu_order,
        )
        status, stdout, stderr = self.run_main(["--dry-run", "--physical-cores-only"])
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.thread_counts(stdout), list(range(2, 21)))
        last = benchmark.create_argument_parser().parse_args(
            self.benchmark_commands(stdout)[-1]
        )
        self.assertEqual(
            dict(last.candidate_environment)["OMP_PLACES"],
            ",".join(f"{{{cpu}}}" for cpu in range(8, 28)),
        )

    def test_effective_cpu_capacity_limits_sweep(self) -> None:
        self.detect.return_value = replace(
            self.topology,
            logical_cpus=2,
            physical_cores=2,
            cpu_order=(8, 10),
            affinity_cpus=(8, 9, 10, 11),
            slurm_cpus_per_task=2,
        )
        status, stdout, stderr = self.run_main(["--dry-run"])
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.thread_counts(stdout), [2])
        command = benchmark.create_argument_parser().parse_args(
            self.benchmark_commands(stdout)[0]
        )
        self.assertEqual(dict(command.candidate_environment)["OMP_PLACES"], "{8},{10}")

    def test_single_available_cpu_cannot_run_parallel_sweep(self) -> None:
        self.detect.return_value = replace(
            self.topology,
            logical_cpus=1,
            physical_cores=1,
            cpu_order=(8,),
            affinity_cpus=(8,),
        )
        for arguments in (["--dry-run"], ["--dry-run", "--threads", "2"]):
            with (
                self.subTest(arguments=arguments),
                mock.patch.object(benchmark, "main") as run,
            ):
                status, _, stderr = self.run_main(arguments)
            self.assertNotEqual(status, 0)
            self.assertIn("error:", stderr)
            run.assert_not_called()

    def test_manual_local_sweeps_cannot_exceed_selected_capacity(self) -> None:
        for arguments in (
            ["--threads", "5"],
            ["--physical-cores-only", "--threads", "3"],
        ):
            with (
                self.subTest(arguments=arguments),
                mock.patch.object(benchmark, "main") as run,
            ):
                status, _, stderr = self.run_main(["--dry-run", *arguments])
            self.assertNotEqual(status, 0)
            self.assertIn("error:", stderr)
            run.assert_not_called()

    def test_custom_candidate_launcher_requires_explicit_thread_counts(self) -> None:
        launcher = shlex.quote(sys.executable)
        for arguments in (
            [],
            ["--physical-cores-only", "--threads", "2"],
        ):
            with (
                self.subTest(arguments=arguments),
                mock.patch.object(benchmark, "main") as run,
            ):
                status, _, stderr = self.run_main(
                    ["--dry-run", "--candidate-launcher", launcher, *arguments]
                )
            self.assertNotEqual(status, 0)
            self.assertIn("launcher", stderr)
            run.assert_not_called()

    def test_custom_candidate_launcher_allows_remote_capacity_without_local_ids(
        self,
    ) -> None:
        launcher = shlex.quote(sys.executable)
        status, stdout, stderr = self.run_main(
            ["--dry-run", "--candidate-launcher", launcher, "--threads", "8", "2"]
        )
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.thread_counts(stdout), [8, 2])
        for command in self.benchmark_commands(stdout):
            arguments = benchmark.create_argument_parser().parse_args(command)
            self.assertEqual(arguments.candidate_launcher, (sys.executable,))
            self.assertEqual(
                dict(arguments.candidate_environment)["OMP_PLACES"], "threads"
            )

    def test_unknown_cpu_ids_use_runtime_placement(self) -> None:
        self.detect.return_value = replace(
            self.topology, physical_cores=None, cpu_order=None, affinity_cpus=None
        )
        self.taskset = None
        status, stdout, stderr = self.run_main(["--dry-run"])
        self.assertEqual(status, 0, stderr)
        self.assertEqual(self.thread_counts(stdout), [2, 3, 4])
        for command in self.benchmark_commands(stdout):
            arguments = benchmark.create_argument_parser().parse_args(command)
            self.assertFalse(arguments.baseline_launcher)
            self.assertEqual(
                dict(arguments.candidate_environment)["OMP_PLACES"], "threads"
            )
            self.assertEqual(
                dict(arguments.baseline_environment)["OMP_PLACES"], "threads"
            )

    def test_physical_core_sweep_requires_known_physical_count(self) -> None:
        self.detect.return_value = replace(self.topology, physical_cores=None)
        with mock.patch.object(benchmark, "main") as run:
            status, _, stderr = self.run_main(["--dry-run", "--physical-cores-only"])
        self.assertNotEqual(status, 0)
        self.assertIn("physical", stderr)
        run.assert_not_called()

    def test_unavailable_taskset_keeps_detected_openmp_places(self) -> None:
        self.taskset = None
        status, stdout, stderr = self.run_main(["--dry-run", "--threads", "2"])
        self.assertEqual(status, 0, stderr)
        arguments = benchmark.create_argument_parser().parse_args(
            self.benchmark_commands(stdout)[0]
        )
        self.assertFalse(arguments.baseline_launcher)
        self.assertEqual(dict(arguments.baseline_environment)["OMP_PLACES"], "{0}")
        self.assertEqual(dict(arguments.candidate_environment)["OMP_PLACES"], "{0},{2}")

    def test_automatic_topology_report_is_written_before_measurements(self) -> None:
        for physical_only in (False, True):
            with (
                self.subTest(physical_only=physical_only),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                directory = Path(temp_directory)
                build = self.create_build(directory)
                output = directory / "reports"
                topology_report = output / "cpu-topology.json"
                captured_documents = []

                def stop_measurements(
                    _: Sequence[str],
                    report_path: Path = topology_report,
                    documents: list[dict[str, object]] = captured_documents,
                ) -> int:
                    documents.append(
                        json.loads(report_path.read_text(encoding="utf-8"))
                    )
                    return 17

                arguments = ["--build-dir", str(build), "--output-dir", str(output)]
                if physical_only:
                    arguments.append("--physical-cores-only")
                with mock.patch.object(
                    benchmark, "main", side_effect=stop_measurements
                ) as run:
                    status, _, _ = self.run_main(arguments)
                self.assertEqual(status, 17)
                run.assert_called_once()
                self.assertFalse((output / "scaling.txt").exists())
                self.assertFalse((output / "scaling.png").exists())
                self.assertFalse((output / "scaling.pdf").exists())
                document = captured_documents[0]
                counts = [2] if physical_only else [2, 3, 4]
                self.assertEqual(document["schema_version"], 1)
                self.assertEqual(document["selection"], "automatic")
                self.assertEqual(
                    document["cpu_kind"], "physical" if physical_only else "logical"
                )
                self.assertEqual(document["candidate_placement"], "detected_cpu_ids")
                self.assertEqual(document["threads"], counts)
                self.assertEqual(
                    document["baseline_launcher"], [sys.executable, "-c", "0"]
                )
                self.assertEqual(
                    document["runs"],
                    [
                        {
                            "threads": count,
                            "cpu_ids": list(self.topology.cpu_order[:count]),
                            "report": str(output / f"omp-{count}.json"),
                        }
                        for count in counts
                    ],
                )

    def test_failed_benchmark_stops_sweep_before_scaling_check(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = self.create_build(directory)
            output = directory / "reports"
            with (
                mock.patch.object(benchmark, "main", side_effect=[0, 17]) as run,
                mock.patch.object(check_parallel_scaling, "evaluate_specs") as evaluate,
                mock.patch.object(check_parallel_scaling, "print_report") as report,
            ):
                status, _, _ = self.run_main(
                    ["--build-dir", str(build), "--output-dir", str(output)]
                )
            self.assertNotEqual(status, 0)
            self.assertEqual(run.call_count, 2)
            evaluate.assert_not_called()
            report.assert_not_called()
            self.assertFalse((output / "scaling.txt").exists())
            self.assertFalse((output / "scaling.png").exists())
            self.assertFalse((output / "scaling.pdf").exists())

    def test_invalid_results_do_not_create_summary_or_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = self.create_build(directory)
            output = directory / "reports"
            with (
                mock.patch.object(benchmark, "main", return_value=0),
                mock.patch.object(
                    check_parallel_scaling,
                    "evaluate_specs",
                    side_effect=check_parallel_scaling.ScalingError("invalid report"),
                ),
                mock.patch.object(paclitaxel_scaling, "save_scaling_plot") as plot,
            ):
                status, _, stderr = self.run_main(
                    ["--build-dir", str(build), "--output-dir", str(output)]
                )
            self.assertEqual(status, 1)
            self.assertIn("invalid report", stderr)
            plot.assert_not_called()
            self.assertFalse((output / "scaling.txt").exists())
            self.assertFalse((output / "scaling.png").exists())
            self.assertFalse((output / "scaling.pdf").exists())

    def test_missing_matplotlib_is_rejected_before_measurements_or_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = self.create_build(directory)
            output = directory / "reports"
            with (
                mock.patch.object(benchmark, "main") as run,
                mock.patch.dict(sys.modules, {"matplotlib.backends.backend_agg": None}),
            ):
                status, _, stderr = self.run_main(
                    ["--build-dir", str(build), "--output-dir", str(output)]
                )
            self.assertEqual(status, 1)
            self.assertIn("Matplotlib", stderr)
            run.assert_not_called()
            self.assertFalse(output.exists())

    def test_existing_reports_are_rejected_before_any_measurement(self) -> None:
        for filename in (
            "omp-4.json",
            "scaling.txt",
            "cpu-topology.json",
            "scaling.png",
            "scaling.pdf",
        ):
            with (
                self.subTest(filename=filename),
                tempfile.TemporaryDirectory() as temp_directory,
            ):
                directory = Path(temp_directory)
                build = self.create_build(directory)
                output = directory / "reports"
                output.mkdir()
                existing = output / filename
                existing.write_text("existing evidence\n", encoding="utf-8")
                with mock.patch.object(benchmark, "main") as run:
                    status, _, stderr = self.run_main(
                        ["--build-dir", str(build), "--output-dir", str(output)]
                    )
                self.assertNotEqual(status, 0)
                self.assertIn(str(existing), stderr)
                run.assert_not_called()
                self.assertEqual(
                    existing.read_text(encoding="utf-8"), "existing evidence\n"
                )
                self.assertEqual(list(output.iterdir()), [existing])

    def test_missing_binary_is_rejected_before_any_measurement(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            build = self.create_build(directory)
            with mock.patch.object(benchmark, "main") as run:
                status, _, stderr = self.run_main(
                    [
                        "--build-dir",
                        str(build),
                        "--output-dir",
                        str(directory / "reports"),
                        "--telemetry",
                    ]
                )
            self.assertNotEqual(status, 0)
            self.assertIn("AssemblyCppOMPTelemetry", stderr)
            run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
