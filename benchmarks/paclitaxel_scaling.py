"""Measure Paclitaxel OpenMP scaling against paired serial calculations."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import shlex
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

    from matplotlib.figure import Figure

if __package__:
    from . import benchmark, check_parallel_scaling, cpu_topology
else:
    import benchmark
    import check_parallel_scaling
    import cpu_topology


DEFAULT_BUILD_DIRECTORY = benchmark.REPOSITORY_ROOT / "build" / "parallel"
DEFAULT_OUTPUT_DIRECTORY = benchmark.REPOSITORY_ROOT / "build" / "paclitaxel-scaling"


@dataclass(frozen=True)
class ScalingRun:
    threads: int
    report: Path
    arguments: tuple[str, ...]


def parallel_threads(value: str) -> int:
    threads = int(value)
    if threads < 2:
        raise argparse.ArgumentTypeError(
            "must be at least two; the single-worker reference uses serial search"
        )
    return threads


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD_DIRECTORY,
        help=(
            "directory containing ParallelAssemblyCpp and ParallelAssemblyCppOMP "
            "(default: build/parallel)"
        ),
    )
    parser.add_argument(
        "--threads",
        type=parallel_threads,
        nargs="+",
        help="explicit thread counts (default: every count from 2 to available CPUs)",
    )
    parser.add_argument(
        "--physical-cores-only",
        action="store_true",
        help="stop at the available physical core count, excluding SMT siblings",
    )
    parser.add_argument(
        "--runs",
        type=benchmark.positive_int,
        default=6,
        help="measured paired rounds per thread count (default: 6)",
    )
    parser.add_argument(
        "--warmup",
        type=benchmark.non_negative_int,
        default=1,
        help="warm-up rounds per thread count (default: 1)",
    )
    parser.add_argument(
        "--timeout",
        type=benchmark.positive_float,
        default=600.0,
        help="per-calculation timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
        help=(
            "report directory (default: build/paclitaxel-scaling); "
            "saves scaling.txt, scaling.png and scaling.pdf; "
            "refuses to overwrite results"
        ),
    )
    for role in ("baseline", "candidate"):
        parser.add_argument(
            f"--{role}-launcher",
            type=benchmark.launcher_prefix,
            metavar="COMMAND",
            help=f"optional {role} placement prefix (for example, 'taskset -c 0-7')",
        )
    parser.add_argument(
        "--telemetry",
        action="store_true",
        help=(
            "collect one extra untimed ParallelAssemblyCppOMPTelemetry run "
            "per thread count"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print commands without requiring executables or creating output",
    )
    return parser


def executable_path(build_directory: Path, name: str) -> Path:
    return build_directory / (name + (".exe" if os.name == "nt" else ""))


def make_runs(
    arguments: argparse.Namespace, topology: cpu_topology.CpuTopology
) -> list[ScalingRun]:
    build_directory = arguments.build_dir.expanduser().resolve()
    output_directory = arguments.output_dir.expanduser().resolve()
    common = [
        "--baseline-executable",
        str(executable_path(build_directory, "ParallelAssemblyCpp")),
        "--baseline-parallel",
        "off",
        "--executable",
        str(executable_path(build_directory, "ParallelAssemblyCppOMP")),
        "--candidate-parallel",
        "on",
        "--suite",
        "profile",
        "--case",
        "paclitaxel",
        "--runs",
        str(arguments.runs),
        "--warmup",
        str(arguments.warmup),
        "--timeout",
        str(arguments.timeout),
    ]
    for role in ("baseline", "candidate"):
        launcher = getattr(arguments, f"{role}_launcher")
        if launcher:
            common.extend((f"--{role}-launcher", shlex.join(launcher)))
        for setting in ("OMP_DYNAMIC=FALSE", "OMP_PROC_BIND=close"):
            common.extend((f"--{role}-env", setting))
    baseline_places = (
        f"{{{topology.cpu_order[0]}}}" if topology.cpu_order else "threads"
    )
    for setting in (
        "OMP_NUM_THREADS=1",
        "OMP_THREAD_LIMIT=1",
        f"OMP_PLACES={baseline_places}",
    ):
        common.extend(("--baseline-env", setting))
    if arguments.telemetry:
        common.extend(
            (
                "--telemetry",
                "--telemetry-executable",
                str(
                    executable_path(build_directory, "ParallelAssemblyCppOMPTelemetry")
                ),
            )
        )

    runs = []
    for threads in arguments.threads:
        report = output_directory / f"omp-{threads}.json"
        places = "threads"
        if topology.cpu_order and not arguments.candidate_launcher:
            places = ",".join(f"{{{cpu}}}" for cpu in topology.cpu_order[:threads])
        run_arguments = [
            *common,
            "--candidate-env",
            f"OMP_NUM_THREADS={threads}",
            "--candidate-env",
            f"OMP_THREAD_LIMIT={threads}",
            "--candidate-env",
            f"OMP_PLACES={places}",
            "--json-output",
            str(report),
        ]
        runs.append(ScalingRun(threads, report, tuple(run_arguments)))
    return runs


def preflight(
    runs: Sequence[ScalingRun],
    summary_path: Path,
    topology_path: Path,
    plot_path: Path,
) -> None:
    """Reject unavailable tools and existing output before expensive calculations."""
    for path in [
        *(run.report for run in runs),
        summary_path,
        topology_path,
        plot_path,
        plot_path.with_suffix(".pdf"),
    ]:
        if path.exists() or path.is_symlink():
            raise benchmark.BenchmarkError(
                f"refusing to overwrite {path}; choose a new --output-dir"
            )
    first = benchmark.create_argument_parser().parse_args(runs[0].arguments)
    benchmark.resolve_requested_cases(first)
    for path in (
        first.executable,
        first.baseline_executable,
        first.telemetry_executable,
    ):
        if path is not None:
            benchmark.resolve_executable(path)
    for role in ("baseline", "candidate"):
        benchmark.create_execution_config(
            getattr(first, f"{role}_launcher"),
            getattr(first, f"{role}_environment"),
            role,
            parallel_mode=getattr(first, f"{role}_parallel"),
        )


def create_scaling_figure() -> Figure:
    """Load the headless plotting dependency before starting measurements."""
    try:
        # Dry runs must work without importing or configuring Matplotlib.
        from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: PLC0415
        from matplotlib.figure import Figure  # noqa: PLC0415
    except ImportError as error:
        raise benchmark.BenchmarkError(
            "Matplotlib is required to save scaling plots; "
            "install the benchmark dependencies before running the sweep"
        ) from error
    figure = Figure(figsize=(10, 4.5), constrained_layout=True)
    FigureCanvasAgg(figure)
    return figure


def save_scaling_plot(
    results: Sequence[check_parallel_scaling.ScalingResult],
    path: Path,
    figure: Figure,
) -> None:
    """Save PNG and PDF plots of validated speedup and parallel efficiency."""
    ordered_results = sorted(results, key=lambda result: result.spec.workers)
    threads = [result.spec.workers for result in ordered_results]
    speedups = [result.suite_speedup for result in ordered_results]
    efficiencies = [
        100 * result.suite_speedup / result.spec.workers for result in ordered_results
    ]
    speedup_axis, efficiency_axis = figure.subplots(1, 2)
    figure.suptitle("Paclitaxel OpenMP scaling")
    speedup_axis.plot(threads, speedups, "o-", label="Measured paired median")
    speedup_axis.plot([1, max(threads)], [1, max(threads)], "--", label="Ideal scaling")
    speedup_axis.axhline(1, color="gray", linewidth=0.8, linestyle=":")
    speedup_axis.set_ylabel("Wall-time speedup (serial / parallel)")
    efficiency_axis.plot(threads, efficiencies, "o-", label="Measured efficiency")
    efficiency_axis.axhline(100, color="gray", linestyle="--", label="Ideal efficiency")
    efficiency_axis.set_ylabel("Parallel efficiency (%)")
    for axis in (speedup_axis, efficiency_axis):
        axis.set_xlabel("OpenMP threads")
        axis.set_xlim(1, max(threads) + 0.5)
        axis.set_ylim(bottom=0)
        if len(threads) <= 16:
            axis.set_xticks(threads)
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize="small")
    figure.savefig(path, format="png", dpi=160)
    figure.savefig(path.with_suffix(".pdf"), format="pdf")


def configure_threads(
    arguments: argparse.Namespace, topology: cpu_topology.CpuTopology
) -> None:
    """Select the full available range unless the user supplies a subset."""
    if arguments.candidate_launcher and arguments.threads is None:
        raise benchmark.BenchmarkError(
            "automatic CPU detection cannot inspect a candidate launcher's allocation; "
            "run this driver inside the launcher, or supply explicit --threads"
        )
    capacity = topology.logical_cpus
    if arguments.physical_cores_only:
        if arguments.candidate_launcher:
            raise benchmark.BenchmarkError(
                "--physical-cores-only cannot inspect a candidate launcher's topology; "
                "run this driver inside the launcher instead"
            )
        if topology.physical_cores is None:
            raise benchmark.BenchmarkError(
                "physical core topology is unavailable; omit --physical-cores-only "
                "and supply explicit --threads if needed"
            )
        capacity = topology.physical_cores
    if arguments.threads is None:
        if capacity < 2:
            raise benchmark.BenchmarkError(
                "parallel scaling needs at least two available CPUs; "
                "the current affinity or allocation allows only one"
            )
        arguments.threads = list(range(2, capacity + 1))
    elif not arguments.candidate_launcher and max(arguments.threads) > capacity:
        raise benchmark.BenchmarkError(
            f"requested {max(arguments.threads)} threads but only {capacity} "
            "CPUs are available for the selected CPU kind"
        )
    if topology.cpu_order and not arguments.baseline_launcher:
        taskset = shutil.which("taskset")
        if taskset:
            arguments.baseline_launcher = (taskset, "-c", str(topology.cpu_order[0]))


def topology_report(
    arguments: argparse.Namespace,
    topology: cpu_topology.CpuTopology,
    runs: Sequence[ScalingRun],
    *,
    automatic: bool,
) -> dict[str, object]:
    detected_placement = bool(topology.cpu_order and not arguments.candidate_launcher)
    return {
        "schema_version": 1,
        "topology": asdict(topology),
        "selection": "automatic" if automatic else "explicit",
        "cpu_kind": "physical" if arguments.physical_cores_only else "logical",
        "threads": arguments.threads,
        "candidate_placement": (
            "launcher"
            if arguments.candidate_launcher
            else "detected_cpu_ids"
            if detected_placement
            else "runtime"
        ),
        "baseline_launcher": list(arguments.baseline_launcher or ()),
        "runs": [
            {
                "threads": run.threads,
                "cpu_ids": (
                    list(topology.cpu_order[: run.threads])
                    if detected_placement and topology.cpu_order
                    else None
                ),
                "report": str(run.report),
            }
            for run in runs
        ],
    }


def describe_topology(
    arguments: argparse.Namespace, topology: cpu_topology.CpuTopology
) -> str:
    physical = (
        "unknown physical core count"
        if topology.physical_cores is None
        else f"{topology.physical_cores} physical cores"
    )
    lines = [
        f"Driver CPU: {topology.model}",
        (
            f"Available: {physical}, {topology.logical_cpus} logical CPUs "
            f"({topology.source})"
        ),
        "Thread counts: " + ", ".join(map(str, arguments.threads)),
    ]
    if arguments.candidate_launcher:
        lines.append(
            "Candidate CPU placement is launcher-managed; "
            "driver limits do not verify it."
        )
    if arguments.baseline_launcher:
        lines.append("Serial launcher: " + shlex.join(arguments.baseline_launcher))
    else:
        lines.append("Serial affinity is inherited from the driver.")
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_argument_parser()
    arguments = parser.parse_args(argv)
    if arguments.threads and len(set(arguments.threads)) != len(arguments.threads):
        parser.error("--threads must not contain duplicates")
    automatic = arguments.threads is None
    try:
        topology = cpu_topology.detect_cpu_topology()
        configure_threads(arguments, topology)
    except (benchmark.BenchmarkError, ValueError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    runs = make_runs(arguments, topology)
    specs = [
        check_parallel_scaling.TopologySpec("omp", run.threads, run.report)
        for run in runs
    ]
    summary_path = arguments.output_dir.expanduser().resolve() / "scaling.txt"
    topology_path = summary_path.with_name("cpu-topology.json")
    plot_path = summary_path.with_name("scaling.png")
    description = describe_topology(arguments, topology)
    calculations = len(runs) * (
        2 * (arguments.runs + arguments.warmup) + int(arguments.telemetry)
    )
    if arguments.dry_run:
        for line in description.splitlines():
            print(f"# {line}")
        print(f"# Calculations including warm-ups: {calculations}")
        for run in runs:
            print(
                shlex.join(
                    (sys.executable, str(Path(benchmark.__file__)), *run.arguments)
                )
            )
        print(
            shlex.join(
                (
                    sys.executable,
                    str(Path(check_parallel_scaling.__file__)),
                    *(f"omp:{run.threads}:{run.report}" for run in runs),
                )
            )
        )
        return 0

    try:
        preflight(runs, summary_path, topology_path, plot_path)
        figure = create_scaling_figure()
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        topology_path.write_text(
            json.dumps(
                topology_report(arguments, topology, runs, automatic=automatic),
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print(description, end="")
        print(
            f"Paclitaxel scaling: {calculations} calculations including warm-ups.",
            flush=True,
        )
        for run in runs:
            print(f"\nOpenMP threads: {run.threads}", flush=True)
            status = benchmark.main(run.arguments)
            if status != 0:
                print(
                    f"Paclitaxel scaling stopped at {run.threads} threads; "
                    "no scaling summary was written.",
                    file=sys.stderr,
                )
                return status
        results = check_parallel_scaling.evaluate_specs(specs)
        summary = io.StringIO()
        summary.write(description + "\n")
        with contextlib.redirect_stdout(summary):
            check_parallel_scaling.print_report(results)
        save_scaling_plot(results, plot_path, figure)
        summary_path.write_text(summary.getvalue(), encoding="utf-8")
        print(f"\n{summary.getvalue()}", end="")
        print(f"\nScaling summary: {summary_path}")
        print(f"Scaling plot: {plot_path}")
        print(f"Scaling plot: {plot_path.with_suffix('.pdf')}")
    except (
        benchmark.BenchmarkError,
        check_parallel_scaling.ScalingError,
        OSError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print(
            "Paclitaxel scaling interrupted; no scaling summary was written.",
            file=sys.stderr,
        )
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
