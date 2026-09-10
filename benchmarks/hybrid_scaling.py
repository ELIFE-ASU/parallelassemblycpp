"""Measure single-node Paclitaxel MPI/OpenMP scaling against paired serial runs."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import re
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
    from . import benchmark, check_parallel_scaling, cpu_topology, paclitaxel_scaling
else:
    import benchmark
    import check_parallel_scaling
    import cpu_topology
    import paclitaxel_scaling


@dataclass(frozen=True)
class Layout:
    ranks: int
    threads: int

    @property
    def workers(self) -> int:
        return self.ranks * self.threads

    def __str__(self) -> str:
        return f"{self.ranks}x{self.threads}"


@dataclass(frozen=True)
class HybridRun:
    layout: Layout
    report: Path
    arguments: tuple[str, ...]

    @property
    def spec(self) -> check_parallel_scaling.TopologySpec:
        return check_parallel_scaling.TopologySpec(
            f"hybrid-{self.layout}", self.layout.workers, self.report
        )


def hybrid_layout(value: str) -> Layout:
    if re.fullmatch(r"[0-9]+x[0-9]+", value) is None:
        raise argparse.ArgumentTypeError("expected RANKSxTHREADS, for example 2x8")
    ranks, threads = map(int, value.split("x"))
    if ranks < 2 or threads < 2:
        raise argparse.ArgumentTypeError(
            "hybrid layouts need at least 2 ranks and 2 threads"
        )
    return Layout(ranks, threads)


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=paclitaxel_scaling.DEFAULT_BUILD_DIRECTORY,
        help="directory containing AssemblyCpp and AssemblyCppHybrid",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=benchmark.REPOSITORY_ROOT / "build" / "hybrid-scaling",
        help="report directory; refuses to overwrite existing results",
    )
    parser.add_argument(
        "--layouts",
        type=hybrid_layout,
        nargs="+",
        required=True,
        metavar="RANKSxTHREADS",
        help="distinct MPI rank / OpenMP thread layouts",
    )
    parser.add_argument(
        "--cpus",
        type=benchmark.positive_int,
        required=True,
        help="allocated physical cores on this node; each layout must fit",
    )
    parser.add_argument("--runs", type=benchmark.positive_int, default=6)
    parser.add_argument("--warmup", type=benchmark.non_negative_int, default=1)
    parser.add_argument("--timeout", type=benchmark.positive_float, default=600.0)
    parser.add_argument(
        "--mpirun",
        type=Path,
        help="Open MPI launcher (default: mpirun on PATH)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print commands without requiring executables or creating output",
    )
    return parser


def mpi_launcher(mpirun: Path, cpus: int, layout: Layout) -> tuple[str, ...]:
    """Launch locally without inheriting a one-task Slurm or parent MPI universe."""
    # Only the MPI child loses scheduler identity; the driver retains metadata
    # and every child remains inside Slurm's physical CPU/memory cgroup.
    inherited_identity = sorted(
        name
        for name in os.environ
        if name.startswith(("SLURM_", "PMI_", "PMIX_", "OMPI_COMM_WORLD_"))
        or name in {"OMPI_UNIVERSE_SIZE", "ORTE_HNP_URI", "PRTE_LAUNCHED"}
    )
    cleanup = tuple(
        argument for name in inherited_identity for argument in ("-u", name)
    )
    return (
        shutil.which("env") or "env",
        *cleanup,
        str(mpirun),
        "--host",
        f"localhost:{cpus}",
        "--map-by",
        f"slot:PE={layout.threads}",
        "--bind-to",
        "core",
        "--nooversubscribe",
        "--report-bindings",
        "-n",
        str(layout.ranks),
    )


def configure_placement(
    arguments: argparse.Namespace, topology: cpu_topology.CpuTopology
) -> None:
    """Check the physical allocation and match the OpenMP driver's serial pinning."""
    if not arguments.dry_run:
        if (
            topology.physical_cores is None
            or not topology.cpu_order
            or not topology.affinity_cpus
        ):
            raise benchmark.BenchmarkError(
                "hybrid scaling requires Linux CPU affinity and physical core topology"
            )
        if arguments.cpus > topology.physical_cores:
            raise benchmark.BenchmarkError(
                f"--cpus declares {arguments.cpus} physical cores but only "
                f"{topology.physical_cores} are available in this process's allocation"
            )
    # The same helper pins the OpenMP sweep's serial baseline to this first CPU.
    arguments.baseline_launcher = None
    arguments.candidate_launcher = None
    arguments.physical_cores_only = False
    arguments.threads = [1]
    paclitaxel_scaling.configure_threads(arguments, topology)
    if not arguments.dry_run and not arguments.baseline_launcher:
        raise benchmark.BenchmarkError("taskset is required to pin the serial baseline")
    if arguments.mpirun is None:
        resolved = shutil.which("mpirun")
        if resolved is None and not arguments.dry_run:
            raise benchmark.BenchmarkError("mpirun was not found on PATH")
        arguments.mpirun = Path(resolved or "mpirun")
    else:
        arguments.mpirun = arguments.mpirun.expanduser().resolve()
    if not arguments.dry_run:
        arguments.mpirun = benchmark.resolve_executable(arguments.mpirun)


def make_runs(
    arguments: argparse.Namespace, topology: cpu_topology.CpuTopology
) -> list[HybridRun]:
    build = arguments.build_dir.expanduser().resolve()
    output = arguments.output_dir.expanduser().resolve()
    common = [
        "--baseline-executable",
        str(paclitaxel_scaling.executable_path(build, "AssemblyCpp")),
        "--baseline-parallel",
        "off",
        "--executable",
        str(paclitaxel_scaling.executable_path(build, "AssemblyCppHybrid")),
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
    if arguments.baseline_launcher:
        common.extend(("--baseline-launcher", shlex.join(arguments.baseline_launcher)))
    runs = []
    for layout in arguments.layouts:
        report = output / f"hybrid-{layout}.json"
        command = [
            *common,
            "--candidate-launcher",
            shlex.join(mpi_launcher(arguments.mpirun, arguments.cpus, layout)),
            "--candidate-env",
            f"OMP_NUM_THREADS={layout.threads}",
            "--candidate-env",
            f"OMP_THREAD_LIMIT={layout.threads}",
            "--candidate-env",
            "OMP_PLACES=cores",
            "--json-output",
            str(report),
        ]
        runs.append(HybridRun(layout, report, tuple(command)))
    return runs


def save_scaling_plot(
    results: Sequence[check_parallel_scaling.ScalingResult],
    runs: Sequence[HybridRun],
    path: Path,
    figure: Figure,
) -> None:
    """Plot a separate curve for each rank count, including equal-worker layouts."""
    by_path = {result.spec.path: result for result in results}
    speedup_axis, efficiency_axis = figure.subplots(1, 2)
    figure.suptitle("Paclitaxel MPI/OpenMP hybrid scaling")
    for ranks in sorted({run.layout.ranks for run in runs}):
        group = sorted(
            (run for run in runs if run.layout.ranks == ranks),
            key=lambda run: run.layout.workers,
        )
        workers = [run.layout.workers for run in group]
        speedups = [by_path[run.report].suite_speedup for run in group]
        label = f"{ranks} MPI ranks"
        speedup_axis.plot(workers, speedups, "o-", label=label)
        efficiency_axis.plot(
            workers,
            [
                100 * speedup / count
                for count, speedup in zip(workers, speedups, strict=True)
            ],
            "o-",
            label=label,
        )
    max_workers = max(run.layout.workers for run in runs)
    speedup_axis.plot([1, max_workers], [1, max_workers], "--", label="Ideal scaling")
    speedup_axis.axhline(1, color="gray", linewidth=0.8, linestyle=":")
    speedup_axis.set_ylabel("Wall-time speedup (serial / hybrid)")
    efficiency_axis.axhline(100, color="gray", linestyle="--", label="Ideal efficiency")
    efficiency_axis.set_ylabel("Parallel efficiency (%)")
    counts = sorted({run.layout.workers for run in runs})
    for axis in (speedup_axis, efficiency_axis):
        axis.set_xlabel("Total workers (MPI ranks x OpenMP threads)")
        axis.set_xlim(1, max_workers + 0.5)
        axis.set_ylim(bottom=0)
        if len(counts) <= 16:
            axis.set_xticks(counts)
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize="small")
    figure.savefig(path, format="png", dpi=160)
    figure.savefig(path.with_suffix(".pdf"), format="pdf")


def topology_report(
    arguments: argparse.Namespace,
    topology: cpu_topology.CpuTopology,
    runs: Sequence[HybridRun],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "topology": asdict(topology),
        "allocated_physical_cores": arguments.cpus,
        "cpu_kind": "physical",
        "selection": "explicit",
        "candidate_placement": "Open MPI core binding",
        "baseline_launcher": list(arguments.baseline_launcher or ()),
        "runs": [
            {
                "ranks": run.layout.ranks,
                "threads": run.layout.threads,
                "workers": run.layout.workers,
                "report": str(run.report),
                "launcher": list(
                    mpi_launcher(arguments.mpirun, arguments.cpus, run.layout)
                ),
            }
            for run in runs
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_argument_parser()
    arguments = parser.parse_args(argv)
    if len(set(arguments.layouts)) != len(arguments.layouts):
        parser.error("--layouts must not contain duplicates")
    if any(layout.workers > arguments.cpus for layout in arguments.layouts):
        parser.error("each layout's ranks x threads must be no greater than --cpus")
    try:
        topology = cpu_topology.detect_cpu_topology()
        configure_placement(arguments, topology)
        runs = make_runs(arguments, topology)
        specs = [run.spec for run in runs]
        output = arguments.output_dir.expanduser().resolve()
        summary_path = output / "scaling.txt"
        topology_path = output / "cpu-topology.json"
        plot_path = output / "scaling.png"
        calculations = len(runs) * 2 * (arguments.runs + arguments.warmup)
        description = (
            f"Driver CPU: {topology.model}\n"
            f"Allocated physical cores: {arguments.cpus}\n"
            f"Available physical cores: {topology.physical_cores} ({topology.source})\n"
            "Hybrid layouts (ranks x threads): "
            f"{', '.join(map(str, arguments.layouts))}\n"
            f"Serial launcher: {shlex.join(arguments.baseline_launcher or ())}\n"
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
                        *(f"{spec.label}:{spec.workers}:{spec.path}" for spec in specs),
                    )
                )
            )
            return 0
        paclitaxel_scaling.preflight(
            [
                paclitaxel_scaling.ScalingRun(
                    run.layout.workers, run.report, run.arguments
                )
                for run in runs
            ],
            summary_path,
            topology_path,
            plot_path,
        )
        figure = paclitaxel_scaling.create_scaling_figure()
        output.mkdir(parents=True, exist_ok=True)
        topology_path.write_text(
            json.dumps(topology_report(arguments, topology, runs), indent=2) + "\n",
            encoding="utf-8",
        )
        print(description, end="")
        print(
            f"Hybrid scaling: {calculations} calculations including warm-ups.",
            flush=True,
        )
        for run in runs:
            print(f"\nHybrid layout: {run.layout}", flush=True)
            status = benchmark.main(run.arguments)
            if status != 0:
                print(
                    f"Hybrid scaling stopped at {run.layout}; completed JSON reports "
                    "were retained, and no scaling summary was written.",
                    file=sys.stderr,
                )
                return status
        results = check_parallel_scaling.evaluate_specs(specs)
        summary = io.StringIO()
        summary.write(description + "\n")
        with contextlib.redirect_stdout(summary):
            check_parallel_scaling.print_report(results)
        save_scaling_plot(results, runs, plot_path, figure)
        summary_path.write_text(summary.getvalue(), encoding="utf-8")
        print(f"\n{summary.getvalue()}", end="")
        print(f"\nScaling summary: {summary_path}")
        print(f"Scaling plot: {plot_path}")
        print(f"Scaling plot: {plot_path.with_suffix('.pdf')}")
    except (
        benchmark.BenchmarkError,
        check_parallel_scaling.ScalingError,
        ValueError,
        OSError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print(
            "Hybrid scaling interrupted; completed JSON reports were retained.",
            file=sys.stderr,
        )
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
