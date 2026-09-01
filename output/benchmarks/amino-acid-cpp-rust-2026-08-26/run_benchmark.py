#!/usr/bin/env python3
"""Benchmark shipped, current C++, and DaymudeLab Rust assembly solvers.

The harness measures end-to-end CLI wall time (spawn through process exit),
validates every assembly index, rotates variant order, and retains raw samples.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE = (
    ROOT.parent / "assemblytheorytools/assemblytheorytools/"
    "precompiled/asscpp_combined_static_linux"
)
DEFAULT_CURRENT_SERIAL = ROOT / "build/parallel/AssemblyCpp"
DEFAULT_CURRENT_OMP = ROOT / "build/parallel/AssemblyCppOMP"
DEFAULT_RUST_REPOSITORY = ROOT / "build/daymudelab-assembly-theory"
DEFAULT_INPUT_DIRECTORY = ROOT / "benchmarks/inputs/scaling"
CPP_RESULT_RE = re.compile(r"has assembly index:\s*(-?\d+)")
RUST_RESULT_RE = re.compile(r"^[0-9]+$")

CASES = (
    ("02c", "amino_acid_scaling_02c_018a.mol", 2, 18, 16, 10),
    ("03c", "amino_acid_scaling_03c_027a.mol", 3, 27, 24, 13),
    ("04c", "amino_acid_scaling_04c_036a.mol", 4, 36, 32, 15),
    ("05c", "amino_acid_scaling_05c_043a.mol", 5, 43, 38, 17),
    ("06c", "amino_acid_scaling_06c_053a.mol", 6, 53, 47, 19),
    ("07c", "amino_acid_scaling_07c_063a.mol", 7, 63, 56, 21),
    ("08c", "amino_acid_scaling_08c_068a.mol", 8, 68, 60, 23),
    ("09c", "amino_acid_scaling_09c_077a.mol", 9, 77, 68, 25),
    ("10c", "amino_acid_scaling_10c_086a.mol", 10, 86, 76, 27),
    ("11c", "amino_acid_scaling_11c_096a.mol", 11, 96, 85, 29),
    ("12c", "amino_acid_scaling_12c_105a.mol", 12, 105, 93, 32),
    ("13c", "amino_acid_scaling_13c_113a.mol", 13, 113, 100, 34),
)


def cpu_mask(workers: int) -> str:
    """Return one hardware thread from each of the requested P-cores."""
    if workers not in (1, 2, 4, 8):
        raise ValueError(f"unsupported worker count: {workers}")
    return ",".join(str(cpu) for cpu in range(0, workers * 2, 2))


def variant(
    label: str,
    implementation: str,
    executable: Path,
    workers: int,
    mode: str,
) -> dict[str, Any]:
    return {
        "label": label,
        "implementation": implementation,
        "executable": executable,
        "workers": workers,
        "cpus": cpu_mask(workers),
        "mode": mode,
    }


VARIANTS = {
    "reference": variant(
        "Shipped reference (serial)",
        "shipped-reference",
        DEFAULT_REFERENCE,
        1,
        "reference",
    ),
    "current_serial": variant(
        "Current C++ serial",
        "current-cpp",
        DEFAULT_CURRENT_SERIAL,
        1,
        "current-serial",
    ),
    "rust_serial": variant(
        "Rust serial",
        "daymudelab-rust",
        DEFAULT_RUST_REPOSITORY / "target/release/assembly-theory",
        1,
        "rust-serial",
    ),
    "current_omp1": variant(
        "Current C++ OpenMP, 1 worker",
        "current-cpp",
        DEFAULT_CURRENT_OMP,
        1,
        "current-openmp",
    ),
    "rust_parallel1": variant(
        "Rust depth-one, 1 worker",
        "daymudelab-rust",
        DEFAULT_RUST_REPOSITORY / "target/release/assembly-theory",
        1,
        "rust-depth-one",
    ),
    "current_omp2": variant(
        "Current C++ OpenMP, 2 workers",
        "current-cpp",
        DEFAULT_CURRENT_OMP,
        2,
        "current-openmp",
    ),
    "rust_parallel2": variant(
        "Rust depth-one, 2 workers",
        "daymudelab-rust",
        DEFAULT_RUST_REPOSITORY / "target/release/assembly-theory",
        2,
        "rust-depth-one",
    ),
    "current_omp4": variant(
        "Current C++ OpenMP, 4 workers",
        "current-cpp",
        DEFAULT_CURRENT_OMP,
        4,
        "current-openmp",
    ),
    "rust_parallel4": variant(
        "Rust depth-one, 4 workers",
        "daymudelab-rust",
        DEFAULT_RUST_REPOSITORY / "target/release/assembly-theory",
        4,
        "rust-depth-one",
    ),
    "current_omp8": variant(
        "Current C++ OpenMP, 8 workers",
        "current-cpp",
        DEFAULT_CURRENT_OMP,
        8,
        "current-openmp",
    ),
    "rust_parallel8": variant(
        "Rust depth-one, 8 workers",
        "daymudelab-rust",
        DEFAULT_RUST_REPOSITORY / "target/release/assembly-theory",
        8,
        "rust-depth-one",
    ),
}


def fingerprint(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(repository: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=repository,
        capture_output=True,
        text=True,
        check=True,
    )
    return completed.stdout.strip()


def controlled_environment(specification: dict[str, Any]) -> dict[str, str]:
    environment = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OMP_DYNAMIC",
        "OMP_PLACES",
        "OMP_PROC_BIND",
        "ASSEMBLYCPP_PARALLEL_MIN_BONDS",
        "RAYON_NUM_THREADS",
    ):
        environment.pop(key, None)
    environment["LC_ALL"] = "C"

    mode = specification["mode"]
    workers = str(specification["workers"])
    if mode == "current-openmp":
        environment.update(
            {
                "OMP_NUM_THREADS": workers,
                "OMP_THREAD_LIMIT": workers,
                "OMP_DYNAMIC": "FALSE",
                "OMP_PLACES": "cores",
                "OMP_PROC_BIND": "close",
                "ASSEMBLYCPP_PARALLEL_MIN_BONDS": "0",
            }
        )
    elif mode in ("rust-serial", "rust-depth-one"):
        environment["RAYON_NUM_THREADS"] = workers
    return environment


def explicit_environment(specification: dict[str, Any]) -> dict[str, str]:
    environment = controlled_environment(specification)
    keys = (
        "LC_ALL",
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OMP_DYNAMIC",
        "OMP_PLACES",
        "OMP_PROC_BIND",
        "ASSEMBLYCPP_PARALLEL_MIN_BONDS",
        "RAYON_NUM_THREADS",
    )
    return {key: environment[key] for key in keys if key in environment}


def command_for(specification: dict[str, Any], input_name: str) -> list[str]:
    executable = str(specification["executable"])
    command = ["taskset", "-c", specification["cpus"], executable]
    mode = specification["mode"]
    if mode == "reference":
        return [
            *command,
            Path(input_name).stem,
            "-pathway=0",
            "-removeHydrogens=1",
            "-compensateDisjoint=0",
        ]
    if mode in ("current-serial", "current-openmp"):
        return [
            *command,
            "--pathway=0",
            "--remove-hydrogens=1",
            "--compensate-disjoint=0",
            "--memory-report=0",
            "--write-intermediate-mas=0",
            "--",
            input_name,
        ]
    if mode == "rust-serial":
        return [*command, "--parallel", "none", input_name]
    if mode == "rust-depth-one":
        return [*command, "--parallel", "depth-one", input_name]
    raise ValueError(f"unknown mode: {mode}")


def validate_result(
    specification: dict[str, Any], directory: Path, stdout: str, expected: int
) -> None:
    mode = specification["mode"]
    if mode.startswith("rust-"):
        rendered = stdout.strip()
        if not RUST_RESULT_RE.fullmatch(rendered):
            raise RuntimeError(f"invalid Rust output: {rendered!r}")
        actual = int(rendered)
    else:
        output_path = directory / "inputOut"
        if not output_path.is_file():
            raise RuntimeError(f"missing C++ output: {output_path}")
        match = CPP_RESULT_RE.search(output_path.read_text(encoding="utf-8"))
        if match is None:
            raise RuntimeError("C++ output did not report an assembly index")
        actual = int(match.group(1))
    if actual != expected:
        raise RuntimeError(f"assembly index {actual}, expected {expected}")


def run_once(name: str, source: Path, expected: int, timeout: float) -> float:
    specification = VARIANTS[name]
    with tempfile.TemporaryDirectory(prefix="cpp-rust-amino-") as raw_directory:
        directory = Path(raw_directory)
        input_path = directory / "input.mol"
        input_path.write_bytes(source.read_bytes())
        command = command_for(specification, input_path.name)

        started_ns = time.perf_counter_ns()
        completed = subprocess.run(
            command,
            cwd=directory,
            env=controlled_environment(specification),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        elapsed_seconds = (time.perf_counter_ns() - started_ns) / 1_000_000_000
        if completed.returncode != 0:
            raise RuntimeError(
                f"{name} exited {completed.returncode}: {completed.stderr.strip()}"
            )
        validate_result(specification, directory, completed.stdout, expected)
        return elapsed_seconds


def summarize(values: list[float]) -> dict[str, Any]:
    median = statistics.median(values)
    ordered = sorted(values)
    p95_index = max(0, min(len(ordered) - 1, round(0.95 * len(ordered) + 0.5) - 1))
    return {
        "samples": values,
        "minimum": min(values),
        "median": median,
        "mean": statistics.mean(values),
        "mad": statistics.median(abs(value - median) for value in values),
        "p95": ordered[p95_index],
        "maximum": max(values),
    }


def write_json_atomic(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def report_metadata(
    selected_variants: list[str],
    runs: int,
    warmup: int,
    first: int,
    last: int,
    rust_repository: Path,
) -> dict[str, Any]:
    variants = {}
    for name in selected_variants:
        specification = VARIANTS[name]
        executable = Path(specification["executable"])
        variants[name] = {
            "label": specification["label"],
            "implementation": specification["implementation"],
            "mode": specification["mode"],
            "workers": specification["workers"],
            "cpus": specification["cpus"],
            "executable": str(executable),
            "sha256": fingerprint(executable),
            "environment": explicit_environment(specification),
        }

    return {
        "schema_version": 1,
        "status": "running",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "metric": "external CLI wall time in seconds via time.perf_counter_ns",
        "corpus": "maintained amino-acid V2000 MOL scaling inputs",
        "case_range": {"first_components": first, "last_components": last},
        "runs": runs,
        "warmup": warmup,
        "schedule": "cyclically rotated variant order per case and measured round",
        "validation": "zero exit status and exact expected assembly index on every run",
        "platform": {
            "description": platform.platform(),
            "machine": platform.machine(),
            "cpu_model": next(
                (
                    line.split(":", 1)[1].strip()
                    for line in Path("/proc/cpuinfo").read_text().splitlines()
                    if line.startswith("model name")
                ),
                "unknown",
            ),
        },
        "source": {
            "current_cpp_commit": git_value(ROOT, "rev-parse", "HEAD"),
            "rust_repository": "https://github.com/DaymudeLab/assembly-theory",
            "rust_commit": git_value(rust_repository, "rev-parse", "HEAD"),
            "rust_describe": git_value(
                rust_repository, "describe", "--always", "--dirty"
            ),
            "rust_build": "cargo build --locked --release",
            "rust_target_tuning": "Cargo release default (no RUSTFLAGS target-cpu override)",
            "current_cpp_build": "CMake parallel Release preset, x86-64-v3, no LTO",
        },
        "variants": variants,
        "cases": [],
    }


def write_csv(path: Path, report: dict[str, Any]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "case",
                "components",
                "atoms",
                "bonds",
                "expected_assembly_index",
                "variant",
                "label",
                "implementation",
                "mode",
                "workers",
                "median_wall_seconds",
                "mad_wall_seconds",
                "minimum_wall_seconds",
                "p95_wall_seconds",
                "maximum_wall_seconds",
            ]
        )
        for case in report["cases"]:
            for name, timing in case["timings"].items():
                specification = report["variants"][name]
                writer.writerow(
                    [
                        case["id"],
                        case["components"],
                        case["atoms"],
                        case["bonds"],
                        case["expected_assembly_index"],
                        name,
                        specification["label"],
                        specification["implementation"],
                        specification["mode"],
                        specification["workers"],
                        timing["median"],
                        timing["mad"],
                        timing["minimum"],
                        timing["p95"],
                        timing["maximum"],
                    ]
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=2)
    parser.add_argument("--last", type=int, default=10)
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--current-serial", type=Path, default=DEFAULT_CURRENT_SERIAL)
    parser.add_argument("--current-openmp", type=Path, default=DEFAULT_CURRENT_OMP)
    parser.add_argument("--rust-repository", type=Path, default=DEFAULT_RUST_REPOSITORY)
    parser.add_argument(
        "--rust-binary",
        type=Path,
        help="Rust executable (defaults to target/release under --rust-repository)",
    )
    parser.add_argument("--input-directory", type=Path, default=DEFAULT_INPUT_DIRECTORY)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=tuple(VARIANTS),
        default=list(VARIANTS),
    )
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("benchmark.json")
    )
    arguments = parser.parse_args()
    if arguments.runs < 1 or arguments.warmup < 0:
        raise SystemExit("runs must be >= 1 and warmup must be >= 0")

    selected_cases = [
        case for case in CASES if arguments.first <= case[2] <= arguments.last
    ]
    if not selected_cases:
        raise SystemExit("no cases selected")

    rust_repository = arguments.rust_repository.resolve()
    rust_binary = (
        arguments.rust_binary.resolve()
        if arguments.rust_binary is not None
        else rust_repository / "target/release/assembly-theory"
    )
    VARIANTS["reference"]["executable"] = arguments.reference.resolve()
    VARIANTS["current_serial"]["executable"] = arguments.current_serial.resolve()
    for name in ("current_omp1", "current_omp2", "current_omp4", "current_omp8"):
        VARIANTS[name]["executable"] = arguments.current_openmp.resolve()
    for name in (
        "rust_serial",
        "rust_parallel1",
        "rust_parallel2",
        "rust_parallel4",
        "rust_parallel8",
    ):
        VARIANTS[name]["executable"] = rust_binary

    for name in arguments.variants:
        executable = Path(VARIANTS[name]["executable"])
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise SystemExit(f"missing executable for {name}: {executable}")

    report = report_metadata(
        arguments.variants,
        arguments.runs,
        arguments.warmup,
        arguments.first,
        arguments.last,
        rust_repository,
    )
    samples: dict[str, dict[str, list[float]]] = {
        case[0]: {name: [] for name in arguments.variants} for case in selected_cases
    }
    variant_order = list(arguments.variants)

    for case_index, (
        case_id,
        filename,
        components,
        atoms,
        bonds,
        expected,
    ) in enumerate(selected_cases):
        source = arguments.input_directory.resolve() / filename
        print(
            f"[{case_index + 1}/{len(selected_cases)}] {case_id}: "
            f"{components} components, {atoms} atoms, {bonds} bonds",
            flush=True,
        )
        for warmup_index in range(arguments.warmup):
            offset = (case_index + warmup_index) % len(variant_order)
            order = variant_order[offset:] + variant_order[:offset]
            for name in order:
                value = run_once(name, source, expected, arguments.timeout)
                print(f"  warm-up {name}: {value:.6f} s", flush=True)

        for run_index in range(arguments.runs):
            offset = (case_index + run_index) % len(variant_order)
            order = variant_order[offset:] + variant_order[:offset]
            for name in order:
                value = run_once(name, source, expected, arguments.timeout)
                samples[case_id][name].append(value)
                print(
                    f"  run {run_index + 1}/{arguments.runs} {name}: {value:.6f} s",
                    flush=True,
                )

        report["cases"].append(
            {
                "id": case_id,
                "input": str(source),
                "input_sha256": fingerprint(source),
                "components": components,
                "atoms": atoms,
                "bonds": bonds,
                "expected_assembly_index": expected,
                "timings": {
                    name: summarize(samples[case_id][name])
                    for name in arguments.variants
                },
            }
        )
        write_json_atomic(arguments.output, report)

    report["status"] = "complete"
    report["completed_at"] = datetime.now(timezone.utc).isoformat()
    write_json_atomic(arguments.output, report)
    csv_path = arguments.output.with_name("summary.csv")
    write_csv(csv_path, report)
    print(f"Wrote {arguments.output}", flush=True)
    print(f"Wrote {csv_path}", flush=True)


if __name__ == "__main__":
    main()
