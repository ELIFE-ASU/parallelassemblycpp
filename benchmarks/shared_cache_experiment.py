"""Compare shared-cache policies with paired timings and separate memory profiles."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

if __package__:
    from . import benchmark
else:
    import benchmark


POLICY_ENV = "PARALLELASSEMBLYCPP_SHARED_CACHE_POLICY"
BYTES_ENV = "PARALLELASSEMBLYCPP_SHARED_CACHE_BYTES"
POLICIES = ("shared", "local", "selective")
DEFAULT_SELECTIVE_BYTES = 256 * 1024 * 1024


@dataclass(frozen=True)
class Variant:
    name: str
    policy: str
    byte_limit: int = 0


def variant_specification(value: str) -> Variant:
    name, separator, settings = value.partition("=")
    policy, limit_separator, limit = settings.partition(":")
    if (
        not separator
        or not benchmark.CASE_NAME_PATTERN.fullmatch(name)
        or name in {"baseline", "profiles"}
        or policy not in POLICIES
        or (limit_separator and not limit)
    ):
        raise argparse.ArgumentTypeError(
            "expected NAME=POLICY[:BYTES]; policies: shared, local, selective; "
            "baseline and profiles are reserved names"
        )
    try:
        byte_limit = benchmark.non_negative_int(limit) if limit else 0
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "BYTES must be a nonnegative integer"
        ) from error
    return Variant(name, policy, byte_limit)


def variant_environment(value: str) -> tuple[str, tuple[str, str]]:
    name, separator, assignment = value.partition(":")
    if not separator:
        raise argparse.ArgumentTypeError("expected NAME:KEY=VALUE")
    return name, benchmark.environment_assignment(assignment)


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-executable", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--baseline-telemetry-executable", type=Path)
    parser.add_argument("--telemetry-executable", type=Path)
    parser.add_argument("--manifest", type=Path, default=benchmark.DEFAULT_MANIFEST)
    parser.add_argument("--suite", choices=benchmark.KNOWN_SUITES)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--threads", type=benchmark.positive_int, default=4)
    parser.add_argument("--parallel", choices=("auto", "on"), default="auto")
    parser.add_argument("--launcher", type=benchmark.launcher_prefix)
    parser.add_argument("--runs", type=benchmark.positive_int, default=6)
    parser.add_argument("--warmup", type=benchmark.non_negative_int, default=1)
    parser.add_argument("--timeout", type=benchmark.positive_float, default=600.0)
    parser.add_argument("--variant", type=variant_specification, action="append")
    parser.add_argument(
        "--env", type=benchmark.environment_assignment, action="append", default=[]
    )
    parser.add_argument(
        "--baseline-env",
        type=benchmark.environment_assignment,
        action="append",
        default=[],
    )
    parser.add_argument(
        "--variant-env", type=variant_environment, action="append", default=[]
    )
    parser.add_argument("--time-executable", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--require-all-faster", action="store_true")
    return parser


def execution_configs(
    arguments: argparse.Namespace, variants: Sequence[Variant]
) -> dict[str, benchmark.ExecutionConfig]:
    if arguments.threads < 2:
        raise benchmark.BenchmarkError(
            "shared-cache experiments need at least two threads"
        )
    names = [variant.name for variant in variants]
    if len(names) != len(set(names)):
        raise benchmark.BenchmarkError("variant names must be distinct")
    overrides: dict[str, list[tuple[str, str]]] = {name: [] for name in names}
    for name, setting in arguments.variant_env:
        if name not in overrides:
            raise benchmark.BenchmarkError(
                f"unknown variant environment target: {name}"
            )
        overrides[name].append(setting)
    defaults = [
        ("OMP_NUM_THREADS", str(arguments.threads)),
        ("OMP_THREAD_LIMIT", str(arguments.threads)),
        ("OMP_DYNAMIC", "FALSE"),
        ("OMP_PROC_BIND", "close"),
        ("OMP_PLACES", "cores"),
        (POLICY_ENV, "shared"),
        (BYTES_ENV, "0"),
    ]
    configs = {}
    for name, extra in [
        ("baseline", arguments.baseline_env),
        *[
            (
                variant.name,
                [
                    (POLICY_ENV, variant.policy),
                    (BYTES_ENV, str(variant.byte_limit)),
                    *overrides[variant.name],
                ],
            )
            for variant in variants
        ],
    ]:
        settings = dict(defaults)
        for layer in (arguments.env, extra):
            keys = [key for key, _ in layer]
            if len(keys) != len(set(keys)):
                raise benchmark.BenchmarkError(f"duplicate {name} environment setting")
            settings.update(layer)
        config = benchmark.create_execution_config(
            arguments.launcher,
            list(settings.items()),
            name,
            parallel_mode=arguments.parallel,
        )
        if config is None:
            raise benchmark.BenchmarkError(
                f"missing execution configuration for {name}"
            )
        configs[name] = config
    return configs


def collect_profiles(
    executable: Path,
    cases: Sequence[benchmark.BenchmarkCase],
    execution: benchmark.ExecutionConfig,
    timeout: float,
    time_executable: Path,
    *,
    telemetry: bool,
) -> dict[str, dict[str, object]]:
    """Collect process peak RSS; OpenMP threads share that address space."""
    profiles = {}
    with tempfile.TemporaryDirectory(
        prefix="parallelassemblycpp-cache-profile-"
    ) as directory:
        for prepared in benchmark.prepare_cases(cases, Path(directory)):
            print(
                f"Untimed memory/telemetry profile [{prepared.case.name}]...",
                flush=True,
            )
            rss_path = prepared.working_directory / "peak-rss-kib.txt"
            profiled_execution = replace(
                execution,
                launcher=(
                    str(time_executable),
                    "-f",
                    "%M",
                    "-o",
                    str(rss_path),
                    *execution.launcher,
                ),
            )
            telemetry_data = None
            if telemetry:
                telemetry_data = benchmark.run_telemetry_once(
                    executable, prepared, timeout, execution=profiled_execution
                )
            else:
                benchmark.run_once(
                    executable, prepared, timeout, execution=profiled_execution
                )
            try:
                peak_rss = int(rss_path.read_text(encoding="utf-8").strip())
                if peak_rss <= 0:
                    raise ValueError("peak RSS must be positive")
            except (OSError, ValueError) as error:
                raise benchmark.BenchmarkError(
                    f"{prepared.case.name}: invalid GNU time peak RSS: {error}"
                ) from error
            profiles[prepared.case.name] = {
                "peak_rss_kib": peak_rss,
                "telemetry": telemetry_data,
            }
    return profiles


def shared_cache_fields(profile: dict[str, object]) -> dict[str, object]:
    telemetry = profile.get("telemetry")
    if not isinstance(telemetry, dict):
        return {}
    parallel = telemetry.get("parallel")
    if not isinstance(parallel, dict):
        return {}
    aggregate = parallel.get("aggregate")
    if not isinstance(aggregate, dict):
        return {}
    cache = aggregate.get("shared_assembly_cache")
    if not isinstance(cache, dict):
        return {}
    return {
        name: value
        for name, value in cache.items()
        if not isinstance(value, (dict, list))
    }


def summary_rows(
    name: str,
    results: Sequence[benchmark.CaseResult],
    profiles: dict[str, dict[str, object]],
    baseline_profiles: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows = []
    for result in results:
        speedup = benchmark.paired_speedup_summary(
            result.measurements, result.baseline_measurements, "wall_seconds"
        )
        clock_speedup = benchmark.paired_speedup_summary(
            result.measurements, result.baseline_measurements, "clock_ticks"
        )
        row: dict[str, object] = {
            "variant": name,
            "case": result.case.name,
            "expected_assembly_index": result.case.expected_assembly_index,
            "wall_seconds_median": benchmark.summarize(
                [measurement.wall_seconds for measurement in result.measurements]
            ).median,
            "baseline_wall_seconds_median": benchmark.summarize(
                [
                    measurement.wall_seconds
                    for measurement in result.baseline_measurements
                ]
            ).median,
            "paired_wall_speedup_median": None if speedup is None else speedup.median,
            "paired_clock_speedup_median": (
                None if clock_speedup is None else clock_speedup.median
            ),
        }
        for role, profile in (
            ("candidate", profiles[result.case.name]),
            ("baseline", baseline_profiles[result.case.name]),
        ):
            row[f"{role}_peak_rss_kib"] = profile["peak_rss_kib"]
            cache = shared_cache_fields(profile)
            for field, value in cache.items():
                row[f"{role}_cache_{field}"] = value
            hits, misses = cache.get("hits"), cache.get("misses")
            if isinstance(hits, int) and isinstance(misses, int):
                row[f"{role}_cache_hit_rate"] = (
                    hits / (hits + misses) if hits + misses else None
                )
        rows.append(row)
    return rows


def write_summary(path: Path, rows: Sequence[dict[str, object]]) -> None:
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)
    try:
        variants = arguments.variant or [
            Variant("shared", "shared"),
            Variant("local", "local"),
            Variant("bounded", "shared", DEFAULT_SELECTIVE_BYTES),
            Variant("selective", "selective", DEFAULT_SELECTIVE_BYTES),
        ]
        manifest, available_cases = benchmark.load_manifest(arguments.manifest)
        cases = benchmark.select_cases(available_cases, arguments.suite, arguments.case)
        configs = execution_configs(arguments, variants)
        if arguments.runs % 2:
            raise benchmark.BenchmarkError(
                "use an even --runs count for balanced AB/BA pairs"
            )
        output = arguments.output_dir.expanduser().resolve()
        print(f"Cases: {len(cases)} unique manifest cases; variants: {len(variants)}")
        count = len(cases) * (
            2 * len(variants) * (arguments.runs + arguments.warmup) + len(variants) + 1
        )
        print(f"Calculations: {count}; profiles are excluded from timing aggregates")
        if arguments.dry_run:
            print(
                json.dumps(
                    {
                        "cases": [case.name for case in cases],
                        "execution": {
                            name: benchmark.execution_config_metadata(config)
                            for name, config in configs.items()
                        },
                        "output_directory": str(output),
                    },
                    indent=2,
                )
            )
            return 0
        if output.exists():
            raise benchmark.BenchmarkError(f"output directory already exists: {output}")
        executable = benchmark.resolve_executable(arguments.executable)
        baseline = benchmark.resolve_executable(arguments.baseline_executable)
        time_executable = benchmark.resolve_executable(arguments.time_executable)
        candidate_profile = (
            benchmark.resolve_executable(arguments.telemetry_executable)
            if arguments.telemetry_executable
            else executable
        )
        baseline_profile = (
            benchmark.resolve_executable(arguments.baseline_telemetry_executable)
            if arguments.baseline_telemetry_executable
            else baseline
        )
        if arguments.telemetry_executable:
            benchmark.ensure_distinct_executables(candidate_profile, executable)
        if arguments.baseline_telemetry_executable:
            benchmark.ensure_distinct_executables(baseline_profile, baseline)
        fingerprints = {
            path: benchmark.executable_metadata(path)
            for path in (executable, baseline, candidate_profile, baseline_profile)
        }
        for variant in variants:
            benchmark.ensure_distinct_execution_identities(
                fingerprints[executable],
                fingerprints[baseline],
                configs[variant.name],
                configs["baseline"],
            )
        corpus = benchmark.benchmark_corpus_metadata(manifest, cases)
        output.mkdir(parents=True)
        baseline_profiles = collect_profiles(
            baseline_profile,
            cases,
            configs["baseline"],
            arguments.timeout,
            time_executable,
            telemetry=arguments.baseline_telemetry_executable is not None,
        )
        profile_report: dict[str, object] = {
            "complete": False,
            "collection": "separate untimed validated run per case and variant",
            "peak_rss_method": (
                "GNU time %M; KiB; instrumented when telemetry executable supplied"
            ),
            "executables": {
                str(path): metadata for path, metadata in fingerprints.items()
            },
            "corpus": corpus,
            "execution": {
                name: benchmark.execution_config_metadata(config)
                for name, config in configs.items()
            },
            "baseline": baseline_profiles,
            "variants": {},
        }
        profiles_by_variant: dict[str, object] = {}
        profile_report["variants"] = profiles_by_variant
        rows = []
        for variant in variants:
            print(f"Comparing {variant.name} against baseline...", flush=True)
            results = benchmark.run_benchmarks(
                executable,
                baseline,
                cases,
                arguments.runs,
                arguments.warmup,
                arguments.timeout,
                candidate_execution=configs[variant.name],
                baseline_execution=configs["baseline"],
            )
            profiles = collect_profiles(
                candidate_profile,
                cases,
                configs[variant.name],
                arguments.timeout,
                time_executable,
                telemetry=arguments.telemetry_executable is not None,
            )
            for path, metadata in fingerprints.items():
                benchmark.verify_executable_unchanged(path, metadata, str(path))
            if benchmark.benchmark_corpus_metadata(manifest, cases) != corpus:
                raise benchmark.BenchmarkError(
                    "benchmark corpus changed during experiment"
                )
            benchmark.write_json_report(
                output / f"{variant.name}.json",
                fingerprints[executable],
                fingerprints[baseline],
                None,
                corpus,
                manifest,
                arguments.suite,
                arguments.runs,
                arguments.warmup,
                arguments.timeout,
                results,
                configs[variant.name],
                configs["baseline"],
            )
            profiles_by_variant[variant.name] = profiles
            (output / "profiles.json").write_text(
                json.dumps(profile_report, indent=2) + "\n", encoding="utf-8"
            )
            rows.extend(
                summary_rows(variant.name, results, profiles, baseline_profiles)
            )
            write_summary(output / "summary.csv", rows)
            benchmark.print_comparison_summary(results)
        regressions = [
            row
            for row in rows
            if any(
                not isinstance(row[metric], (int, float)) or row[metric] <= 1.0
                for metric in (
                    "paired_wall_speedup_median",
                    "paired_clock_speedup_median",
                )
            )
        ]
        profile_report["complete"] = True
        (output / "profiles.json").write_text(
            json.dumps(profile_report, indent=2) + "\n", encoding="utf-8"
        )
        print(
            "Case/variant comparisons at or below 1.0 in wall or clock: "
            f"{len(regressions)}/{len(rows)}"
        )
        for row in regressions:
            print(
                f"  {row['variant']}/{row['case']}: "
                f"wall={row['paired_wall_speedup_median']}, "
                f"clock={row['paired_clock_speedup_median']}"
            )
        print(f"Summary: {output / 'summary.csv'}")
        return 1 if arguments.require_all_faster and regressions else 0
    except (benchmark.BenchmarkError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("error: experiment interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
