"""Generate GCC PGO profiles from a weighted AssemblyCpp corpus."""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

if __package__:
    from . import benchmark
else:
    import benchmark


DEFAULT_WEIGHTS = benchmark.BENCHMARK_DIRECTORY / "pgo-training.tsv"
WEIGHTS_HEADER = ("name", "repetitions")
DIRECTORY_MARKER_FILENAME = ".assemblycpp-pgo-profile-directory"
COMPLETION_FILENAME = "assemblycpp-pgo-complete.json"
COMPLETION_SCHEMA_VERSION = 2


@dataclass(frozen=True)
class WeightedCase:
    case: benchmark.BenchmarkCase
    repetitions: int


def load_training_weights(
    path: Path,
    cases: Sequence[benchmark.BenchmarkCase],
) -> tuple[Path, list[WeightedCase]]:
    """Load a strict, complete weight mapping for the benchmark corpus."""
    weights_path = benchmark.resolve_file(path, "PGO training weights")
    cases_by_name = {case.name: case for case in cases}
    weighted_cases: list[WeightedCase] = []
    seen_names: set[str] = set()

    try:
        with weights_path.open(encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != WEIGHTS_HEADER:
                raise benchmark.BenchmarkError(
                    f"invalid header in {weights_path}: expected {WEIGHTS_HEADER}, "
                    f"got {header}"
                )

            for line_number, row in enumerate(reader, start=2):
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != len(WEIGHTS_HEADER):
                    raise benchmark.BenchmarkError(
                        f"invalid row in {weights_path}:{line_number}: expected "
                        f"{len(WEIGHTS_HEADER)} columns"
                    )

                name, repetitions_text = (value.strip() for value in row)
                if not benchmark.CASE_NAME_PATTERN.fullmatch(name):
                    raise benchmark.BenchmarkError(
                        f"invalid benchmark name in {weights_path}:{line_number}: "
                        f"{name!r}"
                    )
                if name in seen_names:
                    raise benchmark.BenchmarkError(
                        f"duplicate benchmark name in {weights_path}:{line_number}: "
                        f"{name!r}"
                    )
                seen_names.add(name)

                try:
                    repetitions = int(repetitions_text)
                except ValueError as error:
                    raise benchmark.BenchmarkError(
                        f"invalid repetitions in {weights_path}:{line_number}: "
                        f"{repetitions_text!r}"
                    ) from error
                if repetitions < 1:
                    raise benchmark.BenchmarkError(
                        f"invalid repetitions in {weights_path}:{line_number}: "
                        "must be a positive integer"
                    )

                benchmark_case = cases_by_name.get(name)
                if benchmark_case is None:
                    raise benchmark.BenchmarkError(
                        f"unknown benchmark name in {weights_path}:{line_number}: "
                        f"{name!r}"
                    )
                weighted_cases.append(WeightedCase(benchmark_case, repetitions))
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not read PGO training weights {weights_path}: {error}"
        ) from error

    missing_names = [case.name for case in cases if case.name not in seen_names]
    if missing_names:
        raise benchmark.BenchmarkError(
            f"PGO training weights do not cover benchmark case(s): "
            f"{', '.join(missing_names)}"
        )
    if not weighted_cases:
        raise benchmark.BenchmarkError(
            f"PGO training weights are empty: {weights_path}"
        )
    return weights_path, weighted_cases


def prepare_profile_directory(path: Path) -> tuple[Path, int]:
    """Prepare a marked profile directory and remove stale training outputs."""
    profile_directory = path.expanduser().resolve()
    protected_directories = {
        Path(profile_directory.anchor),
        Path.home().resolve(),
        benchmark.REPOSITORY_ROOT.resolve(),
        Path.cwd().resolve(),
        Path(tempfile.gettempdir()).resolve(),
    }
    if profile_directory in protected_directories:
        raise benchmark.BenchmarkError(
            f"protected directory cannot be used as the PGO profile directory: "
            f"{profile_directory}"
        )
    try:
        profile_directory.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not create PGO profile directory {profile_directory}: {error}"
        ) from error
    if not profile_directory.is_dir():
        raise benchmark.BenchmarkError(
            f"PGO profile path is not a directory: {profile_directory}"
        )

    removed_profile_count = 0
    try:
        directory_marker = profile_directory / DIRECTORY_MARKER_FILENAME
        if directory_marker.is_symlink() or directory_marker.is_dir():
            raise benchmark.BenchmarkError(
                f"invalid PGO profile directory marker: {directory_marker}"
            )
        if not directory_marker.is_file():
            if any(profile_directory.iterdir()):
                raise benchmark.BenchmarkError(
                    "refusing to adopt non-empty PGO profile directory without "
                    f"its AssemblyCpp marker: {profile_directory}"
                )
            directory_marker.write_text(
                "Dedicated AssemblyCpp PGO profile directory.\n",
                encoding="utf-8",
            )
        completion_file = profile_directory / COMPLETION_FILENAME
        if completion_file.is_file() or completion_file.is_symlink():
            completion_file.unlink()
        elif completion_file.exists():
            raise benchmark.BenchmarkError(
                f"PGO completion path is not a file: {completion_file}"
            )
        for profile_path in sorted(profile_directory.rglob("*.gcda")):
            if profile_path.is_file() or profile_path.is_symlink():
                profile_path.unlink()
                removed_profile_count += 1
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not clear GCC profile data in {profile_directory}: {error}"
        ) from error
    return profile_directory, removed_profile_count


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of a regular file."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not fingerprint {path}: {error}"
        ) from error
    return digest.hexdigest()


def repository_path(path: Path) -> str:
    """Return a stable repository-relative path when possible."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(benchmark.REPOSITORY_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def corpus_metadata(
    manifest_path: Path,
    cases: Sequence[benchmark.BenchmarkCase],
) -> dict[str, object]:
    """Fingerprint the manifest and every input used by the training corpus."""
    return {
        "manifest": {
            "path": repository_path(manifest_path),
            "sha256": file_sha256(manifest_path),
        },
        "inputs": [
            {
                "name": case.name,
                "path": repository_path(case.source),
                "sha256": file_sha256(case.source),
            }
            for case in cases
        ],
    }


def write_completion_record(
    profile_directory: Path,
    weights_path: Path,
    weights_sha256: str,
    training_corpus: dict[str, object],
    executable_metadata: dict[str, object],
    completed_repetitions: int,
    profile_files: Sequence[Path],
) -> Path:
    """Atomically mark a fully completed training corpus."""
    completion_path = profile_directory / COMPLETION_FILENAME
    temporary_path = profile_directory / f".{COMPLETION_FILENAME}.tmp"
    record: dict[str, object] = {
        "schema_version": COMPLETION_SCHEMA_VERSION,
        "executable": executable_metadata,
        "weights": {
            "path": str(weights_path),
            "sha256": weights_sha256,
        },
        "corpus": training_corpus,
        "completed_repetitions": completed_repetitions,
        "profile_files": [
            str(profile_path.relative_to(profile_directory))
            for profile_path in profile_files
        ],
    }
    try:
        with temporary_path.open("w", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2, sort_keys=True)
            stream.write("\n")
        temporary_path.replace(completion_path)
    except OSError as error:
        with contextlib.suppress(OSError):
            temporary_path.unlink(missing_ok=True)
        raise benchmark.BenchmarkError(
            f"could not write PGO completion record {completion_path}: {error}"
        ) from error
    return completion_path


def find_gcda_files(profile_directory: Path) -> list[Path]:
    """Return the regular GCC profile files beneath a profile directory."""
    try:
        return sorted(
            profile_path
            for profile_path in profile_directory.rglob("*.gcda")
            if profile_path.is_file()
        )
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not inspect GCC profile data in {profile_directory}: {error}"
        ) from error


def train(
    executable: Path,
    weighted_cases: Sequence[WeightedCase],
    timeout: float,
) -> int:
    """Run every weighted case serially in benchmark-style isolation."""
    total_repetitions = sum(
        weighted_case.repetitions for weighted_case in weighted_cases
    )
    completed_repetitions = 0
    scheduled_cases = [
        weighted_case.case
        for weighted_case in weighted_cases
        for _ in range(weighted_case.repetitions)
    ]

    with tempfile.TemporaryDirectory(
        prefix="assemblycpp-pgo-training-"
    ) as temp_directory:
        prepared_cases = iter(
            benchmark.prepare_cases(scheduled_cases, Path(temp_directory))
        )
        for weighted_case in weighted_cases:
            print(
                f"Training {weighted_case.case.name}: "
                f"{weighted_case.repetitions} "
                f"{'repetition' if weighted_case.repetitions == 1 else 'repetitions'}",
                flush=True,
            )
            for repetition in range(1, weighted_case.repetitions + 1):
                prepared_case = next(prepared_cases)
                try:
                    benchmark.run_once(executable, prepared_case, timeout)
                except benchmark.BenchmarkError as error:
                    raise benchmark.BenchmarkError(
                        f"training {weighted_case.case.name} repetition "
                        f"{repetition}/{weighted_case.repetitions}: {error}"
                    ) from error
                completed_repetitions += 1
            print(
                f"  completed {completed_repetitions}/{total_repetitions}",
                flush=True,
            )
    return completed_repetitions


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate GCC PGO profiles with the weighted benchmark corpus. "
            "The executable must be built with -fprofile-generate."
        )
    )
    parser.add_argument(
        "--executable",
        type=Path,
        required=True,
        help="AssemblyCpp executable built with -fprofile-generate",
    )
    parser.add_argument(
        "--profile-dir",
        type=Path,
        required=True,
        help="PGO profile directory; removes existing *.gcda files",
    )
    parser.add_argument(
        "--weights",
        type=Path,
        default=DEFAULT_WEIGHTS,
        help=(
            "weights TSV "
            f"(default: {DEFAULT_WEIGHTS.relative_to(benchmark.REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=benchmark.positive_float,
        default=600.0,
        help="per-calculation timeout in seconds (default: 600)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)
    try:
        executable = benchmark.resolve_executable(arguments.executable)
        manifest_path, benchmark_cases = benchmark.load_manifest(
            benchmark.DEFAULT_MANIFEST
        )
        weights_path, weighted_cases = load_training_weights(
            arguments.weights,
            benchmark_cases,
        )
        weights_sha256 = file_sha256(weights_path)
        training_corpus = corpus_metadata(manifest_path, benchmark_cases)
        profile_directory, removed_profile_count = prepare_profile_directory(
            arguments.profile_dir
        )
        executable_metadata = benchmark.executable_metadata(executable)

        print(f"Executable: {executable}")
        print(f"Weights: {weights_path}")
        print(f"Profile directory: {profile_directory}")
        print(f"Cleared GCC profile files: {removed_profile_count}")
        print(
            f"Training cases: {len(weighted_cases)}, "
            "repetitions: "
            f"{sum(weighted_case.repetitions for weighted_case in weighted_cases)}",
            flush=True,
        )

        completed_repetitions = train(
            executable,
            weighted_cases,
            arguments.timeout,
        )
        benchmark.verify_executable_unchanged(
            executable,
            executable_metadata,
            "training",
        )
        if file_sha256(weights_path) != weights_sha256:
            raise benchmark.BenchmarkError(
                f"PGO training weights changed during training: {weights_path}"
            )
        if corpus_metadata(manifest_path, benchmark_cases) != training_corpus:
            raise benchmark.BenchmarkError(
                "benchmark manifest or corpus inputs changed during PGO training"
            )
        profile_files = find_gcda_files(profile_directory)
        if not profile_files:
            raise benchmark.BenchmarkError(
                f"training produced no .gcda files in {profile_directory}; "
                "ensure the executable was built with -fprofile-generate "
                "targeting this directory"
            )

        completion_path = write_completion_record(
            profile_directory,
            weights_path,
            weights_sha256,
            training_corpus,
            executable_metadata,
            completed_repetitions,
            profile_files,
        )

        print(
            f"Training complete: {completed_repetitions} "
            f"{'repetition' if completed_repetitions == 1 else 'repetitions'}"
        )
        print(f"GCC profile files: {len(profile_files)}")
        print(f"Completion record: {completion_path}")
    except (benchmark.BenchmarkError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
