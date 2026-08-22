"""Train a GCC-instrumented AssemblyCpp executable with a weighted corpus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
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
    entries: list[WeightedCase] = []
    names: set[str] = set()

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
                if name in names:
                    raise benchmark.BenchmarkError(
                        f"duplicate benchmark name in {weights_path}:{line_number}: "
                        f"{name!r}"
                    )
                names.add(name)

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

                case = cases_by_name.get(name)
                if case is None:
                    raise benchmark.BenchmarkError(
                        f"unknown benchmark name in {weights_path}:{line_number}: "
                        f"{name!r}"
                    )
                entries.append(WeightedCase(case, repetitions))
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not read PGO training weights {weights_path}: {error}"
        ) from error

    missing = [case.name for case in cases if case.name not in names]
    if missing:
        raise benchmark.BenchmarkError(
            f"PGO training weights do not cover benchmark case(s): "
            f"{', '.join(missing)}"
        )
    if not entries:
        raise benchmark.BenchmarkError(
            f"PGO training weights are empty: {weights_path}"
        )
    return weights_path, entries


def prepare_profile_directory(path: Path) -> tuple[Path, int]:
    """Create the requested directory and remove only existing GCC data files."""
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

    removed = 0
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
        for profile in sorted(profile_directory.rglob("*.gcda")):
            if profile.is_file() or profile.is_symlink():
                profile.unlink()
                removed += 1
    except OSError as error:
        raise benchmark.BenchmarkError(
            f"could not clear GCC profile data in {profile_directory}: {error}"
        ) from error
    return profile_directory, removed


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
    completed: int,
    profiles: Sequence[Path],
) -> Path:
    """Atomically mark a fully completed training corpus."""
    completion_path = profile_directory / COMPLETION_FILENAME
    temporary_path = profile_directory / f".{COMPLETION_FILENAME}.tmp"
    record = {
        "schema_version": COMPLETION_SCHEMA_VERSION,
        "executable": executable_metadata,
        "weights": {
            "path": str(weights_path),
            "sha256": weights_sha256,
        },
        "corpus": training_corpus,
        "completed_repetitions": completed,
        "profile_files": [
            str(profile.relative_to(profile_directory)) for profile in profiles
        ],
    }
    try:
        with temporary_path.open("w", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary_path, completion_path)
    except OSError as error:
        try:
            temporary_path.unlink(missing_ok=True)
        except OSError:
            pass
        raise benchmark.BenchmarkError(
            f"could not write PGO completion record {completion_path}: {error}"
        ) from error
    return completion_path


def find_gcda_files(profile_directory: Path) -> list[Path]:
    """Return the regular GCC profile files beneath a profile directory."""
    try:
        return sorted(
            path
            for path in profile_directory.rglob("*.gcda")
            if path.is_file()
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
    total_repetitions = sum(entry.repetitions for entry in weighted_cases)
    completed_repetitions = 0
    scheduled_cases = [
        entry.case
        for entry in weighted_cases
        for _ in range(entry.repetitions)
    ]

    with tempfile.TemporaryDirectory(prefix="assemblycpp-pgo-training-") as directory:
        prepared_cases = iter(
            benchmark.prepare_cases(scheduled_cases, Path(directory))
        )
        for entry in weighted_cases:
            print(
                f"Training {entry.case.name}: "
                f"{entry.repetitions} repetition(s)",
                flush=True,
            )
            for repetition in range(1, entry.repetitions + 1):
                prepared = next(prepared_cases)
                try:
                    benchmark.run_once(executable, prepared, timeout)
                except benchmark.BenchmarkError as error:
                    raise benchmark.BenchmarkError(
                        f"training {entry.case.name} repetition "
                        f"{repetition}/{entry.repetitions}: {error}"
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
            "Run the weighted AssemblyCpp corpus to create GCC PGO profile data. "
            "The supplied executable must have been built with "
            "-fprofile-generate."
        )
    )
    parser.add_argument(
        "--executable",
        type=Path,
        required=True,
        help="GCC profile-generation AssemblyCpp executable",
    )
    parser.add_argument(
        "--profile-dir",
        type=Path,
        required=True,
        help="explicit GCC profile directory; existing *.gcda files are removed",
    )
    parser.add_argument(
        "--weights",
        type=Path,
        default=DEFAULT_WEIGHTS,
        help=f"training weights TSV (default: {DEFAULT_WEIGHTS})",
    )
    parser.add_argument(
        "--timeout",
        type=benchmark.positive_float,
        default=300.0,
        help="per-calculation timeout in seconds (default: 300)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)
    try:
        executable = benchmark.resolve_executable(arguments.executable)
        manifest_path, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        weights_path, weighted_cases = load_training_weights(arguments.weights, cases)
        weights_sha256 = file_sha256(weights_path)
        training_corpus = corpus_metadata(manifest_path, cases)
        profile_directory, removed = prepare_profile_directory(arguments.profile_dir)
        metadata = benchmark.executable_metadata(executable)

        print(f"Executable: {executable}")
        print(f"Weights: {weights_path}")
        print(f"Profile directory: {profile_directory}")
        print(f"Cleared GCC profile files: {removed}")
        print(
            f"Training cases: {len(weighted_cases)}, "
            f"repetitions: {sum(entry.repetitions for entry in weighted_cases)}",
            flush=True,
        )

        completed = train(executable, weighted_cases, arguments.timeout)
        benchmark.verify_executable_unchanged(executable, metadata, "training")
        if file_sha256(weights_path) != weights_sha256:
            raise benchmark.BenchmarkError(
                f"PGO training weights changed during training: {weights_path}"
            )
        if corpus_metadata(manifest_path, cases) != training_corpus:
            raise benchmark.BenchmarkError(
                "benchmark manifest or corpus input changed during PGO training"
            )
        profiles = find_gcda_files(profile_directory)
        if not profiles:
            raise benchmark.BenchmarkError(
                f"training produced no .gcda files in {profile_directory}; "
                "verify that the executable and profile directory were built "
                "with matching -fprofile-generate settings"
            )

        completion_path = write_completion_record(
            profile_directory,
            weights_path,
            weights_sha256,
            training_corpus,
            metadata,
            completed,
            profiles,
        )

        print(f"Training complete: {completed} repetition(s)")
        print(f"GCC profile files: {len(profiles)}")
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
