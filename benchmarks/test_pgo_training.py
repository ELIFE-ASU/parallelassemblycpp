from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from benchmarks import benchmark, pgo_training


class PgoTrainingTests(unittest.TestCase):
    def create_case(self, directory: Path, name: str) -> benchmark.BenchmarkCase:
        source = directory / f"{name}.mol"
        source.write_text("fixture\n", encoding="utf-8")
        return benchmark.BenchmarkCase(
            name=name,
            source=source,
            expected_assembly_index=7,
            expectation="reviewed",
            suites=("quick",),
            workload=name,
        )

    def write_weights(self, directory: Path, contents: str) -> Path:
        path = directory / "weights.tsv"
        path.write_text(contents, encoding="utf-8")
        return path

    def test_default_weights_cover_corpus_once_with_tuned_values(self) -> None:
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        _, entries = pgo_training.load_training_weights(
            pgo_training.DEFAULT_WEIGHTS,
            cases,
        )
        weights = {entry.case.name: entry.repetitions for entry in entries}

        self.assertEqual(len(entries), len(cases))
        self.assertEqual(len(weights), len(cases))
        tuned_weights = {
            "sr1001": 512,
            "THC": 256,
            "sucrose": 256,
            "Cefquinome": 128,
            "Cefpirome": 1024,
            "amino-acid-scale-03c": 256,
        }
        for case in cases:
            if case.name in tuned_weights:
                self.assertEqual(weights[case.name], tuned_weights[case.name])
            elif "full" in case.suites:
                self.assertEqual(weights[case.name], 32)
            elif case.name == "phosphatidylcholine":
                self.assertEqual(weights[case.name], 4)
            elif case.name in {"erythromycin", "clarithromycin"}:
                self.assertEqual(weights[case.name], 1)
            elif "scaling" in case.suites:
                self.assertEqual(weights[case.name], 8)
            else:
                self.fail(f"unclassified default training case: {case.name}")
        self.assertEqual(sum(weights.values()), 2854)

    def test_weights_reject_invalid_header_and_column_count(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            case = self.create_case(directory, "one")
            invalid_header = self.write_weights(directory, "case\tweight\none\t1\n")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "invalid header"):
                pgo_training.load_training_weights(invalid_header, [case])

            invalid_row = self.write_weights(directory, "name\trepetitions\none\n")
            with self.assertRaisesRegex(benchmark.BenchmarkError, "expected 2 columns"):
                pgo_training.load_training_weights(invalid_row, [case])

    def test_weights_reject_duplicates_unknown_names_and_missing_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            one = self.create_case(directory, "one")
            two = self.create_case(directory, "two")

            duplicate = self.write_weights(
                directory,
                "name\trepetitions\none\t1\none\t2\ntwo\t1\n",
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "duplicate"):
                pgo_training.load_training_weights(duplicate, [one, two])

            unknown = self.write_weights(
                directory,
                "name\trepetitions\none\t1\nunknown\t1\n",
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "unknown"):
                pgo_training.load_training_weights(unknown, [one, two])

            missing = self.write_weights(
                directory,
                "name\trepetitions\none\t1\n",
            )
            with self.assertRaisesRegex(benchmark.BenchmarkError, "do not cover.*two"):
                pgo_training.load_training_weights(missing, [one, two])

    def test_weights_require_positive_integer_repetitions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            case = self.create_case(directory, "one")
            for invalid in ("0", "-1", "1.5", "many"):
                with self.subTest(invalid=invalid):
                    weights = self.write_weights(
                        directory,
                        f"name\trepetitions\none\t{invalid}\n",
                    )
                    with self.assertRaisesRegex(
                        benchmark.BenchmarkError,
                        "invalid repetitions",
                    ):
                        pgo_training.load_training_weights(weights, [case])

    def test_prepare_profile_directory_removes_only_gcda_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            profiles = directory / "profiles"
            nested = profiles / "objects"
            nested.mkdir(parents=True)
            (profiles / pgo_training.DIRECTORY_MARKER_FILENAME).write_text(
                "marker",
                encoding="utf-8",
            )
            (profiles / "root.gcda").write_text("old", encoding="utf-8")
            (nested / "nested.gcda").write_text("old", encoding="utf-8")
            (profiles / "default.profraw").write_text("keep", encoding="utf-8")
            (nested / "notes.gcno").write_text("keep", encoding="utf-8")
            completion = profiles / pgo_training.COMPLETION_FILENAME
            completion.write_text("stale", encoding="utf-8")
            outside = directory / "outside.gcda"
            outside.write_text("keep", encoding="utf-8")

            resolved, removed = pgo_training.prepare_profile_directory(profiles)

            self.assertEqual(resolved, profiles.resolve())
            self.assertEqual(removed, 2)
            self.assertFalse((profiles / "root.gcda").exists())
            self.assertFalse((nested / "nested.gcda").exists())
            self.assertTrue((profiles / "default.profraw").is_file())
            self.assertTrue((nested / "notes.gcno").is_file())
            self.assertFalse(completion.exists())
            self.assertTrue(outside.is_file())

    def test_prepare_profile_directory_rejects_broad_protected_paths(self) -> None:
        protected = {
            Path.home().resolve(),
            benchmark.REPOSITORY_ROOT.resolve(),
            Path.cwd().resolve(),
            Path(tempfile.gettempdir()).resolve(),
        }
        for path in protected:
            with (
                self.subTest(path=path),
                self.assertRaisesRegex(
                    benchmark.BenchmarkError,
                    "protected directory",
                ),
            ):
                pgo_training.prepare_profile_directory(path)

    def test_prepare_profile_directory_refuses_unmarked_nonempty_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            profiles = Path(temp_directory) / "profiles"
            profiles.mkdir()
            unrelated = profiles / "unrelated.gcda"
            unrelated.write_text("keep", encoding="utf-8")

            with self.assertRaisesRegex(
                benchmark.BenchmarkError,
                "refusing to adopt non-empty",
            ):
                pgo_training.prepare_profile_directory(profiles)

            self.assertTrue(unrelated.is_file())

    def test_train_runs_weighted_cases_serially(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            one = self.create_case(directory, "one")
            two = self.create_case(directory, "two")
            weighted = [
                pgo_training.WeightedCase(one, 2),
                pgo_training.WeightedCase(two, 1),
            ]
            calls: list[str] = []
            working_directories: list[Path] = []

            def fake_run_once(
                executable: Path,
                prepared: benchmark.PreparedCase,
                timeout: float,
            ) -> benchmark.Measurement:
                self.assertEqual(executable, Path("AssemblyCpp"))
                self.assertEqual(timeout, 12.0)
                calls.append(prepared.case.name)
                working_directories.append(prepared.working_directory)
                return benchmark.Measurement(0.1, 10, 7)

            with (
                mock.patch.object(benchmark, "run_once", side_effect=fake_run_once),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                completed = pgo_training.train(
                    Path("AssemblyCpp"), weighted, timeout=12.0
                )

            self.assertEqual(completed, 3)
            self.assertEqual(calls, ["one", "one", "two"])
            self.assertEqual(len(set(working_directories)), 3)

    def test_main_requires_new_profile_data(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = directory / "AssemblyCpp"
            executable.write_bytes(b"instrumented")
            profiles = directory / "profiles"
            stderr = io.StringIO()

            with (
                mock.patch.object(pgo_training, "train", return_value=590),
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(stderr),
            ):
                status = pgo_training.main(
                    [
                        "--executable",
                        str(executable),
                        "--profile-dir",
                        str(profiles),
                    ]
                )

            self.assertEqual(status, 1)
            self.assertIn("produced no .gcda", stderr.getvalue())

    def test_main_accepts_fresh_profile_data(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = directory / "AssemblyCpp"
            executable.write_bytes(b"instrumented")
            profiles = directory / "profiles"
            profiles.mkdir()
            (profiles / pgo_training.DIRECTORY_MARKER_FILENAME).write_text(
                "marker",
                encoding="utf-8",
            )
            stale = profiles / "stale.gcda"
            stale.write_text("stale", encoding="utf-8")
            keep = profiles / "keep.txt"
            keep.write_text("keep", encoding="utf-8")

            def fake_train(
                executable_path: Path,
                weighted_cases: list[pgo_training.WeightedCase],
                timeout: float,
            ) -> int:
                self.assertEqual(executable_path, executable.resolve())
                self.assertEqual(timeout, 9.0)
                self.assertEqual(len(weighted_cases), 31)
                self.assertFalse(stale.exists())
                (profiles / "fresh.gcda").write_text("fresh", encoding="utf-8")
                return sum(entry.repetitions for entry in weighted_cases)

            stdout = io.StringIO()
            with (
                mock.patch.object(pgo_training, "train", side_effect=fake_train),
                contextlib.redirect_stdout(stdout),
            ):
                status = pgo_training.main(
                    [
                        "--executable",
                        str(executable),
                        "--profile-dir",
                        str(profiles),
                        "--timeout",
                        "9",
                    ]
                )

            self.assertEqual(status, 0)
            self.assertTrue((profiles / "fresh.gcda").is_file())
            self.assertTrue(keep.is_file())
            self.assertIn("Training complete", stdout.getvalue())
            completion = profiles / pgo_training.COMPLETION_FILENAME
            record = json.loads(completion.read_text(encoding="utf-8"))
            self.assertEqual(
                record["schema_version"],
                pgo_training.COMPLETION_SCHEMA_VERSION,
            )
            self.assertEqual(record["completed_repetitions"], 2854)
            self.assertEqual(record["profile_files"], ["fresh.gcda"])
            self.assertEqual(
                record["corpus"]["manifest"]["sha256"],
                pgo_training.file_sha256(benchmark.DEFAULT_MANIFEST),
            )
            self.assertEqual(len(record["corpus"]["inputs"]), 31)

    def test_failed_training_leaves_no_completion_record(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = directory / "AssemblyCpp"
            executable.write_bytes(b"instrumented")
            profiles = directory / "profiles"
            profiles.mkdir()
            (profiles / pgo_training.DIRECTORY_MARKER_FILENAME).write_text(
                "marker",
                encoding="utf-8",
            )
            completion = profiles / pgo_training.COMPLETION_FILENAME
            completion.write_text("stale", encoding="utf-8")

            with (
                mock.patch.object(
                    pgo_training,
                    "train",
                    side_effect=benchmark.BenchmarkError("stopped"),
                ),
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                status = pgo_training.main(
                    [
                        "--executable",
                        str(executable),
                        "--profile-dir",
                        str(profiles),
                    ]
                )

            self.assertEqual(status, 1)
            self.assertFalse(completion.exists())

    def test_changed_weights_leave_no_completion_record(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            directory = Path(temp_directory)
            executable = directory / "AssemblyCpp"
            executable.write_bytes(b"instrumented")
            profiles = directory / "profiles"
            weights = directory / "weights.tsv"
            weights.write_bytes(pgo_training.DEFAULT_WEIGHTS.read_bytes())

            def fake_train(
                executable_path: Path,
                weighted_cases: list[pgo_training.WeightedCase],
                timeout: float,
            ) -> int:
                weights.write_text("changed\n", encoding="utf-8")
                (profiles / "fresh.gcda").write_text("fresh", encoding="utf-8")
                return sum(entry.repetitions for entry in weighted_cases)

            stderr = io.StringIO()
            with (
                mock.patch.object(pgo_training, "train", side_effect=fake_train),
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(stderr),
            ):
                status = pgo_training.main(
                    [
                        "--executable",
                        str(executable),
                        "--profile-dir",
                        str(profiles),
                        "--weights",
                        str(weights),
                    ]
                )

            self.assertEqual(status, 1)
            self.assertIn("weights changed during training", stderr.getvalue())
            self.assertFalse(
                (profiles / pgo_training.COMPLETION_FILENAME).exists()
            )


if __name__ == "__main__":
    unittest.main()
