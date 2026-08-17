from __future__ import annotations

import contextlib
import csv
import io
import json
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest import mock

from benchmarks import benchmark


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
                import pathlib
                import sys

                input_name = sys.argv[-1]
                output_name = input_name.removesuffix(".mol") + "Out"
                pathlib.Path(output_name).write_text(
                    f"{{input_name}} has assembly index: {assembly_index}\\n"
                    "time elapsed: {clock_ticks}\\n",
                    encoding="utf-8",
                )
                """
            ),
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return executable

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
        self.assertEqual(len(scaling), 7)
        self.assertTrue(all(case.expectation == "provisional" for case in scaling))
        self.assertEqual(
            [case.expected_assembly_index for case in scaling],
            [10, 13, 15, 17, 19, 21, 23],
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
            ],
        )

    def test_scaling_corpus_is_a_cumulative_size_series(self) -> None:
        _, cases = benchmark.load_manifest(benchmark.DEFAULT_MANIFEST)
        scaling = benchmark.select_cases(cases, "scaling", [])
        expected_sizes = [
            (18, 16, 2),
            (27, 24, 3),
            (36, 32, 4),
            (43, 38, 5),
            (53, 47, 6),
            (63, 56, 7),
            (68, 60, 8),
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

            def root(atom: int) -> int:
                while parents[atom] != atom:
                    parents[atom] = parents[parents[atom]]
                    atom = parents[atom]
                return atom

            for line in bond_lines:
                first = root(int(line[0:3]) - 1)
                second = root(int(line[3:6]) - 1)
                parents[first] = second
            component_count = len({root(atom) for atom in range(atom_count)})

            self.assertEqual(
                (atom_count, bond_count, component_count),
                (expected_atoms, expected_bonds, expected_components),
                case.name,
            )
            self.assertEqual(bond_count, atom_count - component_count, case.name)
            self.assertNotIn("H", {line[31:34].strip() for line in atom_lines})
            self.assertEqual(lines[4 + atom_count + bond_count], "M  END")
            blocks.append((atom_lines, bond_lines))

        for previous, current in zip(blocks, blocks[1:]):
            previous_atoms, previous_bonds = previous
            current_atoms, current_bonds = current
            self.assertEqual(current_atoms[: len(previous_atoms)], previous_atoms)
            self.assertEqual(current_bonds[: len(previous_bonds)], previous_bonds)

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
            report_path = directory / "report.json"

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
                    ]
                )

            self.assertEqual(status, 0)
            self.assertIn("Comparison summary", output.getvalue())
            report = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["schema_version"], 1)
            self.assertEqual(len(report["cases"]), 1)
            case = report["cases"][0]
            self.assertEqual(len(case["candidate"]["measurements"]), 6)
            self.assertEqual(len(case["baseline"]["measurements"]), 6)
            self.assertEqual(case["comparison"]["paired_clock_speedup"]["median"], 2.0)
            self.assertIn("sha256", report["executables"]["candidate"])
            self.assertTrue(report["platform"]["cpu_model"])
            self.assertEqual(report["runs"], 6)
            self.assertEqual(
                report["comparison"]["primary_metric"],
                "paired_round_wall_speedup",
            )

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
