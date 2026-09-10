from __future__ import annotations

import contextlib
import copy
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "output/benchmarks/amino-acid-cpp-rust-2026-08-26/run_benchmark.py"
)
SPEC = importlib.util.spec_from_file_location("archived_amino_benchmark", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
amino_benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(amino_benchmark)


class ArchivedAminoPlotTests(unittest.TestCase):
    def test_completed_run_saves_custom_plot_and_preserves_samples(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output = directory / "reports" / "custom.json"
            for case in amino_benchmark.CASES[:2]:
                (directory / case[1]).write_text("fixture input", encoding="utf-8")
            report = {
                "status": "running",
                "variants": {
                    "current_serial": {**amino_benchmark.VARIANTS["current_serial"]}
                },
                "cases": [],
            }
            report["variants"]["current_serial"].pop("executable")
            arguments = [
                str(SCRIPT),
                "--first",
                "2",
                "--last",
                "3",
                "--runs",
                "2",
                "--warmup",
                "0",
                "--variants",
                "current_serial",
                "--current-serial",
                sys.executable,
                "--input-directory",
                str(directory),
                "--output",
                str(output),
            ]
            with (
                mock.patch.object(sys, "argv", arguments),
                mock.patch.object(
                    amino_benchmark, "VARIANTS", copy.deepcopy(amino_benchmark.VARIANTS)
                ),
                mock.patch.object(
                    amino_benchmark, "report_metadata", return_value=report
                ),
                mock.patch.object(
                    amino_benchmark, "run_once", side_effect=[1.0, 3.0, 2.0, 4.0]
                ) as run,
                contextlib.redirect_stdout(io.StringIO()) as stdout,
            ):
                amino_benchmark.main()
            self.assertEqual(run.call_count, 4)
            saved = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(saved["status"], "complete")
            self.assertEqual(
                [
                    case["timings"]["current_serial"]["samples"]
                    for case in saved["cases"]
                ],
                [[1.0, 3.0], [2.0, 4.0]],
            )
            self.assertTrue(output.with_name("summary.csv").is_file())
            plot = output.with_suffix(".png")
            self.assertTrue(plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertIn(str(plot), stdout.getvalue())
            pdf_plot = output.with_suffix(".pdf")
            self.assertTrue(pdf_plot.read_bytes().startswith(b"%PDF-"))
            self.assertIn(str(pdf_plot), stdout.getvalue())

    def test_plot_renders_all_selected_variants_without_mutating_report(self) -> None:
        report = json.loads(
            SCRIPT.with_name("benchmark.json").read_text(encoding="utf-8")
        )
        original = copy.deepcopy(report)
        with tempfile.TemporaryDirectory() as temporary:
            plot = Path(temporary) / "nested" / "all-variants.png"
            amino_benchmark.write_plot(plot, report)
            self.assertTrue(plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertTrue(plot.with_suffix(".pdf").read_bytes().startswith(b"%PDF-"))
        self.assertEqual(report, original)

    def test_missing_matplotlib_fails_before_measurement(self) -> None:
        with (
            mock.patch.object(sys, "argv", [str(SCRIPT)]),
            mock.patch.dict(sys.modules, {"matplotlib.backends.backend_agg": None}),
            mock.patch.object(amino_benchmark, "run_once") as run,
            mock.patch.object(amino_benchmark, "report_metadata") as metadata,
            self.assertRaisesRegex(SystemExit, "pip install matplotlib"),
        ):
            amino_benchmark.main()
        run.assert_not_called()
        metadata.assert_not_called()

    def test_output_plot_aliases_are_rejected_before_benchmark_work(self) -> None:
        for first_suffix, second_suffix, alias in (
            (".png", ".png", "same-path"),
            (".pdf", ".pdf", "same-path"),
            (".json", ".png", "symlink"),
            (".json", ".png", "hardlink"),
            (".json", ".pdf", "symlink"),
            (".json", ".pdf", "hardlink"),
            (".png", ".pdf", "symlink"),
            (".png", ".pdf", "hardlink"),
        ):
            with (
                self.subTest(first=first_suffix, second=second_suffix, alias=alias),
                tempfile.TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                output = directory / (
                    f"results{first_suffix}" if alias == "same-path" else "results.json"
                )
                output.write_text("existing report", encoding="utf-8")
                first_path = output.with_suffix(first_suffix)
                first_path.write_text("existing report", encoding="utf-8")
                second_path = output.with_suffix(second_suffix)
                if alias == "symlink":
                    second_path.symlink_to(first_path)
                elif alias == "hardlink":
                    second_path.hardlink_to(first_path)
                with (
                    mock.patch.object(
                        sys, "argv", [str(SCRIPT), "--output", str(output)]
                    ),
                    mock.patch.object(amino_benchmark, "require_matplotlib") as require,
                    mock.patch.object(amino_benchmark, "report_metadata") as metadata,
                    mock.patch.object(amino_benchmark, "run_once") as run,
                    self.assertRaisesRegex(SystemExit, "--output.*different files"),
                ):
                    amino_benchmark.main()
                require.assert_not_called()
                metadata.assert_not_called()
                run.assert_not_called()
                self.assertEqual(output.read_text(encoding="utf-8"), "existing report")
                self.assertEqual(
                    first_path.read_text(encoding="utf-8"), "existing report"
                )
                self.assertEqual(
                    second_path.read_text(encoding="utf-8"), "existing report"
                )

    def test_help_does_not_load_matplotlib(self) -> None:
        with (
            mock.patch.object(sys, "argv", [str(SCRIPT), "--help"]),
            mock.patch.object(amino_benchmark, "require_matplotlib") as require,
            contextlib.redirect_stdout(io.StringIO()),
            self.assertRaises(SystemExit) as raised,
        ):
            amino_benchmark.main()
        self.assertEqual(raised.exception.code, 0)
        require.assert_not_called()


if __name__ == "__main__":
    unittest.main()
