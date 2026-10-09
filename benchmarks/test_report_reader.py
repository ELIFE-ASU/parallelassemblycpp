from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from benchmarks import check_parallel_scaling as scaling
from benchmarks import check_speedups as speedups


class ReportReaderTests(unittest.TestCase):
    def test_legacy_execution_and_public_identity_types(self) -> None:
        path = Path("report.json")
        self.assertEqual(
            speedups.execution_identity({}, "candidate", path), ((), (), ())
        )
        self.assertEqual(
            scaling.execution_identity({}, "candidate", path),
            scaling.ExecutionIdentity((), (), ()),
        )
        document = {
            "execution": {
                "candidate": {
                    "launcher": ["launcher", "space in argument"],
                    "arguments": ["--parallel=on"],
                    "environment": {"Z": "last", "A": "first"},
                }
            }
        }
        expected = (
            ("launcher", "space in argument"),
            ("--parallel=on",),
            (("A", "first"), ("Z", "last")),
        )
        self.assertEqual(
            speedups.execution_identity(document, "candidate", path), expected
        )
        self.assertEqual(
            scaling.execution_identity(document, "candidate", path),
            scaling.ExecutionIdentity(*expected),
        )

    def test_corpus_string_normalization_remains_gate_specific(self) -> None:
        corpus = {
            "manifest": {"sha256": "manifest"},
            "inputs": [
                {"name": "sample", "sha256": "first"},
                {"name": " sample ", "sha256": "second"},
            ],
        }
        document = {"corpus": corpus}
        speedups.validate_corpus_identity(document, Path("report.json"), corpus)
        with self.assertRaisesRegex(scaling.ScalingError, "duplicate corpus input"):
            scaling.corpus_identity(document, Path("report.json"))

    def test_empty_input_policy_remains_gate_specific(self) -> None:
        corpus = {"manifest": {"sha256": "manifest"}, "inputs": []}
        document = {"corpus": corpus}
        speedups.validate_corpus_identity(document, Path("report.json"), corpus)
        with self.assertRaisesRegex(scaling.ScalingError, "missing corpus input"):
            scaling.corpus_identity(document, Path("report.json"))

    def test_scaling_wall_samples_do_not_require_cpu_clock_samples(self) -> None:
        case = {
            "candidate": {
                "measurements": [{"round": 1, "wall_seconds": 1.5, "assembly_index": 7}]
            }
        }
        arguments = (case, "candidate", 1, 7, "sample", Path("report.json"))
        self.assertEqual(scaling.parse_wall_samples(*arguments), (1.5,))
        with self.assertRaisesRegex(
            speedups.GateError, "invalid candidate clock sample"
        ):
            speedups.parse_measurements(*arguments)

    def test_load_diagnostics_preserve_unexpanded_path_and_object_wording(self) -> None:
        path = Path("~/expected an object.json")
        with tempfile.TemporaryDirectory() as temporary:
            expanded = Path(temporary) / "expected an object.json"
            for document, detail in (
                ([], "expected a JSON object"),
                ({"schema_version": True}, "expected schema_version 2"),
            ):
                expanded.write_text(json.dumps(document), encoding="utf-8")
                with (
                    self.subTest(document=document),
                    patch.object(Path, "expanduser", return_value=expanded),
                    self.assertRaises(speedups.GateError) as error,
                ):
                    speedups.load_result(path)
                self.assertEqual(
                    str(error.exception), f"invalid benchmark report {path}: {detail}"
                )
            expanded.write_text("{", encoding="utf-8")
            with (
                patch.object(Path, "expanduser", return_value=expanded),
                self.assertRaises(speedups.GateError) as error,
            ):
                speedups.load_result(path)
            self.assertTrue(
                str(error.exception).startswith(
                    f"could not read benchmark report {path}: "
                )
            )


if __name__ == "__main__":
    unittest.main()
