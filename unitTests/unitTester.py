"""Run ParallelAssemblyCpp CLI and regression checks."""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

TEST_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = TEST_DIRECTORY.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "ParallelAssemblyCpp"
DEFAULT_MANIFEST = TEST_DIRECTORY / "regression_cases.tsv"
DEFAULT_PATHWAY_MANIFEST = TEST_DIRECTORY / "pathway_cases.tsv"
MANIFEST_HEADER = ("molecule", "expected_assembly_index")
PATHWAY_MANIFEST_HEADER = ("molecule", "expected_pathway")
PATHWAY_KEYS = {"file_graph", "remnant", "duplicates", "removed_edges"}
ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
MOLFILE_SUFFIXES = (".mol", ".sdf")


class TestConfigurationError(RuntimeError):
    """Raised when the manifest, a fixture, or a build setting is invalid."""


@dataclass(frozen=True)
class TestCase:
    name: str
    source: Path
    expected: int
    expected_pathway: Path | None = None


@dataclass(frozen=True)
class TestResult:
    case: TestCase
    actual: int | None
    duration_seconds: float
    failure: str | None = None
    error: str | None = None

    @property
    def status(self) -> str:
        if self.error is not None:
            return "ERROR"
        if self.failure is not None or self.actual != self.case.expected:
            return "FAIL"
        return "PASS"


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def resolve_test_path(path: Path) -> Path:
    if path.is_absolute() or path.exists():
        return path.resolve()

    test_relative = TEST_DIRECTORY / path
    if test_relative.exists():
        return test_relative.resolve()

    return path.resolve()


def resolve_fixture(name: str, fixture_directory: Path) -> Path:
    requested = Path(name)
    candidate = requested if requested.is_absolute() else fixture_directory / requested

    if candidate.is_file():
        return candidate.resolve()

    mol_candidate = candidate.with_suffix(".mol")
    if mol_candidate.is_file():
        return mol_candidate.resolve()

    raise TestConfigurationError(
        f"fixture {name!r} was not found as {candidate} or {mol_candidate}"
    )


def has_molfile_suffix(path: Path) -> bool:
    return path.name[-4:].lower() in MOLFILE_SUFFIXES


def molecule_output_path(path: Path, output_suffix: str) -> Path:
    input_name = path.name[:-4] if has_molfile_suffix(path) else path.name
    return path.parent / f"{input_name}{output_suffix}"


def load_manifest(path: Path) -> tuple[Path, list[TestCase]]:
    manifest = resolve_test_path(path)
    seen: dict[str, int] = {}
    cases: list[TestCase] = []

    try:
        with manifest.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != MANIFEST_HEADER:
                raise TestConfigurationError(
                    f"invalid header in {manifest}: expected {MANIFEST_HEADER}, "
                    f"got {header}"
                )

            for row in reader:
                line_number = reader.line_num
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != 2:
                    raise TestConfigurationError(
                        f"invalid row in {manifest}:{line_number}: expected 2 columns"
                    )

                name, expected_text = (value.strip() for value in row)
                if not name:
                    raise TestConfigurationError(
                        f"empty fixture name in {manifest}:{line_number}"
                    )
                if name in seen:
                    raise TestConfigurationError(
                        f"duplicate fixture {name!r} in {manifest}:{line_number}; "
                        f"first declared on line {seen[name]}"
                    )

                try:
                    expected = int(expected_text)
                except ValueError as error:
                    raise TestConfigurationError(
                        f"invalid assembly index in {manifest}:{line_number}: "
                        f"{expected_text!r}"
                    ) from error

                seen[name] = line_number
                cases.append(
                    TestCase(
                        name=name,
                        source=resolve_fixture(name, manifest.parent),
                        expected=expected,
                    )
                )
    except OSError as error:
        raise TestConfigurationError(f"cannot read {manifest}: {error}") from error

    if not cases:
        raise TestConfigurationError(f"manifest contains no test cases: {manifest}")

    return manifest, cases


def parse_pathway_document(path: Path) -> dict[str, object]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error

    if not isinstance(document, dict):
        # A wrong JSON shape is malformed file content, not caller misuse.
        raise ValueError(  # noqa: TRY004
            f"pathway document must be a JSON object: {path}"
        )
    if set(document) != PATHWAY_KEYS:
        raise ValueError(
            f"invalid pathway keys in {path}: expected {sorted(PATHWAY_KEYS)}, "
            f"got {sorted(document)}"
        )
    for key in PATHWAY_KEYS:
        if not isinstance(document[key], list):
            # Keep every pathway schema violation under the ValueError contract.
            raise ValueError(  # noqa: TRY004
                f"pathway field {key!r} must be an array in {path}"
            )

    return document


def load_pathway_manifest(path: Path, cases: Sequence[TestCase]) -> list[TestCase]:
    manifest = resolve_test_path(path)
    case_names = {case.name for case in cases}
    pathways: dict[str, Path] = {}

    try:
        with manifest.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = tuple(next(reader, ()))
            if header != PATHWAY_MANIFEST_HEADER:
                raise TestConfigurationError(
                    f"invalid header in {manifest}: expected "
                    f"{PATHWAY_MANIFEST_HEADER}, got {header}"
                )

            for row in reader:
                line_number = reader.line_num
                if not row or all(not value.strip() for value in row):
                    continue
                if len(row) != 2:
                    raise TestConfigurationError(
                        f"invalid row in {manifest}:{line_number}: expected 2 columns"
                    )

                name, relative_path = (value.strip() for value in row)
                if name in pathways:
                    raise TestConfigurationError(
                        f"duplicate pathway case {name!r} in {manifest}:{line_number}"
                    )
                if name not in case_names:
                    raise TestConfigurationError(
                        f"pathway case {name!r} is not in the regression manifest"
                    )

                expected_pathway = (manifest.parent / relative_path).resolve()
                try:
                    parse_pathway_document(expected_pathway)
                except ValueError as error:
                    raise TestConfigurationError(str(error)) from error
                pathways[name] = expected_pathway
    except OSError as error:
        raise TestConfigurationError(f"cannot read {manifest}: {error}") from error

    if not pathways:
        raise TestConfigurationError(
            f"pathway manifest contains no test cases: {manifest}"
        )

    return [
        TestCase(
            name=case.name,
            source=case.source,
            expected=case.expected,
            expected_pathway=pathways.get(case.name),
        )
        for case in cases
    ]


def audit_test_data(manifest: Path, cases: Sequence[TestCase], verbose: bool) -> None:
    fixture_directory = manifest.parent
    mol_fixtures = {
        path.resolve()
        for path in fixture_directory.iterdir()
        if path.is_file() and has_molfile_suffix(path)
    }
    referenced_mol = {case.source for case in cases if has_molfile_suffix(case.source)}
    fixture_only = sorted(mol_fixtures - referenced_mol)

    cases_by_hash: dict[str, list[TestCase]] = defaultdict(list)
    for case in cases:
        digest = hashlib.sha256(case.source.read_bytes()).hexdigest()
        cases_by_hash[digest].append(case)

    shared_content = [group for group in cases_by_hash.values() if len(group) > 1]
    conflicting_content = [
        group for group in shared_content if len({case.expected for case in group}) > 1
    ]
    if conflicting_content:
        details = "; ".join(
            ", ".join(f"{case.name}={case.expected}" for case in group)
            for group in conflicting_content
        )
        raise TestConfigurationError(
            f"byte-identical fixtures have conflicting expectations: {details}"
        )

    graph_cases = sum(not has_molfile_suffix(case.source) for case in cases)
    print(f"Manifest: {manifest}")
    print(
        f"Regression cases: {len(cases)} "
        f"({len(cases) - graph_cases} MOL/SDF, {graph_cases} graph)"
    )
    print(f"Molecule fixtures: {len(mol_fixtures)}")
    print(f"Fixture-only molecules: {len(fixture_only)}")
    print(
        "Shared-content case groups:",
        len(shared_content),
        "(consistent expectations)",
    )
    print(
        "Pathway golden cases: "
        f"{sum(case.expected_pathway is not None for case in cases)}"
    )

    if verbose and fixture_only:
        print("Fixture-only molecule names:")
        for path in fixture_only:
            print(f"  {path.stem}")


def resolve_executable(path: Path) -> Path:
    if path.is_file():
        return path.resolve()

    located = shutil.which(str(path))
    if located is not None:
        return Path(located).resolve()

    raise TestConfigurationError(
        f"executable not found: {path}. Run again with --build or compile it first."
    )


def run_cli_command(
    executable: Path, arguments: Sequence[str], working_directory: Path
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            [str(executable), *arguments],
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise TestConfigurationError(
            f"CLI check timed out after {error.timeout:g} seconds: "
            f"{shlex.join([str(executable), *arguments])}"
        ) from error
    except OSError as error:
        raise TestConfigurationError(f"cannot run CLI check: {error}") from error


def require_cli(
    condition: bool,
    message: str,
    completed: subprocess.CompletedProcess[str] | None = None,
) -> None:
    if condition:
        return

    diagnostics = ""
    if completed is not None:
        diagnostics = format_process_diagnostics(completed)
    suffix = f"; {diagnostics}" if diagnostics else ""
    raise TestConfigurationError(f"CLI check failed: {message}{suffix}")


def read_first_line_assembly_index(path: Path) -> int | None:
    """Return an index only when the first output line is fully numeric."""
    if not path.is_file():
        return None
    lines = path.read_text().splitlines()
    if not lines:
        return None
    match = ASSEMBLY_INDEX_PATTERN.search(lines[0])
    if match is None or match.end() != len(lines[0]):
        return None
    return int(match.group(1))


def read_last_intermediate_index(path: Path) -> int | None:
    """Return the assembly index in the final well-formed intermediate row."""
    if not path.is_file():
        return None
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if not rows or len(rows[-1]) != 2 or not rows[-1][1].lstrip("-").isdigit():
        return None
    return int(rows[-1][1])


def make_mask_capacity_graph(component_sizes: Sequence[int]) -> str:
    """Build a bounded-work native graph whose final atom appears in edge zero."""
    atom_colours: list[str] = []
    component_edges: list[list[tuple[int, int]]] = []
    first_atom = 1
    for component, size in enumerate(component_sizes):
        if size < 3:
            raise ValueError("capacity-test components must be cycles of size >= 3")
        vertices = list(range(first_atom, first_atom + size))
        atom_colours.extend([f"C{component}"] * size)
        component_edges.append(
            list(zip(vertices, vertices[1:] + vertices[:1], strict=True))
        )
        first_atom += size

    # Visit the last component first and start at its wraparound edge. This
    # makes the first connected-subgraph seed include the highest atom index,
    # so a 513-atom case necessarily exercises AtomMask word eight as well as
    # EdgeMask word eight.
    last_edges = component_edges[-1]
    edges = [last_edges[-1], *last_edges[:-1]]
    for cycle_edges in component_edges[:-1]:
        edges.extend(cycle_edges)

    return "\n".join(
        (
            "mask-capacity-boundary",
            str(len(atom_colours)),
            " ".join(f"{first} {second}" for first, second in edges),
            " ".join(atom_colours),
            " ".join("1" for _ in edges),
            "",
        )
    )


def run_string_record_checks(executable: Path) -> int:
    """Preserve record symbols while recognizing only LF and CRLF separators."""
    cases = (
        (b"", []),
        (b"\n", [("", -1)]),
        (b"\n\n", [("", -1), ("", -1)]),
        (b"a\n", [("a", 0)]),
        (b"a\n\n", [("a", 0), ("", -1)]),
        (b" \t\n", [(" \t", 1)]),
        (b'quote"\\\t\x00\n', [('quote"\\\t\x00', 8)]),
        (b"abab\r", [("abab\r", 3)]),
        (b"\r", [("\r", 0)]),
        (b"abab\n", [("abab", 2)]),
        (b"abab\r\n", [("abab", 2)]),
        (b"abab\r\r\n", [("abab\r", 3)]),
        (b"ab\rab\n", [("ab\rab", 3)]),
        # BOM, combining marks, and Unicode line separators remain symbols;
        # indexing counts Unicode scalars, including NUL, rather than bytes.
        ("\ufeffabab\n".encode(), [("\ufeffabab", 3)]),
        ("\u00e9e\u0301\n".encode(), [("\u00e9e\u0301", 2)]),
        ("\u0085\u2028\u2029\n".encode(), [("\u0085\u2028\u2029", 2)]),
        (
            "\x00\U0001f600\x00\U0001f600\n".encode(),
            [("\x00\U0001f600\x00\U0001f600", 2)],
        ),
        (
            b"abab\nabab\r\nabab\r\r\n\r\n\r\r\nabab\r",
            [
                ("abab", 2),
                ("abab", 2),
                ("abab\r", 3),
                ("", -1),
                ("\r", 0),
                ("abab\r", 3),
            ],
        ),
    )
    scenarios = 0
    with tempfile.TemporaryDirectory(prefix="parallelassemblycpp-records-") as name:
        root = Path(name)
        for contents, expected in cases:
            for pathway_enabled in (False, True):
                case_directory = root / f"records-{scenarios}"
                case_directory.mkdir()
                (case_directory / "input").write_bytes(contents)
                completed = run_cli_command(
                    executable,
                    ["input", "--run-strings=1", f"--pathway={int(pathway_enabled)}"],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 0,
                    f"string mode should preserve record bytes {contents!r}",
                    completed,
                )
                # Reading bytes avoids Python's universal newline conversion,
                # which would change literal CR symbols in the reported records.
                output_text = (case_directory / "inputOut").read_bytes().decode()
                actual = [
                    (record, int(index))
                    for record, index in re.findall(
                        r"^(.*) has assembly index: (-?\d+)\r?$",
                        output_text,
                        re.MULTILINE,
                    )
                ]
                require_cli(
                    actual == expected,
                    f"string records {contents!r} changed: {actual!r}, "
                    f"expected {expected!r}",
                    completed,
                )
                pathway_files = list(case_directory.glob("input_*_Pathway"))
                require_cli(
                    len(pathway_files) == (len(expected) if pathway_enabled else 0),
                    f"string records {contents!r} produced unexpected pathway files",
                    completed,
                )
                if pathway_enabled:
                    for index, (record, _) in enumerate(expected):
                        pathway = json.loads(
                            (case_directory / f"input_{index}_Pathway").read_text()
                        )
                        require_cli(
                            pathway["file_graph"][0]["Fragments"] == [record],
                            f"string pathway changed record {contents!r}",
                            completed,
                        )
                scenarios += 1
    return scenarios


def run_string_output_alias_checks(executable: Path) -> int:
    """Reject enabled output aliases without changing the source records."""
    scenarios = 0
    with tempfile.TemporaryDirectory(
        prefix="parallelassemblycpp-output-alias-"
    ) as name:
        root = Path(name)
        for target in ("inputOut", "input_0_Pathway", "input_1_Pathway"):
            for link_mode in ("symlink", "hardlink"):
                for pathway in ((True,) if target == "inputOut" else (False, True)):
                    case_directory = root / f"{target}-{link_mode}-{int(pathway)}"
                    case_directory.mkdir()
                    input_path = case_directory / "input"
                    before = b"abab\nabcabc\n"
                    input_path.write_bytes(before)
                    output_path = case_directory / target
                    if link_mode == "symlink":
                        try:
                            output_path.symlink_to(input_path.name)
                        except OSError as error:
                            # Windows can require a privilege for symlink creation.
                            if os.name == "nt" and error.winerror in (5, 1314):
                                continue
                            raise
                    else:
                        output_path.hardlink_to(input_path)
                    completed = run_cli_command(
                        executable,
                        ["input", "--run-strings=1", f"--pathway={int(pathway)}"],
                        case_directory,
                    )
                    require_cli(
                        input_path.read_bytes() == before
                        and output_path.read_bytes() == before,
                        f"string output corrupted input through {target} {link_mode}",
                        completed,
                    )
                    if target == "inputOut" or pathway:
                        require_cli(
                            completed.returncode == 1
                            and target in completed.stderr
                            and "would overwrite input" in completed.stderr,
                            f"string output must reject {target} {link_mode} aliases",
                            completed,
                        )
                    else:
                        require_cli(
                            completed.returncode == 0,
                            f"disabled pathway must allow {target} {link_mode} aliases",
                            completed,
                        )
                        indices = ASSEMBLY_INDEX_PATTERN.findall(
                            (case_directory / "inputOut").read_text()
                        )
                        require_cli(
                            indices == ["2", "3"],
                            "disabled pathway alias changed the string results",
                            completed,
                        )
                    scenarios += 1
    return scenarios


def run_graph_output_alias_checks(executable: Path, telemetry_supported: bool) -> int:
    """Reject graph output aliases before opening any enabled output file."""
    scenarios = 0
    outputs = [
        ("Out", None),
        ("Pathway", "pathway"),
        ("IntermediateMAs", "write-intermediate-mas"),
    ]
    if telemetry_supported:
        outputs.append(("Telemetry.json", "telemetry"))
    sentinel = b"existing output sentinel\n"
    with tempfile.TemporaryDirectory(
        prefix="parallelassemblycpp-graph-output-alias-"
    ) as name:
        root = Path(name)
        for mode, input_name, argument in (
            ("mol", "input.mol", "input.mol"),
            ("sdf", "input.SDF", "input.SDF"),
            ("mol-fallback", "input.mol", "input"),
            ("native", "input", "input"),
        ):
            for suffix, flag in outputs:
                for link_mode in ("symlink", "hardlink"):
                    for enabled in (True,) if flag is None else (False, True):
                        label = f"{mode}-{suffix}-{link_mode}-{int(enabled)}"
                        case_directory = root / label
                        case_directory.mkdir()
                        input_path = case_directory / input_name
                        if mode == "native":
                            input_path.write_text(
                                "butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n"
                            )
                        else:
                            shutil.copy2(TEST_DIRECTORY / "butane.mol", input_path)
                        before = input_path.read_bytes()
                        target = f"input{suffix}"
                        output_path = case_directory / target
                        if link_mode == "symlink":
                            try:
                                output_path.symlink_to(input_path.name)
                            except OSError as error:
                                # Windows can require a privilege for symlinks.
                                if os.name == "nt" and error.winerror in (5, 1314):
                                    continue
                                raise
                        else:
                            output_path.hardlink_to(input_path)

                        # All other outputs are enabled so that rejection of a
                        # later output must precede truncation of earlier ones.
                        other_outputs = [
                            case_directory / f"input{other_suffix}"
                            for other_suffix, _ in outputs
                            if other_suffix != suffix
                        ]
                        options = [
                            f"--{other_flag}="
                            f"{int(enabled if other_flag == flag else True)}"
                            for _, other_flag in outputs
                            if other_flag is not None
                        ]
                        if sys.platform.startswith("linux"):
                            options.append("--memory-report=1")
                            other_outputs.append(case_directory / "memUsage")
                        for path in other_outputs:
                            path.write_bytes(sentinel)
                        completed = run_cli_command(
                            executable, [argument, *options], case_directory
                        )
                        require_cli(
                            input_path.read_bytes() == before
                            and output_path.read_bytes() == before,
                            f"graph output corrupted input through {label}",
                            completed,
                        )
                        if enabled:
                            require_cli(
                                completed.returncode == 1
                                and target in completed.stderr
                                and "would overwrite input" in completed.stderr,
                                f"graph output must reject {label} aliases",
                                completed,
                            )
                            require_cli(
                                all(
                                    path.read_bytes() == sentinel
                                    for path in other_outputs
                                ),
                                f"graph alias {label} changed another output "
                                "before rejection",
                                completed,
                            )
                        else:
                            require_cli(
                                completed.returncode == 0
                                and read_first_line_assembly_index(
                                    case_directory / "inputOut"
                                )
                                == 2,
                                f"disabled graph output must allow {label} aliases",
                                completed,
                            )
                        scenarios += 1
    return scenarios


def run_flag_matrix_checks(executable: Path, telemetry_supported: bool) -> int:
    """Exercise every flag spelling against actual molecular and string inputs."""
    spellings = {
        "runtime": ("runtime", "runTime"),
        "enum-max": ("enum-max", "enumMax"),
        "pathway": ("pathway",),
        "run-strings": ("run-strings", "runStrings"),
        "accept-palindromes": ("accept-palindromes", "acceptPalindromes", "palindrome"),
        "parallel": ("parallel",),
        "threads": ("threads",),
        "remove-hydrogens": ("remove-hydrogens", "removeHydrogens"),
        "verbose": ("verbose",),
        "compensate-disjoint": (
            "compensate-disjoint",
            "compensateDisjoint",
            "disjointCompensation",
        ),
        "memory-report": ("memory-report", "memTest", "testMemory"),
        "write-intermediate-mas": ("write-intermediate-mas", "writeIntermediateMAs"),
    }
    if telemetry_supported:
        spellings["telemetry"] = ("telemetry",)
    boolean_flags = set(spellings) - {"runtime", "enum-max", "parallel", "threads"}
    scenarios = 0
    with tempfile.TemporaryDirectory(prefix="parallelassemblycpp-flags-") as directory:
        root = Path(directory)

        # Parsing must reject malformed values regardless of input mode, and
        # before opening or overwriting any output files.
        invalid_values = dict.fromkeys(
            sorted(boolean_flags), ("", "2", "-1", "true", "01", "+1", " 1", "1 ")
        )
        invalid_values.update(
            {
                "runtime": ("", "-1", "+1", "1.0", " 1", "1 ", "18446744073709551616"),
                "enum-max": ("", "0", "-1", "+1", "1.0", "2147483648"),
                "threads": ("", "0", "-1", "+1", "1.0", "2147483648", "AUTO"),
                "parallel": ("", "1", "ON", "false", "automatic"),
            }
        )
        for strings in (False, True):
            mode = "string" if strings else "molecular"
            case_directory = root / f"invalid-{mode}"
            case_directory.mkdir()
            input_name = "input" if strings else "input.mol"
            input_path = case_directory / input_name
            if strings:
                input_path.write_text("abcxcba\n")
            else:
                shutil.copy2(TEST_DIRECTORY / "butane.mol", input_path)
            output = case_directory / "inputOut"
            output.write_text("existing result\n")
            for name, values in invalid_values.items():
                mode_arguments = (
                    [] if name == "run-strings" else [f"--run-strings={int(strings)}"]
                )
                for value in values:
                    arguments = [input_name, *mode_arguments, f"--{name}={value}"]
                    completed = run_cli_command(executable, arguments, case_directory)
                    require_cli(
                        completed.returncode == 2 and f"--{name}" in completed.stderr,
                        f"{mode} should reject malformed --{name}={value!r}",
                        completed,
                    )
                    require_cli(
                        output.read_text() == "existing result\n",
                        f"invalid --{name} unexpectedly changed an output",
                        completed,
                    )
                    scenarios += 1
                completed = run_cli_command(
                    executable,
                    [input_name, *mode_arguments, f"--{name}"],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 2
                    and "requires a value" in completed.stderr,
                    f"{mode} should reject --{name} without =VALUE",
                    completed,
                )
                scenarios += 1

        # Every canonical/legacy spelling accepts either dash prefix. Use real
        # calculations so accepting an alias without applying it cannot pass.
        for mode in ("mol", "native", "string"):
            strings = mode == "string"
            for name, names in spellings.items():
                value = {
                    "runtime": "0",
                    "enum-max": "1",
                    "run-strings": str(int(strings)),
                    "accept-palindromes": str(int(strings)),
                    "parallel": "off",
                    "threads": "1",
                    "remove-hydrogens": "0",
                    "compensate-disjoint": "0",
                    "write-intermediate-mas": str(int(not strings)),
                    "telemetry": str(int(not strings)),
                }.get(name, "1")
                for spelling in names:
                    for prefix in ("-", "--"):
                        case_directory = root / f"case-{scenarios}"
                        case_directory.mkdir()
                        input_name = "input.mol" if mode == "mol" else "input"
                        if mode == "mol":
                            shutil.copy2(
                                TEST_DIRECTORY / "butane.mol",
                                case_directory / input_name,
                            )
                        elif mode == "native":
                            (case_directory / input_name).write_text(
                                "butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n"
                            )
                        else:
                            (case_directory / input_name).write_text("abcxcba\n")
                        defaults = {
                            "run-strings": str(int(strings)),
                            "pathway": "0",
                            "parallel": "off",
                        }
                        defaults.pop(name, None)
                        option = f"{prefix}{spelling}={value}"
                        options = [
                            f"--{key}={setting}" for key, setting in defaults.items()
                        ]
                        # Exercise options both before and after INPUT.
                        arguments = [option, input_name, *options]
                        completed = run_cli_command(
                            executable, arguments, case_directory
                        )
                        output = case_directory / "inputOut"
                        incompatible = (
                            strings and name in {"enum-max", "remove-hydrogens"}
                        ) or (
                            name == "memory-report"
                            and not sys.platform.startswith("linux")
                        )
                        if incompatible:
                            require_cli(
                                completed.returncode == 2
                                and f"--{name}" in completed.stderr
                                and not output.exists(),
                                f"{mode} should explain inapplicable {option}",
                                completed,
                            )
                            scenarios += 1
                            continue
                        require_cli(
                            completed.returncode == 0,
                            f"{mode} should execute {option}",
                            completed,
                        )
                        expected_index = (
                            4
                            if strings and name == "accept-palindromes"
                            else (6 if strings else 2)
                        )
                        require_cli(
                            read_first_line_assembly_index(output) == expected_index,
                            f"{mode} {option} returned the wrong index",
                            completed,
                        )
                        pathway = case_directory / (
                            "input_0_Pathway" if strings else "inputPathway"
                        )
                        require_cli(
                            pathway.exists() == (name == "pathway"),
                            f"{mode} {option} did not honor pathway output",
                            completed,
                        )
                        if name == "pathway":
                            document = json.loads(pathway.read_text())
                            require_cli(
                                {"file_graph", "remnant", "duplicates"}
                                <= document.keys(),
                                f"{mode} {option} wrote malformed pathway JSON",
                                completed,
                            )
                        if name == "runtime":
                            require_cli(
                                "status: runtime limit reached" in output.read_text(),
                                f"{mode} {option} did not enforce the runtime cap",
                                completed,
                            )
                        if name == "enum-max":
                            require_cli(
                                "status: enumeration limit reached"
                                in output.read_text(),
                                f"{mode} {option} did not enforce the enumeration cap",
                                completed,
                            )
                        if name == "verbose":
                            require_cli(
                                "Input: input" in completed.stdout,
                                f"{mode} {option} did not print diagnostics",
                                completed,
                            )
                        if name == "memory-report":
                            memory = case_directory / "memUsage"
                            require_cli(
                                memory.exists() == sys.platform.startswith("linux"),
                                f"{mode} {option} did not honor memory output",
                                completed,
                            )
                        if name == "write-intermediate-mas" and not strings:
                            require_cli(
                                read_last_intermediate_index(
                                    case_directory / "inputIntermediateMAs"
                                )
                                == expected_index,
                                f"{mode} {option} did not write the final index",
                                completed,
                            )
                        if name == "telemetry" and not strings:
                            telemetry = case_directory / "inputTelemetry.json"
                            require_cli(
                                telemetry.is_file(),
                                f"{mode} {option} did not write telemetry",
                                completed,
                            )
                            require_cli(
                                isinstance(json.loads(telemetry.read_text()), dict),
                                f"{mode} {option} wrote malformed telemetry",
                                completed,
                            )
                        scenarios += 1

        # Applicability is checked after parsing the complete command, so the
        # position and spelling of the mode selector cannot change the result.
        for name, values, strings in (
            ("enum-max", ("1", "50000000"), True),
            ("remove-hydrogens", ("0", "1"), True),
            ("compensate-disjoint", ("1",), True),
            ("accept-palindromes", ("1",), False),
        ):
            for spelling in spellings[name]:
                for value in values:
                    option = f"-{spelling}={value}"
                    for mode_first in (False, True):
                        mode_option = f"--runStrings={int(strings)}"
                        arguments = (
                            [mode_option, option]
                            if mode_first
                            else [option, mode_option]
                        )
                        completed = run_cli_command(
                            executable, ["missing-input", *arguments], root
                        )
                        require_cli(
                            completed.returncode == 2
                            and f"--{name}" in completed.stderr,
                            f"applicability should reject {arguments!r} "
                            "before reading input",
                            completed,
                        )
                        scenarios += 1
                    completed = run_cli_command(
                        executable, ["--help", *arguments], root
                    )
                    require_cli(
                        completed.returncode == 0,
                        f"help should bypass the applicability of {option}",
                        completed,
                    )
                    scenarios += 1

        for arguments in ([""], ["", "input"], ["input", ""], ["--", "", "input"]):
            completed = run_cli_command(executable, arguments, root)
            require_cli(
                completed.returncode == 2 and "INPUT" in completed.stderr,
                f"an empty input argument should be rejected: {arguments!r}",
                completed,
            )
            scenarios += 1

        # Alias duplicates must refer to the same underlying flag. --help
        # bypasses mode applicability, but must not hide malformed syntax.
        for name, names in spellings.items():
            value = {
                "runtime": "1",
                "enum-max": "1",
                "parallel": "off",
                "threads": "1",
            }.get(name, "0")
            for alias in names:
                completed = run_cli_command(
                    executable,
                    ["--help", f"--{name}={value}", f"-{alias}={value}"],
                    root,
                )
                require_cli(
                    completed.returncode == 2 and "only once" in completed.stderr,
                    f"-{alias} should duplicate --{name}",
                    completed,
                )
                scenarios += 1

        for option in (
            "--runtime=18446744073709551615",
            "--enum-max=2147483647",
            "--threads=2147483647",
        ):
            completed = run_cli_command(executable, ["--help", option], root)
            require_cli(
                completed.returncode == 0,
                f"the exact upper boundary {option} should parse",
                completed,
            )
            scenarios += 1

    return scenarios


def run_input_output_matrix_checks(executable: Path, telemetry_supported: bool) -> int:
    """Check names, optional output isolation, and I/O failures in every mode."""
    scenarios = 0
    with tempfile.TemporaryDirectory(prefix="parallelassemblycpp-inputs-") as directory:
        root = Path(directory)
        for mode in ("mol", "native", "string"):
            strings = mode == "string"
            for requested_name, absolute in (
                ("input with spaces", False),
                ("-h", False),
                ("--runtime=0", False),
                ("nested dir/input=name", False),
                ("absolute input", True),
                # In string mode, MOL/SDF suffixes must not switch the parser
                # or be stripped from output names.
                *(([("input.MOL", False), ("input.sdf", False)]) if strings else []),
            ):
                case_directory = root / f"names-{scenarios}"
                case_directory.mkdir()
                filename = requested_name + (".mol" if mode == "mol" else "")
                input_path = case_directory / filename
                input_path.parent.mkdir(parents=True, exist_ok=True)
                if mode == "mol":
                    shutil.copy2(TEST_DIRECTORY / "butane.mol", input_path)
                elif mode == "native":
                    input_path.write_text("butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n")
                else:
                    input_path.write_text("abab\n")
                output_base = (
                    input_path.with_suffix("") if mode == "mol" else input_path
                )
                output_path = Path(str(output_base) + "Out")
                pathway_path = Path(
                    str(output_base) + ("_0_Pathway" if strings else "Pathway")
                )
                argument = str(input_path) if absolute else filename
                completed = run_cli_command(
                    executable,
                    [f"--run-strings={int(strings)}", "--pathway=1", "--", argument],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 0
                    and read_first_line_assembly_index(output_path) == 2,
                    f"{mode} should preserve input name {argument!r}",
                    completed,
                )
                require_cli(
                    pathway_path.is_file(),
                    f"{mode} should preserve the pathway base for {argument!r}",
                    completed,
                )
                json.loads(pathway_path.read_text())
                scenarios += 1

            for enabled in (False, True):
                case_directory = root / f"outputs-{mode}-{int(enabled)}"
                case_directory.mkdir()
                input_name = "input.mol" if mode == "mol" else "input"
                if mode == "mol":
                    shutil.copy2(
                        TEST_DIRECTORY / "butane.mol", case_directory / input_name
                    )
                elif mode == "native":
                    (case_directory / input_name).write_text(
                        "butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n"
                    )
                else:
                    (case_directory / input_name).write_text("abab\n")
                output_names = [
                    "input_0_Pathway" if strings else "inputPathway",
                    "memUsage",
                ]
                if not strings:
                    output_names.append("inputIntermediateMAs")
                    if telemetry_supported:
                        output_names.append("inputTelemetry.json")
                for name in output_names:
                    (case_directory / name).write_text("existing output sentinel\n")
                memory_enabled = enabled and sys.platform.startswith("linux")
                options = [
                    f"--run-strings={int(strings)}",
                    f"--pathway={int(enabled)}",
                    f"--memory-report={int(memory_enabled)}",
                    "--verbose=0",
                ]
                if not strings:
                    options.append(f"--write-intermediate-mas={int(enabled)}")
                    if telemetry_supported:
                        options.append(f"--telemetry={int(enabled)}")
                completed = run_cli_command(
                    executable, [input_name, *options], case_directory
                )
                require_cli(
                    completed.returncode == 0,
                    f"{mode} output toggles should work",
                    completed,
                )
                require_cli(
                    "Input:" not in completed.stdout
                    and "Graph:" not in completed.stdout,
                    f"{mode} --verbose=0 should suppress diagnostics",
                    completed,
                )
                for name in output_names:
                    should_change = enabled and (
                        name != "memUsage" or sys.platform.startswith("linux")
                    )
                    require_cli(
                        (
                            (case_directory / name).read_text()
                            != "existing output sentinel\n"
                        )
                        == should_change,
                        f"{mode} enabled={enabled} mishandled {name}",
                        completed,
                    )
                scenarios += 1

            for target in (
                "inputOut",
                "input_0_Pathway" if strings else "inputPathway",
            ):
                case_directory = root / f"output-error-{mode}-{target}"
                case_directory.mkdir()
                input_name = "input.mol" if mode == "mol" else "input"
                if mode == "mol":
                    shutil.copy2(
                        TEST_DIRECTORY / "butane.mol", case_directory / input_name
                    )
                elif mode == "native":
                    (case_directory / input_name).write_text(
                        "butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n"
                    )
                else:
                    (case_directory / input_name).write_text("abab\n")
                (case_directory / target).mkdir()
                completed = run_cli_command(
                    executable,
                    [input_name, f"--run-strings={int(strings)}", "--pathway=1"],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 1 and target in completed.stderr,
                    f"{mode} must report an unwritable {target}",
                    completed,
                )
                scenarios += 1

        for contents, byte_offset in (
            (b"\x80", 0),
            (b"\xc3", 0),
            (b"\xc3(", 0),
            (b"\xc0\xaf", 0),
            (b"\xe0\x9f\xbf", 0),
            (b"\xe2\x82", 0),
            (b"\xed\xa0\x80", 0),
            (b"\xf0\x8f\xbf\xbf", 0),
            (b"\xf0\x90\x80", 0),
            (b"\xf4\x90\x80\x80", 0),
            (b"\xf5\x80\x80\x80", 0),
            (b"\xff", 0),
            # Report the failing byte, after a valid multibyte scalar or NUL.
            (b"a\xc3\xa9\xe2(", 3),
            (b"\x00\x80", 1),
        ):
            for pathway, prefix in (
                (0, b""),
                (1, b""),
                (0, b"abab\n"),
                (1, b"abab\n"),
            ):
                case_directory = root / f"invalid-utf8-{scenarios}"
                case_directory.mkdir()
                # A valid prefix record must survive, and processing must stop
                # before either the malformed record or the following record.
                (case_directory / "input").write_bytes(prefix + contents + b"\ncdcd\n")
                completed = run_cli_command(
                    executable,
                    ["input", "--run-strings=1", f"--pathway={pathway}", "--verbose=1"],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 1
                    and f"failed on line {2 if prefix else 1} of 'input'"
                    in completed.stderr
                    and f"string input is not valid UTF-8 at byte {byte_offset}"
                    in completed.stderr,
                    f"string input {contents!r} should fail with pathway={pathway}",
                    completed,
                )
                require_cli(
                    {path.name for path in case_directory.glob("input_*_Pathway")}
                    == ({"input_0_Pathway"} if pathway and prefix else set()),
                    "invalid UTF-8 must preserve only completed pathways",
                    completed,
                )
                require_cli(
                    ASSEMBLY_INDEX_PATTERN.findall(
                        (case_directory / "inputOut").read_text()
                    )
                    == (["2"] if prefix else []),
                    "invalid UTF-8 must stop output after completed records",
                    completed,
                )
                scenarios += 1

        for mode in ("mol", "native", "string"):
            for missing in (True, False):
                case_directory = root / f"input-error-{mode}-{int(missing)}"
                case_directory.mkdir()
                input_name = "input.mol" if mode == "mol" else "input"
                if not missing:
                    (case_directory / input_name).mkdir()
                completed = run_cli_command(
                    executable,
                    [input_name, f"--run-strings={int(mode == 'string')}"],
                    case_directory,
                )
                require_cli(
                    completed.returncode == 1 and "input" in completed.stderr,
                    f"{mode} should reject "
                    f"{'missing' if missing else 'directory'} input",
                    completed,
                )
                scenarios += 1

        if sys.platform.startswith("linux"):
            # Memory reports use a fixed filename in the working directory.
            # A matching input or filesystem alias must never be overwritten.
            for mode in ("native", "string", "mol-fallback"):
                link_modes = (
                    ("symlink", "hardlink")
                    if mode == "mol-fallback"
                    else ("same-path", "symlink", "hardlink")
                )
                for link_mode in link_modes:
                    case_directory = root / f"memory-input-{mode}-{link_mode}"
                    case_directory.mkdir()
                    input_name = "memUsage" if link_mode == "same-path" else "input"
                    input_path = case_directory / input_name
                    if mode == "mol-fallback":
                        input_path = input_path.with_suffix(".mol")
                        shutil.copy2(TEST_DIRECTORY / "butane.mol", input_path)
                    else:
                        input_path.write_text(
                            "abab\n"
                            if mode == "string"
                            else ("butane\n4\n1 2 2 3 3 4\nC C C C\n1 1 1\n")
                        )
                    memory_path = case_directory / "memUsage"
                    if link_mode == "symlink":
                        memory_path.symlink_to(input_path.name)
                    elif link_mode == "hardlink":
                        memory_path.hardlink_to(input_path)
                    before = input_path.read_bytes()
                    completed = run_cli_command(
                        executable,
                        [
                            input_name,
                            f"--run-strings={int(mode == 'string')}",
                            "--memory-report=1",
                        ],
                        case_directory,
                    )
                    require_cli(
                        completed.returncode == 1
                        and "would overwrite input" in completed.stderr,
                        f"memory output must reject {mode} {link_mode} collisions",
                        completed,
                    )
                    require_cli(
                        input_path.read_bytes() == before
                        and memory_path.read_bytes() == before,
                        f"memory output corrupted {mode} {link_mode} input",
                        completed,
                    )
                    require_cli(
                        not (case_directory / f"{input_name}Out").exists(),
                        "memory/input collisions should be detected before calculation",
                        completed,
                    )
                    scenarios += 1
    return scenarios


def run_cli_checks(executable: Path) -> int:
    """Exercise help, validation, aliases, input handling, and output flags."""
    scenarios = 0
    help_tokens = (
        "Usage:",
        "--runtime=<TICKS>",
        "--enum-max=<COUNT>",
        "--pathway=<0|1>",
        "--run-strings=<0|1>",
        "--accept-palindromes=<0|1>",
        "--parallel=<auto|on|off>",
        (
            "Select parallel search automatically, require it, or disable it. "
            "Default: off."
        ),
        "--threads=<auto|N>",
        "--remove-hydrogens=<0|1>",
        "--verbose=<0|1>",
        "--compensate-disjoint=<0|1>",
        "--memory-report=<0|1>",
        "--write-intermediate-mas=<0|1>",
        "Outputs:",
        "Legacy options:",
    )
    telemetry_supported: bool | None = None

    with tempfile.TemporaryDirectory(prefix="parallelassemblycpp-cli-") as directory:
        working_directory = Path(directory)

        for help_option in ("--help", "-h"):
            completed = run_cli_command(executable, [help_option], working_directory)
            require_cli(
                completed.returncode == 0,
                f"{help_option} should exit successfully",
                completed,
            )
            help_has_telemetry = "--telemetry=<0|1>" in completed.stdout
            if telemetry_supported is None:
                telemetry_supported = help_has_telemetry
            require_cli(
                help_has_telemetry == telemetry_supported,
                "help aliases disagree about telemetry support",
                completed,
            )
            expected_help_tokens = help_tokens + (
                ("--telemetry=<0|1>",) if telemetry_supported else ()
            )
            missing_tokens = [
                token for token in expected_help_tokens if token not in completed.stdout
            ]
            require_cli(
                not missing_tokens,
                f"{help_option} output is missing {missing_tokens}",
                completed,
            )
            require_cli(
                "pathwayFolder" not in completed.stdout,
                f"{help_option} still documents the removed pathwayFolder option",
                completed,
            )
            scenarios += 1

        completed = run_cli_command(executable, [], working_directory)
        require_cli(
            completed.returncode == 2,
            "a missing INPUT should be a command-line error",
            completed,
        )
        require_cli(
            "INPUT is required" in completed.stderr,
            "the missing-input error should explain what is required",
            completed,
        )
        scenarios += 1

        completed = run_cli_command(
            executable, ["missing-input-file"], working_directory
        )
        require_cli(
            completed.returncode == 1,
            "a missing input file should fail the calculation",
            completed,
        )
        require_cli(
            "input file not found" in completed.stderr,
            "a missing input file should produce a clear error",
            completed,
        )
        scenarios += 1

        compatibility_options = (
            "-runtime=1000000000",
            "--runTime=1000000000",
            "-runTime=1000000000",
            "-enumMax=1000000",
            "-runStrings=0",
            "--runStrings=0",
            "-acceptPalindromes=0",
            "-palindrome=0",
            "-removeHydrogens=0",
            "-compensateDisjoint=0",
            "-disjointCompensation=0",
            "-memTest=0",
            "-testMemory=0",
            "-writeIntermediateMAs=0",
        )
        for option in compatibility_options:
            completed = run_cli_command(
                executable, ["--help", option], working_directory
            )
            require_cli(
                completed.returncode == 0,
                f"compatibility option {option!r} should remain accepted",
                completed,
            )
            scenarios += 1

        valid_execution_options = (
            "--parallel=auto",
            "--parallel=on",
            "--parallel=off",
            "--threads=auto",
            "--threads=1",
            "--threads=2147483647",
        )
        for option in valid_execution_options:
            completed = run_cli_command(
                executable, ["--help", option], working_directory
            )
            require_cli(
                completed.returncode == 0,
                f"execution option {option!r} should be accepted",
                completed,
            )
            scenarios += 1

        invalid_cases = [
            (["input", "--does-not-exist=1"], "unknown option"),
            (["input", "--pathway"], "requires a value"),
            (["input", "--pathway="], "expected 0 or 1"),
            (["input", "--pathway=2"], "expected 0 or 1"),
            (["input", "--run-strings=2"], "expected 0 or 1"),
            (["input", "--accept-palindromes=yes"], "expected 0 or 1"),
            (["input", "--parallel"], "requires a value"),
            (["input", "--parallel="], "expected auto, on, or off"),
            (["input", "--parallel=ON"], "expected auto, on, or off"),
            (["input", "--parallel=1"], "expected auto, on, or off"),
            (["input", "--threads"], "requires a value"),
            (["input", "--threads="], "expected a non-negative integer"),
            (["input", "--threads=0"], "expected auto or an integer from 1"),
            (["input", "--threads=-1"], "expected a non-negative integer"),
            (["input", "--threads=2junk"], "expected a non-negative integer"),
            (
                ["input", "--threads=2147483648"],
                "expected auto or an integer from 1",
            ),
            (["input", "--remove-hydrogens=yes"], "expected 0 or 1"),
            (["input", "--verbose=2"], "expected 0 or 1"),
            (["input", "--enum-max=0"], "expected an integer from 1"),
            (["input", "--enum-max=12junk"], "non-negative integer"),
            (["input", "--runtime=-1"], "non-negative integer"),
            (["input", f"--runtime={'9' * 100}"], "non-negative integer"),
            (
                ["input", "--pathway=0", "--pathway=1"],
                "may be specified only once",
            ),
            (
                ["input", "--run-strings=0", "-runStrings=1"],
                "may be specified only once",
            ),
            (
                ["input", "--parallel=auto", "--parallel=off"],
                "may be specified only once",
            ),
            (
                ["input", "--threads=auto", "--threads=2"],
                "may be specified only once",
            ),
            (["first-input", "second-input"], "expected one INPUT"),
        ]
        if telemetry_supported:
            invalid_cases.append((["input", "--telemetry=2"], "expected 0 or 1"))
        else:
            invalid_cases.append((["input", "--telemetry=1"], "unknown option"))
        for arguments, error_text in invalid_cases:
            completed = run_cli_command(executable, arguments, working_directory)
            require_cli(
                completed.returncode == 2,
                f"invalid arguments {arguments!r} should be rejected",
                completed,
            )
            require_cli(
                error_text in completed.stderr,
                f"invalid arguments {arguments!r} should report {error_text!r}",
                completed,
            )
            scenarios += 1

        string_directory = working_directory / "string-assembly"
        string_directory.mkdir()
        string_input = string_directory / "strings.txt"
        string_cases = tuple(
            symbol_count * "0" + symbol_count * "1" + symbol_count * "2"
            for symbol_count in (5, 10, 15, 20, 25)
        )
        string_input.write_text("\n".join(string_cases) + "\n")
        completed = run_cli_command(
            executable,
            [
                string_input.name,
                "--run-strings=1",
                "--accept-palindromes=0",
                "--pathway=1",
            ],
            string_directory,
        )
        require_cli(
            completed.returncode == 0,
            "string assembly should process a line-oriented input",
            completed,
        )
        string_output = string_directory / "strings.txtOut"
        require_cli(
            string_output.is_file(),
            "string assembly did not create INPUTOut",
            completed,
        )
        string_output_text = string_output.read_text()
        string_indices = [
            int(match.group(1))
            for match in ASSEMBLY_INDEX_PATTERN.finditer(string_output_text)
        ]
        require_cli(
            string_indices == [11, 14, 17, 17, 20],
            "string assembly disagrees with the upstream reference corpus: "
            f"{string_indices}",
            completed,
        )
        require_cli(
            string_output_text.count("time elapsed:") == len(string_cases),
            "string assembly should record timing for every input line",
            completed,
        )
        for line_index, value in enumerate(string_cases):
            pathway_path = string_directory / f"strings.txt_{line_index}_Pathway"
            require_cli(
                pathway_path.is_file(),
                f"string line {line_index} did not create its pathway",
                completed,
            )
            pathway = json.loads(pathway_path.read_text())
            require_cli(
                pathway["file_graph"][0]["Fragments"] == [value]
                and isinstance(pathway["duplicates"], list),
                f"string line {line_index} pathway has the wrong structure",
                completed,
            )
        scenarios += 1

        line_endings_directory = working_directory / "string-line-endings"
        line_endings_directory.mkdir()
        line_endings_input = line_endings_directory / "input"
        # Keep the bytes explicit: the first two records use CRLF, the second
        # record is empty, and the final record deliberately has no newline.
        line_endings_input.write_bytes(b"abab\r\n\r\nabcdef")
        completed = run_cli_command(
            executable,
            ["input", "--run-strings=1", "--pathway=1"],
            line_endings_directory,
        )
        require_cli(
            completed.returncode == 0,
            "string assembly should accept blank lines, CRLF, and a final "
            "line without a newline",
            completed,
        )
        line_endings_output = line_endings_directory / "inputOut"
        require_cli(
            line_endings_output.is_file(),
            "mixed-line-ending string input did not create INPUTOut",
            completed,
        )
        line_endings_output_text = line_endings_output.read_text()
        line_endings_indices = [
            int(match.group(1))
            for match in ASSEMBLY_INDEX_PATTERN.finditer(line_endings_output_text)
        ]
        require_cli(
            line_endings_indices == [2, -1, 5],
            "mixed-line-ending string records were not preserved: "
            f"{line_endings_indices}",
            completed,
        )
        require_cli(
            line_endings_output_text.count("time elapsed:") == 3,
            "mixed-line-ending string input should produce three timing records",
            completed,
        )
        for line_index, value in enumerate(("abab", "", "abcdef")):
            pathway_path = line_endings_directory / f"input_{line_index}_Pathway"
            require_cli(
                pathway_path.is_file(),
                f"mixed-line-ending string record {line_index} has no pathway",
                completed,
            )
            pathway = json.loads(pathway_path.read_text())
            require_cli(
                pathway["file_graph"][0]["Fragments"] == [value],
                f"mixed-line-ending pathway {line_index} changed its record",
                completed,
            )
        require_cli(
            not (line_endings_directory / "input_3_Pathway").exists(),
            "mixed-line-ending input produced a spurious fourth record",
            completed,
        )
        scenarios += 1

        # These features are unavailable in string mode on this serial target.
        # Keep cases isolated so one rejection cannot mask another branch.
        incompatible_string_options = [
            (
                "parallel",
                "--parallel=on",
                "this executable was built without parallel support",
            ),
            (
                "intermediate-indices",
                "--write-intermediate-mas=1",
                "--write-intermediate-mas is unavailable for string assembly",
            ),
        ]
        if telemetry_supported:
            incompatible_string_options.append(
                (
                    "telemetry",
                    "--telemetry=1",
                    "--telemetry is unavailable for string assembly",
                )
            )
        for case_name, option, error_text in incompatible_string_options:
            incompatible_directory = (
                working_directory / f"string-incompatible-{case_name}"
            )
            incompatible_directory.mkdir()
            (incompatible_directory / "input").write_text("abab\n")
            completed = run_cli_command(
                executable,
                ["input", "--run-strings=1", option],
                incompatible_directory,
            )
            require_cli(
                completed.returncode == 1,
                f"string mode should reject {option}",
                completed,
            )
            require_cli(
                error_text in completed.stderr,
                f"string-mode rejection for {option} should explain the conflict",
                completed,
            )
            require_cli(
                not (incompatible_directory / "inputOut").exists(),
                f"string-mode rejection for {option} should not create output",
                completed,
            )
            scenarios += 1

        no_pathway_directory = working_directory / "string-no-pathway"
        no_pathway_directory.mkdir()
        no_pathway_input = no_pathway_directory / "input"
        no_pathway_input.write_text("abab\n")
        completed = run_cli_command(
            executable,
            ["input", "-runStrings=1", "--pathway=0"],
            no_pathway_directory,
        )
        require_cli(
            completed.returncode == 0,
            "the legacy string flag should run successfully",
            completed,
        )
        require_cli(
            read_first_line_assembly_index(no_pathway_directory / "inputOut") == 2,
            "legacy string mode calculated the wrong index for 'abab'",
            completed,
        )
        require_cli(
            not (no_pathway_directory / "input_0_Pathway").exists(),
            "--pathway=0 should suppress string pathway output",
            completed,
        )
        scenarios += 1

        reversal_directory = working_directory / "string-reversal"
        reversal_directory.mkdir()
        (reversal_directory / "input").write_text("abcxcba\n")
        reversal_indices: list[int | None] = []
        for reversal_option in ("--accept-palindromes=0", "-acceptPalindromes=1"):
            completed = run_cli_command(
                executable,
                [
                    "input",
                    "--run-strings=1",
                    reversal_option,
                    "--pathway=0",
                ],
                reversal_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"string reversal option {reversal_option!r} should succeed",
                completed,
            )
            reversal_indices.append(
                read_first_line_assembly_index(reversal_directory / "inputOut")
            )
        require_cli(
            reversal_indices == [6, 4],
            "--accept-palindromes did not enable reversal equivalence: "
            f"{reversal_indices}",
        )
        scenarios += 1

        unique_string_directory = working_directory / "string-unique-pathway"
        unique_string_directory.mkdir()
        (unique_string_directory / "input").write_text("abcdef\n")
        completed = run_cli_command(
            executable,
            ["input", "--run-strings=1", "--pathway=1"],
            unique_string_directory,
        )
        require_cli(
            completed.returncode == 0,
            "a unique string should still produce a valid pathway document",
            completed,
        )
        unique_pathway = json.loads(
            (unique_string_directory / "input_0_Pathway").read_text()
        )
        require_cli(
            unique_pathway["remnant"][0] == {"Fragments": ["abcdef"], "Positions": [0]}
            and unique_pathway["duplicates"] == [],
            "a no-copy pathway should preserve the whole string as its remnant",
            completed,
        )
        scenarios += 1

        source = TEST_DIRECTORY / "butane.mol"
        input_path = working_directory / "input.mol"
        shutil.copy2(source, input_path)
        canonical_options = [
            "--runtime=1000000000",
            "--enum-max=1000000",
            "--pathway=0",
            "--parallel=off",
            "--threads=2",
            "--remove-hydrogens=0",
            "--verbose=0",
            "--compensate-disjoint=1",
            "--memory-report=0",
            "--write-intermediate-mas=1",
        ]
        if telemetry_supported:
            canonical_options.append("--telemetry=0")
        completed = run_cli_command(
            executable,
            [*canonical_options, input_path.name],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "canonical options before a .mol input should run successfully",
            completed,
        )
        require_cli(
            (working_directory / "inputOut").is_file(),
            "a .mol input should create output without .mol in the output name",
            completed,
        )
        require_cli(
            (working_directory / "inputIntermediateMAs").is_file(),
            "--write-intermediate-mas=1 should create its output",
            completed,
        )
        require_cli(
            not (working_directory / "inputPathway").exists(),
            "--pathway=0 should suppress pathway output",
            completed,
        )
        require_cli(
            not (working_directory / "memUsage").exists(),
            "--memory-report=0 should suppress memory output",
            completed,
        )
        require_cli(
            not (working_directory / "inputTelemetry.json").exists(),
            "--telemetry=0 should suppress telemetry output",
            completed,
        )
        require_cli(
            "Graph:" not in completed.stdout and "  Atom " not in completed.stdout,
            "quiet mode should suppress the parsed graph dump",
            completed,
        )
        scenarios += 1

        extension_directory = working_directory / "molfile-extension-coverage"
        extension_directory.mkdir()
        for case_index, extension in enumerate(
            (".MOL", ".mOl", ".sdf", ".SDF", ".sDf"),
            start=1,
        ):
            case_directory = extension_directory / f"case-{case_index}"
            case_directory.mkdir()
            input_name = f"input{extension}"
            extension_input_path = case_directory / input_name
            if extension.lower() == ".sdf":
                extension_input_path.write_text(
                    source.read_text()
                    + "$$$$\n"
                    + (TEST_DIRECTORY / "alanine.mol").read_text()
                    + "$$$$\n"
                )
            else:
                shutil.copy2(source, extension_input_path)
            completed = run_cli_command(
                executable,
                [
                    input_name,
                    "--pathway=0",
                    "--parallel=off",
                    "--verbose=1",
                    "--remove-hydrogens=0",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"a {extension} input should run successfully",
                completed,
            )
            require_cli(
                f"Input: {input_name}\n" in completed.stdout
                and "Molfile: 4 atoms, 3 bonds" in completed.stdout,
                f"a {extension} input should use MOL parsing",
                completed,
            )
            require_cli(
                (case_directory / "inputOut").is_file(),
                f"a {extension} input should omit its suffix from output names",
                completed,
            )
            require_cli(
                read_first_line_assembly_index(case_directory / "inputOut") == 2,
                f"a {extension} input should use only its first record",
                completed,
            )
            require_cli(
                not (case_directory / f"{input_name}Out").exists(),
                f"a {extension} input should not retain its suffix in output names",
                completed,
            )
            scenarios += 1

        native_suffix_directory = extension_directory / "native-final-suffix"
        native_suffix_directory.mkdir()
        native_suffix_name = "input.sdf.txt"
        shutil.copy2(
            TEST_DIRECTORY / "graphio_test",
            native_suffix_directory / native_suffix_name,
        )
        completed = run_cli_command(
            executable,
            [native_suffix_name, "--pathway=0", "--verbose=1"],
            native_suffix_directory,
        )
        require_cli(
            completed.returncode == 0,
            "a native graph with non-final .sdf text should run successfully",
            completed,
        )
        require_cli(
            "Molfile:" not in completed.stdout and "Graph:" in completed.stdout,
            "only a final MOL/SDF suffix should select MOL parsing",
            completed,
        )
        require_cli(
            (native_suffix_directory / f"{native_suffix_name}Out").is_file(),
            "a native graph should retain its full filename in output names",
            completed,
        )
        scenarios += 1

        format_parity_directory = working_directory / "input-format-parity"
        mol_directory = format_parity_directory / "mol"
        native_directory = format_parity_directory / "native"
        mol_directory.mkdir(parents=True)
        native_directory.mkdir(parents=True)
        shutil.copy2(TEST_DIRECTORY / "alanine.mol", mol_directory / "input.mol")
        (native_directory / "input").write_text(
            "alanine native graph\n6\n1 2 2 3 3 4 3 5 2 6\nN C C O O C\n1 1 2 1 1\n"
        )

        equivalent_results: list[tuple[int | None, dict[str, object]]] = []
        for format_name, input_name, case_directory in (
            ("MOL", "input.mol", mol_directory),
            ("native graph", "input", native_directory),
        ):
            completed = run_cli_command(
                executable,
                [
                    input_name,
                    "--pathway=1",
                    "--parallel=off",
                    "--remove-hydrogens=0",
                    "--memory-report=0",
                    "--write-intermediate-mas=0",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the equivalent {format_name} input should run successfully",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"the equivalent {format_name} input omitted its pathway",
                completed,
            )
            assembly_index = read_first_line_assembly_index(case_directory / "inputOut")
            require_cli(
                assembly_index is not None,
                f"the equivalent {format_name} input omitted its assembly index",
                completed,
            )
            equivalent_results.append(
                (
                    assembly_index,
                    parse_pathway_document(pathway_path),
                )
            )

        require_cli(
            equivalent_results[0][0] == equivalent_results[1][0],
            "equivalent MOL and native graph inputs should produce the same "
            "assembly index",
        )
        require_cli(
            equivalent_results[0][1] == equivalent_results[1][1],
            "equivalent MOL and native graph inputs should produce the same pathway",
        )
        scenarios += 1

        precedence_directory = working_directory / "native-input-precedence"
        precedence_directory.mkdir()
        shutil.copy2(
            TEST_DIRECTORY / "graphio_test",
            precedence_directory / "input",
        )
        shutil.copy2(
            TEST_DIRECTORY / "tridecane.mol",
            precedence_directory / "input.mol",
        )
        completed = run_cli_command(
            executable,
            ["input", "--pathway=0", "--verbose=1"],
            precedence_directory,
        )
        require_cli(
            completed.returncode == 0,
            "an exact native input with a .mol sibling should run successfully",
            completed,
        )
        require_cli(
            read_first_line_assembly_index(precedence_directory / "inputOut") == 5,
            "an exact native input should take precedence over its .mol sibling",
            completed,
        )
        require_cli(
            "Input: input\n" in completed.stdout
            and "Input: input.mol\n" not in completed.stdout,
            "verbose output should identify the exact native input",
            completed,
        )
        scenarios += 1

        completed = run_cli_command(
            executable,
            ["--pathway=0", "--verbose=1", input_path.name],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "--verbose=1 should run successfully",
            completed,
        )
        require_cli(
            "Input: input.mol" in completed.stdout
            and "Molfile: 4 atoms, 3 bonds" in completed.stdout
            and "Graph: 4 atoms, 3 bonds" in completed.stdout
            and "  Atom 1 (C):" in completed.stdout,
            "--verbose=1 should print the input summary and parsed graph",
            completed,
        )
        scenarios += 1

        parity_indices: list[int | None] = []
        for pathway_enabled in (False, True):
            mode = int(pathway_enabled)
            case_directory = working_directory / f"pathway-parity-{mode}"
            case_directory.mkdir()
            shutil.copy2(
                TEST_DIRECTORY / "ketoconazole.mol",
                case_directory / "input.mol",
            )
            parity_options = [
                "input.mol",
                f"--pathway={mode}",
                "--memory-report=0",
                "--write-intermediate-mas=0",
            ]
            if telemetry_supported:
                parity_options.append("--telemetry=0")
            completed = run_cli_command(
                executable,
                parity_options,
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"ketoconazole pathway={mode} parity run should succeed",
                completed,
            )
            parity_indices.append(
                read_first_line_assembly_index(case_directory / "inputOut")
            )
            require_cli(
                (case_directory / "inputPathway").is_file() == pathway_enabled,
                f"ketoconazole pathway={mode} created unexpected pathway output",
                completed,
            )
        require_cli(
            parity_indices == [22, 22],
            "ketoconazole should have assembly index 22 in both pathway modes",
        )
        scenarios += 1

        dash_input = working_directory / "-dash-input.mol"
        shutil.copy2(source, dash_input)
        completed = run_cli_command(
            executable,
            ["--pathway=0", "--", dash_input.name],
            working_directory,
        )
        require_cli(
            completed.returncode == 0,
            "-- should allow an input name beginning with a dash",
            completed,
        )
        require_cli(
            (working_directory / "-dash-inputOut").is_file(),
            "the dash-prefixed input should create the expected output",
            completed,
        )
        scenarios += 1

        cutoff_cases = (
            (
                "runtime-limit",
                TEST_DIRECTORY / "butane.mol",
                "--runtime=0",
                "status: runtime limit reached",
            ),
            # Cyclosporin exercises two active words in both mask domains.
            (
                "wide-runtime-limit",
                TEST_DIRECTORY / "cyclosporin.mol",
                "--runtime=0",
                "status: runtime limit reached",
            ),
            (
                "enumeration-limit",
                TEST_DIRECTORY / "113.mol",
                "--enum-max=1",
                "status: enumeration limit reached",
            ),
        )
        for name, source_path, limit_option, expected_status in cutoff_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            shutil.copy2(source_path, case_directory / "input.mol")
            completed = run_cli_command(
                executable,
                ["input.mol", "--pathway=0", limit_option],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"{name} scenario should return its best result successfully",
                completed,
            )
            output_path = case_directory / "inputOut"
            index = read_first_line_assembly_index(output_path)
            require_cli(
                index is not None and index != 2_147_483_647,
                f"{name} scenario should put a finite assembly index on the first line",
                completed,
            )
            require_cli(
                expected_status in output_path.read_text().splitlines(),
                f"{name} scenario should record {expected_status!r}",
                completed,
            )
            scenarios += 1

        explicit_hydrogen_mol = (
            "Explicit hydrogens\n"
            "ParallelAssemblyCpp CLI test\n"
            "\n"
            "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    0.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    3.0000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  3  1  0  0  0  0\n"
            "  2  3  1  0  0  0  0\n"
            "  3  4  1  0  0  0  0\n"
            "  4  5  1  0  0  0  0\n"
            "  4  6  1  0  0  0  0\n"
            "M  END\n"
        )
        hydrogen_cases = (
            ("hydrogens-default", [], ["C", "C"], 1),
            (
                "hydrogens-on",
                ["--remove-hydrogens=1"],
                ["C", "C"],
                1,
            ),
            (
                "hydrogens-off",
                ["--remove-hydrogens=0"],
                ["H", "H", "C", "C", "H", "H"],
                5,
            ),
        )
        molecular_hydrogen_graphs = {}
        for (
            name,
            hydrogen_options,
            expected_colours,
            expected_edge_count,
        ) in hydrogen_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input.mol").write_text(explicit_hydrogen_mol)
            completed = run_cli_command(
                executable,
                ["input.mol", "--pathway=1", *hydrogen_options],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"explicit-hydrogen scenario {name!r} should succeed",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"explicit-hydrogen scenario {name!r} omitted its pathway",
                completed,
            )
            pathway = parse_pathway_document(pathway_path)
            graph = pathway["file_graph"][0]
            require_cli(
                graph["VertexColours"] == expected_colours
                and len(graph["Edges"]) == expected_edge_count,
                f"explicit-hydrogen scenario {name!r} transformed the wrong "
                "atoms or bonds",
                completed,
            )
            molecular_hydrogen_graphs[tuple(hydrogen_options)] = graph
            scenarios += 1

        native_hydrogen_graph = (
            "native-hydrogens\n6\n1 3 2 3 3 4 4 5 4 6\nH H C C H H\n1 1 1 1 1\n"
        )
        native_hydrogen_cases = (
            ("native-hydrogens-default", [], ["C", "C"], 1),
            (
                "native-hydrogens-on",
                ["--remove-hydrogens=1"],
                ["C", "C"],
                1,
            ),
            (
                "native-hydrogens-off",
                ["--remove-hydrogens=0"],
                ["H", "H", "C", "C", "H", "H"],
                5,
            ),
        )
        for (
            name,
            hydrogen_options,
            expected_colours,
            expected_edge_count,
        ) in native_hydrogen_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input").write_text(native_hydrogen_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=1", *hydrogen_options],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"native-graph hydrogen scenario {name!r} should succeed",
                completed,
            )
            pathway_path = case_directory / "inputPathway"
            require_cli(
                pathway_path.is_file(),
                f"native-graph hydrogen scenario {name!r} omitted its pathway",
                completed,
            )
            graph = parse_pathway_document(pathway_path)["file_graph"][0]
            require_cli(
                graph["VertexColours"] == expected_colours
                and len(graph["Edges"]) == expected_edge_count,
                f"native-graph hydrogen scenario {name!r} transformed the "
                "wrong atoms or bonds",
                completed,
            )
            require_cli(
                graph == molecular_hydrogen_graphs[tuple(hydrogen_options)],
                f"native-graph hydrogen scenario {name!r} should match its "
                "MOL equivalent",
                completed,
            )
            scenarios += 1

        all_hydrogen_mol = (
            "Hydrogen\n"
            "ParallelAssemblyCpp CLI test\n"
            "\n"
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "M  END\n"
        )
        empty_graph_results: list[tuple[int, int]] = []
        for name, compensation_option in (
            ("empty-compensation-off", "--compensate-disjoint=0"),
            ("empty-compensation-on", "--compensate-disjoint=1"),
        ):
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input.mol").write_text(all_hydrogen_mol)
            completed = run_cli_command(
                executable,
                [
                    "input.mol",
                    "--pathway=0",
                    "--remove-hydrogens=1",
                    compensation_option,
                    "--write-intermediate-mas=1",
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"empty-graph compensation scenario {name!r} should succeed",
                completed,
            )
            final_index = read_first_line_assembly_index(case_directory / "inputOut")
            intermediate_index = read_last_intermediate_index(
                case_directory / "inputIntermediateMAs"
            )
            require_cli(
                final_index is not None and intermediate_index == final_index,
                f"empty-graph compensation scenario {name!r} should report "
                "consistent indices",
                completed,
            )
            empty_graph_results.append((final_index, intermediate_index))
            scenarios += 1
        require_cli(
            empty_graph_results[0] == empty_graph_results[1],
            "disjoint compensation should not change an empty processed graph",
        )

        disconnected_graph = "disconnected\n4\n1 2 3 4\nC C C C\n1 1\n"
        compensation_cases = (
            (
                "canonical-off",
                "--compensate-disjoint=0",
                "--write-intermediate-mas=1",
                1,
            ),
            (
                "canonical-on",
                "--compensate-disjoint=1",
                "--write-intermediate-mas=1",
                0,
            ),
            (
                "legacy-parser",
                "-compensateDisjoint=1",
                "-writeIntermediateMAs=1",
                0,
            ),
            (
                "legacy-docs",
                "-disjointCompensation=1",
                "--write-intermediate-mas=1",
                0,
            ),
        )
        for (
            name,
            compensation_option,
            intermediate_option,
            expected_index,
        ) in compensation_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            (case_directory / "input").write_text(disconnected_graph)
            completed = run_cli_command(
                executable,
                [
                    "input",
                    "-pathway=0",
                    compensation_option,
                    intermediate_option,
                ],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"disjoint-compensation scenario {name!r} should succeed",
                completed,
            )

            output_path = case_directory / "inputOut"
            intermediate_path = case_directory / "inputIntermediateMAs"
            require_cli(
                output_path.is_file() and intermediate_path.is_file(),
                f"disjoint-compensation scenario {name!r} omitted an output file",
                completed,
            )
            output_match = ASSEMBLY_INDEX_PATTERN.search(output_path.read_text())
            require_cli(
                output_match is not None
                and int(output_match.group(1)) == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong "
                "final index",
                completed,
            )
            last_intermediate = read_last_intermediate_index(intermediate_path)
            require_cli(
                last_intermediate == expected_index,
                f"disjoint-compensation scenario {name!r} returned the wrong "
                "intermediate index",
                completed,
            )
            require_cli(
                not (case_directory / "inputPathway").exists(),
                "legacy -pathway=0 should suppress pathway output",
                completed,
            )
            scenarios += 1

        for enum_limit, expect_limit_status in ((1, True), (2, False)):
            case_directory = working_directory / f"enum-boundary-{enum_limit}"
            case_directory.mkdir()
            (case_directory / "input").write_text(disconnected_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=0", f"--enum-max={enum_limit}"],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"enum boundary {enum_limit} should return a best result",
                completed,
            )
            output_path = case_directory / "inputOut"
            status_present = (
                output_path.is_file()
                and "status: enumeration limit reached"
                in output_path.read_text().splitlines()
            )
            require_cli(
                status_present == expect_limit_status,
                f"enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            scenarios += 1

        connected_graph = "connected\n3\n1 2 2 3\nC C C\n1 1\n"
        for enum_limit, expect_limit_status in ((2, True), (3, False)):
            case_directory = working_directory / f"enum-connected-boundary-{enum_limit}"
            case_directory.mkdir()
            (case_directory / "input").write_text(connected_graph)
            completed = run_cli_command(
                executable,
                ["input", "--pathway=0", f"--enum-max={enum_limit}"],
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"connected enum boundary {enum_limit} should return a best result",
                completed,
            )
            output_path = case_directory / "inputOut"
            status_present = (
                output_path.is_file()
                and "status: enumeration limit reached"
                in output_path.read_text().splitlines()
            )
            require_cli(
                status_present == expect_limit_status,
                f"connected enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            scenarios += 1

        for enum_limit, expected_index, expect_limit_status in (
            (86, 12, True),
            (87, 8, False),
        ):
            case_directory = working_directory / f"enum-deep-boundary-{enum_limit}"
            case_directory.mkdir()
            shutil.copy2(
                TEST_DIRECTORY / "113.mol",
                case_directory / "input.mol",
            )
            deep_options = [
                "input.mol",
                "--pathway=0",
                f"--enum-max={enum_limit}",
            ]
            if telemetry_supported:
                deep_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                deep_options,
                case_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"deep enum boundary {enum_limit} should return a best result",
                completed,
            )
            output_path = case_directory / "inputOut"
            output_lines = output_path.read_text().splitlines()
            require_cli(
                read_first_line_assembly_index(output_path) == expected_index,
                f"deep enum boundary {enum_limit} returned the wrong index",
                completed,
            )
            require_cli(
                ("status: enumeration limit reached" in output_lines)
                == expect_limit_status,
                f"deep enum boundary {enum_limit} enforced the wrong state cap",
                completed,
            )
            if telemetry_supported:
                phases = json.loads(
                    (case_directory / "inputTelemetry.json").read_text()
                )["memory"]["phases"]
                expected_later_activations = 0 if expect_limit_status else 1
                require_cli(
                    phases["dag_conversion"]["activations"]
                    == expected_later_activations
                    and phases["assembly_search"]["activations"]
                    == expected_later_activations,
                    f"deep enum boundary {enum_limit} reported incorrect phase "
                    "activity",
                    completed,
                )
            scenarios += 1

        # Exercise scalar, wide, adaptive, later-word, and pre-fragment
        # equivalence-quotient paths.
        for edge_count, expected_index, active_words, cache_outcome in (
            (64, 6, 1, "scalar-lookups"),
            (65, 7, 2, "equivalence-quotient"),
            (127, 10, 2, "wide-hits"),
            (128, 7, 2, "adaptive-fallback"),
            (129, 8, 3, "adaptive-fallback"),
        ):
            wide_graph = "\n".join(
                (
                    "wide-path",
                    str(edge_count + 1),
                    " ".join(
                        f"{vertex} {vertex + 1}" for vertex in range(1, edge_count + 1)
                    ),
                    " ".join("C" for _ in range(edge_count + 1)),
                    " ".join("1" for _ in range(edge_count)),
                    "",
                )
            )
            wide_directory = working_directory / f"wide-initial-dag-{edge_count}"
            wide_directory.mkdir()
            (wide_directory / "input").write_text(wide_graph)
            wide_options = ["input", "--pathway=0"]
            if telemetry_supported:
                wide_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                wide_options,
                wide_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the {edge_count}-edge initial DAG scenario should succeed",
                completed,
            )
            require_cli(
                read_first_line_assembly_index(wide_directory / "inputOut")
                == expected_index,
                f"the {edge_count}-edge initial DAG scenario returned the wrong index",
                completed,
            )
            if not telemetry_supported:
                scenarios += 1
                continue
            telemetry = json.loads((wide_directory / "inputTelemetry.json").read_text())
            counters = telemetry["counters"]
            graph = telemetry["processed_graph"]
            residual = telemetry["caches"]["residual_decomposition"]
            canonical = telemetry["caches"]["canonical_mask"]
            require_cli(
                graph["edges"] == edge_count
                and graph["active_mask_words"] == active_words,
                f"the {edge_count}-edge telemetry reported the wrong mask width",
                completed,
            )
            require_cli(
                residual["eligible_for_processed_graph"] is True,
                f"the {edge_count}-edge telemetry reported the wrong cache eligibility",
                completed,
            )
            require_cli(
                counters["retained_mask_attempts"]
                == counters["retained_masks"]
                + counters["duplicate_mask_attempts"]
                + counters["rejected_masks"],
                f"the {edge_count}-edge retained-mask counters are inconsistent",
                completed,
            )
            require_cli(
                counters["canonicalisation_calls"]
                == canonical["hits"] + canonical["misses"],
                f"the {edge_count}-edge canonical counters are inconsistent",
                completed,
            )
            require_cli(
                residual["lookups"] == residual["hits"] + residual["misses"],
                f"the {edge_count}-edge residual-cache counters are inconsistent",
                completed,
            )
            require_cli(
                counters["retained_masks"] > 0
                and counters["matching_visits"] > 0
                and counters["canonicalisation_calls"] > 0,
                f"the {edge_count}-edge telemetry did not exercise search counters",
                completed,
            )
            require_cli(
                residual["requests"] > 0
                and residual["eligible_requests"] == residual["requests"]
                and residual["small_molecule_bypasses"] == 0
                and residual["wide_molecule_bypasses"] == 0
                and residual["eligible_requests"]
                == residual["small_residual_bypasses"]
                + residual["first_occurrence_bypasses"]
                + residual["runtime_disabled_bypasses"]
                + residual["lookups"],
                f"the {edge_count}-edge case did not exercise an eligible cache path",
                completed,
            )
            if cache_outcome == "scalar-lookups":
                require_cli(
                    residual["lookups"] > 0
                    and residual["admissions"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not exercise scalar caching",
                    completed,
                )
            elif cache_outcome == "wide-hits":
                require_cli(
                    residual["lookups"] > 0
                    and residual["hits"] > 0
                    and residual["admissions"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not exercise wide cache hits",
                    completed,
                )
            elif cache_outcome == "equivalence-quotient":
                require_cli(
                    counters["matching_visits"] < 5000
                    and residual["first_occurrence_bypasses"] > 0
                    and residual["runtime_disabled_bypasses"] == 0,
                    f"the {edge_count}-edge case did not reduce equivalent "
                    "matchings before adaptive fallback",
                    completed,
                )
            else:
                require_cli(
                    residual["first_occurrence_bypasses"] > 0
                    and residual["runtime_disabled_bypasses"] > 0,
                    f"the {edge_count}-edge case did not exercise adaptive fallback",
                    completed,
                )
            phases = telemetry["memory"]["phases"]
            require_cli(
                set(phases)
                == {
                    "input_setup",
                    "initial_enumeration",
                    "dag_conversion",
                    "assembly_search",
                    "output",
                },
                f"the {edge_count}-edge telemetry omitted a search phase",
                completed,
            )
            if sys.platform.startswith("linux"):
                require_cli(
                    all(
                        phase["peak_rss_kib"] is None
                        or (
                            phase["peak_rss_kib"] >= phase["start_rss_kib"]
                            and phase["peak_rss_kib"] >= phase["end_rss_kib"]
                        )
                        for phase in phases.values()
                    ),
                    f"the {edge_count}-edge phase RSS peaks are inconsistent",
                    completed,
                )
            scenarios += 1

        # Cross the former fixed 512-bit cap in both mask domains without
        # introducing a large connected-subgraph search. Each component is a
        # small cycle with its own atom colour; the enumeration cap is reached
        # only after every one-edge EdgeMask has been retained and the first
        # high-index AtomMask seed has expanded once.
        for atom_count, edge_count, component_sizes in (
            (512, 512, [4] * 128),
            (513, 513, [4] * 127 + [5]),
        ):
            capacity_directory = (
                working_directory / f"dynamic-mask-capacity-{atom_count}a-{edge_count}e"
            )
            capacity_directory.mkdir()
            capacity_graph = make_mask_capacity_graph(component_sizes)
            (capacity_directory / "input").write_text(capacity_graph)
            capacity_options = [
                "input",
                "--pathway=0",
                f"--enum-max={edge_count + 1}",
            ]
            if telemetry_supported:
                capacity_options.append("--telemetry=1")
            completed = run_cli_command(
                executable,
                capacity_options,
                capacity_directory,
            )
            require_cli(
                completed.returncode == 0,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "should succeed",
                completed,
            )
            output_path = capacity_directory / "inputOut"
            output_lines = output_path.read_text().splitlines()
            require_cli(
                read_first_line_assembly_index(output_path) == edge_count - 1,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "returned the wrong bounded-search index",
                completed,
            )
            require_cli(
                "status: enumeration limit reached" in output_lines,
                f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                "did not reach its deterministic enumeration boundary",
                completed,
            )
            if telemetry_supported:
                telemetry = json.loads(
                    (capacity_directory / "inputTelemetry.json").read_text()
                )
                graph = telemetry["processed_graph"]
                counters = telemetry["counters"]
                require_cli(
                    graph
                    == {
                        "atoms": atom_count,
                        "edges": edge_count,
                        "active_mask_words": (edge_count + 63) // 64,
                    },
                    f"the {atom_count}-atom/{edge_count}-edge capacity telemetry "
                    "reported the wrong dynamic mask dimensions",
                    completed,
                )
                require_cli(
                    counters["retained_masks"] == edge_count + 1
                    and counters["rejected_masks"] == 1,
                    f"the {atom_count}-atom/{edge_count}-edge capacity scenario "
                    "did not traverse every singleton mask and the first child",
                    completed,
                )
            scenarios += 1

        output_failure_cases = [
            (
                "pathway-output-failure",
                "inputPathway",
                ["--pathway=1", "--write-intermediate-mas=0", "--memory-report=0"],
            ),
            (
                "intermediate-output-failure",
                "inputIntermediateMAs",
                ["--pathway=0", "--write-intermediate-mas=1", "--memory-report=0"],
            ),
        ]
        if telemetry_supported:
            output_failure_cases.append(
                (
                    "telemetry-output-failure",
                    "inputTelemetry.json",
                    ["--pathway=0", "--telemetry=1", "--memory-report=0"],
                )
            )
        if sys.platform.startswith("linux"):
            output_failure_cases.append(
                (
                    "memory-output-failure",
                    "memUsage",
                    [
                        "--pathway=0",
                        "--write-intermediate-mas=0",
                        "--memory-report=1",
                    ],
                )
            )

        for name, target_name, output_options in output_failure_cases:
            case_directory = working_directory / name
            case_directory.mkdir()
            shutil.copy2(source, case_directory / "input.mol")
            (case_directory / target_name).mkdir()
            completed = run_cli_command(
                executable,
                ["input.mol", *output_options],
                case_directory,
            )
            require_cli(
                completed.returncode != 0,
                f"{name} should fail when {target_name!r} cannot be opened",
                completed,
            )
            require_cli(
                f"could not open output file '{target_name}'" in completed.stderr,
                f"{name} should identify the output it could not open",
                completed,
            )
            scenarios += 1

        disabled_output_directory = working_directory / "disabled-output-sentinels"
        disabled_output_directory.mkdir()
        shutil.copy2(source, disabled_output_directory / "input.mol")
        sentinels = {
            "inputPathway": "existing pathway sentinel\n",
            "inputIntermediateMAs": "existing intermediate sentinel\n",
            "memUsage": "existing memory sentinel\n",
            "inputTelemetry.json": "existing telemetry sentinel\n",
        }
        for filename, content in sentinels.items():
            (disabled_output_directory / filename).write_text(content)
        disabled_options = [
            "input.mol",
            "--pathway=0",
            "--write-intermediate-mas=0",
            "--memory-report=0",
        ]
        if telemetry_supported:
            disabled_options.append("--telemetry=0")
        completed = run_cli_command(
            executable,
            disabled_options,
            disabled_output_directory,
        )
        require_cli(
            completed.returncode == 0,
            "disabled output flags should not obstruct a successful calculation",
            completed,
        )
        for filename, content in sentinels.items():
            require_cli(
                (disabled_output_directory / filename).read_text() == content,
                f"disabled output flag unexpectedly overwrote {filename!r}",
                completed,
            )
        scenarios += 1

    for memory_option, expected_on_linux in (
        (None, False),
        ("--memory-report=0", False),
        ("--memory-report=1", True),
        ("-memTest=1", True),
        ("-testMemory=1", True),
    ):
        with tempfile.TemporaryDirectory(
            prefix="parallelassemblycpp-memory-"
        ) as directory:
            working_directory = Path(directory)
            shutil.copy2(TEST_DIRECTORY / "butane.mol", working_directory / "input.mol")
            arguments = ["input.mol", "--pathway=0"]
            if memory_option is not None:
                arguments.append(memory_option)
            completed = run_cli_command(executable, arguments, working_directory)
            if expected_on_linux and not sys.platform.startswith("linux"):
                require_cli(
                    completed.returncode == 2 and "only on Linux" in completed.stderr,
                    f"memory scenario {memory_option} should explain platform support",
                    completed,
                )
                scenarios += 1
                continue
            require_cli(
                completed.returncode == 0,
                f"memory scenario {memory_option or 'default'} should succeed",
                completed,
            )

            memory_path = working_directory / "memUsage"
            should_exist = expected_on_linux and sys.platform.startswith("linux")
            require_cli(
                memory_path.exists() == should_exist,
                f"memory scenario {memory_option or 'default'} created an "
                "unexpected report",
                completed,
            )
            if should_exist:
                require_cli(
                    re.fullmatch(r"VmPeak:\s+\d+\s+kB\n?", memory_path.read_text())
                    is not None,
                    "the Linux memory report should contain a VmPeak value in kB",
                    completed,
                )
            scenarios += 1

    scenarios += run_string_record_checks(executable)
    scenarios += run_string_output_alias_checks(executable)
    scenarios += run_graph_output_alias_checks(executable, bool(telemetry_supported))
    scenarios += run_flag_matrix_checks(executable, bool(telemetry_supported))
    scenarios += run_input_output_matrix_checks(executable, bool(telemetry_supported))
    return scenarios


def compiler_command(compiler: str) -> list[str]:
    compiler_command = shlex.split(compiler)
    if not compiler_command:
        raise TestConfigurationError("the compiler command is empty")
    if shutil.which(compiler_command[0]) is None:
        raise TestConfigurationError(f"compiler not found: {compiler_command[0]}")
    return compiler_command


def run_cpp_unit_test(
    compiler: str,
    source_stem: str,
    label: str,
    *,
    stack_limit_bytes: int | None = None,
) -> None:
    """Build and run one standalone C++ unit-test executable."""
    command_prefix = compiler_command(compiler)
    prefix = f"parallelassemblycpp-{label.replace(' ', '-')}-tests-"
    with tempfile.TemporaryDirectory(prefix=prefix) as directory:
        test_executable = Path(directory) / source_stem
        command = [
            *command_prefix,
            str(TEST_DIRECTORY / f"{source_stem}.cpp"),
            "-std=c++20",
            "-O2",
            "-mpopcnt",
            "-march=x86-64-v3",
            "-o",
            str(test_executable),
        ]
        print(f"Building {label} tests: {shlex.join(command)}", flush=True)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise TestConfigurationError(
                f"{label} test build failed with exit code {completed.returncode}"
            )

        run_options: dict[str, object] = {}
        if stack_limit_bytes is not None and sys.platform.startswith("linux"):
            # resource is unavailable on Windows, so keep this import platform-local.
            import resource  # noqa: PLC0415

            def limit_stack() -> None:
                resource.setrlimit(
                    resource.RLIMIT_STACK,
                    (stack_limit_bytes, stack_limit_bytes),
                )

            run_options["preexec_fn"] = limit_stack
        completed = subprocess.run([str(test_executable)], check=False, **run_options)
        if completed.returncode != 0:
            raise TestConfigurationError(
                f"{label} tests failed with exit code {completed.returncode}"
            )
        print(f"{label.capitalize()} tests: passed", flush=True)


def run_mask_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(compiler, "activeWordMaskTester", "mask")


def run_tree_canon_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(
        compiler,
        "treeCanonTester",
        "tree canon",
        stack_limit_bytes=64 * 1024,
    )


def run_cyclic_canon_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(compiler, "cyclicCanonTester", "cyclic canon")


def run_string_assembly_unit_tests(compiler: str) -> None:
    run_cpp_unit_test(compiler, "stringAssemblyTester", "string assembly")


def build_executable(executable: Path, compiler: str) -> Path:
    executable = executable.resolve()
    executable.parent.mkdir(parents=True, exist_ok=True)

    run_mask_unit_tests(compiler)
    run_tree_canon_unit_tests(compiler)
    run_cyclic_canon_unit_tests(compiler)
    run_string_assembly_unit_tests(compiler)
    command_prefix = compiler_command(compiler)

    command = [
        *command_prefix,
        str(REPOSITORY_ROOT / "src" / "main.cpp"),
        "-std=c++20",
        "-O3",
        "-DNDEBUG",
        "-mpopcnt",
        "-march=x86-64-v3",
        "-DASSEMBLY_ENABLE_TELEMETRY",
    ]
    command.extend(["-o", str(executable)])

    print(f"Building: {shlex.join(command)}", flush=True)
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise TestConfigurationError(
            f"build failed with exit code {completed.returncode}"
        )

    return executable


def format_process_diagnostics(completed: subprocess.CompletedProcess[str]) -> str:
    diagnostics = []
    if completed.returncode != 0:
        diagnostics.append(f"process exited with code {completed.returncode}")
    if completed.stdout.strip():
        diagnostics.append(f"stdout: {completed.stdout.strip()}")
    if completed.stderr.strip():
        diagnostics.append(f"stderr: {completed.stderr.strip()}")
    return "; ".join(diagnostics)


def compare_pathway_output(actual_path: Path, expected_path: Path) -> str | None:
    if not actual_path.is_file():
        return f"expected pathway output was not created: {actual_path}"

    try:
        actual = parse_pathway_document(actual_path)
    except ValueError as error:
        return str(error)
    expected = parse_pathway_document(expected_path)
    if actual == expected:
        return None

    expected_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
    actual_lines = json.dumps(actual, indent=2, sort_keys=True).splitlines()
    difference = "\n".join(
        difflib.unified_diff(
            expected_lines,
            actual_lines,
            fromfile=str(expected_path),
            tofile=actual_path.name,
            lineterm="",
        )
    )
    return f"pathway output differs from its golden file:\n{difference}"


def run_test_case(executable: Path, case: TestCase, timeout: float) -> TestResult:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "-", Path(case.name).name)
    started = time.perf_counter()

    try:
        with tempfile.TemporaryDirectory(
            prefix=f"parallelassemblycpp-{safe_name}-"
        ) as directory:
            working_directory = Path(directory)
            copied_input = working_directory / case.source.name
            shutil.copy2(case.source, copied_input)
            input_argument = (
                copied_input.with_suffix("")
                if copied_input.name.endswith(".mol")
                else copied_input
            )

            completed = subprocess.run(
                [
                    str(executable),
                    str(input_argument),
                    f"--pathway={int(case.expected_pathway is not None)}",
                ],
                cwd=working_directory,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
            duration = time.perf_counter() - started
            output_path = molecule_output_path(copied_input, "Out")

            if completed.returncode != 0:
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=format_process_diagnostics(completed),
                )
            if not output_path.is_file():
                diagnostics = format_process_diagnostics(completed)
                suffix = f"; {diagnostics}" if diagnostics else ""
                missing_output_error = (
                    f"expected output file was not created: {output_path}{suffix}"
                )
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=missing_output_error,
                )

            output = output_path.read_text()
            match = ASSEMBLY_INDEX_PATTERN.search(output)
            if match is None:
                return TestResult(
                    case=case,
                    actual=None,
                    duration_seconds=duration,
                    error=f"assembly index was not found in {output_path.name}",
                )

            actual_index = int(match.group(1))
            pathway_failure = None
            if case.expected_pathway is not None:
                pathway_failure = compare_pathway_output(
                    molecule_output_path(copied_input, "Pathway"),
                    case.expected_pathway,
                )

            return TestResult(
                case=case,
                actual=actual_index,
                duration_seconds=duration,
                failure=pathway_failure,
            )
    except subprocess.TimeoutExpired:
        return TestResult(
            case=case,
            actual=None,
            duration_seconds=time.perf_counter() - started,
            error=f"timed out after {timeout:g} seconds",
        )
    except OSError as error:
        return TestResult(
            case=case,
            actual=None,
            duration_seconds=time.perf_counter() - started,
            error=str(error),
        )


def run_test_cases(
    executable: Path,
    cases: Sequence[TestCase],
    timeout: float,
    jobs: int,
    on_result: Callable[[TestResult], None] | None = None,
) -> list[TestResult]:
    results: list[TestResult] = []

    def record(result: TestResult) -> None:
        results.append(result)
        if on_result is not None:
            on_result(result)

    if jobs == 1:
        for case in cases:
            record(run_test_case(executable, case, timeout))
        return results

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for result in executor.map(
            lambda case: run_test_case(executable, case, timeout), cases
        ):
            record(result)

    return results


def print_result_header() -> None:
    print()
    print(f"{'Molecule':<32} {'Expected':>8} {'Actual':>8} {'Time (s)':>10}  Status")
    print("-" * 78)


def print_result(result: TestResult) -> None:
    actual = "-" if result.actual is None else str(result.actual)
    print(
        f"{result.case.name:<32.32} {result.case.expected:>8} "
        f"{actual:>8} {result.duration_seconds:>10.3f}  {result.status}"
    )
    if result.error:
        print(f"  {result.error}")
    if result.failure:
        print(f"  {result.failure}")


def print_summary(results: Sequence[TestResult]) -> None:
    passed = sum(result.status == "PASS" for result in results)
    failed = sum(result.status == "FAIL" for result in results)
    errors = sum(result.status == "ERROR" for result in results)
    duration = sum(result.duration_seconds for result in results)

    print(
        f"Summary: {passed} passed, {failed} failed, {errors} errors, "
        f"{len(results)} total ({duration:.3f}s cumulative)"
    )


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run ParallelAssemblyCpp CLI and regression checks."
    )
    parser.add_argument(
        "executable",
        nargs="?",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=(
            "ParallelAssemblyCpp executable "
            f"(default: {DEFAULT_EXECUTABLE.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=(
            "regression manifest (TSV; default: "
            f"{DEFAULT_MANIFEST.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--pathway-manifest",
        type=Path,
        default=DEFAULT_PATHWAY_MANIFEST,
        help=(
            "pathway manifest (TSV; default: "
            f"{DEFAULT_PATHWAY_MANIFEST.relative_to(REPOSITORY_ROOT)})"
        ),
    )
    parser.add_argument(
        "--audit",
        action="store_true",
        help="check manifests and fixture coverage without running tests",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="build x86-64-v3 test executables before testing",
    )
    parser.add_argument(
        "--compiler",
        default=os.environ.get("CXX") or "c++",
        help="compiler used by --build (uses $CXX, then c++)",
    )
    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=1,
        help="parallel test processes (default: 1)",
    )
    parser.add_argument(
        "--limit",
        type=positive_int,
        help="run the first N regression cases",
    )
    parser.add_argument(
        "--pathways-only",
        action="store_true",
        help="select regression cases with pathway golden files",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=300.0,
        help="timeout per case in seconds (default: 300)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="list passing checks; in audit mode, list fixture-only molecules",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = create_argument_parser().parse_args(argv)

    try:
        manifest, all_cases = load_manifest(arguments.manifest)
        all_cases = load_pathway_manifest(arguments.pathway_manifest, all_cases)
        if arguments.audit:
            audit_test_data(manifest, all_cases, arguments.verbose)
            return 0

        cases = (
            [case for case in all_cases if case.expected_pathway is not None]
            if arguments.pathways_only
            else all_cases
        )
        cases = cases[: arguments.limit]
        executable_path = arguments.executable
        if arguments.build:
            executable_path = build_executable(executable_path, arguments.compiler)
        executable = resolve_executable(executable_path)
        cli_scenarios = run_cli_checks(executable)
    except TestConfigurationError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"CLI checks: {cli_scenarios} passed", flush=True)
    print(
        f"Running {len(cases)} tests with {arguments.jobs} worker(s); "
        f"timeout={arguments.timeout:g}s",
        flush=True,
    )
    if arguments.verbose:
        print_result_header()
    results = run_test_cases(
        executable=executable,
        cases=cases,
        timeout=arguments.timeout,
        jobs=arguments.jobs,
        on_result=print_result if arguments.verbose else None,
    )

    failures = [result for result in results if result.status != "PASS"]
    if failures and not arguments.verbose:
        print_result_header()
        for result in failures:
            print_result(result)
    print_summary(results)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
