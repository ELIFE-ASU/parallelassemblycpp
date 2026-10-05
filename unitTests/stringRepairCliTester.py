"""Check string Re-Pair selection, Unicode certificates, and CLI error policy."""

# Assertions implement this executable integration harness.
# ruff: noqa: S101, S603

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import subprocess
import tempfile
from typing import Any

BOUND = re.compile(r"has assembly upper bound: (-?\d+)\n")
FLAGS = ["--run-strings=1", "--algorithm=re-pair"]
LAUNCHER: list[str] = []
STATUS = "status: heuristic upper bound (minimum not proven)\n"


def replay(document: dict[str, Any], original: str, reversed_mode: bool) -> int:
    """Reconstruct scalars from the exported grammar without any solver code."""
    assert document["schema"] == "string-repair-assembly-v1"
    assert document["input"] == original
    assert document["length"] == len(original), "certificate uses UTF-8 byte length"
    assert document["accept_reversed"] is reversed_mode
    assert document["trivial_upper_bound"] == len(original) - 1
    assert document["rule_count"] == len(document["rules"])
    assert document["remaining_fragments"] == len(document["residual"])
    assert (
        document["upper_bound"]
        == len(document["rules"]) + len(document["residual"]) - 1
    )
    assert document["upper_bound"] <= len(original) - 1
    expansions: dict[int, str] = {}
    terminals = document["terminals"]
    assert [chr(t["code_point"]) for t in terminals] == list(dict.fromkeys(original))
    for terminal in terminals:
        assert terminal["id"] not in expansions
        expansions[terminal["id"]] = chr(terminal["code_point"])
    for rule in document["rules"]:
        assert rule["id"] not in expansions
        assert rule["left"] in expansions and rule["right"] in expansions
        left = expansions[rule["left"]]
        right = expansions[rule["right"]]
        assert reversed_mode or not (rule["left_reversed"] or rule["right_reversed"])
        if rule["left_reversed"]:
            left = left[::-1]
        if rule["right_reversed"]:
            right = right[::-1]
        expansion = left + right
        assert rule["length"] == len(expansion) > 1
        assert not reversed_mode or expansion <= expansion[::-1]
        expansions[rule["id"]] = expansion
    recovered = ""
    for fragment in document["residual"]:
        assert fragment["symbol"] in expansions
        assert reversed_mode or not fragment["reversed"]
        expansion = expansions[fragment["symbol"]]
        if fragment["reversed"]:
            expansion = expansion[::-1]
        assert fragment["offset"] == len(recovered)
        assert fragment["length"] == len(expansion)
        assert original[len(recovered) : len(recovered) + len(expansion)] == expansion
        recovered += expansion
    assert recovered == original
    return int(document["upper_bound"])


def test_checker_rejects_mutations(document: dict[str, Any]) -> None:
    """Check the independent verifier rejects corrupted but parseable certificates."""
    original = str(document["input"])
    assert replay(document, original, False) >= 0
    mutations: list[tuple[tuple[str | int, ...], object]] = [
        (("length",), len(original.encode("utf-8")) + 1),
        (("upper_bound",), -10),
        (("rule_count",), len(document["rules"]) + 1),
        (("terminals", 0, "code_point"), ord("z")),
        (("rules", 0, "id"), document["terminals"][0]["id"]),
        (("rules", 0, "left"), document["rules"][0]["id"]),
        (("rules", 0, "left"), document["rules"][0]["right"]),
        (("rules", 0, "length"), 999),
        (("rules", 0, "left_reversed"), True),
        (("residual", 0, "symbol"), -1),
        (("residual", 0, "offset"), 1),
        (("residual", 0, "length"), 999),
        (("residual", 0, "reversed"), True),
    ]
    for path, replacement in mutations:
        damaged = json.loads(json.dumps(document))
        parent = damaged
        for component in path[:-1]:
            parent = parent[component]
        parent[path[-1]] = replacement
        try:
            replay(damaged, original, False)
        except (AssertionError, KeyError):
            continue
        raise AssertionError(f"certificate verifier accepted mutation at {path}")


def invoke(
    executable: pathlib.Path,
    directory: pathlib.Path,
    flags: list[str],
    *,
    name: str = "input.txt",
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [*LAUNCHER, str(executable), name, *flags],
        cwd=directory,
        capture_output=True,
        text=True,
        encoding="utf-8",
        timeout=30,
        check=False,
    )


def success(
    executable: pathlib.Path,
    root: pathlib.Path,
    name: str,
    records: list[str],
    *,
    data: bytes | None = None,
    extra: tuple[str, ...] = (),
    reversed_mode: bool = False,
    pathway: bool = True,
    expected: list[int] | None = None,
) -> tuple[str, list[dict[str, Any]]]:
    directory = root / name
    directory.mkdir()
    source = directory / "input.txt"
    if data is None:
        data = ("\n".join(records) + "\n").encode("utf-8") if records else b""
    source.write_bytes(data)
    flags = [
        *FLAGS,
        f"--accept-palindromes={int(reversed_mode)}",
        f"--pathway={int(pathway)}",
        *extra,
    ]
    completed = invoke(executable, directory, flags)
    assert completed.returncode == 0, f"{name}: {completed.stderr}"
    assert source.read_bytes() == data, f"{name}: input was modified"
    # newline='' preserves a final CR that is an input symbol, rather than CRLF.
    with (directory / "input.txtOut").open(encoding="utf-8", newline="") as stream:
        output = stream.read().replace("\r\n", "\n")
    bounds = [int(match) for match in BOUND.findall(output)]
    assert len(bounds) == len(records), (
        f"{name}: missing or spurious records: {output!r}"
    )
    if expected is not None:
        assert bounds == expected, f"{name}: wrong bounds {bounds} != {expected}"
    assert "has assembly index:" not in output, (
        "a heuristic was presented as an exact result"
    )
    assert output.count(STATUS) == len(records)
    assert output.count("time elapsed: ") == len(records)
    assert output.count("string-repair rules: ") == len(records)
    assert output.count("remaining fragments: ") == len(records)
    assert "runtime limit reached" not in output
    expected_files = {"input.txt", "input.txtOut"}
    if pathway:
        expected_files |= {f"input.txt_{i}_Pathway" for i in range(len(records))}
    assert {p.name for p in directory.iterdir()} == expected_files
    certificates = []
    if pathway:
        for index, original in enumerate(records):
            certificate = json.loads(
                (directory / f"input.txt_{index}_Pathway").read_text(encoding="ascii")
            )
            assert replay(certificate, original, reversed_mode) == bounds[index]
            certificates.append(certificate)
    if "--parallel=auto" in extra:
        assert (
            "parallel: serial fallback: string-repair upper bound uses serial execution"
            in completed.stderr
        )
    stable = "\n".join(
        line for line in output.split("\n") if not line.startswith("time elapsed: ")
    )
    return stable, certificates


def invalid_options(
    executable: pathlib.Path,
    root: pathlib.Path,
    name: str,
    flags: list[str],
    diagnostic: str,
) -> None:
    directory = root / name
    directory.mkdir()
    source = directory / "input.txt"
    source.write_bytes(b"abab\n")
    completed = invoke(executable, directory, flags)
    assert completed.returncode == 2, f"{name}: {completed.stderr}"
    assert diagnostic in completed.stderr, f"{name}: {completed.stderr}"
    assert source.read_bytes() == b"abab\n"
    assert {p.name for p in directory.iterdir()} == {"input.txt"}, (
        "rejected options created output"
    )


def io_failure(executable: pathlib.Path, root: pathlib.Path, fixture: str) -> None:
    directory = root / ("io-" + fixture)
    directory.mkdir()
    source = directory / "input.txt"
    original = b"abab\n"
    if fixture != "missing-input":
        source.write_bytes(original)
    if fixture == "output-directory":
        (directory / "input.txtOut").mkdir()
    elif fixture == "pathway-directory":
        (directory / "input.txt_0_Pathway").mkdir()
    elif fixture in ("output-hardlink", "pathway-hardlink"):
        target = (
            "input.txtOut" if fixture == "output-hardlink" else "input.txt_0_Pathway"
        )
        (directory / target).hardlink_to(source)
    elif fixture in ("output-symlink", "pathway-symlink"):
        target = (
            "input.txtOut" if fixture == "output-symlink" else "input.txt_0_Pathway"
        )
        try:
            (directory / target).symlink_to(source.name)
        except OSError as error:
            # Match the existing CLI suite: Windows may require symlink privilege.
            if os.name == "nt" and error.winerror in (5, 1314):
                return
            raise
    completed = invoke(executable, directory, FLAGS)
    assert completed.returncode == 1, f"{fixture}: {completed.stderr}"
    assert "error:" in completed.stderr
    if fixture != "missing-input":
        assert source.read_bytes() == original, f"{fixture}: input was overwritten"
    if "link" in fixture:
        assert "would overwrite input file" in completed.stderr


def test_pathway_disabled_preserves_file(
    executable: pathlib.Path, root: pathlib.Path
) -> None:
    directory = root / "preserve-pathway"
    directory.mkdir()
    (directory / "input.txt").write_bytes(b"abab\n")
    certificate = directory / "input.txt_0_Pathway"
    before = b"existing certificate\0must remain unchanged\n"
    certificate.write_bytes(before)
    completed = invoke(executable, directory, [*FLAGS, "--pathway=0"])
    assert completed.returncode == 0, completed.stderr
    assert certificate.read_bytes() == before
    assert "has assembly upper bound: 2" in (directory / "input.txtOut").read_text()


def test_exact_selection(executable: pathlib.Path, root: pathlib.Path) -> None:
    results = []
    for name, options in (("default", []), ("full", ["--algorithm=full"])):
        directory = root / ("exact-" + name)
        directory.mkdir()
        (directory / "input.txt").write_bytes(b"aaaaaaaa\n")
        completed = invoke(executable, directory, ["--run-strings=1", *options])
        assert completed.returncode == 0, completed.stderr
        output = (directory / "input.txtOut").read_text(encoding="utf-8")
        assert "has assembly index: 3\n" in output
        assert "upper bound" not in output
        assert "Initial Re-Pair" not in completed.stdout
        certificate = json.loads(
            (directory / "input.txt_0_Pathway").read_text(encoding="ascii")
        )
        assert set(certificate) == {"file_graph", "remnant", "duplicates"}
        results.append(
            (
                [
                    line
                    for line in output.splitlines()
                    if not line.startswith("time elapsed: ")
                ],
                certificate,
            )
        )
        limited_directory = root / ("exact-limited-" + name)
        limited_directory.mkdir()
        (limited_directory / "input.txt").write_bytes(b"aaaaaaaa\n")
        limited = invoke(
            executable,
            limited_directory,
            ["--run-strings=1", *options, "--runtime=0", "--pathway=0"],
        )
        assert limited.returncode == 0, limited.stderr
        limited_output = (limited_directory / "input.txtOut").read_text(
            encoding="utf-8"
        )
        assert "has assembly index: 7\n" in limited_output
        assert "status: runtime limit reached\n" in limited_output
        assert "heuristic upper bound" not in limited_output
    assert results[0] == results[1], (
        "explicit full changed the default string calculation"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument("--telemetry-executable", type=pathlib.Path)
    parser.add_argument("--mpiexec", type=pathlib.Path)
    parser.add_argument("--mpiexec-numproc-flag", default="-n")
    options = parser.parse_args()
    executable = options.executable.resolve()
    if options.mpiexec is not None:
        LAUNCHER.extend(
            [str(options.mpiexec.resolve()), options.mpiexec_numproc_flag, "2"]
        )
    with tempfile.TemporaryDirectory(prefix="string-repair-cli-") as temporary:
        root = pathlib.Path(temporary)
        test_pathway_disabled_preserves_file(executable, root)
        test_exact_selection(executable, root)
        records = [
            "",
            "a",
            "abcdef",
            "aaa",
            "aaaa",
            "abab",
            "abcabc",
            "aaaaaaaa",
            "abba",
            "abcxcba",
        ]
        _, certificates = success(
            executable, root, "basic", records, expected=[-1, 0, 5, 2, 2, 2, 3, 3, 3, 6]
        )
        test_checker_rejects_mutations(certificates[6])
        success(
            executable,
            root,
            "reverse",
            records,
            reversed_mode=True,
            expected=[-1, 0, 5, 2, 2, 2, 3, 3, 2, 4],
        )
        special = [
            "éaéa",
            "😀é中😀é中",
            "\u03b1\u03b2\u03b3x\u03b3\u03b2\u03b1",
            "e\u0301é",
            "\0😀\0😀",
            'a"\\\b\f\t\x01z',
            "\x7f\u0080\u07ff\u0800\ud7ff\ue000\uffff\U00010000\U0010ffff",
            "\u0085\u2028\u2029",
            "\ufeffabab",
            " \t \t",
        ]
        for reversed_mode in (False, True):
            success(
                executable,
                root,
                f"unicode-{reversed_mode}",
                special,
                reversed_mode=reversed_mode,
            )
        success(executable, root, "empty-file", [], data=b"", expected=[])
        success(executable, root, "empty-line", [""], data=b"\n", expected=[-1])
        success(
            executable,
            root,
            "line-endings",
            ["abab", "", "x", "abab\r"],
            data=b"abab\r\n\r\nx\r\nabab\r",
            expected=[2, -1, 0, 3],
        )
        success(
            executable, root, "no-final-newline", ["abab"], data=b"abab", expected=[2]
        )
        success(executable, root, "no-pathway", special, pathway=False)
        serial = success(executable, root, "serial", records, extra=("--parallel=off",))
        automatic = success(
            executable,
            root,
            "automatic",
            records,
            extra=("--parallel=auto", "--threads=2"),
        )
        assert serial == automatic, "auto fallback changed bounds or certificates"
        repeated = success(executable, root, "repeat", records)
        assert serial == repeated, "repeated CLI run changed deterministic output"
        success(
            executable,
            root,
            "disabled-options",
            ["abab"],
            extra=("--write-intermediate-mas=0",),
            expected=[2],
        )
        invalid = {
            "parallel-on": ([*FLAGS, "--parallel=on"], "--parallel=on"),
            "runtime-zero": ([*FLAGS, "--runtime=0"], "--runtime"),
            "runtime-positive": ([*FLAGS, "--runtime=1000000"], "--runtime"),
            "enum-limit": ([*FLAGS, "--enum-max=1000"], "--enum-max"),
            "intermediate": (
                [*FLAGS, "--write-intermediate-mas=1"],
                "--write-intermediate-mas=1",
            ),
            "legacy-graph": (
                ["--run-strings=1", "--upper-bound=graph-repair"],
                "string assembly",
            ),
            "conflicting-selectors": (
                [*FLAGS, "--upper-bound=graph-repair"],
                "cannot be combined",
            ),
            "duplicate-selector": ([*FLAGS, "--algorithm=re-pair"], "only once"),
            "hydrogens": ([*FLAGS, "--remove-hydrogens=0"], "--remove-hydrogens"),
            "disjoint": ([*FLAGS, "--compensate-disjoint=1"], "--compensate-disjoint"),
        }
        for name, (flags, diagnostic) in invalid.items():
            invalid_options(executable, root, name, flags, diagnostic)
        for fixture in (
            "missing-input",
            "output-directory",
            "pathway-directory",
            "output-hardlink",
            "pathway-hardlink",
            "output-symlink",
            "pathway-symlink",
        ):
            io_failure(executable, root, fixture)
        for index, malformed in enumerate(
            (b"\x80", b"\xc0\xaf", b"\xed\xa0\x80", b"\xf4\x90\x80\x80", b"\xe2\x82")
        ):
            directory = root / f"malformed-{index}"
            directory.mkdir()
            source = directory / "input.txt"
            original = b"abab\n" + malformed + b"\n"
            source.write_bytes(original)
            completed = invoke(executable, directory, FLAGS)
            assert completed.returncode == 1, completed.stderr
            assert "string calculation failed on line 2" in completed.stderr
            assert source.read_bytes() == original
            assert (directory / "input.txt_0_Pathway").is_file()
            assert not (directory / "input.txt_1_Pathway").exists()
        if options.telemetry_executable is not None:
            telemetry = options.telemetry_executable.resolve()
            invalid_options(
                telemetry,
                root,
                "telemetry-on",
                [*FLAGS, "--telemetry=1"],
                "--telemetry=1",
            )
            success(
                telemetry,
                root,
                "telemetry-off",
                ["abab"],
                extra=("--telemetry=0",),
                expected=[2],
            )
    print("String Re-Pair CLI checks passed")


if __name__ == "__main__":
    main()
