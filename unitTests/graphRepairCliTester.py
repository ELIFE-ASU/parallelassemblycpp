"""Check full/graph-repair CLI selection, output semantics, and option rejection."""

# Assertions are the checks in this harness; the executable is selected explicitly.
# ruff: noqa: S101, S603

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile


def run_case(
    executable: pathlib.Path,
    root: pathlib.Path,
    name: str,
    flags: list[str],
    *,
    succeeds: bool = True,
    pathway: bool = True,
    fallback: bool = False,
    upper_bound: bool = True,
    strings: bool = False,
    diagnostics: tuple[str, ...] = (),
) -> tuple[list[str], dict[str, object] | None] | None:
    directory = root / name
    directory.mkdir()
    source = directory / ("chain.txt" if strings else "chain.graph")
    source.write_text(
        "AAAAAAAA\n"
        if strings
        else (
            "eight-bond chain\n9\n"
            "1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9\n"
            "C C C C C C C C C\n1 1 1 1 1 1 1 1\n"
        ),
        encoding="utf-8",
    )
    before = source.read_bytes()
    completed = subprocess.run(
        [str(executable), str(source), *flags],
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert source.read_bytes() == before, f"{name}: input changed"
    if not succeeds:
        assert completed.returncode == 2, (
            f"{name}: expected option rejection, got {completed.returncode}: "
            f"{completed.stderr}"
        )
        assert "error:" in completed.stderr, f"{name}: no diagnostic"
        for diagnostic in diagnostics:
            assert diagnostic in completed.stderr, (
                f"{name}: missing {diagnostic!r} diagnostic: {completed.stderr}"
            )
        assert sorted(p.name for p in directory.iterdir()) == [source.name], (
            f"{name}: invalid options created outputs"
        )
        return None
    assert completed.returncode == 0, f"{name}: {completed.stderr}"
    output = pathlib.Path(f"{source}Out").read_text(encoding="utf-8")
    assert "has assembly index: 3\n" in output, f"{name}: wrong index: {output}"
    if upper_bound:
        assert "status: heuristic upper bound (minimum not proven)\n" in output
        assert "trivial upper bound: 7\n" in output
    else:
        assert "heuristic upper bound" not in output, f"{name}: ran bound only"
        assert "trivial upper bound:" not in output, f"{name}: wrong output format"
    assert "time elapsed: " in output
    certificate_path = pathlib.Path(f"{source}{'_0_' if strings else ''}Pathway")
    assert certificate_path.exists() == pathway, f"{name}: wrong pathway policy"
    certificate = None
    if pathway:
        certificate = json.loads(certificate_path.read_text(encoding="ascii"))
        if upper_bound:
            assert certificate["schema"] == "graph-repair-assembly-v1"
        else:
            expected_keys = {"file_graph", "remnant", "duplicates"}
            if not strings:
                expected_keys.add("removed_edges")
            assert set(certificate) == expected_keys, (
                f"{name}: full calculation did not preserve the exact pathway schema"
            )
    if fallback:
        assert "serial fallback" in completed.stderr, f"{name}: no fallback notice"
    stable_output = [
        line.replace(str(source), "INPUT")
        for line in output.splitlines()
        if not line.startswith("time elapsed: ")
    ]
    return stable_output, certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument("--telemetry-executable", type=pathlib.Path)
    options = parser.parse_args()
    executable = options.executable.resolve()
    help_result = subprocess.run(
        [str(executable), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert help_result.returncode == 0, help_result.stderr
    assert "--algorithm=<full|re-pair>" in help_result.stdout
    with tempfile.TemporaryDirectory(prefix="graph-repair-cli-") as temporary:
        root = pathlib.Path(temporary)
        bound = ["--upper-bound=graph-repair"]
        algorithm = ["--algorithm=re-pair"]
        legacy_result = run_case(executable, root, "legacy-bound", bound)
        selected_result = run_case(executable, root, "selected-bound", algorithm)
        assert selected_result == legacy_result, "Re-Pair selectors differ"
        for strings in (False, True):
            input_mode = "strings" if strings else "graph"
            flags = ["--run-strings=1"] if strings else []
            default_result = run_case(
                executable,
                root,
                f"default-{input_mode}",
                flags,
                upper_bound=False,
                strings=strings,
            )
            full_result = run_case(
                executable,
                root,
                f"full-{input_mode}",
                [*flags, "--algorithm=full"],
                upper_bound=False,
                strings=strings,
            )
            assert full_result == default_result, (
                f"explicit full calculation differs from default for {input_mode}"
            )
        run_case(executable, root, "no-pathway", [*bound, "--pathway=0"], pathway=False)
        run_case(
            executable,
            root,
            "selected-no-pathway",
            [*algorithm, "--pathway=0"],
            pathway=False,
        )
        run_case(
            executable, root, "automatic", [*bound, "--parallel=auto"], fallback=True
        )
        run_case(
            executable,
            root,
            "disabled-search-options",
            [*bound, "--parallel=off", "--write-intermediate-mas=0"],
        )
        invalid_options = {
            "invalid-name": ["--upper-bound=unknown"],
            "strings": [*bound, "--run-strings=1"],
            "runtime": [*bound, "--runtime=0"],
            "enum-limit": [*bound, "--enum-max=50000000"],
            "forced-parallel": [*bound, "--parallel=on"],
            "intermediate": [*bound, "--write-intermediate-mas=1"],
            "duplicate": [*bound, *bound],
        }
        for name, flags in invalid_options.items():
            run_case(executable, root, name, flags, succeeds=False)
        invalid_algorithms = {
            "unknown": ["--algorithm=unknown"],
            "empty": ["--algorithm="],
            "missing": ["--algorithm"],
            "duplicate": [*algorithm, *algorithm],
            "duplicate-full": ["--algorithm=full", "--algorithm=full"],
            "duplicate-full-then-re-pair": ["--algorithm=full", *algorithm],
            "duplicate-re-pair-then-full": [*algorithm, "--algorithm=full"],
        }
        for name, flags in invalid_algorithms.items():
            run_case(
                executable,
                root,
                f"algorithm-{name}",
                flags,
                succeeds=False,
                diagnostics=("--algorithm",),
            )
        for value in ("full", "re-pair"):
            selected = [f"--algorithm={value}"]
            for order, flags in enumerate(([*selected, *bound], [*bound, *selected])):
                run_case(
                    executable,
                    root,
                    f"conflicting-selectors-{value}-{order}",
                    flags,
                    succeeds=False,
                    diagnostics=("--algorithm", "--upper-bound"),
                )
        unsupported = {
            "strings": ("--run-strings=1", "string assembly"),
            "runtime": ("--runtime=0", "--runtime"),
            "enum-limit": ("--enum-max=50000000", "--enum-max"),
            "forced-parallel": ("--parallel=on", "--parallel=on"),
            "intermediate": (
                "--write-intermediate-mas=1",
                "--write-intermediate-mas=1",
            ),
        }
        for name, (flag, diagnostic) in unsupported.items():
            run_case(
                executable,
                root,
                f"algorithm-unsupported-{name}",
                [*algorithm, flag],
                succeeds=False,
                diagnostics=(diagnostic,),
            )
        if options.telemetry_executable is not None:
            telemetry = options.telemetry_executable.resolve()
            run_case(
                telemetry,
                root,
                "telemetry-enabled",
                [*bound, "--telemetry=1"],
                succeeds=False,
            )
            run_case(telemetry, root, "telemetry-disabled", [*bound, "--telemetry=0"])
            run_case(
                telemetry,
                root,
                "selected-telemetry-enabled",
                [*algorithm, "--telemetry=1"],
                succeeds=False,
                diagnostics=("--telemetry=1",),
            )
            run_case(
                telemetry,
                root,
                "selected-telemetry-disabled",
                [*algorithm, "--telemetry=0"],
            )
    print("Graph-repair CLI checks passed")


if __name__ == "__main__":
    main()
