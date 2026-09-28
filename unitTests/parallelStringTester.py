"""Check string assembly parity, deterministic pathways, and parallel CLI policy."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from parallelSolverTester import run_command

ASSEMBLY_INDEX_PATTERN = re.compile(r"has assembly index:\s*(-?\d+)")
FINITE_RUNTIME = "--runtime=1000000000"
INPUTS = (
    "",
    "x",
    "abcdef",
    "abab",
    "abcxcba",
    "ababcdcd",
    "abcababc",
    "abcabzzxyyxabcab",
    "0" * 25 + "1" * 25 + "2" * 25,
    "\u00e9",
    "\U0001f600\U0001f601\U0001f602",
    "\U0001f600\U0001f601\U0001f600\U0001f601",
    "\u03b1\u03b2\u03b3x\u03b3\u03b2\u03b1",
)


def run_case(
    executable: Path,
    *,
    launcher: list[str],
    mode: str,
    threads: int,
    accept_reversed: bool,
    thread_option: str | None = None,
    pathway: bool = True,
    input_exists: bool = True,
    output_blocked: bool = False,
    extra: tuple[str, ...] = (),
    expected_error: str | None = None,
    expected_fallback: str | None = None,
    timeout: float,
) -> tuple[list[int], list[object]]:
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": str(threads),
            "OMP_THREAD_LIMIT": str(threads),
            "OMP_DYNAMIC": "FALSE",
        }
    )
    with tempfile.TemporaryDirectory(
        prefix="parallelassemblycpp-strings-"
    ) as temporary:
        directory = Path(temporary)
        # Include empty lines, CRLF, and a final unterminated record in every run.
        if input_exists:
            (directory / "input").write_bytes("\r\n".join(INPUTS).encode("utf-8"))
        if output_blocked:
            (directory / "inputOut").mkdir()
        arguments = [
            *launcher,
            str(executable),
            "input",
            "--run-strings=1",
            f"--parallel={mode}",
            f"--threads={threads if thread_option is None else thread_option}",
            f"--pathway={int(pathway)}",
            f"--accept-palindromes={int(accept_reversed)}",
            *extra,
        ]
        completed = run_command(arguments, directory, environment, timeout)
        details = (
            f"{executable.name} {mode}: exit {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
        if expected_error is not None:
            if completed.returncode == 0 or expected_error not in completed.stderr:
                raise AssertionError(
                    f"expected rejection {expected_error!r}\n{details}"
                )
            if not output_blocked and (directory / "inputOut").exists():
                raise AssertionError(
                    f"rejected options created an output file\n{details}"
                )
            return [], []
        if completed.returncode != 0:
            raise AssertionError(details)
        if (
            expected_fallback is not None
            and ("parallel: serial fallback: " + expected_fallback)
            not in completed.stderr
        ):
            raise AssertionError(f"missing expected serial fallback\n{details}")
        output = (directory / "inputOut").read_text(encoding="utf-8")
        indices = [
            int(match.group(1)) for match in ASSEMBLY_INDEX_PATTERN.finditer(output)
        ]
        if len(indices) != len(INPUTS) or output.count("time elapsed:") != len(INPUTS):
            raise AssertionError(
                f"missing or repeated string records: {indices}\n{details}"
            )
        pathways = []
        if pathway:
            pathways = [
                json.loads(
                    (directory / f"input_{index}_Pathway").read_text(encoding="utf-8")
                )
                for index in range(len(INPUTS))
            ]
        elif list(directory.glob("input_*_Pathway")):
            raise AssertionError(f"--pathway=0 created a pathway file\n{details}")
        if (
            FINITE_RUNTIME in extra
            and mode == "auto"
            and "serial fallback" not in completed.stderr
        ):
            raise AssertionError(
                f"finite runtime did not explain its fallback\n{details}"
            )
        return indices, pathways


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", type=Path, required=True)
    parser.add_argument("--openmp", type=Path)
    parser.add_argument("--mpi", type=Path)
    parser.add_argument("--hybrid", type=Path)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--mpiexec-numproc-flag", default="-n")
    parser.add_argument("--timeout", type=float, default=30.0)
    arguments = parser.parse_args()
    if not any((arguments.openmp, arguments.mpi, arguments.hybrid)):
        parser.error("at least one parallel executable is required")
    if (arguments.mpi or arguments.hybrid) and arguments.mpiexec is None:
        parser.error("--mpiexec is required for distributed targets")

    baselines = {}
    for accept_reversed in (False, True):
        baselines[accept_reversed] = run_case(
            arguments.serial.resolve(),
            launcher=[],
            mode="off",
            threads=1,
            accept_reversed=accept_reversed,
            timeout=arguments.timeout,
        )
    if baselines[False][0] != [-1, 0, 5, 2, 6, 5, 4, 10, 20, 0, 2, 2, 6]:
        raise AssertionError(
            f"serial fixture indices are incorrect: {baselines[False][0]}"
        )
    if baselines[True][0][4] != 4 or baselines[True][0][-1] != 4:
        raise AssertionError("serial reversal fixture has an incorrect index")

    runs = 2
    for backend in ("openmp", "mpi", "hybrid"):
        executable = getattr(arguments, backend)
        if executable is None:
            continue
        launcher = []
        if backend in {"mpi", "hybrid"}:
            launcher = [
                str(arguments.mpiexec.resolve()),
                arguments.mpiexec_numproc_flag,
                "2",
            ]
        threads = 1 if backend == "mpi" else 2
        for accept_reversed in (False, True):
            for mode in ("on", "on", "auto", "off"):
                actual = run_case(
                    executable.resolve(),
                    launcher=launcher,
                    mode=mode,
                    threads=threads,
                    accept_reversed=accept_reversed,
                    timeout=arguments.timeout,
                )
                if actual != baselines[accept_reversed]:
                    raise AssertionError(
                        f"{backend} {mode}, reversed={accept_reversed}: "
                        "indices or ordered pathway JSON differ from serial"
                    )
                runs += 1
        for mode in ("on", "auto"):
            actual = run_case(
                executable.resolve(),
                launcher=launcher,
                mode=mode,
                threads=threads,
                accept_reversed=False,
                extra=(FINITE_RUNTIME,),
                expected_error="finite --runtime" if mode == "on" else None,
                timeout=arguments.timeout,
            )
            if mode == "auto" and actual != baselines[False]:
                raise AssertionError(f"{backend}: runtime fallback changed results")
            runs += 1
        for mode in ("on", "auto"):
            actual = run_case(
                executable.resolve(),
                launcher=launcher,
                mode=mode,
                threads=threads,
                thread_option="auto",
                accept_reversed=False,
                timeout=arguments.timeout,
            )
            if actual != baselines[False]:
                raise AssertionError(f"{backend}: automatic threads changed results")
            runs += 1
        invalid_threads_reason = (
            "--threads=2 requires an OpenMP-enabled executable"
            if backend == "mpi"
            else "--threads=3 exceeds the OpenMP thread limit 2"
        )
        for mode in ("on", "auto"):
            actual = run_case(
                executable.resolve(),
                launcher=launcher,
                mode=mode,
                threads=threads,
                thread_option=str(threads + 1),
                accept_reversed=False,
                expected_error=invalid_threads_reason if mode == "on" else None,
                expected_fallback=invalid_threads_reason if mode == "auto" else None,
                timeout=arguments.timeout,
            )
            if mode == "auto" and actual != baselines[False]:
                raise AssertionError(
                    f"{backend}: thread-limit fallback changed results"
                )
            runs += 1
        without_pathways = run_case(
            executable.resolve(),
            launcher=launcher,
            mode="on",
            threads=threads,
            accept_reversed=False,
            pathway=False,
            timeout=arguments.timeout,
        )
        if without_pathways != (baselines[False][0], []):
            raise AssertionError(f"{backend}: disabling pathways changed indices")
        runs += 1
        if launcher:
            for output_blocked in (False, True):
                run_case(
                    executable.resolve(),
                    launcher=launcher,
                    mode="on",
                    threads=threads,
                    accept_reversed=False,
                    input_exists=output_blocked,
                    output_blocked=output_blocked,
                    expected_error=(
                        "could not open output file"
                        if output_blocked
                        else "could not open input file"
                    ),
                    timeout=arguments.timeout,
                )
                runs += 1
        run_case(
            executable.resolve(),
            launcher=launcher,
            mode="on",
            threads=threads,
            accept_reversed=False,
            extra=("--write-intermediate-mas=1",),
            expected_error=(
                "--write-intermediate-mas is unavailable for string assembly"
            ),
            timeout=arguments.timeout,
        )
        runs += 1
        print(f"PASS {backend}: string indices, deterministic pathways, and policy")
    print(f"PASS: {runs} string CLI runs")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, subprocess.TimeoutExpired, ValueError) as error:
        print(f"parallel string test failed: {error}", file=sys.stderr)
        sys.exit(1)
