"""Run child commands with timeout and interruption cleanup."""

from __future__ import annotations

import os
import signal
import subprocess
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from pathlib import Path


def terminate_command(process: subprocess.Popen[str]) -> None:
    """Terminate a configured launch, including its POSIX process group."""
    killed_group = False
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
            killed_group = True
        except ProcessLookupError:
            killed_group = True
        except OSError:
            pass
    if not killed_group:
        process.kill()


def run_command(
    command: Sequence[str],
    working_directory: Path,
    timeout: float,
    environment: Mapping[str, str] | None,
) -> subprocess.CompletedProcess[str]:
    """Run a command and clean up its process group on timeout or interruption."""
    popen_arguments: dict[str, object] = {
        "cwd": working_directory,
        "env": environment,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        popen_arguments["start_new_session"] = True

    process = subprocess.Popen(command, **popen_arguments)  # noqa: S603 - explicit caller command.
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as error:
        terminate_command(process)
        stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(
            command,
            error.timeout,
            output=stdout,
            stderr=stderr,
        ) from error
    except BaseException:
        # start_new_session deliberately keeps launcher workers out of the
        # driver's foreground process group, so an interrupt must stop them
        # explicitly before it propagates to the caller.
        terminate_command(process)
        process.communicate()
        raise

    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
