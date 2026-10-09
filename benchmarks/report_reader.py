"""Shared syntax and identities for schema-v2 benchmark evidence.

Acceptance policies belong to the promotion and scaling gates.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from pathlib import Path

ENVIRONMENT_KEY_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
PAIRED_COMPARISON_ORDER = (
    "baseline/candidate on odd rounds, candidate/baseline on even rounds"
)


@dataclass(frozen=True)
class ExecutionIdentity:
    launcher: tuple[str, ...]
    arguments: tuple[str, ...]
    environment: tuple[tuple[str, str], ...]


def execution_identity(
    document: dict[str, object],
    role: str,
    path: Path,
    error_type: type[RuntimeError],
) -> ExecutionIdentity:
    execution = document.get("execution")
    if execution is None:
        # Early schema-v2 reports always launched both roles directly.
        return ExecutionIdentity((), (), ())
    if not isinstance(execution, dict):
        raise error_type(f"invalid execution configurations in {path}")
    config = execution.get(role)
    if not isinstance(config, dict):
        raise error_type(f"missing {role} execution configuration in {path}")

    launcher = config.get("launcher")
    if not isinstance(launcher, list) or any(
        not isinstance(argument, str) or not argument for argument in launcher
    ):
        raise error_type(f"invalid {role} launcher configuration in {path}")
    arguments = config.get("arguments", [])
    if not isinstance(arguments, list) or any(
        not isinstance(argument, str) or not argument or "\x00" in argument
        for argument in arguments
    ):
        raise error_type(f"invalid {role} arguments configuration in {path}")
    environment = config.get("environment")
    if not isinstance(environment, dict):
        raise error_type(f"invalid {role} environment configuration in {path}")

    normalized_environment: list[tuple[str, str]] = []
    for key, value in environment.items():
        if (
            not isinstance(key, str)
            or ENVIRONMENT_KEY_PATTERN.fullmatch(key) is None
            or not isinstance(value, str)
            or "\x00" in value
        ):
            raise error_type(f"invalid {role} environment configuration in {path}")
        normalized_environment.append((key, value))
    return ExecutionIdentity(
        tuple(launcher),
        tuple(arguments),
        tuple(sorted(normalized_environment)),
    )


def load_result(path: Path, error_type: type[RuntimeError]) -> dict[str, object]:
    try:
        with path.open(encoding="utf-8") as stream:
            document: object = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise error_type(f"could not read benchmark report {path}: {error}") from error
    if not isinstance(document, dict):
        raise error_type(f"invalid benchmark report {path}: expected an object")
    schema_version = document.get("schema_version")
    if type(schema_version) is not int or schema_version != 2:
        raise error_type(f"invalid benchmark report {path}: expected schema_version 2")
    return document


def corpus_input_fingerprints(
    corpus: dict[str, object],
    path: Path,
    read_string: Callable[[object, Sequence[str], str], str],
    error_type: type[RuntimeError],
) -> dict[str, str]:
    """Read unique input fingerprints using the caller's string normalization."""
    inputs = corpus.get("inputs")
    if not isinstance(inputs, list):
        raise error_type(f"missing corpus input fingerprints in {path}")
    fingerprints: dict[str, str] = {}
    for index, entry in enumerate(inputs):
        if not isinstance(entry, dict):
            raise error_type(f"invalid corpus input fingerprint {index} in {path}")
        name = read_string(entry, ("name",), f"corpus input {index} name in {path}")
        sha256 = read_string(
            entry, ("sha256",), f"corpus input {name!r} SHA-256 in {path}"
        )
        if name in fingerprints:
            raise error_type(f"duplicate corpus input fingerprint {name!r} in {path}")
        fingerprints[name] = sha256
    return fingerprints
