"""Validate repository-wide text-file invariants."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BINARY_SUFFIXES = frozenset({".pdf"})
TRAILING_WHITESPACE_EXEMPT_SUFFIXES = frozenset({".mol", ".sdf", ".tsv"})
IGNORED_DIRECTORIES = frozenset(
    {
        ".agents",
        ".codex",
        ".git",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        ".venv",
        "__pycache__",
        "build",
        "venv",
    }
)


def source_archive_paths() -> list[Path]:
    """Return source-archive files while pruning generated and metadata trees."""
    paths: list[Path] = []
    for directory, child_directories, filenames in os.walk(REPOSITORY_ROOT):
        child_directories[:] = sorted(
            name for name in child_directories if name not in IGNORED_DIRECTORIES
        )
        paths.extend(
            Path(directory) / filename
            for filename in sorted(filenames)
            if not (Path(directory) / filename).is_symlink()
        )
    return paths


def repository_paths() -> list[Path]:
    """Return tracked and untracked, non-ignored source files."""
    if not (REPOSITORY_ROOT / ".git").exists():
        return source_archive_paths()

    git_executable = shutil.which("git")
    if git_executable is None:
        raise SystemExit("git is required to check a repository working tree")
    completed = subprocess.run(
        [
            git_executable,
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
        ],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        check=True,
    )
    return [
        REPOSITORY_ROOT / raw_path.decode("utf-8")
        for raw_path in completed.stdout.split(b"\0")
        if raw_path and not (REPOSITORY_ROOT / raw_path.decode("utf-8")).is_symlink()
    ]


def text_policy_issues(path: Path) -> list[str]:
    """Return policy failures for one non-binary repository file."""
    relative_path = path.relative_to(REPOSITORY_ROOT)
    raw_content = path.read_bytes()
    try:
        content = raw_content.decode("utf-8")
    except UnicodeDecodeError as error:
        return [f"{relative_path}: is not valid UTF-8 ({error})"]

    issues: list[str] = []
    if "\r" in content:
        issues.append(f"{relative_path}: contains a carriage return; use LF endings")
    if raw_content and not raw_content.endswith(b"\n"):
        issues.append(f"{relative_path}: is missing its final newline")

    if path.name[-4:].lower() not in TRAILING_WHITESPACE_EXEMPT_SUFFIXES:
        for line_number, line in enumerate(content.splitlines(), start=1):
            if line.endswith((" ", "\t")):
                issues.append(f"{relative_path}:{line_number}: trailing whitespace")
    return issues


def main() -> None:
    """Check every repository text file and fail with actionable diagnostics."""
    paths = repository_paths()
    issues = [
        issue
        for path in paths
        if path.suffix not in BINARY_SUFFIXES
        for issue in text_policy_issues(path)
    ]
    if issues:
        raise SystemExit("\n".join(issues))
    print(f"Repository text policy passed for {len(paths)} files.")


if __name__ == "__main__":
    main()
