"""Discover the CPU capacity and physical-first placement available to this process."""

from __future__ import annotations

import contextlib
import os
import platform
from dataclasses import dataclass
from pathlib import Path

SYSFS_CPU_ROOT = Path("/sys/devices/system/cpu")


@dataclass(frozen=True)
class CpuTopology:
    logical_cpus: int
    physical_cores: int | None
    cpu_order: tuple[int, ...] | None
    model: str
    affinity_cpus: tuple[int, ...] | None
    slurm_cpus_per_task: int | None
    source: str


def cpu_model() -> str:
    """Read the processor model without invoking a platform command."""
    try:
        with Path("/proc/cpuinfo").open(encoding="utf-8") as stream:
            for line in stream:
                key, separator, value = line.partition(":")
                if (
                    separator
                    and key.strip() in {"model name", "Processor"}
                    and (description := value.strip())
                ):
                    return description
    except OSError:
        pass
    return platform.machine() or "unknown"


def slurm_cpu_limit() -> int | None:
    value = os.environ.get("SLURM_CPUS_PER_TASK")
    if value is None:
        return None
    if not value.isascii() or not value.isdigit() or int(value) < 1:
        raise ValueError(
            "SLURM_CPUS_PER_TASK must be a positive integer; "
            "correct or unset it before benchmarking"
        )
    return int(value)


def physical_cpu_groups(
    allowed_cpus: tuple[int, ...], sysfs_root: Path
) -> tuple[tuple[int, ...], ...] | None:
    """Group permitted Linux CPU IDs, requiring complete core identity evidence."""
    groups: dict[tuple[int, int | None, int], list[int]] = {}
    known_dies: set[bool] = set()
    for cpu in allowed_cpus:
        topology = sysfs_root / f"cpu{cpu}" / "topology"
        try:
            package = int((topology / "physical_package_id").read_text())
            core = int((topology / "core_id").read_text())
            try:
                die = int((topology / "die_id").read_text())
            except FileNotFoundError:
                die = None
        except (OSError, ValueError):
            return None
        if package < 0 or core < 0:
            return None
        if die is not None and die < 0:
            die = None
        known_dies.add(die is not None)
        groups.setdefault((package, die, core), []).append(cpu)
    if len(known_dies) > 1:
        return None
    return tuple(tuple(group) for group in groups.values())


def available_cpu_count() -> tuple[int, str]:
    """Prefer a process-aware count when the runtime provides one."""
    process_count = getattr(os, "process_cpu_count", None)
    if callable(process_count):
        try:
            count = process_count()
        except (OSError, NotImplementedError):
            count = None
        if count is not None:
            if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                raise ValueError("os.process_cpu_count returned no usable CPUs")
            return count, "os.process_cpu_count"
    count = os.cpu_count()
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("could not determine an available CPU count")
    return count, "os.cpu_count"


def detect_cpu_topology() -> CpuTopology:
    """Apply affinity and scheduler limits before reporting capacity or placement."""
    slurm_limit = slurm_cpu_limit()
    get_affinity = getattr(os, "sched_getaffinity", None)
    allowed = None
    if callable(get_affinity):
        with contextlib.suppress(OSError, NotImplementedError):
            allowed = tuple(sorted(get_affinity(0)))
    physical_cores = None
    order = None
    if allowed is not None:
        if not allowed:
            raise ValueError("process CPU affinity contains no available CPUs")
        logical_cpus = len(allowed)
        groups = physical_cpu_groups(allowed, SYSFS_CPU_ROOT)
        if groups is None:
            order = allowed
            source = "os.sched_getaffinity (physical topology unavailable)"
        else:
            physical_cores = len(groups)
            order = tuple(
                group[sibling]
                for sibling in range(max(map(len, groups)))
                for group in groups
                if sibling < len(group)
            )
            source = "os.sched_getaffinity + Linux sysfs"
    else:
        logical_cpus, source = available_cpu_count()
    if slurm_limit is not None:
        logical_cpus = min(logical_cpus, slurm_limit)
        if order is not None:
            order = order[:logical_cpus]
        if physical_cores is not None:
            physical_cores = min(physical_cores, logical_cpus)
        source += "; SLURM_CPUS_PER_TASK"
    return CpuTopology(
        logical_cpus=logical_cpus,
        physical_cores=physical_cores,
        cpu_order=order,
        model=cpu_model(),
        affinity_cpus=allowed,
        slurm_cpus_per_task=slurm_limit,
        source=source,
    )
