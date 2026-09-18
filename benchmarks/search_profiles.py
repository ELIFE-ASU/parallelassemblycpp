"""Validate and summarize additive search-work and cache-capacity telemetry."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path


STATE_COUNTERS = (
    "states_expanded",
    "states_pruned",
    "states_bound_pruned",
    "duplicate_classes_pruned",
    "occurrence_pairs_pruned",
    "fragment_pair_blocks_pruned",
)
CACHE_COMPONENTS = (
    "canonical_mask_retained_bytes",
    "canonical_graph_retained_bytes",
    "canonical_tree_retained_bytes",
    "assembly_state_retained_bytes",
    "residual_decomposition_retained_bytes",
)
CACHE_MEASUREMENT = "retained_capacity_estimate_allocator_overhead_excluded"
TRAJECTORY_CLOCK = "steady_wall_since_search_start"


class SearchProfileError(ValueError):
    """Raised for inconsistent additive telemetry."""


def nonnegative_integer(value: object) -> bool:
    return type(value) is int and value >= 0


def validate_search_profile(telemetry: dict[str, object]) -> None:
    """Keep old schema-v1 reports valid, rejecting partial or inconsistent groups."""
    parallel = telemetry.get("parallel")
    workers = parallel["workers"] if isinstance(parallel, dict) else []
    aggregate = parallel["aggregate"] if isinstance(parallel, dict) else None
    records = [telemetry, *workers]
    counter_records = [*records, *([aggregate] if aggregate is not None else [])]
    if any(
        set(STATE_COUNTERS).intersection(record["counters"])
        for record in counter_records
    ):
        for record in counter_records:
            counters = record["counters"]
            if any(
                not nonnegative_integer(counters.get(name)) for name in STATE_COUNTERS
            ):
                raise SearchProfileError(
                    "invalid or missing search-state counter group"
                )
            if counters["states_bound_pruned"] > counters["states_pruned"]:
                raise SearchProfileError("bound-pruned states exceed pruned states")
        if aggregate is not None:
            for name in STATE_COUNTERS:
                total = sum(worker["counters"][name] for worker in workers)
                if (
                    telemetry["counters"][name] != total
                    or aggregate["counters"][name] != total
                ):
                    raise SearchProfileError(
                        f"search-state counter {name} does not match workers"
                    )

    cache_records = [*records]
    if aggregate is not None and "local_caches" in aggregate:
        cache_records.append(aggregate)
    if any("local_caches" in record for record in cache_records):
        for record in cache_records:
            cache = record.get("local_caches")
            if (
                not isinstance(cache, dict)
                or cache.get("measurement") != CACHE_MEASUREMENT
                or any(
                    not nonnegative_integer(cache.get(name))
                    for name in (*CACHE_COMPONENTS, "total_retained_bytes")
                )
            ):
                raise SearchProfileError(
                    "invalid or missing local-cache retained bytes"
                )
            if cache["total_retained_bytes"] != sum(
                cache[name] for name in CACHE_COMPONENTS
            ):
                raise SearchProfileError(
                    "local-cache retained-byte total does not match components"
                )
        if workers:
            for name in (*CACHE_COMPONENTS, "total_retained_bytes"):
                total = sum(worker["local_caches"][name] for worker in workers)
                if telemetry["local_caches"][name] != total or (
                    aggregate is not None
                    and "local_caches" in aggregate
                    and aggregate["local_caches"][name] != total
                ):
                    raise SearchProfileError(
                        f"local-cache {name} does not match workers"
                    )

    if (
        any("incumbent_trajectory" in record for record in records)
        or "trajectory_clock" in telemetry
    ):
        if telemetry.get("trajectory_clock") != TRAJECTORY_CLOCK:
            raise SearchProfileError("invalid incumbent trajectory clock")
        identities = {
            (worker["rank"], worker["global_worker_index"]) for worker in workers
        } or {(0, 0)}
        for record in records:
            events = record.get("incumbent_trajectory")
            if not isinstance(events, list):
                raise SearchProfileError("invalid or missing incumbent trajectory")
            previous_time, previous_index = -1, None
            for event in events:
                if not isinstance(event, dict) or any(
                    not nonnegative_integer(event.get(name))
                    for name in (
                        "elapsed_nanoseconds",
                        "assembly_index",
                        "rank",
                        "global_worker_index",
                    )
                ):
                    raise SearchProfileError("invalid incumbent trajectory event")
                elapsed, index = event["elapsed_nanoseconds"], event["assembly_index"]
                identity = (event["rank"], event["global_worker_index"])
                if identity not in identities or (
                    record is not telemetry
                    and identity != (record["rank"], record["global_worker_index"])
                ):
                    raise SearchProfileError(
                        "incumbent trajectory event has invalid worker identity"
                    )
                if elapsed < previous_time or (
                    previous_index is not None and index >= previous_index
                ):
                    raise SearchProfileError(
                        "incumbent trajectory must be time-ordered strict improvements"
                    )
                previous_time, previous_index = elapsed, index


def summary_fields(telemetry: object) -> dict[str, object]:
    """Return scalar columns; the original report keeps every trajectory event."""
    if not isinstance(telemetry, dict):
        return {}
    result = {}
    counters = telemetry.get("counters", {})
    for name in (*STATE_COUNTERS, "matching_visits"):
        if name in counters:
            result[name] = counters[name]
    caches = telemetry.get("caches", {})
    canonical = caches.get("canonical_mask", {})
    for name in ("hits", "misses", "hit_rate"):
        if name in canonical:
            result[f"canonical_mask_{name}"] = canonical[name]
    local_caches = telemetry.get("local_caches", {})
    for name in (*CACHE_COMPONENTS, "total_retained_bytes"):
        if name in local_caches:
            result[name] = local_caches[name]
    events = telemetry.get("incumbent_trajectory")
    if isinstance(events, list):
        result["incumbent_events"] = len(events)
        if events:
            result["initial_incumbent_index"] = events[0]["assembly_index"]
            result["final_incumbent_index"] = events[-1]["assembly_index"]
            result["first_incumbent_nanoseconds"] = events[0]["elapsed_nanoseconds"]
            result["final_incumbent_nanoseconds"] = events[-1]["elapsed_nanoseconds"]
        for name in ("trajectory_clock", "trajectory_index", "trajectory_scope"):
            if name in telemetry:
                result[name] = telemetry[name]
    return result


def placement_profiles(reports: Sequence[tuple[str, Path]]) -> list[dict[str, object]]:
    """Collect per-placement profiles without mixing them into timed samples."""
    profiles = []
    for placement, path in reports:
        report = json.loads(path.read_text(encoding="utf-8"))
        for case in report["cases"]:
            telemetry = case["candidate"].get("telemetry")
            if not isinstance(telemetry, dict):
                continue
            profiles.append(
                {
                    "placement": placement,
                    "case": case["name"],
                    "report": str(path),
                    "execution": report["execution"]["telemetry"],
                    **summary_fields(telemetry),
                    "incumbent_trajectory": telemetry.get("incumbent_trajectory"),
                    "local_caches": telemetry.get("local_caches"),
                }
            )
    return profiles


def write_placement_profiles(path: Path, reports: Sequence[tuple[str, Path]]) -> str:
    profiles = placement_profiles(reports)
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "collection": "one untimed telemetry run per placement",
                "placements": profiles,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    lines = ["\nUntimed search profiles by placement"]
    for profile in profiles:
        fields = [
            f"{name}={profile[name]}"
            for name in (
                "states_expanded",
                "states_pruned",
                "canonical_mask_misses",
                "final_incumbent_nanoseconds",
                "total_retained_bytes",
            )
            if name in profile
        ]
        lines.append(
            f"  {profile['placement']} {profile['case']}: " + ", ".join(fields)
        )
    lines.append(f"  Full trajectories and retained bytes: {path}")
    return "\n".join(lines) + "\n"
