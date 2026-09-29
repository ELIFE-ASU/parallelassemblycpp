"""Plot measured calculation times; run after graph_repair_speed.py completes."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parent
report = json.loads((ROOT / "speed.json").read_text())
if "finished_utc" not in report:
    raise RuntimeError("Benchmark has not finished")
cases = [
    case
    for case in report["cases"]
    if case["completed_comparison"] and "scaling" not in case["suites"].split(",")
]
cases.sort(key=lambda case: case["methods"]["exact"]["median_algorithm_seconds"])
fig, ax = plt.subplots(figsize=(11, 8))
fig.patch.set_facecolor("white")
for row, case in enumerate(cases):
    exact_ms = case["methods"]["exact"]["median_algorithm_seconds"] * 1000
    graph_ms = case["methods"]["graph"]["median_algorithm_seconds"] * 1000
    ax.plot([graph_ms, exact_ms], [row, row], color="#c3d0d8", linewidth=2, zorder=1)
    ax.scatter(
        exact_ms,
        row,
        color="#203e60",
        s=44,
        zorder=2,
        label="Exact search" if row == 0 else None,
    )
    ax.scatter(
        graph_ms,
        row,
        color="#078c78",
        s=44,
        zorder=3,
        label="Graph-pair upper bound" if row == 0 else None,
    )
ax.set_yticks(
    range(len(cases)),
    [
        f"{case['name']}   "
        f"[{case['exact_assembly_index']} → {case['graph_upper_bound']}]"
        for case in cases
    ],
)
ax.set_xscale("log")
ax.set_xlabel(
    "Calculation time in milliseconds (log scale; lower is faster)", labelpad=10
)
ax.set_title(
    "Graph-pair bound versus exact molecular assembly search",
    loc="left",
    pad=28,
    fontsize=15,
    fontweight="bold",
)
ax.text(
    0,
    1.015,
    "Labels show exact index → upper bound. Six samples per method; one CPU core.",
    transform=ax.transAxes,
    fontsize=10,
    color="#4a5966",
)
ax.grid(axis="x", which="major", color="#e0e6ea")
ax.set_axisbelow(True)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color("#aab7c2")
ax.tick_params(axis="y", length=0)
ax.legend(loc="lower right", frameon=False)
fig.text(
    0.02,
    0.012,
    "Fresh processes, same Release build. Parsing/startup excluded. "
    "Paclitaxel timed out and is omitted here.",
    fontsize=9,
    color="#4a5966",
)
fig.tight_layout(rect=(0, 0.035, 1, 1))
preview_dir = ROOT.parents[1] / "build" / "graph-repair"
preview_dir.mkdir(parents=True, exist_ok=True)
fig.savefig(preview_dir / "speed.png", dpi=180)
svg_path = ROOT / "speed.svg"
fig.savefig(svg_path)
svg_path.write_text(
    "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n"
)
