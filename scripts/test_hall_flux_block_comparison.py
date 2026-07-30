#!/usr/bin/env python3
"""Bounded regression for adjacent Hall flux block comparison."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMPARATOR = ROOT / "scripts" / "compare_hall_flux_blocks.py"
METRICS = (
    "left_electron_charge_rate_a",
    "left_ion_charge_rate_a",
    "left_net_charge_rate_a",
    "reverse_step_fraction",
    "reverse_mean_demand_macroparticles",
)


def write_report(
    path: Path, start: int, end: int, means: float, cv: float
) -> None:
    value = {
        "schema_version": 1,
        "start_step": start,
        "end_step": end,
        "window_steps": 10,
        "summary": {
            metric: {
                "mean": means,
                "minimum": means,
                "maximum": means,
                "coefficient_of_variation": cv,
            }
            for metric in METRICS
        },
    }
    path.write_text(json.dumps(value), encoding="utf-8")


def run(
    previous: Path, current: Path, report: Path
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable, str(COMPARATOR), str(previous), str(current),
            "--max-mean-relative-change", "0.1",
            "--max-window-cv", "0.1",
            "--report", str(report),
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_flux_blocks_"
    ) as temporary:
        work = Path(temporary)
        previous = work / "previous.json"
        passing = work / "passing.json"
        failing = work / "failing.json"
        write_report(previous, 0, 100, 10.0, 0.05)
        write_report(passing, 100, 200, 10.5, 0.08)
        write_report(failing, 100, 200, 12.0, 0.12)

        passing_output = work / "passing-output.json"
        result = run(previous, passing, passing_output)
        value = json.loads(passing_output.read_text(encoding="utf-8"))
        require(
            result.returncode == 0
            and value["stationarity_screen_passed"]
            and all(item["passed"] for item in value["metrics"].values())
            and value["criteria"]["max_mean_absolute_relative_change"] == 0.1
            and value["physics_claim"] == "none",
            f"passing block comparison failed: {result.stderr}",
        )

        failing_output = work / "failing-output.json"
        result = run(previous, failing, failing_output)
        value = json.loads(failing_output.read_text(encoding="utf-8"))
        require(
            result.returncode == 1
            and not value["stationarity_screen_passed"]
            and not any(item["passed"] for item in value["metrics"].values()),
            f"failing block comparison was accepted: {result.stderr}",
        )

        nonadjacent = work / "nonadjacent.json"
        write_report(nonadjacent, 210, 310, 10.0, 0.01)
        result = run(
            previous, nonadjacent, work / "nonadjacent-output.json"
        )
        require(
            result.returncode == 2
            and "adjacent blocks" in result.stderr,
            "non-adjacent reports were accepted",
        )

    print("Hall flux block comparison passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
