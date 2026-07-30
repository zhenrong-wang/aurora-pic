#!/usr/bin/env python3
"""Bounded regression for Hall boundary-flux stationarity analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_hall_flux_stationarity.py"


def write_csv(
    path: Path, fields: list[str], values: list[dict[str, object]]
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(values)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_flux_stationarity_"
    ) as temporary:
        work = Path(temporary)
        output = work / "output"
        output.mkdir()
        flux_rows: list[dict[str, object]] = []
        charge_rates = {
            ("electrons", "left"): -2.0,
            ("ions", "left"): 0.5,
            ("electrons", "right"): -0.1,
            ("ions", "right"): 0.7,
        }
        for step in range(1, 5):
            for species in ("electrons", "ions"):
                for boundary in ("left", "right"):
                    rate = charge_rates[(species, boundary)]
                    flux_rows.append({
                        "step": step,
                        "time": step,
                        "window_start_step": step - 1,
                        "window_duration": 1,
                        "species": species,
                        "boundary": boundary,
                        "represented_charge": rate,
                        "charge_rate": rate,
                    })
        write_csv(
            output / "boundary_flux.csv",
            [
                "step", "time", "window_start_step", "window_duration",
                "species", "boundary", "represented_charge",
                "charge_rate",
            ],
            flux_rows,
        )
        write_csv(
            output / "current_source.csv",
            [
                "step", "reverse_distribution_steps",
                "reverse_one_macro_steps", "reverse_two_macro_steps",
                "reverse_multi_macro_steps",
                "distributed_reverse_demand_macroparticles",
            ],
            [
                {
                    "step": 0, "reverse_distribution_steps": 0,
                    "reverse_one_macro_steps": 0,
                    "reverse_two_macro_steps": 0,
                    "reverse_multi_macro_steps": 0,
                    "distributed_reverse_demand_macroparticles": 0,
                },
                {
                    "step": 2, "reverse_distribution_steps": 1,
                    "reverse_one_macro_steps": 1,
                    "reverse_two_macro_steps": 0,
                    "reverse_multi_macro_steps": 0,
                    "distributed_reverse_demand_macroparticles": 1,
                },
                {
                    "step": 4, "reverse_distribution_steps": 3,
                    "reverse_one_macro_steps": 2,
                    "reverse_two_macro_steps": 1,
                    "reverse_multi_macro_steps": 0,
                    "distributed_reverse_demand_macroparticles": 4,
                },
            ],
        )
        report_path = work / "report.json"
        result = subprocess.run(
            [
                sys.executable, str(ANALYZER), str(output),
                "--start-step", "0", "--end-step", "4",
                "--window-steps", "2", "--report", str(report_path),
            ],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        report = json.loads(report_path.read_text(encoding="utf-8"))
        require(
            result.returncode == 0
            and len(report["windows"]) == 2
            and report["windows"][0]["boundaries"]["left"][
                "net_charge_rate_a"
            ] == -1.5
            and report["windows"][1]["reverse"]["steps"] == 2
            and report["summary"]["reverse_step_fraction"]["mean"]
                == 0.75
            and report["physics_claim"] == "none",
            f"stationarity analysis failed: {result.stderr}",
        )
    print("Hall flux stationarity analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
