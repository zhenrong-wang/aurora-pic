#!/usr/bin/env python3
"""Regression tests for the AuroraPIC eduPIC horizon screen."""

from __future__ import annotations

import csv
from pathlib import Path
import tempfile

from extend_aurorapic_edupic_horizon import (
    HorizonError,
    endpoint,
    report_end_cycle,
    solver_command,
    stationarity,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    require(
        solver_command(Path("aurorapic_cli"), Path("stage.cfg")) == [
            "aurorapic_cli", "--allow-large-run",
            "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", "stage.cfg",
        ],
        "horizon solver launch lost its explicit CLI acknowledgement",
    )
    require(
        report_end_cycle({
            "completed_through_cycle": 4, "all_gates_passed": True
        }) == 4
        and report_end_cycle({
            "block": {"end_cycle": 8, "hard_safety_gates_passed": True}
        }) == 8,
        "safe report chaining lost a pilot or horizon report",
    )
    try:
        report_end_cycle({
            "block": {"end_cycle": 8, "hard_safety_gates_passed": False}
        })
    except HorizonError:
        pass
    else:
        raise RuntimeError("unsafe horizon report was accepted for chaining")
    with tempfile.TemporaryDirectory(prefix="aurorapic-edupic-horizon-") as text:
        output = Path(text)
        write_csv(
            output / "scalars.csv",
            [
                "step", "live_particles", "live_particles_electrons",
                "live_particles_ions", "field_energy",
            ],
            [
                {"step": 16000, "live_particles": 7989,
                 "live_particles_electrons": 2483,
                 "live_particles_ions": 5506, "field_energy": 2.8e-5},
                {"step": 20000, "live_particles": 8000,
                 "live_particles_electrons": 2490,
                 "live_particles_ions": 5510, "field_energy": 2.9e-5},
            ],
        )
        key = "cumulative_collisions_electron_mcc.ionization"
        write_csv(
            output / "collisions.csv", ["step", key],
            [{"step": 16000, key: 4719}, {"step": 20000, key: 5519}],
        )
        write_csv(
            output / "fields_20000.csv", ["x", "E"],
            [{"x": 0, "E": -30000}, {"x": 0.025, "E": 35000}],
        )
        item = endpoint(output, 5)
        require(
            item["total_particles"] == 8000
            and item["ionization_pairs_in_cycle"] == 800
            and item["maximum_sampled_absolute_field_V_m"] == 35000,
            "horizon endpoint extraction lost population, collision, or field data",
        )

    stable = [
        {
            "cycle": cycle,
            "total_particles": 10000 + (cycle - 4) * 20,
            "field_energy_J_m-2": 3e-5 + (cycle - 4) * 1e-8,
            "maximum_sampled_absolute_field_V_m": 35000 + (cycle - 4) * 20,
            "ionization_pairs_in_cycle": 800 + ((cycle % 2) * 4),
        }
        for cycle in range(4, 9)
    ]
    stable_result = stationarity(stable)
    require(stable_result["passed"], "stable synthetic horizon failed")
    trending = [dict(item) for item in stable]
    for index, item in enumerate(trending):
        item["total_particles"] = 10000 * (1.2 ** index)
        item["field_energy_J_m-2"] = 3e-5 * (1.1 ** index)
        item["maximum_sampled_absolute_field_V_m"] = 35000 * (1.05 ** index)
        item["ionization_pairs_in_cycle"] = 700 + 100 * index
    trending_result = stationarity(trending)
    require(
        not trending_result["passed"]
        and not trending_result["gates"]["total_population_slope"]
        and not trending_result["gates"]["field_energy_slope"],
        "trending synthetic horizon passed stationarity",
    )
    print("AuroraPIC eduPIC horizon stationarity regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
