#!/usr/bin/env python3
"""Regression tests for the bounded AuroraPIC eduPIC pilot runner."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile

from run_aurorapic_edupic_pilot import PilotError, analyze_stage, stage_deck


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    base = """steps = 4000
output_interval = 100
output_dir = old
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 1
spatial_average_end_step = 4000
spatial_average_rf_frequency = 13560000
spatial_average_rf_cycles = 1
spatial_average_phase_bins = 16
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 1000000
checkpoint_output = true
checkpoint_interval = 4000
[species.electrons]
mass = 1
"""
    deck = stage_deck(
        base, 1, Path("/tmp/pilot-cycle-2"), Path("/tmp/checkpoint_4000.apc")
    )
    require(
        "steps = 8000\n" in deck
        and "spatial_average_reset_on_restart = true\n" in deck
        and "spatial_average_start_step = 4001\n" in deck
        and "spatial_average_end_step = 8000\n" in deck
        and "restart_path = /tmp/checkpoint_4000.apc\n" in deck,
        "continuation deck lost its absolute-step or averaging-reset contract",
    )

    with tempfile.TemporaryDirectory(prefix="aurorapic-edupic-pilot-") as text:
        output = Path(text)
        scalar_fields = [
            "step", "live_particles_electrons", "live_particles_ions",
            "kinetic_energy", "field_energy", "total_energy", "charge_l1",
        ]
        write_csv(output / "scalars.csv", scalar_fields, [
            {"step": 4000, "live_particles_electrons": 863,
             "live_particles_ions": 3369, "kinetic_energy": 1,
             "field_energy": 2, "total_energy": 3, "charge_l1": 4},
            {"step": 8000, "live_particles_electrons": 1000,
             "live_particles_ions": 5000, "kinetic_energy": 2,
             "field_energy": 3, "total_energy": 5, "charge_l1": 6},
        ])
        collision_key = "cumulative_collisions_electron_mcc.ionization"
        write_csv(output / "collisions.csv", ["step", collision_key], [
            {"step": 4000, collision_key: 2381},
            {"step": 8000, collision_key: 4381},
        ])
        boundary_fields = ["step"] + [
            f"absorbed_{side}_count_{species}"
            for species in ("electrons", "ions")
            for side in ("left", "right")
        ]
        write_csv(output / "boundary_losses.csv", boundary_fields, [
            {
                "step": 4000,
                "absorbed_left_count_electrons": 1486,
                "absorbed_right_count_electrons": 1032,
                "absorbed_left_count_ions": 8,
                "absorbed_right_count_ions": 4,
            },
            {
                "step": 8000,
                "absorbed_left_count_electrons": 3349,
                "absorbed_right_count_electrons": 1032,
                "absorbed_left_count_ions": 377,
                "absorbed_right_count_ions": 4,
            },
        ])
        write_csv(output / "fields_8000.csv", ["x", "E"], [
            {"x": 0, "E": -10000}, {"x": 0.025, "E": 20000}
        ])
        (output / "energy-budget.json").write_text(json.dumps({
            "relative_closure_residual": 1e-15, "passes": True
        }), encoding="utf-8")
        (output / "spatial-collision.json").write_text(json.dumps({
            "closure": {
                "maximum_spatial_global_residual_J_m-2": 1e-20,
                "maximum_phase_spatial_residual_J_m-2": 2e-20,
                "passes": True,
            }
        }), encoding="utf-8")
        result = analyze_stage(
            output, 1, {"electrons": 863, "ions": 3369}, 2_000_000
        )
        require(
            result["passes"]
            and result["population"]["ionization_pairs"] == 2000
            and result["population"]["electron_wall_losses"] == 1863
            and result["population"]["ion_wall_losses"] == 369
            and result["collision_totals"][
                "collisions_electron_mcc.ionization"
            ] == 2000
            and result["state"]["final_total_energy_J_m-2"] == 5,
            "valid continuation diagnostics did not pass exact balance gates",
        )
        bad_rows = [
            {"step": 4000, "live_particles_electrons": 863,
             "live_particles_ions": 3369, "kinetic_energy": 1,
             "field_energy": 2, "total_energy": 3, "charge_l1": 4},
            {"step": 8000, "live_particles_electrons": 1001,
             "live_particles_ions": 5000, "kinetic_energy": 2,
             "field_energy": 3, "total_energy": 5, "charge_l1": 6},
        ]
        write_csv(output / "scalars.csv", scalar_fields, bad_rows)
        try:
            analyze_stage(
                output, 1, {"electrons": 863, "ions": 3369}, 2_000_000
            )
        except PilotError as error:
            require("population balance" in str(error), "unexpected rejection")
        else:
            raise RuntimeError("pilot accepted a non-closing population balance")
    print("AuroraPIC eduPIC multi-cycle pilot regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
