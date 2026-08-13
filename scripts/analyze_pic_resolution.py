#!/usr/bin/env python3
"""Screen explicit 1D PIC resolution using fresh phase-space diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


EPSILON_0 = 8.8541878128e-12
ELEMENTARY_CHARGE = 1.602176634e-19
ELECTRON_MASS = 9.1093837015e-31


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"empty diagnostic table: {path}")
    return result


def sections(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {"global": {}}
    current = result["global"]
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            name = line[1:-1].strip()
            if not name or name in result:
                raise ValueError("configuration contains duplicate/empty section")
            current = result.setdefault(name, {})
            continue
        key, separator, value = line.partition("=")
        if not separator or not key.strip() or key.strip() in current:
            raise ValueError("configuration contains malformed/duplicate key")
        current[key.strip()] = value.strip()
    return result


def positive(value: str, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return parsed


def analyze(output: Path, config: Path,
            populated_density_fraction: float = 0.01) -> dict[str, object]:
    if not 0.0 < populated_density_fraction < 1.0:
        raise ValueError("populated density fraction must lie between zero and one")
    cfg = sections(config)
    global_cfg = cfg["global"]
    if (global_cfg.get("units") != "si" or
            int(global_cfg.get("dimension", "0")) != 1 or
            int(global_cfg.get("velocity_dimensions", "0")) != 3):
        raise ValueError("resolution audit requires an SI 1D3V configuration")
    nodes = int(global_cfg["nx"])
    if nodes < 2:
        raise ValueError("nx must contain at least two nodes")
    cells = nodes - 1
    length = positive(global_cfg["length"], "length")
    dt = positive(global_cfg["dt"], "dt")
    dx = length / cells

    moments = [row for row in table(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    maximum_density = max(float(row["number_density_mean_m-3"])
                          for row in moments)
    populated = []
    for row in moments:
        density = float(row["number_density_mean_m-3"])
        temperature = float(row["drift_separated_temperature_eV"])
        if (density >= populated_density_fraction * maximum_density and
                temperature > 0.0 and math.isfinite(temperature)):
            debye = math.sqrt(
                EPSILON_0 * temperature / (density * ELEMENTARY_CHARGE))
            populated.append((debye, density, temperature))
    if not populated:
        raise ValueError("no populated finite-temperature electron samples")
    minimum_debye, debye_density, debye_temperature = min(populated)
    maximum_plasma_frequency = math.sqrt(
        maximum_density * ELEMENTARY_CHARGE**2 / (EPSILON_0 * ELECTRON_MASS))

    scalar_rows = table(output / "scalars.csv")
    mean_ppc = {}
    for species in ("electrons", "ions"):
        field = f"live_particles_{species}"
        mean_ppc[species] = math.fsum(float(row[field]) for row in scalar_rows) / (
            len(scalar_rows) * cells)

    species_multiplier = {
        name.removeprefix("species."): int(values.get("timestep_multiplier", "1"))
        for name, values in cfg.items() if name.startswith("species.")}
    collision_numbers = {}
    for name, values in cfg.items():
        if not name.startswith("collisions.") or ".channel." in name:
            continue
        species = values.get("species")
        if species not in species_multiplier:
            raise ValueError(f"collision section {name} has unknown species")
        maximum_frequency = positive(values["max_frequency"], "max_frequency")
        collision_numbers[name.removeprefix("collisions.")] = (
            maximum_frequency * dt * species_multiplier[species])

    frequency = positive(global_cfg["spatial_average_rf_frequency"],
                         "spatial_average_rf_frequency")
    rf_steps = 1.0 / (frequency * dt)
    metrics = {
        "cell_width_m": dx,
        "minimum_populated_electron_debye_length_m": minimum_debye,
        "cell_width_to_debye_length": dx / minimum_debye,
        "debye_limiting_density_m-3": debye_density,
        "debye_limiting_temperature_eV": debye_temperature,
        "maximum_phase_electron_density_m-3": maximum_density,
        "maximum_electron_plasma_angular_frequency_s-1": maximum_plasma_frequency,
        "maximum_omega_pe_dt": maximum_plasma_frequency * dt,
        "rf_steps_per_cycle": rf_steps,
        "mean_macro_particles_per_cell": mean_ppc,
        "maximum_null_collision_numbers": collision_numbers,
    }
    thresholds = {
        "maximum_cell_width_to_debye_length": 1.0,
        "maximum_omega_pe_dt": 0.2,
        "minimum_rf_steps_per_cycle": 1000.0,
        "minimum_mean_macro_particles_per_cell_per_species": 100.0,
        "maximum_null_collision_number": 0.1,
    }
    gates = {
        "debye_length": metrics["cell_width_to_debye_length"] <= 1.0,
        "electron_plasma_frequency": metrics["maximum_omega_pe_dt"] <= 0.2,
        "rf_drive": rf_steps >= 1000.0,
        "particles_per_cell": min(mean_ppc.values()) >= 100.0,
        "null_collision_bounds": max(collision_numbers.values(), default=0.0) <= 0.1,
    }
    return {
        "schema_version": 1,
        "scope": "post_measurement_explicit_pic_resolution_adequacy_screen",
        "config_sha256": sha256(config),
        "input_hashes": {
            name: sha256(output / name) for name in
            ("spatial_phase_moments.csv", "scalars.csv")},
        "populated_electron_density_fraction": populated_density_fraction,
        "metrics": metrics,
        "engineering_thresholds": thresholds,
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "claim_boundary": (
            "These conservative project engineering screens identify obvious "
            "under-resolution. Passing is necessary evidence but is not a "
            "mesh, timestep, particle-count, or statistical convergence proof."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--populated-density-fraction", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.candidate_output.resolve(), args.config.resolve(),
                     args.populated_density_fraction)
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
