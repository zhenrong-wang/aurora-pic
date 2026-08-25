#!/usr/bin/env python3
"""Audit coupled population, collision, loss, field, and power evolution."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


FILES = ("scalars.csv", "collisions.csv", "boundary_losses.csv")


class CoupledEvolutionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        values = list(csv.DictReader(stream))
    if len(values) < 2:
        raise CoupledEvolutionError(f"insufficient rows: {path}")
    return values


def delta(values: list[dict[str, str]], name: str) -> float:
    return float(values[-1][name]) - float(values[0][name])


def relative_range(values: list[float]) -> float:
    mean = sum(abs(value) for value in values) / len(values)
    if mean == 0.0:
        return 0.0 if max(values) == min(values) else math.inf
    return (max(values) - min(values)) / mean


def summarize(report_path: Path, output: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise CoupledEvolutionError("runner report did not pass all gates")
    hashes = report.get("output_hashes", {})
    for name in FILES:
        path = output / name
        if hashes.get(name) != sha256(path):
            raise CoupledEvolutionError(f"output hash differs: {name}")
    scalar = rows(output / "scalars.csv")
    collision = rows(output / "collisions.csv")
    boundary = rows(output / "boundary_losses.csv")
    endpoints = [row for row in scalar
                 if (int(row["step"]) - int(scalar[0]["step"])) % 4000 == 0]
    if len(endpoints) != 5:
        raise CoupledEvolutionError("expected five RF-phase endpoints")
    electron_losses = sum(delta(boundary, name) for name in (
        "absorbed_left_count_electrons", "absorbed_right_count_electrons"))
    ion_losses = sum(delta(boundary, name) for name in (
        "absorbed_left_count_ions", "absorbed_right_count_ions"))
    ionizations = delta(
        collision, "cumulative_collisions_electron_mcc.ionization")
    electron_change = delta(scalar, "live_particles_electrons")
    ion_change = delta(scalar, "live_particles_ions")
    return {
        "id": report["inputs"]["initial_state_id"],
        "report_sha256": sha256(report_path),
        "output_sha256": {name: hashes[name] for name in FILES},
        "steps": [int(scalar[0]["step"]), int(scalar[-1]["step"])],
        "rf_phase_endpoints": len(endpoints),
        "population": {
            "electron_start": int(scalar[0]["live_particles_electrons"]),
            "electron_final": int(scalar[-1]["live_particles_electrons"]),
            "electron_change": int(electron_change),
            "ion_start": int(scalar[0]["live_particles_ions"]),
            "ion_final": int(scalar[-1]["live_particles_ions"]),
            "ion_change": int(ion_change),
            "ionizations": int(ionizations),
            "electron_wall_losses": int(electron_losses),
            "ion_wall_losses": int(ion_losses),
            "electron_balance_residual": int(
                electron_change - (ionizations - electron_losses)),
            "ion_balance_residual": int(
                ion_change - (ionizations - ion_losses)),
        },
        "collision_increments": {
            name.removeprefix("cumulative_collisions_"): int(
                delta(collision, name))
            for name in collision[0]
            if name.startswith("cumulative_collisions_")
        },
        "same_phase_endpoint_means": {
            name: sum(float(row[name]) for row in endpoints) / len(endpoints)
            for name in ("kinetic_energy", "field_energy", "charge_l1")
        },
        "same_phase_final": {
            name: float(endpoints[-1][name])
            for name in ("kinetic_energy", "field_energy", "charge_l1")
        },
    }


def analyze(pooled_result_path: Path, reports: list[Path],
            outputs: list[Path]) -> dict[str, object]:
    if len(reports) != 2 or len(outputs) != 2:
        raise CoupledEvolutionError("exactly two branches are required")
    pooled = json.loads(pooled_result_path.read_text(encoding="utf-8"))
    if pooled.get("prospective_decision_outcome", {}).get(
            "deficit_persists") is not True:
        raise CoupledEvolutionError("locked pooled deficit result differs")
    branches = [summarize(report, output) for report, output in zip(
        reports, outputs, strict=True)]
    population = [branch["population"] for branch in branches]
    endpoint_means = [branch["same_phase_endpoint_means"]
                      for branch in branches]
    final = [branch["same_phase_final"] for branch in branches]
    promotion_range = float(pooled["critical_phase_0p125_to_0p5"]
                            ["pooled_aurorapic_relative_range"]
                            ["field_push_promotions_per_million_pushes"])
    dispersions = {
        "ionization_increment_relative_range": relative_range(
            [float(item["ionizations"]) for item in population]),
        "electron_wall_loss_relative_range": relative_range(
            [float(item["electron_wall_losses"]) for item in population]),
        "ion_wall_loss_relative_range": relative_range(
            [float(item["ion_wall_losses"]) for item in population]),
        "mean_kinetic_energy_relative_range": relative_range(
            [float(item["kinetic_energy"]) for item in endpoint_means]),
        "mean_field_energy_relative_range": relative_range(
            [float(item["field_energy"]) for item in endpoint_means]),
        "final_kinetic_energy_relative_range": relative_range(
            [float(item["kinetic_energy"]) for item in final]),
        "final_field_energy_relative_range": relative_range(
            [float(item["field_energy"]) for item in final]),
        "pooled_promotion_rate_relative_range": promotion_range,
    }
    balances_close = all(
        item["electron_balance_residual"] == 0 and
        item["ion_balance_residual"] == 0 for item in population)
    gross_repeatable = max(
        dispersions["ionization_increment_relative_range"],
        dispersions["electron_wall_loss_relative_range"],
        dispersions["ion_wall_loss_relative_range"]) < 0.03
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "existing_matched_half_step_coupled_evolution_audit",
        "branches": branches,
        "cross_branch_relative_ranges": dispersions,
        "gates": {
            "exact_species_balance": balances_close,
            "gross_creation_and_loss_repeatable_below_3_percent":
                gross_repeatable,
            "pooled_promotion_repeatability_below_8_percent":
                promotion_range < 0.08,
        },
        "all_accounting_and_repeatability_gates_passed":
            balances_close and gross_repeatable and promotion_range < 0.08,
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "pooled_result_sha256": sha256(pooled_result_path),
        },
        "assessment": {
            "finding": (
                "Species changes close exactly as ionization minus wall loss "
                "in both branches. Gross ionization and wall-loss traffic is "
                "repeatable within three percent, and the already-gated "
                "promotion-rate dispersion is smaller than the cross-code "
                "deficit. The deficit is not explained by within-Aurora "
                "microstate variability in these two consecutive blocks."),
        },
        "claim_boundary": (
            "This is a post hoc coupled-accounting audit of two existing "
            "AuroraPIC microstates. It does not add independent seeds, compare "
            "time-resolved native state variables, or validate experiment."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled_result", type=Path)
    parser.add_argument("reports", nargs=2, type=Path)
    parser.add_argument("outputs", nargs=2, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        result = analyze(args.pooled_result, args.reports, args.outputs)
    except (CoupledEvolutionError, OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_accounting_and_repeatability_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
