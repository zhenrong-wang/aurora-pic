#!/usr/bin/env python3
"""Focused tests for the coupled-evolution audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile

from analyze_edupic_coupled_evolution import (
    CoupledEvolutionError, analyze, sha256,
)


def write_csv(path: Path, names: list[str], values: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(names)
        writer.writerows(values)


def make_branch(root: Path, state: str, ionizations: int,
                electron_losses: int, ion_losses: int) -> tuple[Path, Path]:
    output = root / state
    output.mkdir()
    scalar_names = [
        "step", "live_particles_electrons", "live_particles_ions",
        "kinetic_energy", "field_energy", "charge_l1"]
    electron_change = ionizations - electron_losses
    ion_change = ionizations - ion_losses
    scalar_rows = [[
        step, 1000 + electron_change * index // 4,
        1100 + ion_change * index // 4,
        10.0 + 0.01 * index, 20.0 + 0.02 * index,
        1.0 + 0.001 * index]
        for index, step in enumerate(range(0, 16001, 4000))]
    write_csv(output / "scalars.csv", scalar_names, scalar_rows)
    collision_names = [
        "step", "cumulative_collisions_electron_mcc.elastic",
        "cumulative_collisions_electron_mcc.excitation",
        "cumulative_collisions_electron_mcc.ionization",
        "cumulative_collisions_ion_mcc.isotropic",
        "cumulative_collisions_ion_mcc.backward"]
    write_csv(output / "collisions.csv", collision_names, [
        [0, 0, 0, 0, 0, 0],
        [16000, 10000, 100, ionizations, 200, 190]])
    boundary_names = [
        "step", "absorbed_left_count_electrons",
        "absorbed_right_count_electrons", "absorbed_left_count_ions",
        "absorbed_right_count_ions"]
    write_csv(output / "boundary_losses.csv", boundary_names, [
        [0, 0, 0, 0, 0],
        [16000, electron_losses // 2,
         electron_losses - electron_losses // 2,
         ion_losses // 2, ion_losses - ion_losses // 2]])
    report = root / f"{state}.json"
    report.write_text(json.dumps({
        "all_gates_passed": True,
        "inputs": {"initial_state_id": state},
        "output_hashes": {
            name: sha256(output / name)
            for name in ("scalars.csv", "collisions.csv",
                         "boundary_losses.csv")},
    }) + "\n", encoding="utf-8")
    return report, output


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        pooled = root / "pooled.json"
        pooled.write_text(json.dumps({
            "prospective_decision_outcome": {"deficit_persists": True},
            "critical_phase_0p125_to_0p5": {
                "pooled_aurorapic_relative_range": {
                    "field_push_promotions_per_million_pushes": 0.05}},
        }) + "\n", encoding="utf-8")
        first_report, first_output = make_branch(
            root, "first", 900, 1000, 1010)
        second_report, second_output = make_branch(
            root, "second", 910, 1010, 1020)
        result = analyze(
            pooled, [first_report, second_report],
            [first_output, second_output])
        assert result["all_accounting_and_repeatability_gates_passed"]
        assert all(branch["population"]["electron_balance_residual"] == 0
                   for branch in result["branches"])
        (first_output / "scalars.csv").write_text(
            "corrupt\n", encoding="utf-8")
        try:
            analyze(pooled, [first_report, second_report],
                    [first_output, second_output])
        except CoupledEvolutionError:
            pass
        else:
            raise AssertionError("modified bound output was accepted")
    print("coupled-evolution analyzer tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
