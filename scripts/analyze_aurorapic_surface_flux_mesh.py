#!/usr/bin/env python3
"""Evaluate paired 400/799-node CCP surface-transport refinement."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_aurorapic_surface_flux_timestep import evaluate
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "4bea77b968db89ca6a2e066a599d3e85b99c480de2f0cb6e56b12bdaeb891f54")
ELEMENTARY_CHARGE_C = 1.602176634e-19
RF_FREQUENCY_HZ = 13.56e6
PHASES = 200


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def branch(output: Path, name: str, expected_hash: str
           ) -> tuple[dict[str, object], Path]:
    path = output.parent / "branch-report.json"
    if sha256(path) != expected_hash:
        raise ValueError(f"{name} branch report differs")
    report = json.loads(path.read_text(encoding="utf-8"))
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("branch") != name or
            report.get("all_gates_passed") is not True):
        raise ValueError(f"invalid {name} branch report")
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"{name} output differs: {filename}")
    return report, path


def transport(output: Path) -> dict[str, object]:
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    cycles = int(metadata["rf_cycles"])
    surfaces = [float(value) for value in
                metadata["phase_surface_flux_positions"]]
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    nodes = len(moments) // PHASES
    fields = rows(output / "spatial_phase_fields.csv")
    if (nodes * PHASES != len(moments) or len(fields) != len(moments) or
            surfaces != [0.005, 0.015]):
        raise ValueError("mesh surface-transport shape differs")
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    energy_density = [
        value * float(row["mean_kinetic_energy_eV"]) * ELEMENTARY_CHARGE_C
        for value, row in zip(density, moments)]
    electric_power = [
        -ELEMENTARY_CHARGE_C * value * float(moment["mean_velocity_x"]) *
        float(field["electric_field_mean_V_m"])
        for value, moment, field in zip(density, moments, fields)]
    phase_dt = 1.0 / RF_FREQUENCY_HZ / PHASES
    storage = [
        (energy_density[((phase + 1) % PHASES) * nodes + node] -
         energy_density[((phase - 1) % PHASES) * nodes + node]) /
        (2.0 * phase_dt)
        for phase in range(PHASES) for node in range(nodes)]
    collision = [0.0] * (PHASES * nodes)
    for row in rows(output / "spatial_phase_collision_power.csv"):
        if row["channel"].startswith("electron_mcc."):
            collision[int(row["phase_bin"]) * nodes + int(row["node"])] += (
                float(row["mean_power_density_W_m-3"]))
    inferred = [power + loss - change for power, loss, change in
                zip(electric_power, collision, storage)]
    positions = [float(row["x_m"]) for row in moments[:nodes]]
    dx = positions[1] - positions[0]
    interior = [node for node, position in enumerate(positions)
                if surfaces[0] <= position <= surfaces[1]]
    summary = {}
    overflow = 0.0
    for row in rows(output / "phase_surface_flux_summary.csv"):
        summary[(int(row["phase_bin"]), int(row["surface_id"]),
                 row["direction"])] = float(
                     row["kinetic_energy_flux_W_m-2"])
        overflow = max(overflow, float(row["overflow_fraction"]))
    tail: dict[tuple[int, int, str], float] = {}
    duration = cycles / RF_FREQUENCY_HZ / PHASES
    for row in rows(output / "phase_surface_flux.csv"):
        energy = float(row["energy_eV"])
        if energy < 15.8:
            continue
        key = (int(row["phase_bin"]), int(row["surface_id"]), row["direction"])
        tail[key] = tail.get(key, 0.0) + (
            float(row["represented_crossings"]) * energy *
            ELEMENTARY_CHARGE_C / duration)

    def window(first: int, end: int) -> dict[str, float]:
        inferred_value = math.fsum(
            inferred[phase * nodes + node] for phase in range(first, end)
            for node in interior) * dx / (end - first)
        net, tail_net = [], []
        for surface in range(2):
            net.append(math.fsum(
                summary[(phase, surface, "left_to_right")] -
                summary[(phase, surface, "right_to_left")]
                for phase in range(first, end)) / (end - first))
            tail_net.append(math.fsum(
                tail.get((phase, surface, "left_to_right"), 0.0) -
                tail.get((phase, surface, "right_to_left"), 0.0)
                for phase in range(first, end)) / (end - first))
        direct = net[1] - net[0]
        return {
            "direct_outward_energy_flux_divergence_W_m-2": direct,
            "inferred_outward_energy_flux_divergence_W_m-2": inferred_value,
            "relative_closure_error": abs(direct - inferred_value) /
                max(abs(direct), abs(inferred_value)),
            "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2":
                tail_net[1] - tail_net[0],
        }

    return {
        "nodes": nodes,
        "critical_phase_0p125_to_0p5": window(25, 100),
        "exceptional_phase_0p375_to_0p5": window(75, 100),
        "maximum_surface_overflow_fraction": overflow,
    }


def analyze(baseline_output: Path, refined_output: Path, rule_path: Path,
            baseline_hash: str, refined_hash: str) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux mesh rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_report, baseline_path = branch(
        baseline_output, "baseline_grid", baseline_hash)
    refined_report, refined_path = branch(
        refined_output, "refined_grid", refined_hash)
    baseline = transport(baseline_output)
    refined = transport(refined_output)
    if baseline["nodes"] != 400 or refined["nodes"] != 799:
        raise ValueError("mesh refinement ratio differs")
    limits = dict(rule["prospective_acceptance"])
    limits.pop("all_gates_required")
    metrics, gates = evaluate(baseline, refined, limits)
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_surface_transport_mesh_refinement_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "baseline_branch_report_sha256": sha256(baseline_path),
            "refined_branch_report_sha256": sha256(refined_path),
        },
        "mesh_refinement": {"baseline_nodes": baseline["nodes"],
                            "refined_nodes": refined["nodes"],
                            "fixed_total_particles": True},
        "baseline": baseline,
        "refined_grid": refined,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "assessment": {
            "exceptional_transport_mesh_stable": gates[
                "exceptional_octant_direct_flux"],
            "critical_window_has_larger_mesh_sensitivity": metrics[
                "critical_phase_direct_flux_relative_change"] > metrics[
                    "exceptional_octant_direct_flux_relative_change"],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("refined_output", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--baseline-report-sha256", required=True)
    parser.add_argument("--refined-report-sha256", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.baseline_output.resolve(), args.refined_output.resolve(),
        args.rule.resolve(), args.baseline_report_sha256.lower(),
        args.refined_report_sha256.lower())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
