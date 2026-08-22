#!/usr/bin/env python3
"""Evaluate the paired 2:1 timestep refinement of CCP surface transport."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "b617c59d7c9e2837bcec1c82ec86028f7fdf0772327f438474e3549c4678318e")
ELEMENTARY_CHARGE_C = 1.602176634e-19
RF_FREQUENCY_HZ = 13.56e6
PHASES = 200
NODES = 400


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def relative_change(first: float, second: float) -> float:
    return abs(second - first) / max(abs(first), abs(second), 1e-300)


def branch(output: Path, expected_name: str,
           expected_report_sha256: str) -> tuple[dict[str, object], Path]:
    report_path = output.parent.parent / "branch-report.json"
    if sha256(report_path) != expected_report_sha256:
        raise ValueError(f"{expected_name} branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("branch") != expected_name or
            report.get("all_gates_passed") is not True):
        raise ValueError(f"invalid {expected_name} branch report")
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"{expected_name} output differs: {name}")
    return report, report_path


def transport(output: Path) -> dict[str, object]:
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    cycles = int(metadata["rf_cycles"])
    surfaces = [float(value) for value in
                metadata["phase_surface_flux_positions"]]
    if metadata["phase_bins"] != PHASES or surfaces != [0.005, 0.015]:
        raise ValueError("surface transport contract differs")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    fields = rows(output / "spatial_phase_fields.csv")
    if len(moments) != PHASES * NODES or len(fields) != len(moments):
        raise ValueError("surface transport phase-space shape differs")
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
        (energy_density[((phase + 1) % PHASES) * NODES + node] -
         energy_density[((phase - 1) % PHASES) * NODES + node]) /
        (2.0 * phase_dt)
        for phase in range(PHASES) for node in range(NODES)]
    collision = [0.0] * (PHASES * NODES)
    for row in rows(output / "spatial_phase_collision_power.csv"):
        if row["channel"].startswith("electron_mcc."):
            collision[int(row["phase_bin"]) * NODES + int(row["node"])] += (
                float(row["mean_power_density_W_m-3"]))
    inferred = [power + loss - change for power, loss, change in
                zip(electric_power, collision, storage)]
    positions = [float(row["x_m"]) for row in moments[:NODES]]
    dx = positions[1] - positions[0]
    interior = [node for node, position in enumerate(positions)
                if surfaces[0] <= position <= surfaces[1]]

    summary = {}
    maximum_overflow = 0.0
    for row in rows(output / "phase_surface_flux_summary.csv"):
        key = (int(row["phase_bin"]), int(row["surface_id"]), row["direction"])
        summary[key] = float(row["kinetic_energy_flux_W_m-2"])
        maximum_overflow = max(maximum_overflow, float(row["overflow_fraction"]))
    tail: dict[tuple[int, int, str], float] = {}
    phase_duration = cycles / RF_FREQUENCY_HZ / PHASES
    for row in rows(output / "phase_surface_flux.csv"):
        energy = float(row["energy_eV"])
        if energy < 15.8:
            continue
        key = (int(row["phase_bin"]), int(row["surface_id"]), row["direction"])
        tail[key] = tail.get(key, 0.0) + (
            float(row["represented_crossings"]) * energy *
            ELEMENTARY_CHARGE_C / phase_duration)

    def window(first: int, end: int) -> dict[str, float]:
        inferred_value = math.fsum(
            inferred[phase * NODES + node] for phase in range(first, end)
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
        "critical_phase_0p125_to_0p5": window(25, 100),
        "exceptional_phase_0p375_to_0p5": window(75, 100),
        "maximum_surface_overflow_fraction": maximum_overflow,
    }


def evaluate(baseline: dict[str, object], refined: dict[str, object],
             limits: dict[str, object]) -> tuple[dict[str, float], dict[str, bool]]:
    base_critical = baseline["critical_phase_0p125_to_0p5"]
    fine_critical = refined["critical_phase_0p125_to_0p5"]
    base_exceptional = baseline["exceptional_phase_0p375_to_0p5"]
    fine_exceptional = refined["exceptional_phase_0p375_to_0p5"]
    metrics = {
        "critical_phase_direct_flux_relative_change": relative_change(
            base_critical["direct_outward_energy_flux_divergence_W_m-2"],
            fine_critical["direct_outward_energy_flux_divergence_W_m-2"]),
        "exceptional_octant_direct_flux_relative_change": relative_change(
            base_exceptional["direct_outward_energy_flux_divergence_W_m-2"],
            fine_exceptional["direct_outward_energy_flux_divergence_W_m-2"]),
        "exceptional_octant_tail_flux_relative_change": relative_change(
            base_exceptional[
                "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2"],
            fine_exceptional[
                "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2"]),
        "baseline_critical_phase_closure_error":
            base_critical["relative_closure_error"],
        "refined_critical_phase_closure_error":
            fine_critical["relative_closure_error"],
        "baseline_exceptional_octant_closure_error":
            base_exceptional["relative_closure_error"],
        "refined_exceptional_octant_closure_error":
            fine_exceptional["relative_closure_error"],
    }
    gates = {
        "critical_phase_direct_flux": metrics[
            "critical_phase_direct_flux_relative_change"] <= limits[
                "maximum_critical_phase_direct_flux_relative_change"],
        "exceptional_octant_direct_flux": metrics[
            "exceptional_octant_direct_flux_relative_change"] <= limits[
                "maximum_exceptional_octant_direct_flux_relative_change"],
        "exceptional_octant_tail_flux": metrics[
            "exceptional_octant_tail_flux_relative_change"] <= limits[
                "maximum_exceptional_octant_tail_flux_relative_change"],
        "critical_phase_closure": max(
            metrics["baseline_critical_phase_closure_error"],
            metrics["refined_critical_phase_closure_error"]) <= limits[
                "maximum_critical_phase_closure_error_each_branch"],
        "exceptional_octant_closure": max(
            metrics["baseline_exceptional_octant_closure_error"],
            metrics["refined_exceptional_octant_closure_error"]) <= limits[
                "maximum_exceptional_octant_closure_error_each_branch"],
    }
    return metrics, gates


def analyze(baseline_output: Path, refined_output: Path, rule_path: Path,
            baseline_report_sha256: str,
            refined_report_sha256: str) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux timestep rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_report, baseline_path = branch(
        baseline_output, "baseline_dt", baseline_report_sha256)
    refined_report, refined_path = branch(
        refined_output, "half_dt", refined_report_sha256)
    if (float(baseline_report["numerics"]["timestep_s"]) /
            float(refined_report["numerics"]["timestep_s"]) != 2.0):
        raise ValueError("branch timestep ratio differs")
    baseline = transport(baseline_output)
    refined = transport(refined_output)
    limits = dict(rule["prospective_acceptance"])
    limits.pop("all_gates_required")
    metrics, gates = evaluate(baseline, refined, limits)
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_surface_transport_timestep_refinement_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "baseline_branch_report_sha256": sha256(baseline_path),
            "refined_branch_report_sha256": sha256(refined_path),
        },
        "timestep_refinement_ratio": 2.0,
        "baseline": baseline,
        "half_dt": refined,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "claim_boundary": "A passing paired 2:1 refinement rules out material ordinary timestep sensitivity at the declared tolerances. It does not establish formal temporal order, asymptotic convergence, cross-code flux agreement, or experimental validation."
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
