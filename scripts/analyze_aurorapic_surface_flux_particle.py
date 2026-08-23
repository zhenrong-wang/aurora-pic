#!/usr/bin/env python3
"""Evaluate paired 1x/2x-particle CCP surface-transport refinement."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analyze_aurorapic_surface_flux_mesh import transport
from analyze_aurorapic_surface_flux_timestep import evaluate
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "7fd0cf0eeb432b12a9a63d446b12c6e9e24ca0b03bd40f9a8cbc28e08841d2d7")
MESH_RULE_SHA256 = (
    "4bea77b968db89ca6a2e066a599d3e85b99c480de2f0cb6e56b12bdaeb891f54")


def branch(output: Path, expected_name: str, expected_report_sha256: str,
           expected_rule_sha256: str) -> tuple[dict[str, object], Path]:
    report_path = output.parent / "branch-report.json"
    if sha256(report_path) != expected_report_sha256:
        raise ValueError(f"{expected_name} branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (report.get("rule_sha256") != expected_rule_sha256 or
            report.get("branch") != expected_name or
            report.get("all_gates_passed") is not True):
        raise ValueError(f"invalid {expected_name} branch report")
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"{expected_name} output differs: {filename}")
    return report, report_path


def analyze(baseline_output: Path, doubled_output: Path, rule_path: Path,
            baseline_hash: str, doubled_hash: str) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux particle rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_contract = rule["paired_baseline"]
    doubled_contract = rule["branches"]["double_particles"]
    if baseline_hash != baseline_contract["surface_branch_report_sha256"]:
        raise ValueError("baseline report is not the predeclared branch")
    baseline_report, baseline_path = branch(
        baseline_output, "refined_grid", baseline_hash,
        MESH_RULE_SHA256)
    doubled_report, doubled_path = branch(
        doubled_output, "double_particles", doubled_hash, RULE_SHA256)
    for key in ("input_checkpoint_sha256", "prior_branch_report_sha256"):
        if doubled_report["inputs"][key] != doubled_contract[key]:
            raise ValueError(f"double_particles {key} differs")
    for key in ("start_step", "end_step"):
        if doubled_report["numerics"][key] != rule["fixed_inputs"][key]:
            raise ValueError(f"double_particles {key} differs")
    baseline = transport(baseline_output)
    doubled = transport(doubled_output)
    expected_nodes = int(rule["fixed_inputs"]["nodes"])
    if (baseline["nodes"] != expected_nodes or
            doubled["nodes"] != expected_nodes or
            int(baseline_report["numerics"]["nodes"]) != expected_nodes or
            int(doubled_report["numerics"]["nodes"]) != expected_nodes):
        raise ValueError("particle-refinement grid differs")
    ratio = (float(baseline_contract["macro_weight"]) /
             float(doubled_contract["macro_weight"]))
    if ratio != float(doubled_contract["particle_refinement_ratio"]):
        raise ValueError("particle-refinement ratio differs")
    limits = dict(rule["prospective_acceptance"])
    limits.pop("all_gates_required")
    metrics, gates = evaluate(baseline, doubled, limits)
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_surface_transport_particle_refinement_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "baseline_branch_report_sha256": sha256(baseline_path),
            "double_particles_branch_report_sha256": sha256(doubled_path),
        },
        "particle_refinement": {
            "ratio": ratio,
            "nodes_each_branch": expected_nodes,
            "baseline_macro_weight": baseline_contract["macro_weight"],
            "double_particles_macro_weight": doubled_contract["macro_weight"],
        },
        "baseline": baseline,
        "double_particles": doubled,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "assessment": {
            "ordinary_particle_count_sensitivity_constrained": passed,
            "exceptional_transport_particle_stable": gates[
                "exceptional_octant_direct_flux"],
            "ionizing_tail_particle_stable": gates[
                "exceptional_octant_tail_flux"],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("doubled_output", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--baseline-report-sha256", required=True)
    parser.add_argument("--doubled-report-sha256", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.baseline_output.resolve(), args.doubled_output.resolve(),
        args.rule.resolve(), args.baseline_report_sha256.lower(),
        args.doubled_report_sha256.lower())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
