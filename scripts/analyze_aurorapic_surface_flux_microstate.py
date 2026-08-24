#!/usr/bin/env python3
"""Evaluate the constrained-microstate CCP surface-transport ensemble."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from analyze_aurorapic_surface_flux_mesh import transport
from analyze_aurorapic_surface_flux_seed import relative_range
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "95b72b70c95d95980d8270156bdeeb796e2de8cedde389057173d8705787ce63")
REPORT_SHA256S = {
    "locked_source_microstate":
        "6f99b7a0c76eeaee2bf838fd388708bceda7d901be3a17ab96e7aa439590a637",
    "microstate_51949":
        "20b528d9c75518d91a34a3806df3eea5f5465b179c68c5782ef78d22100efa31",
    "microstate_63059":
        "1ee0e01f7301335433c8dc4ad0e7ac3dd05b340b10485ba5a98e9a93364a0429",
}


def load_member(name: str, root: Path, rule: dict[str, object]
                ) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != REPORT_SHA256S[name]:
        raise ValueError(f"{name} surface branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    contract = rule["branches"][name]
    fixed = rule["fixed_inputs"]
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("branch") != name or
            report.get("all_gates_passed") is not True or
            report["inputs"]["input_checkpoint_sha256"] !=
            contract["input_checkpoint_sha256"] or
            report["inputs"]["prior_branch_report_sha256"] !=
            contract["prior_branch_report_sha256"] or
            report["numerics"]["start_step"] != fixed["start_step"] or
            report["numerics"]["end_step"] != fixed["end_step"]):
        raise ValueError(f"{name} surface branch contract differs")
    output = root / "output"
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"{name} surface output differs: {filename}")
    result = transport(output)
    if result["nodes"] != int(fixed["nodes"]):
        raise ValueError(f"{name} surface grid differs")
    return {
        "name": name,
        "particle_state_sha256": contract["particle_state_sha256"],
        "branch_report_sha256": REPORT_SHA256S[name],
        "peak_resident_set_kib": report["resources"]["peak_resident_set_kib"],
        "wall_seconds": report["resources"]["wall_seconds"],
        **result,
    }


def analyze(rule_path: Path, roots: dict[str, Path]) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux microstate rule differs")
    if set(roots) != set(REPORT_SHA256S):
        raise ValueError("surface-flux microstate members differ")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    order = rule["microstate_contract"]["members"]
    members = [load_member(name, roots[name], rule) for name in order]
    critical = [member["critical_phase_0p125_to_0p5"] for member in members]
    exceptional = [member["exceptional_phase_0p375_to_0p5"] for member in members]
    direct = "direct_outward_energy_flux_divergence_W_m-2"
    tail = "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2"
    metrics = {
        "critical_phase_direct_flux_relative_range": relative_range(
            [item[direct] for item in critical]),
        "exceptional_octant_direct_flux_relative_range": relative_range(
            [item[direct] for item in exceptional]),
        "exceptional_octant_tail_flux_relative_range": relative_range(
            [item[tail] for item in exceptional]),
        "maximum_critical_phase_closure_error": max(
            item["relative_closure_error"] for item in critical),
        "maximum_exceptional_octant_closure_error": max(
            item["relative_closure_error"] for item in exceptional),
        "mean_critical_phase_direct_flux_W_m-2": statistics.fmean(
            item[direct] for item in critical),
        "mean_exceptional_octant_direct_flux_W_m-2": statistics.fmean(
            item[direct] for item in exceptional),
        "mean_exceptional_octant_tail_flux_W_m-2": statistics.fmean(
            item[tail] for item in exceptional),
    }
    declared = rule["prospective_acceptance"]
    thresholds = {key: value for key, value in declared.items()
                  if key not in ("relative_range_definition", "all_gates_required")}
    gates = {
        "critical_phase_direct_flux": metrics[
            "critical_phase_direct_flux_relative_range"] <= thresholds[
                "maximum_critical_phase_direct_flux_relative_range"],
        "exceptional_octant_direct_flux": metrics[
            "exceptional_octant_direct_flux_relative_range"] <= thresholds[
                "maximum_exceptional_octant_direct_flux_relative_range"],
        "exceptional_octant_tail_flux": metrics[
            "exceptional_octant_tail_flux_relative_range"] <= thresholds[
                "maximum_exceptional_octant_tail_flux_relative_range"],
        "critical_phase_closure": metrics[
            "maximum_critical_phase_closure_error"] <= thresholds[
                "maximum_critical_phase_closure_error_each_member"],
        "exceptional_octant_closure": metrics[
            "maximum_exceptional_octant_closure_error"] <= thresholds[
                "maximum_exceptional_octant_closure_error_each_member"],
    }
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_surface_transport_microstate_ensemble_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "branch_report_sha256s": REPORT_SHA256S,
        },
        "measurement_cycles_each_member": rule["fixed_inputs"][
            "measurement_cycles"],
        "microstate_contract": rule["microstate_contract"],
        "members": members,
        "metrics": metrics,
        "thresholds": thresholds,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "assessment": {
            "critical_transport_microstate_robust": gates[
                "critical_phase_direct_flux"],
            "exceptional_transport_microstate_robust": gates[
                "exceptional_octant_direct_flux"],
            "ionizing_tail_microstate_robust": gates[
                "exceptional_octant_tail_flux"],
            "conservation_closure_accepted": gates["critical_phase_closure"] and
                gates["exceptional_octant_closure"],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("locked_source_root", type=Path)
    parser.add_argument("microstate_51949_root", type=Path)
    parser.add_argument("microstate_63059_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule.resolve(), {
        "locked_source_microstate": args.locked_source_root.resolve(),
        "microstate_51949": args.microstate_51949_root.resolve(),
        "microstate_63059": args.microstate_63059_root.resolve(),
    })
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
