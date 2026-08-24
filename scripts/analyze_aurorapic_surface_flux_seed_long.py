#!/usr/bin/env python3
"""Evaluate the four-cycle CCP surface-transport continuation-seed ensemble."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from analyze_aurorapic_surface_flux_mesh import transport
from analyze_aurorapic_surface_flux_seed import relative_range
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "8c52619b22fc4a9f8d8a6d412691660ddd4d1c4f9d9a24fbb9d4a64304cab626")
REPORT_SHA256S = {
    13507: "20522a5613c80eadae6e6801ce95714eedb2c8d06529f633364004e098c5d58c",
    24601: "cee5e10586d043a3230b66e40cd6b5fd49aacdd6ac1db4dfce7989877a90f571",
    35713: "6b190903e218c902eda156edcbc421b1d2702e03adc0ce57d44094843c5a2557",
}


def load_member(seed: int, root: Path, rule: dict[str, object]
                ) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != REPORT_SHA256S[seed]:
        raise ValueError(f"seed-{seed} long-window branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    contract = rule["branches"][f"seed_{seed}"]
    fixed = rule["fixed_inputs"]
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("branch") != f"seed_{seed}" or
            report.get("all_gates_passed") is not True or
            report["inputs"]["input_checkpoint_sha256"] !=
            contract["input_checkpoint_sha256"] or
            report["inputs"]["prior_branch_report_sha256"] !=
            contract["prior_branch_report_sha256"] or
            report["numerics"]["start_step"] != fixed["start_step"] or
            report["numerics"]["end_step"] != fixed["end_step"]):
        raise ValueError(f"seed-{seed} long-window branch contract differs")
    output = root / "output"
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"seed-{seed} long-window output differs: {filename}")
    result = transport(output)
    if result["nodes"] != int(fixed["nodes"]):
        raise ValueError(f"seed-{seed} long-window grid differs")
    return {
        "seed": seed,
        "branch_report_sha256": REPORT_SHA256S[seed],
        "peak_resident_set_kib": report["resources"]["peak_resident_set_kib"],
        "wall_seconds": report["resources"]["wall_seconds"],
        **result,
    }


def analyze(rule_path: Path, roots: dict[int, Path]) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux long-window seed rule differs")
    if set(roots) != set(REPORT_SHA256S):
        raise ValueError("surface-flux long-window seed members differ")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    members = [load_member(seed, roots[seed], rule) for seed in sorted(roots)]
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
        "scope": "prospective_long_window_surface_transport_seed_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "branch_report_sha256s": {
                str(seed): REPORT_SHA256S[seed] for seed in sorted(REPORT_SHA256S)},
        },
        "measurement_cycles_each_member": rule["fixed_inputs"][
            "measurement_cycles"],
        "members": members,
        "metrics": metrics,
        "thresholds": thresholds,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "assessment": {
            "critical_window_resolved_at_four_cycles": gates[
                "critical_phase_direct_flux"],
            "exceptional_transport_seed_robust": gates[
                "exceptional_octant_direct_flux"],
            "ionizing_tail_seed_robust": gates["exceptional_octant_tail_flux"],
            "conservation_closure_accepted": gates["critical_phase_closure"] and
                gates["exceptional_octant_closure"],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("seed_13507_root", type=Path)
    parser.add_argument("seed_24601_root", type=Path)
    parser.add_argument("seed_35713_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule.resolve(), {
        13507: args.seed_13507_root.resolve(),
        24601: args.seed_24601_root.resolve(),
        35713: args.seed_35713_root.resolve(),
    })
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
