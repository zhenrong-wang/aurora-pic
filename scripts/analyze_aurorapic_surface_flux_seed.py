#!/usr/bin/env python3
"""Evaluate the three-member CCP surface-transport continuation-seed ensemble."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from analyze_aurorapic_surface_flux_mesh import transport
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "f0d94adedf0bb6f74cb23472b6515a58d2c7c39ddbf3a7e363e6a3ce83ac904d")
REPORT_SHA256S = {
    13507: "f6c6a6a97258f4e8827df29ff6fd5fbfb9521b84d969cfe8df943b125eccc172",
    24601: "d6f510b7ac59725253523fce0577c6ceedbc913d1a62165f516f730bc5386e04",
    35713: "ecf61ec2f3c4eef8509fb462d4aef70a067f33f0f963846c6c5980a4f795e3ba",
}


def relative_range(values: list[float]) -> float:
    mean = statistics.fmean(values)
    if mean <= 0.0:
        raise ValueError("ensemble range requires a positive mean")
    return (max(values) - min(values)) / mean


def load_member(seed: int, root: Path, rule: dict[str, object]
                ) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != REPORT_SHA256S[seed]:
        raise ValueError(f"seed-{seed} branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    contract = rule["branches"][f"seed_{seed}"]
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("branch") != f"seed_{seed}" or
            report.get("all_gates_passed") is not True or
            report["inputs"]["input_checkpoint_sha256"] !=
            contract["input_checkpoint_sha256"] or
            report["inputs"]["prior_branch_report_sha256"] !=
            contract["prior_branch_report_sha256"]):
        raise ValueError(f"seed-{seed} branch contract differs")
    output = root / "output"
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"seed-{seed} output differs: {filename}")
    result = transport(output)
    if result["nodes"] != int(rule["fixed_inputs"]["nodes"]):
        raise ValueError(f"seed-{seed} grid differs")
    return {
        "seed": seed,
        "branch_report_sha256": REPORT_SHA256S[seed],
        "peak_resident_set_kib": report["resources"]["peak_resident_set_kib"],
        "wall_seconds": report["resources"]["wall_seconds"],
        **result,
    }


def analyze(rule_path: Path, roots: dict[int, Path]) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("surface-flux seed rule differs")
    if set(roots) != set(REPORT_SHA256S):
        raise ValueError("surface-flux seed members differ")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    members = [load_member(seed, roots[seed], rule) for seed in sorted(roots)]
    critical = [member["critical_phase_0p125_to_0p5"] for member in members]
    exceptional = [member["exceptional_phase_0p375_to_0p5"] for member in members]
    key = "direct_outward_energy_flux_divergence_W_m-2"
    tail_key = "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2"
    metrics = {
        "critical_phase_direct_flux_relative_range": relative_range(
            [item[key] for item in critical]),
        "exceptional_octant_direct_flux_relative_range": relative_range(
            [item[key] for item in exceptional]),
        "exceptional_octant_tail_flux_relative_range": relative_range(
            [item[tail_key] for item in exceptional]),
        "maximum_critical_phase_closure_error": max(
            item["relative_closure_error"] for item in critical),
        "maximum_exceptional_octant_closure_error": max(
            item["relative_closure_error"] for item in exceptional),
        "mean_exceptional_octant_direct_flux_W_m-2": statistics.fmean(
            item[key] for item in exceptional),
        "mean_exceptional_octant_tail_flux_W_m-2": statistics.fmean(
            item[tail_key] for item in exceptional),
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
        "scope": "prospective_surface_transport_seed_ensemble_result",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "branch_report_sha256s": {
                str(seed): REPORT_SHA256S[seed] for seed in sorted(REPORT_SHA256S)},
        },
        "members": members,
        "metrics": metrics,
        "thresholds": thresholds,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "assessment": {
            "exceptional_transport_seed_robust": gates[
                "exceptional_octant_direct_flux"],
            "ionizing_tail_seed_robust": gates["exceptional_octant_tail_flux"],
            "critical_window_requires_more_sampling": not gates[
                "critical_phase_direct_flux"],
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
