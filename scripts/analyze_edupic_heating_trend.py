#!/usr/bin/env python3
"""Analyze the prospectively declared cycle-80 to cycle-116 heating trend."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = "f8d4a05fc5fb359cc5808fb35d08d39af3392b8d5108a2474716ee57a2c845be"
BASELINE_SHA256 = "f2d698b183c8bf48e10b4c0a62d984d009f2c6be70d2725ab7c8cdf3e21dd728"
FOLLOWUP_SHA256 = "58163e9027db925ff8b9adbd82b3382db9d465a7924eb765e3a2a6fb5cdb2151"


def relative_change(baseline: float, followup: float) -> float:
    if baseline == 0.0:
        raise ValueError("trend baseline must be nonzero")
    return followup / baseline - 1.0


def classify(changes: dict[str, float]) -> str:
    expected = {
        "electron_density_reference_ratio": 1,
        "electron_rf_power_per_particle_reference_ratio": -1,
        "electron_mean_energy_relative_l2": -1,
        "effective_ionization_frequency_reference_ratio": -1,
    }
    matches = {
        name: changes[name] * direction > 0.0
        for name, direction in expected.items()
    }
    if all(matches.values()):
        return "all_predeclared_directions_support_transient_heating_hypothesis"
    if not any(matches.values()):
        return "no_predeclared_direction_supports_transient_heating_hypothesis"
    return "mixed_predeclared_directional_outcome"


def value(report: dict[str, object], name: str) -> float:
    derived = report["derived_diagnostics"]
    comparisons = report["comparisons"]
    if name == "electron_density_reference_ratio":
        return float(derived["electron_rf_power_per_particle"]
                     ["candidate_to_reference_average_number_density_ratio"])
    if name == "electron_rf_power_per_particle_reference_ratio":
        return float(derived["electron_rf_power_per_particle"]
                     ["candidate_to_reference_power_per_particle_ratio"])
    if name == "electron_mean_energy_relative_l2":
        return float(comparisons["electron_mean_energy"]["relative_l2"])
    if name == "effective_ionization_frequency_reference_ratio":
        return float(derived["ionization_per_electron"]
                     ["candidate_to_reference_effective_ionization_frequency_ratio"])
    if name == "ionization_rate_relative_l2":
        return float(comparisons["ionization_rate"]["relative_l2"])
    if name == "electron_density_relative_l2":
        return float(comparisons["electron_density"]["relative_l2"])
    if name == "eedf_collision_ledger_ratio":
        return float(derived["ionization_eedf_ledger_closure"]
                     ["measured_to_predicted_average_frequency_ratio"])
    raise KeyError(name)


def analyze(baseline_path: Path, followup_path: Path,
            rule_path: Path) -> dict[str, object]:
    for path, expected, label in (
            (baseline_path, BASELINE_SHA256, "baseline"),
            (followup_path, FOLLOWUP_SHA256, "followup"),
            (rule_path, RULE_SHA256, "rule")):
        if sha256(path) != expected:
            raise ValueError(f"locked {label} evidence differs")
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    followup = json.loads(followup_path.read_text(encoding="utf-8"))
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    if (baseline["candidate_measurement_context"]["window"] != {
            "start_cycle": 76, "end_cycle": 80, "measurement_cycles": 4,
            "equilibration_statistics_excluded": True} or
            followup["candidate_measurement_context"]["window"] != {
                "start_cycle": 112, "end_cycle": 116,
                "measurement_cycles": 4,
                "equilibration_statistics_excluded": True}):
        raise ValueError("comparison windows differ from the trend contract")
    followup_inputs = followup["candidate_measurement_context"]["inputs"]
    locked = rule["locked_inputs"]
    if (followup_inputs["checkpoint_sha256"] !=
            locked["cycle112_checkpoint_sha256"] or
            followup_inputs["rule_sha256"] != RULE_SHA256):
        raise ValueError("followup inputs differ from the prospective rule")

    names = (
        "electron_density_reference_ratio",
        "electron_rf_power_per_particle_reference_ratio",
        "electron_mean_energy_relative_l2",
        "effective_ionization_frequency_reference_ratio",
        "ionization_rate_relative_l2",
        "electron_density_relative_l2",
        "eedf_collision_ledger_ratio",
    )
    observables = {}
    changes = {}
    for name in names:
        before, after = value(baseline, name), value(followup, name)
        change = relative_change(before, after)
        observables[name] = {
            "cycle80": before, "cycle116": after,
            "relative_change": change,
        }
        changes[name] = change
    classification = classify(changes)
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "prospective_cycle80_to116_heating_trend_diagnostic",
        "evidence_sha256": {
            "prospective_rule": RULE_SHA256,
            "cycle80_phase_space": BASELINE_SHA256,
            "cycle116_phase_space": FOLLOWUP_SHA256,
        },
        "observables": observables,
        "classification": classification,
        "directional_hypothesis_result": {
            "transient_heating_supported": classification.startswith("all_"),
            "thresholded_acceptance_performed": False,
            "stationarity_established": False,
            "validation_established": False,
        },
        "interpretation": (
            "Increasing density coincides with reduced RF power per electron, "
            "electron-energy error, and ionization frequency. This supports a "
            "population-filling contribution to the hot EEDF but does not exclude "
            "a residual heating-model discrepancy."),
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": "none_directional_transient_discriminator_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("followup", type=Path)
    parser.add_argument("rule", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.baseline.resolve(), args.followup.resolve(),
                     args.rule.resolve())
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
