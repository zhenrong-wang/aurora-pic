#!/usr/bin/env python3
"""Analyze the preregistered collision-enabled common-state ensemble."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics


COLLISION_THRESHOLDS = {
    "electron_elastic": 0.05,
    "electron_excitation": 0.05,
    "electron_ionization": 0.10,
    "ion_isotropic": 0.10,
    "ion_backward": 0.10,
}
WALL_KEYS = (
    "electron_absorbed_left", "electron_absorbed_right",
    "ion_absorbed_left", "ion_absorbed_right")
POPULATION_KEYS = ("electrons", "ions")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def symmetric_relative(a: float, b: float) -> float:
    return abs(a - b) / max(0.5 * (abs(a) + abs(b)), 1.0)


def summary(values: list[float]) -> dict[str, float]:
    return {
        "mean": statistics.fmean(values),
        "sample_standard_deviation": statistics.stdev(values),
        "minimum": min(values), "maximum": max(values),
    }


def read_field(path: Path, implementation: str) -> tuple[list[float], list[float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    x_name = "x_m" if implementation == "edupic" else "x"
    e_name = "electric_field_V_m" if implementation == "edupic" else "E"
    return ([float(row[x_name]) for row in rows],
            [float(row[e_name]) for row in rows])


def analyze(rule_path: Path, lock_path: Path, report_path: Path,
            output_root: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    integrity: dict[str, bool] = {
        "rule_hash_matches_report": report["rule_sha256"] == sha256(rule_path),
        "rule_hash_matches_lock": lock["rule_sha256"] == sha256(rule_path),
        "lock_hash_matches_report": (
            report["execution_lock_sha256"] == sha256(lock_path)),
        "edupic_binary_hash_matches": (
            report["edupic_binary_sha256"] == lock["edupic_binary_sha256"]),
        "aurorapic_binary_hash_matches": (
            report["aurorapic_binary_sha256"] ==
            lock["aurorapic_binary_sha256"]),
        "all_ten_members_complete": bool(report["all_ten_members_complete"]),
        "all_resource_gates_passed": bool(report["all_resource_gates_passed"]),
    }
    by_code = {"edupic": [], "aurorapic": []}
    expected_seeds = rule["ensemble_contract"]["seeds_each_implementation"]
    for member in report["members"]:
        by_code[member["implementation"]].append(member)
    for implementation, members in by_code.items():
        integrity[f"{implementation}_five_unique_locked_seeds"] = (
            sorted(item["seed"] for item in members) == sorted(expected_seeds))
        integrity[f"{implementation}_endpoint_alignment"] = all(
            item["endpoint"]["pre_push_step"] == 4000 for item in members)

    observable_results: dict[str, object] = {}
    for key, threshold in COLLISION_THRESHOLDS.items():
        left = [float(item["endpoint"][key]) for item in by_code["edupic"]]
        right = [float(item["endpoint"][key]) for item in by_code["aurorapic"]]
        left_summary, right_summary = summary(left), summary(right)
        relative = symmetric_relative(left_summary["mean"], right_summary["mean"])
        nonzero = left_summary["mean"] > 0.0 and right_summary["mean"] > 0.0
        observable_results[key] = {
            "edupic": left_summary, "aurorapic": right_summary,
            "symmetric_relative_mean_difference": relative,
            "threshold": threshold, "nonzero_both": nonzero,
            "passed": nonzero and relative <= threshold,
        }

    wall_results: dict[str, object] = {}
    for key in WALL_KEYS:
        left_summary = summary([
            float(item["endpoint"][key]) for item in by_code["edupic"]])
        right_summary = summary([
            float(item["endpoint"][key]) for item in by_code["aurorapic"]])
        relative = symmetric_relative(left_summary["mean"], right_summary["mean"])
        low_count = left_summary["mean"] < 20.0 and right_summary["mean"] < 20.0
        status = "low_count_inconclusive" if low_count else (
            "passed" if relative <= 0.10 else "failed")
        wall_results[key] = {
            "edupic": left_summary, "aurorapic": right_summary,
            "symmetric_relative_mean_difference": relative,
            "threshold": 0.10, "status": status,
        }

    population_results: dict[str, object] = {}
    population_passes = []
    for key in POPULATION_KEYS:
        initial = float(rule["locked_inputs"]["initial_populations"][key])
        left_summary = summary([
            float(item["endpoint"][key]) for item in by_code["edupic"]])
        right_summary = summary([
            float(item["endpoint"][key]) for item in by_code["aurorapic"]])
        difference_fraction = abs(
            left_summary["mean"] - right_summary["mean"]) / initial
        left_change = left_summary["mean"] - initial
        right_change = right_summary["mean"] - initial
        same_sign = ((left_change == 0.0 and right_change == 0.0) or
                     (left_change * right_change > 0.0))
        passed = difference_fraction <= 0.005 and same_sign
        population_passes.append(passed)
        population_results[key] = {
            "initial": initial, "edupic": left_summary,
            "aurorapic": right_summary,
            "mean_difference_fraction_of_initial": difference_fraction,
            "mean_changes": {"edupic": left_change, "aurorapic": right_change},
            "change_signs_agree": same_sign, "threshold": 0.005,
            "passed": passed,
        }

    field_sets: dict[str, list[list[float]]] = {"edupic": [], "aurorapic": []}
    coordinates: dict[str, list[float]] = {}
    field_hashes_pass = True
    finite_fields = True
    for implementation, members in by_code.items():
        for member in members:
            if implementation == "edupic":
                path = (output_root / implementation / f"seed-{member['seed']}" /
                        "edupic_collision_endpoint_field.csv")
            else:
                path = (output_root / implementation / f"seed-{member['seed']}" /
                        "output" / "fields_3999.csv")
            field_hashes_pass &= sha256(path) == member["field_sha256"]
            x, electric = read_field(path, implementation)
            if implementation not in coordinates:
                coordinates[implementation] = x
            integrity[f"{implementation}_field_shape_{member['seed']}"] = (
                len(x) == 400 and x == coordinates[implementation])
            finite_fields &= all(math.isfinite(value) for value in x + electric)
            field_sets[implementation].append(electric)
    integrity["field_hashes_match"] = field_hashes_pass
    integrity["field_values_finite"] = finite_fields
    integrity["crosscode_coordinates_match"] = (
        coordinates.get("edupic") == coordinates.get("aurorapic"))
    mean_fields = {
        implementation: [statistics.fmean(values) for values in zip(*fields)]
        for implementation, fields in field_sets.items()
    }
    numerator = sum((a - b) ** 2 for a, b in zip(
        mean_fields["edupic"], mean_fields["aurorapic"]))
    denominator = sum(value ** 2 for value in mean_fields["edupic"])
    relative_rms = math.sqrt(numerator / denominator)
    member_energy = {
        implementation: [sum(value * value for value in field) for field in fields]
        for implementation, fields in field_sets.items()
    }
    energy_summary = {
        implementation: summary(values)
        for implementation, values in member_energy.items()
    }
    energy_ratio = (energy_summary["aurorapic"]["mean"] /
                    energy_summary["edupic"]["mean"])
    field_passed = relative_rms <= 0.02 and 0.95 <= energy_ratio <= 1.05
    field_result = {
        "ensemble_mean_profile_relative_rms": relative_rms,
        "relative_rms_reference": "edupic_ensemble_mean_profile",
        "relative_rms_threshold": 0.02,
        "field_energy_proxy": energy_summary,
        "aurorapic_to_edupic_mean_field_energy_proxy_ratio": energy_ratio,
        "energy_ratio_interval": [0.95, 1.05],
        "passed": field_passed,
    }

    integrity_passed = all(integrity.values())
    collision_passed = all(item["passed"] for item in observable_results.values())
    wall_has_failure = any(item["status"] == "failed" for item in wall_results.values())
    population_passed = all(population_passes)
    pilot_supported = (integrity_passed and collision_passed and
                       not wall_has_failure and population_passed and field_passed)
    if not integrity_passed:
        outcome = "inconclusive_integrity_failure"
    elif pilot_supported:
        outcome = "one_period_collision_enabled_stochastic_consistency_supported"
    else:
        outcome = "localized_collision_enabled_common_state_discrepancy"
    return {
        "schema_version": 1,
        "scope": "collision_enabled_common_state_ensemble_result",
        "rule_sha256": sha256(rule_path), "execution_lock_sha256": sha256(lock_path),
        "execution_report_sha256": sha256(report_path),
        "integrity": integrity, "integrity_passed": integrity_passed,
        "collision_channels": observable_results,
        "collision_traffic_compatible": collision_passed,
        "wall_absorption": wall_results, "wall_traffic_has_failure": wall_has_failure,
        "populations": population_results, "population_compatible": population_passed,
        "field": field_result, "field_compatible": field_passed,
        "pilot_support_gate_passed": pilot_supported,
        "formal_outcome": outcome,
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    value = analyze(args.rule.resolve(), args.lock.resolve(), args.report.resolve(),
                    args.output_root.resolve())
    args.result.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(value, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
