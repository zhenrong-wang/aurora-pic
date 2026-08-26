#!/usr/bin/env python3
"""Analyze the corrected-cadence common-state phase-EEDF ensemble."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

from analyze_edupic_phase_eedf_crosscode import (
    ENERGY_BINS, load_histogram, relative_range, summarize)
from audit_edupic_ionization_path import read_cross_section
from run_aurorapic_corrected_phase_eedf import sha256


SCALAR_NAMES = (
    "histogram_mean_energy_eV", "fraction_11p5_to_15p8_eV",
    "fraction_15p8_to_30_eV", "fraction_at_or_above_30_eV",
    "fraction_at_or_above_15p8_eV",
    "eedf_folded_ionization_frequency_s-1",
)
SCOPES = {
    "critical_x020_to_x060_phase_0p125_to_0p5":
        (("x020_040", "x040_060"), 25, 100),
    "upstream_x010_to_x020_phase_0p125_to_0p5":
        (("x010_020",), 25, 100),
    "x020_to_x040_phase_0p125_to_0p5": (("x020_040",), 25, 100),
    "x040_to_x060_phase_0p125_to_0p5": (("x040_060",), 25, 100),
    "critical_phase_0p125_to_0p25":
        (("x020_040", "x040_060"), 25, 50),
    "critical_phase_0p25_to_0p375":
        (("x020_040", "x040_060"), 50, 75),
    "critical_phase_0p375_to_0p5":
        (("x020_040", "x040_060"), 75, 100),
}


def ensemble_mean(values: list[dict[str, object]]) -> dict[str, object]:
    return {
        **{name: statistics.fmean(float(item[name]) for item in values)
           for name in SCALAR_NAMES},
        "probability_by_energy_bin": [statistics.fmean(
            float(item["probability_by_energy_bin"][energy_bin])
            for item in values) for energy_bin in range(ENERGY_BINS)],
    }


def compare_ensembles(candidates: list[dict[str, object]],
                      natives: list[dict[str, object]]) -> dict[str, object]:
    candidate = ensemble_mean(candidates)
    native = ensemble_mean(natives)
    ratios = {name: float(candidate[name]) / float(native[name])
              if float(native[name]) else math.inf for name in SCALAR_NAMES}
    total_variation = 0.5 * math.fsum(abs(float(a) - float(b)) for a, b in zip(
        candidate["probability_by_energy_bin"], native["probability_by_energy_bin"]))
    return {
        "aurorapic_ensemble_mean": {
            name: candidate[name] for name in SCALAR_NAMES},
        "native_edupic_ensemble_mean": {
            name: native[name] for name in SCALAR_NAMES},
        "aurorapic_to_native_edupic_ratio": ratios,
        "total_variation_distance": total_variation,
    }


def moment_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def analyze(rule_path: Path, lock_path: Path, report_path: Path,
            candidate_root: Path, native_root: Path,
            ionization_table: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    integrity = {
        "rule_hash_matches_lock": lock["rule_sha256"] == sha256(rule_path),
        "rule_hash_matches_report": report["rule_sha256"] == sha256(rule_path),
        "lock_hash_matches_report": (
            report["execution_lock_sha256"] == sha256(lock_path)),
        "analyzer_hash_matches_lock": (
            lock["analyzer_sha256"] == sha256(Path(__file__).resolve())),
        "binary_hash_matches": (
            report["aurorapic_binary_sha256"] ==
            lock["aurorapic_binary_sha256"]),
        "particle_state_hash_matches": (
            report["particle_state_sha256"] ==
            rule["locked_inputs"]["particle_state_sha256"]),
        "all_three_members_complete": bool(report["all_three_members_complete"]),
        "all_resource_gates_passed": bool(report["all_resource_gates_passed"]),
    }
    if sha256(ionization_table) != rule["locked_inputs"]["ionization_table_sha256"]:
        integrity["ionization_table_hash_matches"] = False
    else:
        integrity["ionization_table_hash_matches"] = True
    energies, cross_sections = read_cross_section(ionization_table)

    candidate_histograms = []
    candidate_members = []
    expected_seeds = rule["diagnostic_contract"]["candidate_seeds"]
    integrity["candidate_seed_set_matches"] = (
        sorted(item["seed"] for item in report["members"]) == sorted(expected_seeds))
    for item in report["members"]:
        seed = int(item["seed"])
        root = candidate_root / f"seed-{seed}" / "output"
        histogram_path = root / "phase_eedf.csv"
        moments_path = root / "phase_eedf_moments.csv"
        integrity[f"candidate_{seed}_histogram_hash_matches"] = (
            sha256(histogram_path) == item["phase_eedf_sha256"])
        integrity[f"candidate_{seed}_moments_hash_matches"] = (
            sha256(moments_path) == item["phase_eedf_moments_sha256"])
        histogram = load_histogram(histogram_path, "represented_count")
        candidate_histograms.append(histogram)
        rows = moment_rows(moments_path)
        critical = set(rule["diagnostic_contract"]["critical_regions"])
        selected = [row for row in rows if row["region"] in critical]
        candidate_members.append({
            "seed": seed,
            "critical_minimum_macro_observations_per_region_phase_bin": min(
                int(row["macro_observations"]) for row in selected),
            "maximum_histogram_overflow_fraction": max(
                float(row["overflow_fraction"]) for row in rows),
        })

    native_histograms = []
    native_members = []
    for seed in rule["diagnostic_contract"]["native_reference_seeds"]:
        expected = rule["locked_inputs"]["native_reference_members"][str(seed)]
        root = native_root / f"seed-{seed}"
        histogram_path = root / "edupic_phase_eedf.csv"
        moments_path = root / "edupic_phase_eedf_moments.csv"
        integrity[f"native_{seed}_histogram_hash_matches"] = (
            sha256(histogram_path) == expected["phase_eedf_sha256"])
        integrity[f"native_{seed}_moments_hash_matches"] = (
            sha256(moments_path) == expected["moments_sha256"])
        native_histograms.append(load_histogram(
            histogram_path, "represented_count"))
        native_members.append({"seed": seed})

    comparisons = {}
    candidate_scope_values = {}
    for name, (regions, first, end) in SCOPES.items():
        candidates = [summarize(histogram, regions, first, end,
                                energies, cross_sections)
                      for histogram in candidate_histograms]
        natives = [summarize(histogram, regions, first, end,
                             energies, cross_sections)
                   for histogram in native_histograms]
        comparisons[name] = compare_ensembles(candidates, natives)
        candidate_scope_values[name] = candidates

    critical_values = candidate_scope_values[
        "critical_x020_to_x060_phase_0p125_to_0p5"]
    repeatability = {
        "candidate_critical_folded_ionization_relative_range": relative_range([
            float(item["eedf_folded_ionization_frequency_s-1"])
            for item in critical_values]),
        "candidate_critical_tail_relative_range": relative_range([
            float(item["fraction_at_or_above_15p8_eV"])
            for item in critical_values]),
    }
    measurement = {
        "candidate_observation_population": all(
            item["critical_minimum_macro_observations_per_region_phase_bin"] >=
            100000 for item in candidate_members),
        "candidate_histogram_overflow": all(
            item["maximum_histogram_overflow_fraction"] <= 1.0e-6
            for item in candidate_members),
        "candidate_folded_ionization_repeatability": repeatability[
            "candidate_critical_folded_ionization_relative_range"] <= 0.15,
        "candidate_tail_repeatability": repeatability[
            "candidate_critical_tail_relative_range"] <= 0.15,
    }
    measurement_passed = all(measurement.values())
    critical = comparisons["critical_x020_to_x060_phase_0p125_to_0p5"]
    ratios = critical["aurorapic_to_native_edupic_ratio"]
    ionization_ratio = ratios["eedf_folded_ionization_frequency_s-1"]
    tail_ratio = ratios["fraction_at_or_above_15p8_eV"]
    phase_ratios = [comparisons[name]["aurorapic_to_native_edupic_ratio"][
        "eedf_folded_ionization_frequency_s-1"] for name in (
            "critical_phase_0p125_to_0p25", "critical_phase_0p25_to_0p375",
            "critical_phase_0p375_to_0p5")]
    strong = (0.95 <= ionization_ratio <= 1.05 and
              0.95 <= tail_ratio <= 1.05 and
              all(0.90 <= value <= 1.10 for value in phase_ratios) and
              critical["total_variation_distance"] <= 0.06)
    old_ratio = float(rule["basis"]["previous_critical_folded_ionization_ratio"])
    error_fraction = (abs(math.log(ionization_ratio)) / abs(math.log(old_ratio)))
    major = (not strong and 0.94 <= ionization_ratio <= 1.06 and
             error_fraction <= 0.40 and 0.94 <= tail_ratio <= 1.06)
    integrity_passed = all(integrity.values())
    if not integrity_passed or not measurement_passed:
        outcome = "inconclusive_integrity_or_measurement_failure"
    elif strong:
        outcome = "strong_corrected_cadence_phase_eedf_closure"
    elif major:
        outcome = "major_corrected_cadence_phase_eedf_support"
    else:
        outcome = "corrected_cadence_phase_eedf_not_supported"
    return {
        "schema_version": 1,
        "scope": "corrected_cadence_common_state_phase_eedf_result",
        "rule_sha256": sha256(rule_path),
        "execution_lock_sha256": sha256(lock_path),
        "execution_report_sha256": sha256(report_path),
        "integrity": integrity, "integrity_passed": integrity_passed,
        "measurement_gates": measurement,
        "measurement_gates_passed": measurement_passed,
        "candidate_members": candidate_members,
        "native_members": native_members,
        "candidate_repeatability": repeatability,
        "comparisons": comparisons,
        "primary_reduction": {
            "critical_folded_ionization_ratio": ionization_ratio,
            "critical_ionizing_tail_ratio": tail_ratio,
            "critical_phase_slice_folded_ionization_ratios": phase_ratios,
            "critical_total_variation_distance": critical[
                "total_variation_distance"],
            "absolute_log_error_fraction_of_previous": error_fraction,
        },
        "strong_closure_gate_passed": strong,
        "major_cadence_support_gate_passed": major,
        "formal_outcome": outcome,
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--ionization-table", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(
        args.rule.resolve(), args.lock.resolve(), args.report.resolve(),
        args.candidate_root.resolve(), args.native_root.resolve(),
        args.ionization_table.resolve())
    args.result.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
