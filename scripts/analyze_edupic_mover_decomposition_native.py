#!/usr/bin/env python3
"""Seal the native half of the locked mover-decomposition comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


METRICS = (
    "mean_origin_energy_eV",
    "origin_longitudinal_energy_fraction",
    "mean_positive_linear_work_eV",
    "mean_quadratic_work_eV",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path, region: str, lower: float,
         upper: float) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            row for row in csv.DictReader(stream)
            if row["region"] == region
            and lower <= float(row["phase_fraction"]) < upper
        ]


def total(values: list[dict[str, str]], field: str) -> float:
    return sum(float(row[field]) for row in values)


def relative_range(members: list[dict[str, float]], metric: str) -> float:
    values = [float(member[metric]) for member in members]
    mean = sum(values) / len(values)
    return (max(values) - min(values)) / max(abs(mean), 1.0e-300)


def resource(path: Path, label: str) -> float:
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"{re.escape(label)}:\s*([^\n]+)", text)
    if not match:
        raise RuntimeError(f"missing resource field {label} in {path}")
    value = match.group(1).strip()
    if label == "Elapsed (wall clock) time (h:mm:ss or m:ss)":
        parts = [float(part) for part in value.split(":")]
        return sum(part * 60 ** power
                   for power, part in enumerate(reversed(parts)))
    return float(value)


def member(directory: Path, seed: int, region: str,
           lower: float, upper: float) -> dict[str, float | int | str]:
    decomposition_path = directory / "edupic_phase_eedf_mover_decomposition.csv"
    work_path = directory / "edupic_phase_eedf_promotion_band_work.csv"
    decomposition = rows(decomposition_path, region, lower, upper)
    work = rows(work_path, region, lower, upper)
    if len(decomposition) != len(work) or not decomposition:
        raise RuntimeError(f"decomposition/work shape mismatch in {directory}")
    observations = int(total(
        decomposition, "field_push_promotion_band_observations"))
    work_observations = int(total(
        work, "field_push_promotion_band_observations"))
    if observations != work_observations:
        raise RuntimeError(f"decomposition/work population mismatch in {directory}")
    promotions = int(total(
        work, "field_push_promotion_band_promotions"))
    origin = total(
        decomposition,
        "field_push_promotion_band_origin_macro_energy_sum_eV")
    longitudinal = total(
        decomposition,
        "field_push_promotion_band_origin_longitudinal_macro_energy_sum_eV")
    linear = total(
        decomposition,
        "field_push_promotion_band_linear_macro_work_sum_eV")
    positive_linear = total(
        decomposition,
        "field_push_promotion_band_positive_linear_macro_work_sum_eV")
    negative_linear = total(
        decomposition,
        "field_push_promotion_band_negative_linear_macro_work_sum_eV")
    quadratic = total(
        decomposition,
        "field_push_promotion_band_quadratic_macro_work_sum_eV")
    signed = total(
        work, "field_push_promotion_band_signed_macro_work_sum_eV")
    return {
        "seed": seed,
        "final_checkpoint_sha256": sha256(directory / "picdata.bin"),
        "promotion_band_work_sha256": sha256(work_path),
        "mover_decomposition_sha256": sha256(decomposition_path),
        "promotion_band_observations": observations,
        "promotion_band_promotions": promotions,
        "origin_energy_sum_eV": origin,
        "origin_longitudinal_energy_sum_eV": longitudinal,
        "linear_work_sum_eV": linear,
        "positive_linear_work_sum_eV": positive_linear,
        "negative_linear_work_sum_eV": negative_linear,
        "quadratic_work_sum_eV": quadratic,
        "mean_origin_energy_eV": origin / observations,
        "origin_longitudinal_energy_fraction": longitudinal / origin,
        "mean_positive_linear_work_eV": positive_linear / observations,
        "mean_quadratic_work_eV": quadratic / observations,
        "linear_work_closure_residual_eV":
            linear - (positive_linear - negative_linear),
        "total_work_decomposition_closure_residual_eV":
            signed - (linear + quadratic),
        "wall_seconds": resource(
            directory / "stderr.txt",
            "Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "peak_resident_set_kib": int(resource(
            directory / "stderr.txt", "Maximum resident set size (kbytes)")),
    }


def analyze(rule_path: Path, prior_path: Path,
            directories: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    prior = json.loads(prior_path.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    seeds = rule["locked_inputs"]["native_seeds"]
    members = [
        member(directory, seed, diagnostic["critical_region"], lower, upper)
        for directory, seed in zip(directories, seeds, strict=True)
    ]
    prior_by_seed = {entry["seed"]: entry for entry in prior["members"]}
    checkpoint_hashes = rule["locked_inputs"][
        "expected_passive_native_checkpoint_sha256"]
    repeatability = {
        metric: relative_range(members, metric) for metric in METRICS
    }
    closure_scale = [
        max(abs(float(entry["positive_linear_work_sum_eV"])) +
            abs(float(entry["negative_linear_work_sum_eV"])), 1.0)
        for entry in members
    ]
    total_scale = [
        max(abs(float(entry["linear_work_sum_eV"])) +
            abs(float(entry["quadratic_work_sum_eV"])), 1.0)
        for entry in members
    ]
    tolerance = float(diagnostic["work_closure_relative_tolerance"])
    gates = {
        "member_count": len(members) == len(seeds) == 3,
        "passive_final_checkpoints": all(
            entry["final_checkpoint_sha256"] == checkpoint_hashes[str(seed)]
            for entry, seed in zip(members, seeds, strict=True)),
        "prior_work_ledgers_exact": all(
            entry["promotion_band_work_sha256"] ==
                prior_by_seed[seed]["work_ledger_sha256"]
            for entry, seed in zip(members, seeds, strict=True)),
        "prior_populations_exact": all(
            entry["promotion_band_observations"] ==
                prior_by_seed[seed]["promotion_band_observations"]
            and entry["promotion_band_promotions"] ==
                prior_by_seed[seed]["promotion_band_promotions"]
            for entry, seed in zip(members, seeds, strict=True)),
        "population_floors": all(
            int(entry["promotion_band_observations"]) >=
                int(diagnostic["minimum_band_observations_per_member"])
            and int(entry["promotion_band_promotions"]) >=
                int(diagnostic["minimum_band_promotions_per_member"])
            for entry in members),
        "finite_metrics": all(
            math.isfinite(float(entry[metric]))
            for entry in members for metric in METRICS),
        "origin_energy_partition": all(
            0.0 <= float(entry["origin_longitudinal_energy_fraction"]) <= 1.0
            for entry in members),
        "linear_work_closure": all(
            abs(float(entry["linear_work_closure_residual_eV"])) <=
                tolerance * scale
            for entry, scale in zip(members, closure_scale, strict=True)),
        "total_work_decomposition_closure": all(
            abs(float(entry[
                "total_work_decomposition_closure_residual_eV"])) <=
                tolerance * scale
            for entry, scale in zip(members, total_scale, strict=True)),
        "native_repeatability": all(value <= 0.15
                                    for value in repeatability.values()),
        "resident_memory": all(
            int(entry["peak_resident_set_kib"]) <=
                int(rule["execution_contract"][
                    "maximum_peak_resident_set_kib"])
            for entry in members),
    }
    means = {
        metric: sum(float(entry[metric]) for entry in members) / len(members)
        for metric in METRICS
    }
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_near_threshold_mover_decomposition_ensemble_result",
        "rule_sha256": sha256(rule_path),
        "prior_promotion_band_result_sha256": sha256(prior_path),
        "gates": gates,
        "all_native_integrity_population_repeatability_and_closure_gates_passed":
            all(gates.values()),
        "critical_scope": {
            "region": diagnostic["critical_region"],
            "phase_fraction": diagnostic["critical_phase_fraction"],
            "promotion_band_eV": diagnostic["promotion_band_eV"],
        },
        "members": members,
        "ensemble_means": means,
        "relative_ranges": repeatability,
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "mover_decomposition_sha256": [
                entry["mover_decomposition_sha256"] for entry in members
            ],
        },
        "physics_claim": "none_until_both_locked_aurorapic_members_and_all_joint_gates_complete",
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("prior_result", type=Path)
    parser.add_argument("native_directories", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.prior_result, args.native_directories)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_native_integrity_population_repeatability_and_closure_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
