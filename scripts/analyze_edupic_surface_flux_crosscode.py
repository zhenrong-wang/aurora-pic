#!/usr/bin/env python3
"""Compare native eduPIC and AuroraPIC internal electron surface transport."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

from run_aurorapic_edupic_pilot import atomic_json, sha256


ELEMENTARY_CHARGE_C = 1.602176634e-19
DIRECTIONS = ("left_to_right", "right_to_left")
SURFACES_M = (0.005, 0.015)
PHASE_BINS = 200
ENERGY_BINS = 320
IONIZATION_THRESHOLD_EV = 15.8


def relative_range(values: list[float]) -> float:
    mean = statistics.fmean(values)
    return (max(values) - min(values)) / abs(mean) if mean else math.inf


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def summarize_run(root: Path, cycles: int, frequency_hz: float) -> dict[str, object]:
    summary_path = root / "edupic_phase_surface_flux_summary.csv"
    histogram_path = root / "edupic_phase_surface_flux.csv"
    summary_rows = read_rows(summary_path)
    histogram_rows = read_rows(histogram_path)
    if len(summary_rows) != PHASE_BINS * 2 * 2:
        raise ValueError(f"{root.name}: surface summary shape differs")
    if len(histogram_rows) != PHASE_BINS * 2 * 2 * ENERGY_BINS:
        raise ValueError(f"{root.name}: surface histogram shape differs")

    summary: dict[tuple[int, int, str], dict[str, str]] = {}
    macro_by_surface = [0, 0]
    maximum_overflow = 0.0
    for row in summary_rows:
        phase = int(row["phase_bin"])
        surface = int(row["surface_id"])
        direction = row["direction"]
        key = (phase, surface, direction)
        if (not 0 <= phase < PHASE_BINS or not 0 <= surface < 2 or
                direction not in DIRECTIONS or key in summary):
            raise ValueError(f"{root.name}: invalid or duplicate summary key")
        values = [float(value) for name, value in row.items()
                  if name not in ("direction",) and value is not None]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{root.name}: non-finite summary value")
        if not math.isclose(float(row["position_m"]), SURFACES_M[surface],
                            rel_tol=0.0, abs_tol=1e-15):
            raise ValueError(f"{root.name}: surface position differs")
        expected_phase = (phase + 0.5) / PHASE_BINS
        if not math.isclose(float(row["phase_fraction"]), expected_phase,
                            rel_tol=0.0, abs_tol=1e-15):
            raise ValueError(f"{root.name}: phase coordinate differs")
        macro_by_surface[surface] += int(row["macro_crossings"])
        maximum_overflow = max(maximum_overflow,
                               float(row["overflow_fraction"]))
        summary[key] = row

    phase_duration = cycles / frequency_hz / PHASE_BINS
    tail_energy: dict[tuple[int, int, str], float] = {}
    seen: set[tuple[int, int, str, int]] = set()
    for row in histogram_rows:
        phase = int(row["phase_bin"])
        surface = int(row["surface_id"])
        direction = row["direction"]
        energy_bin = int(row["energy_bin"])
        key = (phase, surface, direction, energy_bin)
        if (not 0 <= phase < PHASE_BINS or not 0 <= surface < 2 or
                direction not in DIRECTIONS or
                not 0 <= energy_bin < ENERGY_BINS or key in seen):
            raise ValueError(f"{root.name}: invalid or duplicate histogram key")
        seen.add(key)
        energy = float(row["energy_eV"])
        represented = float(row["represented_crossings_m-2"])
        if not math.isfinite(energy) or not math.isfinite(represented):
            raise ValueError(f"{root.name}: non-finite histogram value")
        expected_energy = (energy_bin + 0.5) * 80.0 / ENERGY_BINS
        if not math.isclose(energy, expected_energy, rel_tol=0.0,
                            abs_tol=1e-14):
            raise ValueError(f"{root.name}: histogram energy coordinate differs")
        if energy >= IONIZATION_THRESHOLD_EV:
            scalar = (phase, surface, direction)
            tail_energy[scalar] = tail_energy.get(scalar, 0.0) + (
                represented * energy * ELEMENTARY_CHARGE_C / phase_duration)

    def window(first: int, end: int) -> dict[str, float]:
        direct_net: list[float] = []
        tail_net: list[float] = []
        for surface in range(2):
            direct_net.append(math.fsum(
                float(summary[(phase, surface, "left_to_right")]
                      ["kinetic_energy_flux_W_m-2"]) -
                float(summary[(phase, surface, "right_to_left")]
                      ["kinetic_energy_flux_W_m-2"])
                for phase in range(first, end)) / (end - first))
            tail_net.append(math.fsum(
                tail_energy.get((phase, surface, "left_to_right"), 0.0) -
                tail_energy.get((phase, surface, "right_to_left"), 0.0)
                for phase in range(first, end)) / (end - first))
        return {
            "direct_outward_energy_flux_divergence_W_m-2":
                direct_net[1] - direct_net[0],
            "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2":
                tail_net[1] - tail_net[0],
        }

    return {
        "output_sha256": {
            summary_path.name: sha256(summary_path),
            histogram_path.name: sha256(histogram_path),
        },
        "macro_crossings_by_surface": macro_by_surface,
        "maximum_overflow_fraction": maximum_overflow,
        "critical_phase_0p125_to_0p5": window(25, 100),
        "exceptional_phase_0p375_to_0p5": window(75, 100),
    }


def evaluate(rule: dict[str, object], aurora: dict[str, object],
             members: list[dict[str, object]]) -> tuple[dict[str, object],
                                                        dict[str, bool]]:
    direct = "direct_outward_energy_flux_divergence_W_m-2"
    tail = "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2"
    critical = [float(member["critical_phase_0p125_to_0p5"][direct])
                for member in members]
    exceptional = [float(member["exceptional_phase_0p375_to_0p5"][direct])
                   for member in members]
    exceptional_tail = [float(member["exceptional_phase_0p375_to_0p5"][tail])
                        for member in members]
    aurora_metrics = aurora["metrics"]
    density_ratio = float(rule["comparison"][
        "aurorapic_to_edupic_electron_density_ratio"])
    aurora_normalized = {
        "critical_direct_W_m-2": float(aurora_metrics[
            "mean_critical_phase_direct_flux_W_m-2"]) / density_ratio,
        "exceptional_direct_W_m-2": float(aurora_metrics[
            "mean_exceptional_octant_direct_flux_W_m-2"]) / density_ratio,
        "exceptional_tail_W_m-2": float(aurora_metrics[
            "mean_exceptional_octant_tail_flux_W_m-2"]) / density_ratio,
    }
    means = {
        "critical_direct_W_m-2": statistics.fmean(critical),
        "exceptional_direct_W_m-2": statistics.fmean(exceptional),
        "exceptional_tail_W_m-2": statistics.fmean(exceptional_tail),
    }
    metrics: dict[str, object] = {
        "edupic_ensemble_mean": means,
        "edupic_relative_range": {
            "critical_direct": relative_range(critical),
            "exceptional_direct": relative_range(exceptional),
            "exceptional_tail": relative_range(exceptional_tail),
        },
        "aurorapic_density_normalized": aurora_normalized,
        "density_normalized_aurorapic_to_edupic_ratio": {
            key: aurora_normalized[key] / means[key]
            if means[key] != 0.0 else math.inf for key in means
        },
        "absolute_aurorapic_to_edupic_ratio": {
            "critical_direct": float(aurora_metrics[
                "mean_critical_phase_direct_flux_W_m-2"]) / means[
                    "critical_direct_W_m-2"],
            "exceptional_direct": float(aurora_metrics[
                "mean_exceptional_octant_direct_flux_W_m-2"]) / means[
                    "exceptional_direct_W_m-2"],
            "exceptional_tail": float(aurora_metrics[
                "mean_exceptional_octant_tail_flux_W_m-2"]) / means[
                    "exceptional_tail_W_m-2"],
        },
    }
    limits = rule["prospective_acceptance"]
    ranges = metrics["edupic_relative_range"]
    gates = {
        "critical_repeatability": ranges["critical_direct"] <= limits[
            "maximum_critical_direct_relative_range"],
        "exceptional_repeatability": ranges["exceptional_direct"] <= limits[
            "maximum_exceptional_direct_relative_range"],
        "tail_repeatability": ranges["exceptional_tail"] <= limits[
            "maximum_exceptional_tail_relative_range"],
        "crossing_population": all(
            min(member["macro_crossings_by_surface"]) >= limits[
                "minimum_macro_crossings_each_surface"] for member in members),
        "histogram_overflow": all(
            member["maximum_overflow_fraction"] <= limits[
                "maximum_histogram_overflow_fraction"] for member in members),
    }
    return metrics, gates


def analyze(rule_path: Path, aurora_path: Path,
            roots: dict[str, Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    aurora = json.loads(aurora_path.read_text(encoding="utf-8"))
    provenance = rule["provenance"]
    if sha256(Path(__file__).resolve()) != provenance["analyzer_sha256"]:
        raise ValueError("cross-code analyzer differs from prospective rule")
    if sha256(aurora_path) != provenance["aurorapic_result_sha256"]:
        raise ValueError("AuroraPIC comparison result differs")
    if set(roots) != set(rule["execution"]["seeds"]):
        raise ValueError("native eduPIC replicate set differs")
    members = []
    for seed in rule["execution"]["seeds"]:
        member = summarize_run(roots[seed], int(rule["execution"][
            "measurement_cycles"]), float(rule["case"]["rf_frequency_hz"]))
        member["seed"] = int(seed)
        members.append(member)
    metrics, gates = evaluate(rule, aurora, members)
    ratios = metrics["density_normalized_aurorapic_to_edupic_ratio"]
    hypothesis = rule["prospective_hypothesis"]
    hypothesis_outcome = {
        "exceptional_direct_materially_larger": ratios[
            "exceptional_direct_W_m-2"] >= hypothesis[
                "minimum_exceptional_direct_ratio"],
        "exceptional_tail_materially_larger": ratios[
            "exceptional_tail_W_m-2"] >= hypothesis[
                "minimum_exceptional_tail_ratio"],
    }
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_internal_surface_transport",
        "provenance": {
            "rule_sha256": sha256(rule_path),
            "analyzer_sha256": provenance["analyzer_sha256"],
            "aurorapic_result_sha256": provenance["aurorapic_result_sha256"],
        },
        "members": members,
        "metrics": metrics,
        "gates": gates,
        "all_measurement_gates_passed": all(gates.values()),
        "prospective_hypothesis_outcome": hypothesis_outcome,
        "hypothesis_supported": any(hypothesis_outcome.values()),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("aurorapic_result", type=Path)
    parser.add_argument("runs_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rule = json.loads(args.rule.read_text(encoding="utf-8"))
    roots = {str(seed): args.runs_root / f"seed-{seed}"
             for seed in rule["execution"]["seeds"]}
    result = analyze(args.rule.resolve(), args.aurorapic_result.resolve(), roots)
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_measurement_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
