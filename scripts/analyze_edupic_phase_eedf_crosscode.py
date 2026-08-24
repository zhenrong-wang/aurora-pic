#!/usr/bin/env python3
"""Compare regional phase EEDFs from native eduPIC and AuroraPIC."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

from analyze_aurorapic_ionizing_tail import kernel_frequency
from audit_edupic_ionization_path import read_cross_section
from run_aurorapic_edupic_pilot import atomic_json, sha256


PHASES = 200
REGIONS = ("x000_010", "x010_020", "x020_040", "x040_060",
           "x060_080", "x080_090", "x090_100")
ENERGY_BINS = 320
BIN_WIDTH_EV = 0.25


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def relative_range(values: list[float]) -> float:
    mean = statistics.fmean(values)
    return (max(values) - min(values)) / abs(mean) if mean else math.inf


def load_histogram(path: Path, represented_field: str
                   ) -> dict[tuple[int, str], list[float]]:
    output: dict[tuple[int, str], list[float]] = {}
    source = rows(path)
    if len(source) != PHASES * len(REGIONS) * ENERGY_BINS:
        raise ValueError(f"phase EEDF shape differs: {path}")
    seen: set[tuple[int, str, int]] = set()
    for row in source:
        phase, region, energy_bin = (int(row["phase_bin"]), row["region"],
                                     int(row["energy_bin"]))
        key = (phase, region, energy_bin)
        if (not 0 <= phase < PHASES or region not in REGIONS or
                not 0 <= energy_bin < ENERGY_BINS or key in seen):
            raise ValueError(f"invalid phase EEDF key: {path}")
        seen.add(key)
        energy = float(row["energy_eV"])
        count = float(row[represented_field])
        if (not math.isfinite(energy) or not math.isfinite(count) or
                not math.isclose(energy, (energy_bin + 0.5) * BIN_WIDTH_EV,
                                 rel_tol=0.0, abs_tol=1e-14)):
            raise ValueError(f"invalid phase EEDF value: {path}")
        output.setdefault((phase, region), [0.0] * ENERGY_BINS)[energy_bin] = count
    return output


def summarize(histogram: dict[tuple[int, str], list[float]],
              names: tuple[str, ...], first: int, end: int,
              energies: list[float], cross_sections: list[float]
              ) -> dict[str, object]:
    bins = [math.fsum(histogram[(phase, region)][energy_bin]
                      for phase in range(first, end) for region in names)
            for energy_bin in range(ENERGY_BINS)]
    total = math.fsum(bins)
    if total <= 0.0:
        raise ValueError("empty comparison phase EEDF scope")
    probabilities = [value / total for value in bins]
    centers = [(energy_bin + 0.5) * BIN_WIDTH_EV
               for energy_bin in range(ENERGY_BINS)]
    def fraction(low: float, high: float = math.inf) -> float:
        return math.fsum(probability for probability, energy in
                         zip(probabilities, centers) if low <= energy < high)
    return {
        "represented_histogram_observations": total,
        "histogram_mean_energy_eV": math.fsum(
            probability * energy for probability, energy in
            zip(probabilities, centers)),
        "fraction_11p5_to_15p8_eV": fraction(11.5, 15.8),
        "fraction_15p8_to_30_eV": fraction(15.8, 30.0),
        "fraction_at_or_above_30_eV": fraction(30.0),
        "fraction_at_or_above_15p8_eV": fraction(15.8),
        "eedf_folded_ionization_frequency_s-1": math.fsum(
            probability * kernel_frequency(energy, energies, cross_sections)
            for probability, energy in zip(probabilities, centers)),
        "probability_by_energy_bin": probabilities,
    }


def compare(candidate: dict[str, object], natives: list[dict[str, object]]) -> dict[str, object]:
    scalar_names = (
        "histogram_mean_energy_eV", "fraction_11p5_to_15p8_eV",
        "fraction_15p8_to_30_eV", "fraction_at_or_above_30_eV",
        "fraction_at_or_above_15p8_eV",
        "eedf_folded_ionization_frequency_s-1",
    )
    native_mean = {name: statistics.fmean(float(value[name]) for value in natives)
                   for name in scalar_names}
    ratios = {name: float(candidate[name]) / native_mean[name]
              if native_mean[name] else math.inf for name in scalar_names}
    native_probability = [statistics.fmean(
        float(value["probability_by_energy_bin"][energy_bin])
        for value in natives) for energy_bin in range(ENERGY_BINS)]
    total_variation = 0.5 * math.fsum(abs(float(a) - b) for a, b in zip(
        candidate["probability_by_energy_bin"], native_probability))
    return {
        "aurorapic": {name: candidate[name] for name in scalar_names},
        "native_edupic_ensemble_mean": native_mean,
        "aurorapic_to_native_edupic_ratio": ratios,
        "total_variation_distance": total_variation,
    }


def analyze(rule_path: Path, candidate_root: Path, ionization_table: Path,
            native_roots: dict[str, Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    provenance = rule["provenance"]
    if sha256(Path(__file__).resolve()) != provenance["analyzer_sha256"]:
        raise ValueError("phase EEDF analyzer differs from prospective rule")
    if sha256(ionization_table) != provenance["ionization_table_sha256"]:
        raise ValueError("ionization table differs")
    candidate_report = candidate_root / "block-report.json"
    if sha256(candidate_report) != provenance["aurorapic_block_report_sha256"]:
        raise ValueError("AuroraPIC phase EEDF report differs")
    report = json.loads(candidate_report.read_text(encoding="utf-8"))
    candidate_path = candidate_root / "output" / "phase_eedf.csv"
    if sha256(candidate_path) != provenance["aurorapic_phase_eedf_sha256"]:
        raise ValueError("AuroraPIC phase EEDF differs")
    if set(native_roots) != set(rule["execution"]["seeds"]):
        raise ValueError("native phase EEDF replicate set differs")
    energies, cross_sections = read_cross_section(ionization_table)
    candidate_histogram = load_histogram(candidate_path, "represented_count")
    native_histograms = []
    members = []
    for seed in rule["execution"]["seeds"]:
        root = native_roots[seed]
        histogram_path = root / "edupic_phase_eedf.csv"
        moments_path = root / "edupic_phase_eedf_moments.csv"
        histogram = load_histogram(histogram_path, "represented_count")
        moment_rows = rows(moments_path)
        if len(moment_rows) != PHASES * len(REGIONS):
            raise ValueError(f"native phase EEDF moments shape differs: {seed}")
        minimum = min(int(row["macro_observations"]) for row in moment_rows
                      if row["region"] in ("x020_040", "x040_060"))
        overflow = max(float(row["overflow_fraction"]) for row in moment_rows)
        native_histograms.append(histogram)
        members.append({
            "seed": int(seed),
            "output_sha256": {
                histogram_path.name: sha256(histogram_path),
                moments_path.name: sha256(moments_path),
            },
            "critical_minimum_macro_observations_per_region_phase_bin": minimum,
            "maximum_overflow_fraction": overflow,
        })

    scopes = {
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
    comparisons = {}
    native_values_by_scope = {}
    for name, (regions, first, end) in scopes.items():
        candidate = summarize(candidate_histogram, regions, first, end,
                              energies, cross_sections)
        natives = [summarize(histogram, regions, first, end,
                             energies, cross_sections)
                   for histogram in native_histograms]
        comparisons[name] = compare(candidate, natives)
        native_values_by_scope[name] = natives

    critical_native = native_values_by_scope[
        "critical_x020_to_x060_phase_0p125_to_0p5"]
    limits = rule["prospective_acceptance"]
    metrics = {
        "critical_folded_ionization_relative_range": relative_range([
            float(value["eedf_folded_ionization_frequency_s-1"])
            for value in critical_native]),
        "critical_ionizing_tail_relative_range": relative_range([
            float(value["fraction_at_or_above_15p8_eV"])
            for value in critical_native]),
    }
    gates = {
        "critical_folded_ionization_repeatability": metrics[
            "critical_folded_ionization_relative_range"] <= limits[
                "maximum_critical_folded_ionization_relative_range"],
        "critical_tail_repeatability": metrics[
            "critical_ionizing_tail_relative_range"] <= limits[
                "maximum_critical_tail_relative_range"],
        "critical_observation_population": all(member[
            "critical_minimum_macro_observations_per_region_phase_bin"] >=
            limits["minimum_macro_observations_per_critical_region_phase_bin"]
            for member in members),
        "histogram_overflow": all(member["maximum_overflow_fraction"] <=
            limits["maximum_histogram_overflow_fraction"] for member in members),
    }
    critical_ratio = comparisons[
        "critical_x020_to_x060_phase_0p125_to_0p5"][
            "aurorapic_to_native_edupic_ratio"]
    upstream_ratio = comparisons[
        "upstream_x010_to_x020_phase_0p125_to_0p5"][
            "aurorapic_to_native_edupic_ratio"]
    hypothesis = rule["prospective_hypothesis"]
    outcome = {
        "critical_ionization_kernel_deficit": critical_ratio[
            "eedf_folded_ionization_frequency_s-1"] <= hypothesis[
                "maximum_critical_folded_ionization_ratio"],
        "critical_tail_deficit": critical_ratio[
            "fraction_at_or_above_15p8_eV"] <= hypothesis[
                "maximum_critical_tail_ratio"],
        "deficit_develops_beyond_upstream_region": upstream_ratio[
            "eedf_folded_ionization_frequency_s-1"] > critical_ratio[
                "eedf_folded_ionization_frequency_s-1"],
    }
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_regional_phase_eedf",
        "provenance": {
            "rule_sha256": sha256(rule_path),
            "analyzer_sha256": provenance["analyzer_sha256"],
            "aurorapic_block_report_sha256": sha256(candidate_report),
            "aurorapic_phase_eedf_sha256": sha256(candidate_path),
            "ionization_table_sha256": sha256(ionization_table),
        },
        "aurorapic_global_gate_context": report["gates"],
        "members": members,
        "comparisons": comparisons,
        "repeatability_metrics": metrics,
        "gates": gates,
        "all_measurement_gates_passed": all(gates.values()),
        "prospective_hypothesis_outcome": outcome,
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("aurorapic_root", type=Path)
    parser.add_argument("ionization_table", type=Path)
    parser.add_argument("native_runs_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rule = json.loads(args.rule.read_text(encoding="utf-8"))
    roots = {str(seed): args.native_runs_root / f"seed-{seed}"
             for seed in rule["execution"]["seeds"]}
    result = analyze(args.rule.resolve(), args.aurorapic_root.resolve(),
                     args.ionization_table.resolve(), roots)
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_measurement_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
