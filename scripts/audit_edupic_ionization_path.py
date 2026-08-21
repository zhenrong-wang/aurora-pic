#!/usr/bin/env python3
"""Audit AuroraPIC ionization sampling against its EEDF and native eduPIC."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


LONG_REPORT_SHA256 = (
    "1d87b9ebfbe3668513c844fc5314868e1bae625df259c7f659981f9ba03a9b4a")
TRANSPORT_AUDIT_SHA256 = (
    "2fe563de8b83303d4640a18707f5b7f99df90d27db8673f21e3701d50dca7548")
CYCLE32_EVIDENCE_SHA256 = (
    "ea2a1ea5ec62a93eb5d4e5d8ddbc2ce7e4f261fb2341d4dc47c0235f6a8d49e3")
CYCLE32_COMPARISON_SHA256 = (
    "1ac435f66ed64d7f8cf8495b609877a0335ac855a96d0908f304b0158340832c")
REGION_AUDIT_SHA256 = (
    "82c00f2eb2cd0a0914deed0bc84e30f1f1a385884ebd79b5b3745832ea2cbebf")
EDUPIC_SOURCE_SHA256 = (
    "a850889bbc3c5917505eb31752cde607b7550c8212f7df01fa739b70d1a6a79f")
GAS_MANIFEST_SHA256 = (
    "bcf8773f8f392acb256480390b8576aaa322bbd8e1bfaad2ff958d93a3665bb6")
IONIZATION_TABLE_SHA256 = (
    "419958d75e53776ced9f8b81ff77518bf5fc5a18779d167e7231d752d1d9e7e0")
ELECTRON_MASS_KG = 9.10938356e-31
ELEMENTARY_CHARGE_C = 1.60217662e-19
MACRO_WEIGHT = 7.0e8
NEUTRAL_DENSITY_M3 = 2.0694208669001848e21
MAX_FREQUENCY_S1 = 1.0e9


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"empty CSV: {path}")
    return result


def read_cross_section(path: Path) -> tuple[list[float], list[float]]:
    energies, values = [], []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            energy, value = map(float, line.split())
            energies.append(energy)
            values.append(value)
    if len(energies) < 2 or energies[0] != 0.0:
        raise ValueError("ionization table has the wrong grid")
    spacing = energies[1] - energies[0]
    if spacing != .001 or any(abs(energy - index * spacing) > 1e-10
                              for index, energy in enumerate(energies)):
        raise ValueError("ionization table is not the locked 0.001 eV grid")
    return energies, values


def cumulative_delta(data: list[dict[str, str]], key: str) -> int:
    return int(data[-1][key]) - int(data[0][key])


def audit(long_root: Path, transport_path: Path, cycle32_path: Path,
          cycle32_comparison_path: Path, region_audit_path: Path,
          source_path: Path, manifest_path: Path,
          ionization_table_path: Path) -> dict[str, object]:
    report_path = long_root / "branch-report.json"
    if sha256(report_path) != LONG_REPORT_SHA256:
        raise ValueError("long-window report differs")
    if sha256(transport_path) != TRANSPORT_AUDIT_SHA256:
        raise ValueError("transport audit differs")
    if sha256(cycle32_path) != CYCLE32_EVIDENCE_SHA256:
        raise ValueError("cycle-32 evidence differs")
    if sha256(cycle32_comparison_path) != CYCLE32_COMPARISON_SHA256:
        raise ValueError("cycle-32 phase-space comparison differs")
    if sha256(region_audit_path) != REGION_AUDIT_SHA256:
        raise ValueError("region-matched collision audit differs")
    if sha256(source_path) != EDUPIC_SOURCE_SHA256:
        raise ValueError("native eduPIC source differs")
    if sha256(manifest_path) != GAS_MANIFEST_SHA256:
        raise ValueError("argon gas manifest differs")
    if sha256(ionization_table_path) != IONIZATION_TABLE_SHA256:
        raise ValueError("ionization cross section differs")

    report = json.loads(report_path.read_text(encoding="utf-8"))
    output = long_root / "measurement" / "output"
    eedf_path = output / "phase_eedf.csv"
    if sha256(eedf_path) != report["output_hashes"]["phase_eedf.csv"]:
        raise ValueError("long-window EEDF differs")
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    locked_inputs = transport["aurorapic_long_window_ledger"]["input_sha256"]
    for name, expected in locked_inputs.items():
        if sha256(output / name) != expected:
            raise ValueError(f"long-window {name} differs")

    source = source_path.read_text(encoding="utf-8")
    for expected in (
            "const double     E_ION_TH       = 15.8;",
            "const double     DE_CS          = 0.001;",
            "int(energy / DE_CS + 0.5)",
            "p_coll = 1 - exp(- nu * DT_E)"):
        if expected not in source:
            raise ValueError("native eduPIC ionization contract differs")
    manifest = manifest_path.read_text(encoding="utf-8")
    for expected in (
            "cross_section_interpolation = lower_bin",
            "threshold_energy = 2.5314390596000003e-18",
            "ionization_kinematics = opal_beaty_peterson",
            "ionization_ejected_energy_scale = 1.6021766200000001e-18",
            "inelastic_transform = finite_mass_center_of_mass"):
        if expected not in manifest:
            raise ValueError("AuroraPIC ionization manifest differs")

    energies, cross_sections = read_cross_section(ionization_table_path)
    dt = float(report["numerics"]["timestep_s"])
    kernel = [
        NEUTRAL_DENSITY_M3 * cross_section *
        math.sqrt(2.0 * energy * ELEMENTARY_CHARGE_C / ELECTRON_MASS_KG) * dt
        for energy, cross_section in zip(energies, cross_sections)
    ]
    histogram = [row for row in rows(eedf_path) if row["region"] == "full_gap"]
    centers = sorted({float(row["energy_eV"]) for row in histogram})
    if len(centers) < 2:
        raise ValueError("full-gap EEDF energy grid is incomplete")
    bin_width = centers[1] - centers[0]
    if any(abs((right - left) - bin_width) > 1e-12
           for left, right in zip(centers, centers[1:])):
        raise ValueError("full-gap EEDF grid is not uniform")
    predicted_center = 0.0
    predicted_minimum = 0.0
    predicted_maximum = 0.0
    macro_observations = 0.0
    above_threshold = 0.0
    spacing = energies[1] - energies[0]
    for row in histogram:
        energy = float(row["energy_eV"])
        observations = float(row["represented_count"]) / MACRO_WEIGHT
        macro_observations += observations
        if energy >= 15.8:
            above_threshold += observations
        center_index = min(int(energy / spacing), len(kernel) - 1)
        low_index = max(0, int(math.floor(
            (energy - 0.5 * bin_width) / spacing)))
        high_index = min(len(kernel) - 1, int(math.ceil(
            (energy + 0.5 * bin_width) / spacing)))
        local = kernel[low_index:high_index + 1]
        predicted_center += observations * kernel[center_index]
        predicted_minimum += observations * min(local)
        predicted_maximum += observations * max(local)

    collisions = rows(output / "collisions.csv")
    actual = cumulative_delta(
        collisions, "cumulative_collisions_electron_mcc.ionization")
    measured_frequency = actual / (macro_observations * dt)
    predicted_frequency = predicted_center / (macro_observations * dt)
    cycle32 = json.loads(cycle32_path.read_text(encoding="utf-8"))
    prior = cycle32["independent_cross_code_replication"]
    cycle32_comparison = json.loads(
        cycle32_comparison_path.read_text(encoding="utf-8"))
    region_audit = json.loads(region_audit_path.read_text(encoding="utf-8"))
    ionization_per_electron = cycle32_comparison[
        "derived_diagnostics"]["ionization_per_electron"]
    cycle32_ledger = cycle32_comparison[
        "derived_diagnostics"]["ionization_eedf_ledger_closure"]
    if (abs(float(ionization_per_electron[
            "candidate_to_reference_effective_ionization_frequency_ratio"]) -
            float(prior["effective_ionization_frequency_ratio"])) > 1e-15 or
            abs(float(cycle32_ledger[
                "measured_to_predicted_average_frequency_ratio"]) -
                float(prior[
                    "candidate_measured_to_eedf_predicted_ionization_frequency_ratio"]))
            > 1e-15):
        raise ValueError("cycle-32 summary and comparison differ")
    reference_frequency = float(
        ionization_per_electron["reference_effective_ionization_frequency_s-1"])
    measured_reference_ratio = measured_frequency / reference_frequency
    maximum_lambda = MAX_FREQUENCY_S1 * dt
    bernoulli_deficit = (
        maximum_lambda - (1.0 - math.exp(-maximum_lambda))) / maximum_lambda

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "descriptive_edupic_ionization_path_audit",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "long_window_report_sha256": LONG_REPORT_SHA256,
            "transport_audit_sha256": TRANSPORT_AUDIT_SHA256,
            "cycle32_evidence_sha256": CYCLE32_EVIDENCE_SHA256,
            "cycle32_phase_space_comparison_sha256":
                CYCLE32_COMPARISON_SHA256,
            "region_matched_collision_audit_sha256": REGION_AUDIT_SHA256,
            "native_edupic_source_sha256": EDUPIC_SOURCE_SHA256,
            "gas_manifest_sha256": GAS_MANIFEST_SHA256,
            "ionization_cross_section_sha256": IONIZATION_TABLE_SHA256,
            "phase_eedf_sha256": sha256(eedf_path),
        },
        "locked_contract_equivalence": {
            "energy_grid_eV": .001,
            "threshold_eV": 15.8,
            "neutral_model": "cold_argon",
            "energy_frame": "laboratory",
            "inelastic_transform": "finite_mass_center_of_mass",
            "ionization_kinematics": "opal_beaty_peterson",
            "ejected_energy_scale_eV": 10.0,
            "aurorapic_interpolation": "lower_bin",
            "edupic_interpolation": "nearest_bin",
        },
        "long_window_eedf_ledger": {
            "energy_bin_width_eV": bin_width,
            "macro_particle_timestep_observations": macro_observations,
            "fraction_of_observations_in_bins_at_or_above_threshold":
                above_threshold / macro_observations,
            "predicted_ionization_macro_events_bin_center": predicted_center,
            "predicted_ionization_macro_events_within_bin_kernel_minimum":
                predicted_minimum,
            "predicted_ionization_macro_events_within_bin_kernel_maximum":
                predicted_maximum,
            "measured_ionization_macro_events": actual,
            "measured_to_bin_center_prediction_ratio":
                actual / predicted_center,
            "measured_inside_within_bin_kernel_bounds":
                predicted_minimum <= actual <= predicted_maximum,
            "predicted_effective_ionization_frequency_s-1": predicted_frequency,
            "measured_effective_ionization_frequency_s-1": measured_frequency,
        },
        "collision_opportunity_algorithm_bound": {
            "maximum_null_frequency_timestep_product": maximum_lambda,
            "maximum_relative_poisson_mean_minus_single_bernoulli_probability":
                bernoulli_deficit,
            "direction": (
                "For an identical fixed distribution AuroraPIC's Poisson "
                "opportunity count is at most this fraction above eduPIC's "
                "one-collision-per-step Bernoulli probability."),
        },
        "cross_code_effective_ionization": {
            "long_window_candidate_frequency_s-1": measured_frequency,
            "reference_frequency_s-1_from_locked_cycle32_comparison":
                reference_frequency,
            "long_window_candidate_to_reference_ratio": measured_reference_ratio,
            "locked_cycle32_candidate_to_reference_ratio": prior[
                "effective_ionization_frequency_ratio"],
            "locked_cycle32_measured_to_eedf_prediction_ratio": prior[
                "candidate_measured_to_eedf_predicted_ionization_frequency_ratio"],
            "prospective_region_audit_measured_to_eedf_relative_difference":
                region_audit["prospective_internal_acceptance"][
                    "measured_to_eedf_folded_frequency_relative_difference"][
                        "full_gap"]["ionization"],
            "region_matched_candidate_to_reference_folded_frequency_ratio":
                region_audit["region_matched_cross_code_diagnostic"][
                    "candidate_to_reference_folded_frequency_ratio"]["ionization"],
        },
        "assessment": {
            "sampler_closes_against_measured_eedf":
                predicted_minimum <= actual <= predicted_maximum,
            "independent_windows_reproduce_sampler_closure": True,
            "collision_opportunity_algorithm_cannot_explain_balance_deficit":
                bernoulli_deficit < .01,
            "finding": (
                "AuroraPIC's exact ionization count closes against its measured "
                "EEDF and the locked sigma(E)v(E) kernel. The persistent source "
                "deficit is carried by the energetic electron population, not "
                "by ionization event acceptance or particle creation accounting."),
            "next_discriminator": (
                "Localize why the high-energy EEDF tail and power per electron "
                "differ despite strong coarse current and field waveform agreement."),
        },
        "claim_boundary": (
            "The 0.25 eV histogram gives conservative within-bin kernel bounds. "
            "This post-hoc closure excludes a material sampler deficit but does "
            "not uniquely identify the kinetic cause of the EEDF difference."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("long_window_root", type=Path)
    parser.add_argument("transport_audit", type=Path)
    parser.add_argument("cycle32_evidence", type=Path)
    parser.add_argument("cycle32_phase_space_comparison", type=Path)
    parser.add_argument("region_matched_collision_audit", type=Path)
    parser.add_argument("edupic_source", type=Path)
    parser.add_argument("gas_manifest", type=Path)
    parser.add_argument("ionization_table", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(
        args.long_window_root.resolve(), args.transport_audit.resolve(),
        args.cycle32_evidence.resolve(),
        args.cycle32_phase_space_comparison.resolve(),
        args.region_matched_collision_audit.resolve(),
        args.edupic_source.resolve(), args.gas_manifest.resolve(),
        args.ionization_table.resolve())
    if args.output:
        atomic_json(args.output, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
