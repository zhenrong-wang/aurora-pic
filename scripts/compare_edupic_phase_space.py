#!/usr/bin/env python3
"""Compare AuroraPIC phase-space diagnostics with locked eduPIC matrices."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


NODES = 400
REFERENCE_PHASES = 200
CANDIDATE_PHASES = 16
LENGTH_M = 0.025
ELEMENTARY_CHARGE_C = 1.60217662e-19
EDUPIC_IMPLEMENTATION_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41")

REFERENCE_FILES = {
    "electron_density": ("Figure11a)_ne_xt.dat",
                         "43e70f92717daf5ccb9a7fb491b9f0ba42b952bb24bc851dcab84401f99619e8"),
    "ion_density": ("Figure11b)_ni_xt.dat",
                    "ba1529cf97c340748f78611c3e0f05bbcfd6d7b1fdd2e582c321e98369f454e4"),
    "potential": ("Figure11c)_pot_xt.dat",
                  "e2eaad4698f3113c922910066288e3282848fd5810de09592fa3f0a5a528d812"),
    "electric_field": ("Figure11d)_efield_xt.dat",
                       "41fa8afa93501473ff69cdbe4ac5f1e29377823a1518e517c16ec973360a6186"),
    "electron_current_density": ("Figure11e)_je_xt.dat",
                                 "87ea3b30a52d540546e21317f1aff9a9f1705b1a7c9b4107a271657bd1bf2ea1"),
    "ion_current_density": ("Figure11f)_ji_xt.dat",
                            "77fea7696ebb2d8895625d5dbcf9c72054b73c32ac3ec60f9fdecac6dc0abdec"),
    "electron_ohmic_power_density": ("Figure11g)_heate_xt.dat",
                                     "8b8c6c453d5c1c0865f0697938794377cd8b695154aca97714d304060e7ca335"),
    "ion_ohmic_power_density": ("Figure11h)_heati_xt.dat",
                                "cff2f8d20171e5d87748103acdf991dc80fafacb99de1aeeb147e49c6510a7ac"),
    "electron_mean_energy": ("Figure11i)_meanee_xt.dat",
                             "3e334607f0c23fee92c36e2868714de99f92fa841096b51b9636955f9dcf9d05"),
    "ion_mean_energy": ("Figure11j)_meanei_xt.dat",
                        "479fcb42b8970fc68ba8520ea8a404c8850667cbbfbdf139e605e9d70b9f13fa"),
    "ionization_rate": ("Figure11k)_ioniz_xt.dat",
                        "59baed449b6972eddbd719ad6019fe9c9e6eb981c2ecf952870f0b3925856be4"),
}


def read_matrix(path: Path, rows: int = NODES,
                columns: int = REFERENCE_PHASES) -> list[list[float]]:
    matrix = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            values = [float(value) for value in line.split()]
            if len(values) != columns:
                raise ValueError(
                    f"{path} row {line_number} has {len(values)} columns; "
                    f"expected {columns}")
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path} row {line_number} is not finite")
            matrix.append(values)
    if len(matrix) != rows:
        raise ValueError(f"{path} has {len(matrix)} rows; expected {rows}")
    return matrix


def periodic_overlap_average(values: list[float], target_bins: int) -> list[float]:
    """Conservatively average uniform source bins into uniform target bins."""
    source_bins = len(values)
    if source_bins == 0 or target_bins <= 0:
        raise ValueError("phase-bin counts must be positive")
    result = []
    for target in range(target_bins):
        low = target / target_bins
        high = (target + 1) / target_bins
        weighted = 0.0
        for source in range(source_bins):
            source_low = source / source_bins
            source_high = (source + 1) / source_bins
            overlap = max(0.0, min(high, source_high) - max(low, source_low))
            if overlap:
                weighted += overlap * values[source]
        result.append(weighted * target_bins)
    return result


def resample_matrix(matrix: list[list[float]],
                    target_bins: int) -> list[list[float]]:
    return [periodic_overlap_average(row, target_bins) for row in matrix]


def flatten_phase_major(matrix: list[list[float]]) -> list[float]:
    """Convert a space-major matrix into candidate-compatible phase-major order."""
    return [matrix[node][phase] for phase in range(len(matrix[0]))
            for node in range(len(matrix))]


def metrics(candidate: list[float], reference: list[float]) -> dict[str, float | None]:
    if len(candidate) != len(reference) or not candidate:
        raise ValueError("metric vectors must have equal positive length")
    if not all(math.isfinite(value) for value in candidate + reference):
        raise ValueError("metric vectors must be finite")
    reference_norm = math.fsum(value * value for value in reference)
    if reference_norm == 0.0:
        raise ValueError("reference metric norm is zero")
    relative_l2 = math.sqrt(
        math.fsum((a - b) ** 2 for a, b in zip(candidate, reference)) /
        reference_norm)
    rms = math.sqrt(reference_norm / len(reference))
    mean_bias_over_reference_rms = (
        math.fsum(a - b for a, b in zip(candidate, reference)) /
        len(reference) / rms)
    candidate_mean = math.fsum(candidate) / len(candidate)
    reference_mean = math.fsum(reference) / len(reference)
    covariance = math.fsum((a - candidate_mean) * (b - reference_mean)
                           for a, b in zip(candidate, reference))
    variance_a = math.fsum((a - candidate_mean) ** 2 for a in candidate)
    variance_b = math.fsum((b - reference_mean) ** 2 for b in reference)
    correlation = (covariance / math.sqrt(variance_a * variance_b)
                   if variance_a > 0.0 and variance_b > 0.0 else None)
    reference_peak = max(abs(value) for value in reference)
    return {
        "relative_l2": relative_l2,
        "mean_bias_over_reference_rms": mean_bias_over_reference_rms,
        "pearson_correlation": correlation,
        "candidate_peak_absolute": max(abs(value) for value in candidate),
        "reference_peak_absolute": reference_peak,
        "candidate_to_reference_peak_absolute_ratio":
            max(abs(value) for value in candidate) / reference_peak,
    }


def phase_space_metrics(candidate: list[float], reference: list[float],
                        phases: int = CANDIDATE_PHASES,
                        nodes: int = NODES) -> dict[str, object]:
    result: dict[str, object] = metrics(candidate, reference)
    if len(candidate) != phases * nodes:
        raise ValueError("phase-space vector size differs from the comparison grid")
    phase_l2 = []
    for phase in range(phases):
        start = phase * nodes
        stop = start + nodes
        phase_l2.append(metrics(candidate[start:stop], reference[start:stop])[
            "relative_l2"])
    worst = max(range(phases), key=lambda index: phase_l2[index])
    candidate_average = [
        math.fsum(candidate[phase * nodes + node] for phase in range(phases)) /
        phases for node in range(nodes)]
    reference_average = [
        math.fsum(reference[phase * nodes + node] for phase in range(phases)) /
        phases for node in range(nodes)]
    result.update({
        "phase_profile_relative_l2": phase_l2,
        "maximum_phase_profile_relative_l2": phase_l2[worst],
        "maximum_phase_profile_relative_l2_bin": worst,
        "cycle_average_spatial_profile_relative_l2": metrics(
            candidate_average, reference_average)["relative_l2"],
    })
    return result


def spatial_phase_average(values: list[float], phases: int = CANDIDATE_PHASES,
                          nodes: int = NODES) -> float:
    if len(values) != phases * nodes or nodes < 2:
        raise ValueError("spatial-phase average vector has the wrong size")
    weighted = 0.0
    for phase in range(phases):
        for node in range(nodes):
            weight = 0.5 if node in (0, nodes - 1) else 1.0
            weighted += weight * values[phase * nodes + node]
    return weighted / (phases * (nodes - 1))


def phase_effective_frequency(rate: list[float], density: list[float],
                              phases: int = CANDIDATE_PHASES,
                              nodes: int = NODES) -> list[float]:
    if len(rate) != phases * nodes or len(density) != len(rate):
        raise ValueError("effective-frequency vectors have the wrong size")
    result = []
    for phase in range(phases):
        start = phase * nodes
        weights = [0.5] + [1.0] * (nodes - 2) + [0.5]
        numerator = math.fsum(
            weight * value for weight, value in
            zip(weights, rate[start:start + nodes]))
        denominator = math.fsum(
            weight * value for weight, value in
            zip(weights, density[start:start + nodes]))
        if denominator <= 0.0:
            raise ValueError("phase-integrated electron density is not positive")
        result.append(numerator / denominator)
    return result


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"empty candidate table: {path}")
    return result


def validate_candidate_grid(rows: list[dict[str, str]], expected_rows: int,
                            species: str | None = None,
                            count_field: str = "samples") -> None:
    if len(rows) != expected_rows:
        raise ValueError(f"candidate has {len(rows)} rows; expected {expected_rows}")
    samples = set()
    for index, row in enumerate(rows):
        phase, node = divmod(index, NODES)
        if int(row["phase_bin"]) != phase or int(row["node"]) != node:
            raise ValueError("candidate rows are not ordered phase-major/node-minor")
        if abs(float(row["phase_fraction"]) - (phase + 0.5) / CANDIDATE_PHASES) > 1e-12:
            raise ValueError("candidate phase centers differ from the comparison contract")
        if abs(float(row["x_m"]) - node * LENGTH_M / (NODES - 1)) > 1e-12:
            raise ValueError("candidate spatial coordinates differ from the comparison contract")
        count = int(row[count_field])
        if count <= 0:
            raise ValueError("candidate sample counts must be positive")
        samples.add(count)
        if species is not None and row["species"] != species:
            raise ValueError(f"unexpected candidate species: {row['species']}")
    if len(samples) != 1:
        raise ValueError("candidate phase-space sample counts are not uniform")


def analyze(candidate: Path, reference: Path) -> dict[str, object]:
    report_path = candidate.parent / "measurement-report.json"
    pilot = json.loads(report_path.read_text(encoding="utf-8"))
    if (pilot.get("scope") != "fresh_window_aurorapic_measurement_pilot" or
            pilot.get("all_gates_passed") is not True or
            pilot.get("window", {}).get("equilibration_statistics_excluded") is not True):
        raise ValueError("candidate is not a passing fresh-window pilot")

    reference_matrices = {}
    reference_hashes = {}
    for observable, (filename, expected_hash) in REFERENCE_FILES.items():
        path = reference / filename
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"locked eduPIC matrix differs: {filename}")
        reference_hashes[filename] = actual_hash
        reference_matrices[observable] = resample_matrix(read_matrix(path),
                                                          CANDIDATE_PHASES)

    fields_path = candidate / "spatial_phase_fields.csv"
    moments_path = candidate / "spatial_phase_moments.csv"
    reported_hashes = pilot.get("output_hashes", {})
    for name, path in (("spatial_phase_fields.csv", fields_path),
                       ("spatial_phase_moments.csv", moments_path)):
        if name in reported_hashes and sha256(path) != reported_hashes[name]:
            raise ValueError(
                f"candidate {name} differs from its measurement report")
    fields = csv_rows(fields_path)
    validate_candidate_grid(fields, CANDIDATE_PHASES * NODES)
    all_moments = csv_rows(moments_path)
    moments = {}
    for species in ("electrons", "ions"):
        selected = [row for row in all_moments if row["species"] == species]
        validate_candidate_grid(selected, CANDIDATE_PHASES * NODES, species)
        moments[species] = selected
    if len(all_moments) != 2 * CANDIDATE_PHASES * NODES:
        raise ValueError("candidate moment table has unexpected species rows")

    candidate_values = {
        "potential": [float(row["potential_mean_V"]) for row in fields],
        "electric_field": [float(row["electric_field_mean_V_m"]) for row in fields],
        "electron_density": [float(row["number_density_mean_m-3"])
                             for row in moments["electrons"]],
        "ion_density": [float(row["number_density_mean_m-3"])
                        for row in moments["ions"]],
        "electron_mean_energy": [float(row["mean_kinetic_energy_eV"])
                                 for row in moments["electrons"]],
        "ion_mean_energy": [float(row["mean_kinetic_energy_eV"])
                            for row in moments["ions"]],
    }
    candidate_values["electron_current_density"] = [
        -ELEMENTARY_CHARGE_C * density * float(row["mean_velocity_x"])
        for density, row in zip(candidate_values["electron_density"],
                                moments["electrons"])]
    candidate_values["ion_current_density"] = [
        ELEMENTARY_CHARGE_C * density * float(row["mean_velocity_x"])
        for density, row in zip(candidate_values["ion_density"], moments["ions"])]
    candidate_values["electron_ohmic_power_density"] = [
        current * field for current, field in zip(
            candidate_values["electron_current_density"],
            candidate_values["electric_field"])]
    candidate_values["ion_ohmic_power_density"] = [
        current * field for current, field in zip(
            candidate_values["ion_current_density"],
            candidate_values["electric_field"])]

    comparisons = {
        name: phase_space_metrics(
            values, flatten_phase_major(reference_matrices[name]))
        for name, values in candidate_values.items()
    }
    comparisons["powered_electrode_potential"] = metrics(
        [candidate_values["potential"][phase * NODES]
         for phase in range(CANDIDATE_PHASES)],
        reference_matrices["potential"][0])

    rate_path = candidate / "spatial_phase_collision_rate.csv"
    unavailable_comparisons = {}
    derived_diagnostics = {}
    candidate_hashes = {
        "measurement_report_sha256": sha256(report_path),
        "spatial_phase_fields_sha256": sha256(fields_path),
        "spatial_phase_moments_sha256": sha256(moments_path),
    }
    if rate_path.exists():
        if ("spatial_phase_collision_rate.csv" in reported_hashes and
                sha256(rate_path) !=
                reported_hashes["spatial_phase_collision_rate.csv"]):
            raise ValueError(
                "candidate collision rate differs from its measurement report")
        rate_rows = csv_rows(rate_path)
        ionization_rows = [
            row for row in rate_rows
            if row["channel"] == "ionization" or
            row["channel"].endswith(".ionization")]
        if not ionization_rows:
            raise ValueError(
                "candidate collision-rate output has no ionization channel")
        rate = [0.0] * (CANDIDATE_PHASES * NODES)
        seen_channels = sorted({row["channel"] for row in ionization_rows})
        for channel in seen_channels:
            selected = [row for row in ionization_rows
                        if row["channel"] == channel]
            validate_candidate_grid(
                selected, CANDIDATE_PHASES * NODES,
                count_field="timesteps")
            for index, row in enumerate(selected):
                rate[index] += float(row["mean_event_rate_m-3_s-1"])
        comparisons["ionization_rate"] = phase_space_metrics(
            rate, flatten_phase_major(reference_matrices["ionization_rate"]))
        reference_rate = flatten_phase_major(
            reference_matrices["ionization_rate"])
        reference_density = flatten_phase_major(
            reference_matrices["electron_density"])
        candidate_rate_average = spatial_phase_average(rate)
        reference_rate_average = spatial_phase_average(reference_rate)
        candidate_density_average = spatial_phase_average(
            candidate_values["electron_density"])
        reference_density_average = spatial_phase_average(reference_density)
        candidate_frequency = candidate_rate_average / candidate_density_average
        reference_frequency = reference_rate_average / reference_density_average
        candidate_phase_frequency = phase_effective_frequency(
            rate, candidate_values["electron_density"])
        reference_phase_frequency = phase_effective_frequency(
            reference_rate, reference_density)
        derived_diagnostics["ionization_per_electron"] = {
            "candidate_volume_phase_average_rate_m-3_s-1":
                candidate_rate_average,
            "reference_volume_phase_average_rate_m-3_s-1":
                reference_rate_average,
            "candidate_to_reference_average_rate_ratio":
                candidate_rate_average / reference_rate_average,
            "candidate_volume_phase_average_electron_density_m-3":
                candidate_density_average,
            "reference_volume_phase_average_electron_density_m-3":
                reference_density_average,
            "candidate_to_reference_average_electron_density_ratio":
                candidate_density_average / reference_density_average,
            "candidate_effective_ionization_frequency_s-1": candidate_frequency,
            "reference_effective_ionization_frequency_s-1": reference_frequency,
            "candidate_to_reference_effective_ionization_frequency_ratio":
                candidate_frequency / reference_frequency,
            "candidate_to_reference_phase_effective_frequency_ratio": [
                candidate_value / reference_value
                for candidate_value, reference_value in
                zip(candidate_phase_frequency, reference_phase_frequency)],
        }
        candidate_hashes["spatial_phase_collision_rate_sha256"] = sha256(
            rate_path)
    else:
        unavailable_comparisons["ionization_rate"] = (
            "This checkpoint predates AuroraPIC's phase-resolved volumetric "
            "collision-event-rate diagnostic; rerun the fresh measurement window.")

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "descriptive_transient_position_rf_phase_cross_code_comparison",
        "comparison_contract": {
            "spatial_nodes": NODES,
            "reference_phase_bins": REFERENCE_PHASES,
            "candidate_phase_bins": CANDIDATE_PHASES,
            "gap_length_m": LENGTH_M,
            "phase_mapping": "periodic_interval_overlap_average",
            "orientation": "direct_powered_left_no_reflection",
            "phase_alignment": "direct_no_fitted_shift",
            "current_derivation": "q_number_density_mean_velocity_x",
            "ohmic_power_derivation": "current_density_times_electric_field",
        },
        "candidate": candidate_hashes,
        "candidate_measurement_context": {
            "window": pilot.get("window"),
            "resources": pilot.get("resources"),
            "inputs": pilot.get("inputs"),
            "all_gates_passed": pilot.get("all_gates_passed"),
        },
        "reference": {
            "repository": "https://github.com/donkozoltan/eduPIC",
            "commit": "32050728c961a317d6d6acd6bc86d026da403326",
            "implementation_file": "C/eduPIC.cc",
            "implementation_sha256_at_commit": EDUPIC_IMPLEMENTATION_SHA256,
            "raw_matrix_sha256": reference_hashes,
        },
        "comparisons": comparisons,
        "derived_diagnostics": derived_diagnostics,
        "unavailable_comparisons": unavailable_comparisons,
        "acceptance": {"thresholds_declared": False, "passes": None},
        "candidate_state_boundary": (
            "AuroraPIC contributes a four-cycle fresh measurement window from "
            "a transient, non-stationary source-loss state."),
        "reference_boundary": (
            "The locked matrices are public eduPIC reference-case outputs; this "
            "comparison does not treat them as experimental measurements."),
        "claim_boundary": (
            "The result localizes phase-space discrepancies. It does not establish "
            "converged cross-code agreement, physical validation, or predictive accuracy."),
        "physics_claim": "none_descriptive_transient_cross_code_comparison_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.candidate_output.resolve(),
                     args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
