#!/usr/bin/env python3
"""Analyze immutable eduPIC measurement blocks without false pooling claims."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import statistics
import sys

from advance_edupic_measurement import inspect_stage
from run_edupic_measurement_stage import N_EEPF, N_GRID, N_IFED
from run_edupic_stage import atomic_json, checkpoint_state, sha256


class BlockAnalysisError(RuntimeError):
    pass


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def read_table(path: Path, rows: int, columns: int) -> list[list[float]]:
    result: list[list[float]] = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise BlockAnalysisError(f"cannot read {path}: {error}") from error
    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        if len(fields) != columns:
            raise BlockAnalysisError(
                f"{path.name} line {line_number} has {len(fields)} columns; "
                f"expected {columns}")
        try:
            values = [float(value) for value in fields]
        except ValueError as error:
            raise BlockAnalysisError(
                f"{path.name} line {line_number} is non-numeric") from error
        if not all(math.isfinite(value) for value in values):
            raise BlockAnalysisError(f"{path.name} contains a non-finite value")
        result.append(values)
    if len(result) != rows:
        raise BlockAnalysisError(
            f"{path.name} has {len(result)} rows; expected {rows}")
    return result


def require_common_axis(tables: list[list[list[float]]], name: str) -> list[float]:
    axis = [row[0] for row in tables[0]]
    for block, table in enumerate(tables[1:], 2):
        if [row[0] for row in table] != axis:
            raise BlockAnalysisError(
                f"{name} axis differs in measurement block {block}")
    return axis


def scalar_statistics(values: list[float]) -> dict:
    mean = statistics.fmean(values)
    sample_stddev = statistics.stdev(values) if len(values) >= 2 else None
    standard_error = (
        sample_stddev / math.sqrt(len(values))
        if sample_stddev is not None else None)
    relative_range = (
        (max(values) - min(values)) / abs(mean) if mean != 0.0 else None)
    return {
        "values_by_block": values,
        "mean": mean,
        "sample_standard_deviation": sample_stddev,
        "naive_block_standard_error": standard_error,
        "minimum": min(values),
        "maximum": max(values),
        "range_relative_to_absolute_mean": relative_range,
    }


def column_statistics(values: list[float]) -> tuple[float, float | None, float | None]:
    mean = statistics.fmean(values)
    sample_stddev = statistics.stdev(values) if len(values) >= 2 else None
    standard_error = (
        sample_stddev / math.sqrt(len(values))
        if sample_stddev is not None else None)
    return mean, sample_stddev, standard_error


def relative_l2(candidate: list[float], reference: list[float]) -> float | None:
    denominator = sum(value * value for value in reference)
    if denominator == 0.0:
        return None
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(candidate, reference)) /
                     denominator)


def trapezoid(coordinates: list[float], values: list[float]) -> float:
    if len(coordinates) != len(values) or len(coordinates) < 2:
        raise BlockAnalysisError("profile integration requires at least two points")
    return math.fsum(
        0.5 * (left_value + right_value) * (right_x - left_x)
        for left_x, right_x, left_value, right_value in zip(
            coordinates, coordinates[1:], values, values[1:]))


def linear_slope(values: list[float]) -> float:
    center = (len(values) - 1) / 2.0
    denominator = math.fsum((index - center) ** 2
                            for index in range(len(values)))
    if denominator == 0.0:
        return 0.0
    mean = statistics.fmean(values)
    return math.fsum(
        (index - center) * (value - mean)
        for index, value in enumerate(values)) / denominator


def lag_one(values: list[float]) -> float | None:
    if len(values) < 3:
        return None
    mean = statistics.fmean(values)
    centered = [value - mean for value in values]
    denominator = math.sqrt(
        math.fsum(value * value for value in centered[:-1]) *
        math.fsum(value * value for value in centered[1:]))
    if denominator == 0.0:
        return None
    return math.fsum(left * right for left, right in zip(
        centered, centered[1:])) / denominator


def ar1_effective_count(count: int, correlation: float | None) -> float | None:
    if correlation is None:
        return None
    bounded = min(0.99, max(-0.99, correlation))
    return min(float(count), max(
        1.0, count * (1.0 - bounded) / (1.0 + bounded)))


def batch_mean_diagnostics(values: list[float]) -> list[dict]:
    """Report non-overlapping power-of-two batches with at least four means."""
    result: list[dict] = []
    batch_size = 1
    while len(values) // batch_size >= 4:
        batch_count = len(values) // batch_size
        used = batch_count * batch_size
        means = [statistics.fmean(values[start:start + batch_size])
                 for start in range(0, used, batch_size)]
        mean = statistics.fmean(means)
        sample_stddev = statistics.stdev(means)
        standard_error = sample_stddev / math.sqrt(batch_count)
        correlation = lag_one(means)
        result.append({
            "batch_size_blocks": batch_size,
            "batch_count": batch_count,
            "used_blocks": used,
            "dropped_trailing_blocks": len(values) - used,
            "lag_one_batch_mean_correlation": correlation,
            "sample_standard_deviation": sample_stddev,
            "standard_error": standard_error,
            "standard_error_relative_to_absolute_mean": (
                standard_error / abs(mean) if mean != 0.0 else None),
        })
        batch_size *= 2
    return result


def profile_series_statistics(coordinates: list[float],
                              profiles: list[list[float]]) -> dict:
    integrals = [trapezoid(coordinates, profile) for profile in profiles]
    mean_integral = statistics.fmean(integrals)
    if mean_integral == 0.0:
        raise BlockAnalysisError("mean line-integrated density is zero")
    correlation = lag_one(integrals)
    effective_blocks = ar1_effective_count(len(integrals), correlation)
    sample_stddev = statistics.stdev(integrals)
    naive_standard_error = sample_stddev / math.sqrt(len(integrals))
    ar1_standard_error = (
        sample_stddev / math.sqrt(effective_blocks)
        if effective_blocks is not None else None)
    half = len(integrals) // 2
    split_half_change: float | None = None
    if half:
        split_half_change = (
            statistics.fmean(integrals[-half:]) -
            statistics.fmean(integrals[:half])) / mean_integral
    projected_drift = (
        linear_slope(integrals) * (len(integrals) - 1) / mean_integral)
    adjacent = [relative_l2(left, right) for left, right in zip(
        profiles, profiles[1:])]
    return {
        "line_integrated_density_m-2_by_block": integrals,
        "mean_line_integrated_density_m-2": mean_integral,
        "projected_fractional_drift_across_series": projected_drift,
        "split_half_fractional_change": split_half_change,
        "lag_one_integrated_density_correlation": correlation,
        "ar1_effective_blocks": effective_blocks,
        "sample_standard_deviation_m-2": sample_stddev,
        "naive_standard_error_m-2": naive_standard_error,
        "ar1_corrected_standard_error_m-2": ar1_standard_error,
        "naive_standard_error_relative_to_absolute_mean":
            naive_standard_error / abs(mean_integral),
        "ar1_corrected_standard_error_relative_to_absolute_mean": (
            ar1_standard_error / abs(mean_integral)
            if ar1_standard_error is not None else None),
        "ar1_variance_inflation_relative_to_independent_blocks": (
            len(integrals) / effective_blocks
            if effective_blocks is not None else None),
        "adjacent_profile_relative_l2": adjacent,
        "maximum_adjacent_profile_relative_l2": (
            max(value for value in adjacent if value is not None)
            if any(value is not None for value in adjacent) else None),
        "nonoverlapping_batch_mean_diagnostics":
            batch_mean_diagnostics(integrals),
    }


def atomic_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    os.replace(temporary, path)


def duration_mixture(tables: list[list[list[float]]], column: int,
                     weights: list[int]) -> list[float]:
    total = float(sum(weights))
    return [
        sum(weight * table[row][column]
            for weight, table in zip(weights, tables)) / total
        for row in range(len(tables[0]))
    ]


def load_blocks(campaign_dir: Path) -> tuple[dict, list[Path], list[dict]]:
    manifest = campaign_dir / "campaign-report.json"
    try:
        campaign = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BlockAnalysisError(f"cannot read campaign manifest: {error}") from error
    if campaign.get("scope") != "bounded_native_measurement_campaign":
        raise BlockAnalysisError("input is not a native measurement campaign")
    binary_sha256 = campaign.get("source_binary_sha256")
    initial = campaign.get("initial_state", {})
    if not isinstance(binary_sha256, str) or len(binary_sha256) != 64:
        raise BlockAnalysisError("campaign binary identity is invalid")
    current = initial
    stage_dirs: list[Path] = []
    stages: list[dict] = []
    for recorded in campaign.get("stages", []):
        if recorded.get("start_cycle") != current.get("cycles"):
            raise BlockAnalysisError("campaign stage chain is discontinuous")
        stage_dir = campaign_dir / (
            f"stage-{recorded['start_cycle']:06d}-{recorded['end_cycle']:06d}")
        try:
            stage, final = inspect_stage(stage_dir, current, binary_sha256)
        except Exception as error:
            raise BlockAnalysisError(str(error)) from error
        if sha256(stage_dir / "stage-report.json") != recorded.get(
                "stage_report_sha256"):
            raise BlockAnalysisError("campaign stage-report hash differs")
        stage_dirs.append(stage_dir)
        stages.append(stage)
        current = final
    if not stages:
        raise BlockAnalysisError("measurement campaign has no completed blocks")
    on_disk = checkpoint_state(stage_dirs[-1] / "picdata.bin")
    if campaign.get("latest_state") != on_disk:
        raise BlockAnalysisError("campaign latest state differs from checkpoint chain")
    return campaign, stage_dirs, stages


def analyze(campaign_dirs: list[Path], minimum_blocks: int,
            density_csv: Path | None, eepf_csv: Path | None,
            ifed_csv: Path | None) -> dict:
    if not campaign_dirs:
        raise BlockAnalysisError("at least one campaign is required")
    campaigns: list[dict] = []
    stage_dirs: list[Path] = []
    stages: list[dict] = []
    previous_latest: dict | None = None
    binary_sha256: str | None = None
    for index, campaign_dir in enumerate(campaign_dirs, 1):
        campaign, campaign_stage_dirs, campaign_stages = load_blocks(
            campaign_dir)
        if binary_sha256 is None:
            binary_sha256 = campaign["source_binary_sha256"]
        elif campaign.get("source_binary_sha256") != binary_sha256:
            raise BlockAnalysisError(
                f"continuation campaign {index} changes the source binary")
        if (previous_latest is not None and
                campaign.get("initial_state") != previous_latest):
            raise BlockAnalysisError(
                f"continuation campaign {index} does not begin at the prior "
                "checkpoint")
        campaigns.append(campaign)
        stage_dirs.extend(campaign_stage_dirs)
        stages.extend(campaign_stages)
        previous_latest = campaign.get("latest_state")
    block_cycles = [stage["stage"]["requested_cycles"] for stage in stages]
    equal_block_cycles = len(set(block_cycles)) == 1
    target_reached = all(campaign.get("target_reached") is True
                         for campaign in campaigns)
    eligible = (len(stages) >= minimum_blocks and equal_block_cycles and
                target_reached)

    density_tables = [read_table(path / "density.dat", N_GRID, 3)
                      for path in stage_dirs]
    x = require_common_axis(density_tables, "density.dat")
    electron_profiles = [[row[1] for row in table]
                         for table in density_tables]
    ion_profiles = [[row[2] for row in table] for table in density_tables]
    electron_density = duration_mixture(density_tables, 1, block_cycles)
    ion_density = duration_mixture(density_tables, 2, block_cycles)
    density_rows: list[list[float]] = []
    for index, coordinate in enumerate(x):
        electron_values = [table[index][1] for table in density_tables]
        ion_values = [table[index][2] for table in density_tables]
        e_mean, e_sd, e_se = column_statistics(electron_values)
        i_mean, i_sd, i_se = column_statistics(ion_values)
        density_rows.append([
            coordinate, electron_density[index], ion_density[index],
            e_sd if e_sd is not None else math.nan,
            e_se if e_se is not None else math.nan,
            i_sd if i_sd is not None else math.nan,
            i_se if i_se is not None else math.nan,
        ])

    eepf_tables = [read_table(path / "eepf.dat", N_EEPF, 2)
                   for path in stage_dirs]
    eepf_energy = require_common_axis(eepf_tables, "eepf.dat")
    eepf_mixture = duration_mixture(eepf_tables, 1, block_cycles)
    eepf_rows: list[list[float]] = []
    for index, energy in enumerate(eepf_energy):
        values = [table[index][1] for table in eepf_tables]
        _, sample_stddev, standard_error = column_statistics(values)
        eepf_rows.append([
            energy, eepf_mixture[index],
            sample_stddev if sample_stddev is not None else math.nan,
            standard_error if standard_error is not None else math.nan,
        ])
    de_eepf = eepf_energy[1] - eepf_energy[0]
    eepf_tv = [
        0.5 * sum(abs(row[1] - mixture) * math.sqrt(row[0]) * de_eepf
                  for row, mixture in zip(table, eepf_mixture))
        for table in eepf_tables
    ]

    ifed_tables = [read_table(path / "ifed.dat", N_IFED, 3)
                   for path in stage_dirs]
    ifed_energy = require_common_axis(ifed_tables, "ifed.dat")
    ifed_powered = duration_mixture(ifed_tables, 1, block_cycles)
    ifed_grounded = duration_mixture(ifed_tables, 2, block_cycles)
    ifed_rows: list[list[float]] = []
    for index, energy in enumerate(ifed_energy):
        powered_values = [table[index][1] for table in ifed_tables]
        grounded_values = [table[index][2] for table in ifed_tables]
        _, powered_sd, powered_se = column_statistics(powered_values)
        _, grounded_sd, grounded_se = column_statistics(grounded_values)
        ifed_rows.append([
            energy, ifed_powered[index], ifed_grounded[index],
            powered_sd if powered_sd is not None else math.nan,
            powered_se if powered_se is not None else math.nan,
            grounded_sd if grounded_sd is not None else math.nan,
            grounded_se if grounded_se is not None else math.nan,
        ])
    de_ifed = ifed_energy[1] - ifed_energy[0]
    ifed_powered_tv = [
        0.5 * sum(abs(row[1] - mixture) * de_ifed
                  for row, mixture in zip(table, ifed_powered))
        for table in ifed_tables
    ]
    ifed_grounded_tv = [
        0.5 * sum(abs(row[2] - mixture) * de_ifed
                  for row, mixture in zip(table, ifed_grounded))
        for table in ifed_tables
    ]

    output_hashes: dict[str, str] = {}
    if density_csv is not None:
        atomic_csv(density_csv, [
            "x_m", "electron_density_duration_mean_m3",
            "ion_density_duration_mean_m3",
            "electron_density_block_sample_stddev_m3",
            "electron_density_naive_block_standard_error_m3",
            "ion_density_block_sample_stddev_m3",
            "ion_density_naive_block_standard_error_m3",
        ], density_rows)
        output_hashes["density_csv_sha256"] = sha256(density_csv)
    if eepf_csv is not None:
        atomic_csv(eepf_csv, [
            "energy_ev", "equal_time_block_mixture_eepf_ev-1p5",
            "block_sample_stddev_ev-1p5",
            "naive_block_standard_error_ev-1p5",
        ], eepf_rows)
        output_hashes["eepf_csv_sha256"] = sha256(eepf_csv)
    if ifed_csv is not None:
        atomic_csv(ifed_csv, [
            "energy_ev", "equal_time_block_mixture_powered_ev-1",
            "equal_time_block_mixture_grounded_ev-1",
            "powered_block_sample_stddev_ev-1",
            "powered_naive_block_standard_error_ev-1",
            "grounded_block_sample_stddev_ev-1",
            "grounded_naive_block_standard_error_ev-1",
        ], ifed_rows)
        output_hashes["ifed_csv_sha256"] = sha256(ifed_csv)

    observable_names = sorted(stages[0]["reported_observables"])
    scalar_block_statistics = {
        name: scalar_statistics([
            float(stage["reported_observables"][name]) for stage in stages])
        for name in observable_names
        if name != "measurement_cycles"
    }
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "native_measurement_block_analysis",
        "physics_claim": "none",
        "claim_boundary": (
            "Block analysis quantifies short-window consistency and preserves "
            "upstream normalization boundaries. It does not establish "
            "cross-code agreement or independent physical validation."),
        "campaign": {
            "path": str(campaign_dirs[0].resolve()),
            "campaign_report_sha256": sha256(
                campaign_dirs[0] / "campaign-report.json"),
            "continuations": [
                {
                    "path": str(path.resolve()),
                    "campaign_report_sha256": sha256(
                        path / "campaign-report.json"),
                }
                for path in campaign_dirs[1:]
            ],
            "target_reached": target_reached,
            "completed_measurement_cycles":
                sum(int(campaign.get("completed_measurement_cycles", 0))
                    for campaign in campaigns),
            "block_count": len(stages),
            "block_cycles": block_cycles,
            "equal_block_cycles": equal_block_cycles,
            "first_cycle": stages[0]["stage"]["start_cycle"] + 1,
            "last_cycle": stages[-1]["stage"]["end_cycle"],
        },
        "eligibility": {
            "minimum_blocks": minimum_blocks,
            "passes_minimum_blocks": len(stages) >= minimum_blocks,
            "requires_equal_block_cycles": True,
            "passes_equal_block_cycles": equal_block_cycles,
            "requires_campaign_target_reached": True,
            "passes_campaign_target_reached": target_reached,
            "eligible": eligible,
        },
        "aggregation_contract": {
            "exact_duration_weighted": [
                "density.dat", "pot_xt.dat", "efield_xt.dat", "ne_xt.dat",
                "ni_xt.dat", "je_xt.dat", "ji_xt.dat", "ioniz_xt.dat",
            ],
            "equal_time_normalized_mixture_not_native_pooled": [
                "eepf.dat", "ifed.dat",
            ],
            "blockwise_only_missing_raw_weights_or_nonlinear": [
                "meanee_xt.dat", "meanei_xt.dat", "powere_xt.dat",
                "poweri_xt.dat", "mean ion energies", "plasma characteristics",
                "collision frequencies",
            ],
            "rounded_scalar_duration_means": [
                "electrode fluxes",
            ],
            "uncertainty_boundary": (
                "Reported standard errors treat contiguous blocks as independent "
                "and are not corrected for autocorrelation."),
        },
        "density_profile": {
            "aggregation": "exact_duration_weighted",
            "electron_block_relative_l2_to_aggregate": [
                relative_l2(profile, electron_density)
                for profile in electron_profiles],
            "ion_block_relative_l2_to_aggregate": [
                relative_l2(profile, ion_density) for profile in ion_profiles],
            "electron_series": profile_series_statistics(
                x, electron_profiles),
            "ion_series": profile_series_statistics(x, ion_profiles),
        },
        "eepf": {
            "aggregation": "equal_time_normalized_mixture_not_native_pooled",
            "mixture_weighted_normalization": sum(
                value * math.sqrt(energy) * de_eepf
                for value, energy in zip(eepf_mixture, eepf_energy)),
            "block_total_variation_to_mixture": eepf_tv,
        },
        "ifed": {
            "aggregation": "equal_time_normalized_mixture_not_native_pooled",
            "powered_mixture_normalization": sum(ifed_powered) * de_ifed,
            "grounded_mixture_normalization": sum(ifed_grounded) * de_ifed,
            "powered_block_total_variation_to_mixture": ifed_powered_tv,
            "grounded_block_total_variation_to_mixture": ifed_grounded_tv,
        },
        "scalar_block_statistics": scalar_block_statistics,
        "outputs": output_hashes,
        "analysis_eligible": eligible,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_dir", type=Path)
    parser.add_argument("--continuation-campaign-dir", type=Path,
                        action="append", default=[])
    parser.add_argument("--minimum-blocks", type=positive_integer, default=4)
    parser.add_argument("--density-csv", type=Path)
    parser.add_argument("--eepf-csv", type=Path)
    parser.add_argument("--ifed-csv", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-eligible", action="store_true")
    args = parser.parse_args()
    try:
        report = analyze([args.campaign_dir] + args.continuation_campaign_dir,
                         args.minimum_blocks,
                         args.density_csv, args.eepf_csv, args.ifed_csv)
    except (BlockAnalysisError, OSError, json.JSONDecodeError) as error:
        print(f"eduPIC measurement-block analysis failed: {error}", file=sys.stderr)
        return 2
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(args.output, report)
    print(rendered, end="")
    return 1 if args.require_eligible and not report["analysis_eligible"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
