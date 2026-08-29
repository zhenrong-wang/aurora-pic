#!/usr/bin/env python3
"""Audit consecutive post-benchmark Turner ion-density blocks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


STEPS_PER_RF_CYCLE = {1: 400, 2: 800, 3: 1600, 4: 3200}
STATIONARITY_TARGET_BLOCKS = 16
MINIMUM_EFFECTIVE_BLOCKS = 8.0
MAXIMUM_ABSOLUTE_PROJECTED_DRIFT = 0.01
MAXIMUM_ABSOLUTE_SPLIT_HALF_CHANGE = 0.01
MAXIMUM_ADJACENT_PROFILE_RELATIVE_L2 = 0.025


class BlockError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BlockError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BlockError(f"cannot read {description} {path}: {error}") from error
    require(isinstance(value, dict), f"{description} {path} must be an object")
    return value


def checked_path(raw: object, expected_hash: object,
                 description: str) -> Path:
    require(isinstance(raw, str) and raw, f"{description} path is missing")
    require(
        isinstance(expected_hash, str) and len(expected_hash) == 64,
        f"{description} SHA-256 is missing",
    )
    path = Path(raw)
    require(path.is_file(), f"{description} does not exist: {path}")
    require(
        sha256(path) == expected_hash,
        f"{description} SHA-256 differs from its comparison report",
    )
    return path


def load_profile(path: Path, species: str) -> tuple[list[float], list[float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, f"{path} has no CSV header")
            fields = set(reader.fieldnames)
            generic = {
                "species", "x_m", "number_density_mean_m-3"
            }.issubset(fields)
            direct = {"x_m", "ion_density_mean_m-3"}.issubset(fields)
            require(generic or direct, f"{path} has no supported density columns")
            coordinates: list[float] = []
            densities: list[float] = []
            for line, row in enumerate(reader, 2):
                if generic and row["species"] != species:
                    continue
                try:
                    coordinate = float(row["x_m"])
                    density = float(
                        row["number_density_mean_m-3"]
                        if generic else row["ion_density_mean_m-3"]
                    )
                except (TypeError, ValueError) as error:
                    raise BlockError(f"{path}:{line}: invalid numeric value") from error
                require(
                    math.isfinite(coordinate)
                    and math.isfinite(density)
                    and density >= 0.0,
                    f"{path}:{line}: invalid density sample",
                )
                coordinates.append(coordinate)
                densities.append(density)
    except OSError as error:
        raise BlockError(f"cannot read density profile {path}: {error}") from error
    require(len(coordinates) >= 3, f"{path} has too few {species!r} samples")
    require(
        all(right > left for left, right in zip(coordinates, coordinates[1:])),
        f"{path} coordinates are not strictly increasing",
    )
    return coordinates, densities


def trapezoid(coordinates: list[float], values: list[float]) -> float:
    return math.fsum(
        0.5 * (left_value + right_value) * (right_x - left_x)
        for left_x, right_x, left_value, right_value in zip(
            coordinates, coordinates[1:], values, values[1:]
        )
    )


def relative_l2(left: list[float], right: list[float]) -> float:
    denominator = math.fsum(value * value for value in left)
    require(denominator > 0.0, "density profile has zero norm")
    return math.sqrt(
        math.fsum(
            (right_value - left_value) ** 2
            for left_value, right_value in zip(left, right)
        )
        / denominator
    )


def root_mean_square(values: list[float]) -> float | None:
    if not values:
        return None
    return math.sqrt(math.fsum(value * value for value in values) / len(values))


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def linear_slope(values: list[float]) -> float:
    center = 0.5 * (len(values) - 1)
    denominator = math.fsum((index - center) ** 2 for index in range(len(values)))
    if denominator == 0.0:
        return 0.0
    mean = math.fsum(values) / len(values)
    return math.fsum(
        (index - center) * (value - mean)
        for index, value in enumerate(values)
    ) / denominator


def lag_one(values: list[float]) -> float | None:
    if len(values) < 3:
        return None
    mean = math.fsum(values) / len(values)
    denominator = math.fsum((value - mean) ** 2 for value in values)
    if denominator == 0.0:
        return None
    return math.fsum(
        (left - mean) * (right - mean)
        for left, right in zip(values, values[1:])
    ) / denominator


def analyze(report_paths: list[Path], minimum_blocks: int) -> dict[str, object]:
    require(minimum_blocks >= 4, "minimum blocks must be at least four")
    require(report_paths, "at least one comparison report is required")

    blocks: list[dict[str, object]] = []
    coordinates: list[float] | None = None
    profiles: list[list[float]] = []
    case: int | None = None
    species: str | None = None
    previous_end: int | None = None

    for index, report_path in enumerate(report_paths):
        report = load_json(report_path, "comparison report")
        require(
            report.get("turner_comparison_version") == 1,
            f"comparison report {index + 1} has the wrong version",
        )
        require(
            report.get("comparison_scope")
            == "post_benchmark_density_diagnostic_only",
            f"comparison report {index + 1} is not a post-benchmark block",
        )
        report_case = report.get("case")
        require(
            isinstance(report_case, int) and report_case in range(1, 5),
            f"comparison report {index + 1} has an invalid Turner case",
        )
        case = report_case if case is None else case
        require(report_case == case, "comparison reports mix Turner cases")

        candidate = report.get("candidate")
        require(isinstance(candidate, dict), "comparison candidate is missing")
        report_species = candidate.get("species")
        require(
            isinstance(report_species, str) and report_species,
            "comparison candidate species is missing",
        )
        species = report_species if species is None else species
        require(report_species == species, "comparison reports mix species")
        profile_path = checked_path(
            candidate.get("path"), candidate.get("sha256"),
            f"block {index + 1} density profile",
        )
        metadata_path = checked_path(
            candidate.get("averaging_metadata_path"),
            candidate.get("averaging_metadata_sha256"),
            f"block {index + 1} averaging metadata",
        )
        metadata = load_json(metadata_path, "averaging metadata")
        start = metadata.get("start_step")
        end = metadata.get("end_step")
        samples = metadata.get("samples")
        required_samples = 32 * STEPS_PER_RF_CYCLE[report_case]
        require(
            metadata.get("complete") is True
            and metadata.get("reset_on_restart") is True
            and metadata.get("interval") == 1
            and isinstance(start, int)
            and isinstance(end, int)
            and isinstance(samples, int)
            and samples == required_samples
            and metadata.get("expected_samples") == required_samples
            and end - start + 1 == samples
            and metadata.get("final_step") == end
            and metadata.get("rf_cycles") == 32,
            f"block {index + 1} is not a complete reset 32-cycle window",
        )
        if previous_end is not None:
            require(
                start == previous_end + 1,
                f"block {index + 1} is not contiguous with the previous block",
            )
        previous_end = end

        block_coordinates, density = load_profile(profile_path, report_species)
        if coordinates is None:
            coordinates = block_coordinates
        else:
            require(
                len(block_coordinates) == len(coordinates)
                and all(
                    math.isclose(left, right, rel_tol=0.0, abs_tol=1e-15)
                    for left, right in zip(block_coordinates, coordinates)
                ),
                f"block {index + 1} uses a different spatial grid",
            )
        profiles.append(density)
        statistic = report.get("statistic")
        secondary = report.get("secondary_metrics")
        x_squared = (
            statistic.get("x_squared") if isinstance(statistic, dict) else None
        )
        reference_l2 = (
            secondary.get("relative_l2") if isinstance(secondary, dict) else None
        )
        require(
            isinstance(statistic, dict)
            and statistic.get("published_acceptance_applicable") is False
            and isinstance(x_squared, (int, float))
            and math.isfinite(float(x_squared))
            and float(x_squared) >= 0.0
            and isinstance(secondary, dict)
            and isinstance(reference_l2, (int, float))
            and math.isfinite(float(reference_l2))
            and float(reference_l2) >= 0.0,
            f"comparison report {index + 1} lacks comparison metrics",
        )
        blocks.append({
            "comparison_report": str(report_path.resolve()),
            "comparison_report_sha256": sha256(report_path),
            "start_step": start,
            "end_step": end,
            "x_squared": float(x_squared),
            "reference_relative_l2": float(reference_l2),
            "line_integrated_density_m-2":
                trapezoid(block_coordinates, density),
            "midplane_density_m-3": density[len(density) // 2],
        })

    assert coordinates is not None
    integrals = [
        float(block["line_integrated_density_m-2"]) for block in blocks
    ]
    mean_integral = math.fsum(integrals) / len(integrals)
    require(mean_integral > 0.0, "mean line-integrated density is zero")
    slope = linear_slope(integrals)
    correlation = lag_one(integrals)
    effective_blocks: float | None = None
    if correlation is not None:
        bounded = min(0.99, max(-0.99, correlation))
        effective_blocks = min(
            float(len(integrals)),
            max(1.0, len(integrals) * (1.0 - bounded) / (1.0 + bounded)),
        )

    movements = [
        relative_l2(left, right)
        for left, right in zip(profiles, profiles[1:])
    ]
    normalized_profiles = [
        [value / integral for value in profile]
        for profile, integral in zip(profiles, integrals)
    ]
    normalized_movements = [
        relative_l2(left, right)
        for left, right in zip(normalized_profiles, normalized_profiles[1:])
    ]
    changes = [
        (right - left) / left
        for left, right in zip(integrals, integrals[1:])
    ]
    half = len(integrals) // 2
    split_half_change: float | None = None
    if half > 0:
        first_mean = math.fsum(integrals[:half]) / half
        second_mean = math.fsum(integrals[-half:]) / half
        split_half_change = (second_mean - first_mean) / mean_integral
    projected_drift = slope * (len(integrals) - 1) / mean_integral
    maximum_movement = max(movements) if movements else None
    enough = len(blocks) >= minimum_blocks
    stationarity_horizon_complete = (
        len(blocks) >= STATIONARITY_TARGET_BLOCKS
    )
    stationarity_gates = {
        "minimum_total_blocks": {
            "threshold": STATIONARITY_TARGET_BLOCKS,
            "value": len(blocks),
            "passed": stationarity_horizon_complete,
        },
        "minimum_ar1_effective_blocks": {
            "threshold": MINIMUM_EFFECTIVE_BLOCKS,
            "value": effective_blocks,
            "passed": (
                effective_blocks is not None
                and effective_blocks >= MINIMUM_EFFECTIVE_BLOCKS
            ),
        },
        "maximum_absolute_projected_fractional_drift": {
            "threshold": MAXIMUM_ABSOLUTE_PROJECTED_DRIFT,
            "value": abs(projected_drift),
            "passed": abs(projected_drift)
                <= MAXIMUM_ABSOLUTE_PROJECTED_DRIFT,
        },
        "maximum_absolute_split_half_fractional_change": {
            "threshold": MAXIMUM_ABSOLUTE_SPLIT_HALF_CHANGE,
            "value": (
                abs(split_half_change)
                if split_half_change is not None else None
            ),
            "passed": (
                split_half_change is not None
                and abs(split_half_change)
                    <= MAXIMUM_ABSOLUTE_SPLIT_HALF_CHANGE
            ),
        },
        "maximum_adjacent_profile_relative_l2": {
            "threshold": MAXIMUM_ADJACENT_PROFILE_RELATIVE_L2,
            "value": maximum_movement,
            "passed": (
                maximum_movement is not None
                and maximum_movement
                    <= MAXIMUM_ADJACENT_PROFILE_RELATIVE_L2
            ),
        },
    }
    stationarity_passed = (
        stationarity_horizon_complete
        and all(
            bool(gate["passed"]) for gate in stationarity_gates.values()
        )
    )
    if not enough:
        classification = "insufficient_consecutive_blocks"
    elif not stationarity_horizon_complete:
        classification = "stationarity_horizon_incomplete"
    elif stationarity_passed:
        classification = "internal_stationarity_screen_passed"
    else:
        classification = "internal_stationarity_screen_failed"
    return {
        "turner_density_block_analysis_version": 3,
        "case": case,
        "species": species,
        "block_contract": {
            "required_rf_cycles_per_block": 32,
            "contiguous": True,
            "reset_on_restart": True,
            "minimum_blocks_for_diagnostic_series": minimum_blocks,
            "stationarity_target_blocks": STATIONARITY_TARGET_BLOCKS,
        },
        "blocks": blocks,
        "adjacent_profile_relative_l2": movements,
        "adjacent_integral_normalized_profile_relative_l2":
            normalized_movements,
        "adjacent_integrated_density_fractional_change": changes,
        "series_metrics": {
            "mean_line_integrated_density_m-2": mean_integral,
            "linear_slope_per_block_m-2": slope,
            "projected_fractional_drift_across_series": projected_drift,
            "split_half_fractional_change": split_half_change,
            "lag_one_integrated_density_correlation": correlation,
            "ar1_effective_blocks": effective_blocks,
            "maximum_adjacent_profile_relative_l2": maximum_movement,
            "rms_adjacent_profile_relative_l2":
                root_mean_square(movements),
            "p95_adjacent_profile_relative_l2": quantile(movements, 0.95),
            "maximum_adjacent_integral_normalized_profile_relative_l2":
                max(normalized_movements) if normalized_movements else None,
            "rms_adjacent_integral_normalized_profile_relative_l2":
                root_mean_square(normalized_movements),
            "p95_adjacent_integral_normalized_profile_relative_l2":
                quantile(normalized_movements, 0.95),
        },
        "diagnostic_series_ready": enough,
        "stationarity_screen": {
            "scope": "internal_post_benchmark_diagnostic_only",
            "threshold_basis":
                "predeclared_before_blocks_9_through_16; one-percent "
                "density-equivalence scale follows the approximately "
                "one-percent median Turner reference population scatter",
            "horizon_complete": stationarity_horizon_complete,
            "gates": stationarity_gates,
            "passed": stationarity_passed
                if stationarity_horizon_complete else None,
        },
        "classification": classification,
        "published_acceptance_applicable": False,
        "physics_claim": "none_post_benchmark_diagnostic_only",
    }


def write_report(path: Path, report: dict[str, object]) -> None:
    require(not path.exists(), f"refusing to overwrite existing report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--minimum-blocks", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = analyze(args.reports, args.minimum_blocks)
        write_report(args.output, report)
    except (BlockError, OSError) as error:
        print(f"Turner density-block error: {error}", file=sys.stderr)
        return 2
    print(
        "Turner density-block audit written to "
        f"{args.output} ({report['classification']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
