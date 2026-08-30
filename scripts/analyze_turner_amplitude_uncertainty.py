#!/usr/bin/env python3
"""Explore Turner density-amplitude drift under a fitted stationary AR(1) null."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import sys
import tempfile


class UncertaintyError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise UncertaintyError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_analysis(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise UncertaintyError(f"cannot read block analysis {path}: {error}") from error
    require(isinstance(value, dict), "block analysis must be an object")
    require(
        value.get("turner_density_block_analysis_version") == 3,
        "block analysis must use Turner density-block analysis version 3",
    )
    require(
        value.get("published_acceptance_applicable") is False,
        "block analysis must retain its non-acceptance claim boundary",
    )
    return value


def mean(values: list[float]) -> float:
    require(bool(values), "cannot average an empty series")
    return math.fsum(values) / len(values)


def linear_slope(values: list[float]) -> float:
    center = 0.5 * (len(values) - 1)
    denominator = math.fsum((index - center) ** 2 for index in range(len(values)))
    require(denominator > 0.0, "at least two values are required for a slope")
    average = mean(values)
    return math.fsum(
        (index - center) * (value - average)
        for index, value in enumerate(values)
    ) / denominator


def projected_fractional_drift(values: list[float]) -> float:
    average = mean(values)
    require(average > 0.0, "density series mean must be positive")
    return linear_slope(values) * (len(values) - 1) / average


def lag(values: list[float], offset: int) -> float | None:
    require(offset >= 1, "lag offset must be positive")
    if len(values) <= offset + 1:
        return None
    average = mean(values)
    denominator = math.fsum((value - average) ** 2 for value in values)
    if denominator == 0.0:
        return None
    return math.fsum(
        (values[index] - average) * (values[index + offset] - average)
        for index in range(len(values) - offset)
    ) / denominator


def quantile(values: list[float], probability: float) -> float:
    require(bool(values), "cannot take a quantile of an empty series")
    require(0.0 <= probability <= 1.0, "quantile probability is invalid")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def fit_stationary_ar1(values: list[float]) -> dict[str, object]:
    average = mean(values)
    centered = [value / average - 1.0 for value in values]
    trend_slope = linear_slope(centered)
    trend_center = 0.5 * (len(centered) - 1)
    detrended = [
        value - trend_slope * (index - trend_center)
        for index, value in enumerate(centered)
    ]
    denominator = math.fsum(value * value for value in detrended[:-1])
    require(denominator > 0.0, "density series has no amplitude variance")
    unconstrained_phi = math.fsum(
        left * right for left, right in zip(detrended, detrended[1:])
    ) / denominator
    phi = min(0.95, max(-0.95, unconstrained_phi))
    innovations = [
        right - phi * left for left, right in zip(detrended, detrended[1:])
    ]
    innovation_mean = mean(innovations)
    residuals = [value - innovation_mean for value in innovations]
    require(len(residuals) >= 2, "too few AR(1) residuals")
    variance = math.fsum(value * value for value in residuals) / (
        len(residuals) - 1
    )
    require(variance > 0.0 and math.isfinite(variance), "invalid AR(1) variance")
    innovation_sigma = math.sqrt(variance)
    stationary_sigma = innovation_sigma / math.sqrt(1.0 - phi * phi)
    return {
        "mean_line_integrated_density_m-2": average,
        "fit_series": "linear_detrended_fractional_amplitude_residuals",
        "removed_projected_fractional_drift": trend_slope * (
            len(centered) - 1
        ),
        "unconstrained_phi": unconstrained_phi,
        "phi": phi,
        "phi_clipped_to_stationary_bound": phi != unconstrained_phi,
        "innovation_fractional_mean_before_centering": innovation_mean,
        "innovation_fractional_standard_deviation": innovation_sigma,
        "stationary_fractional_standard_deviation": stationary_sigma,
        "residual_lag_one_correlation": lag(residuals, 1),
        "residual_lag_two_correlation": lag(residuals, 2),
    }


def analyze(
    analysis_path: Path,
    replicates: int,
    random_seed: int,
    alpha: float,
) -> dict[str, object]:
    require(replicates >= 1000, "at least 1000 null replicates are required")
    require(0.0 < alpha < 0.5, "alpha must lie between zero and one half")
    source = load_analysis(analysis_path)
    blocks = source.get("blocks")
    require(isinstance(blocks, list) and len(blocks) >= 16,
            "at least 16 consecutive blocks are required")
    values: list[float] = []
    for index, block in enumerate(blocks, 1):
        require(isinstance(block, dict), f"block {index} must be an object")
        value = block.get("line_integrated_density_m-2")
        require(
            isinstance(value, (int, float))
            and math.isfinite(float(value))
            and float(value) > 0.0,
            f"block {index} has invalid integrated density",
        )
        values.append(float(value))

    fit = fit_stationary_ar1(values)
    phi = float(fit["phi"])
    innovation_sigma = float(fit["innovation_fractional_standard_deviation"])
    stationary_sigma = float(fit["stationary_fractional_standard_deviation"])
    observed_drift = projected_fractional_drift(values)
    observed_absolute = abs(observed_drift)
    generator = random.Random(random_seed)
    null_absolute_drifts: list[float] = []
    for _ in range(replicates):
        state = generator.gauss(0.0, stationary_sigma)
        simulated = [1.0 + state]
        for _ in range(1, len(values)):
            state = phi * state + generator.gauss(0.0, innovation_sigma)
            simulated.append(1.0 + state)
        null_absolute_drifts.append(abs(projected_fractional_drift(simulated)))

    exceedances = sum(
        value >= observed_absolute for value in null_absolute_drifts
    )
    p_value = (exceedances + 1.0) / (replicates + 1.0)
    monte_carlo_standard_error = math.sqrt(
        p_value * (1.0 - p_value) / (replicates + 1.0)
    )
    rejected = p_value < alpha
    return {
        "turner_amplitude_uncertainty_version": 1,
        "scope": "post_protocol_exploratory_stationary_null_diagnostic",
        "source": {
            "path": str(analysis_path.resolve()),
            "sha256": sha256(analysis_path),
            "blocks": len(values),
            "case": source.get("case"),
            "species": source.get("species"),
        },
        "observed": {
            "projected_fractional_drift": observed_drift,
            "absolute_projected_fractional_drift": observed_absolute,
            "exceeds_historical_one_percent_scale": observed_absolute > 0.01,
        },
        "stationary_ar1_fit": fit,
        "parametric_stationary_null": {
            "replicates": replicates,
            "random_seed": random_seed,
            "two_sided_alpha": alpha,
            "absolute_drift_quantile_at_one_minus_alpha": quantile(
                null_absolute_drifts, 1.0 - alpha
            ),
            "absolute_drift_median": quantile(null_absolute_drifts, 0.5),
            "exceedances": exceedances,
            "p_value": p_value,
            "monte_carlo_standard_error": monte_carlo_standard_error,
            "stationary_null_rejected": rejected,
        },
        "classification": (
            "exploratory_stationary_ar1_null_rejected"
            if rejected
            else "exploratory_stationary_ar1_null_not_rejected"
        ),
        "model_boundary": (
            "The null assumes Gaussian stationary AR(1) errors fitted to the "
            "linear-detrended amplitude series. Residual correlation must be "
            "inspected; failure to reject is not proof of stationarity."
        ),
        "decision_boundary": (
            "This post-protocol diagnostic cannot replace, rescue, or reinterpret "
            "a preregistered gate and is not eligible for Turner acceptance."
        ),
        "published_acceptance_applicable": False,
        "physics_claim": "none_exploratory_uncertainty_diagnostic_only",
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
    parser.add_argument("analysis", type=Path)
    parser.add_argument("--replicates", type=int, default=20000)
    parser.add_argument("--random-seed", type=int, default=20260830)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = analyze(
            args.analysis, args.replicates, args.random_seed, args.alpha
        )
        write_report(args.output, report)
    except (UncertaintyError, OSError) as error:
        print(f"Turner amplitude uncertainty error: {error}", file=sys.stderr)
        return 2
    print(
        "Turner amplitude uncertainty report written to "
        f"{args.output} ({report['classification']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
