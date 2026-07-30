#!/usr/bin/env python3
"""Analyze a completed Hall population/duration convergence campaign."""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile


class ConvergenceError(RuntimeError):
    pass


FIELD_OBSERVABLES = ("potential", "electric_x", "charge_density")
SPECIES_OBSERVABLES = (
    "number_density",
    "temperature_ev",
    "mean_velocity_x",
    "mean_velocity_y",
    "current_density_x",
    "current_density_y",
)
MODE_OBSERVABLES = (
    ("electric_x", ""),
    ("number_density", "electrons"),
)
ELEMENTARY_CHARGE_C = 1.602176634e-19


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ConvergenceError(f"cannot hash {path}: {error}") from error


def load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ConvergenceError(f"cannot read {label}: {error}") from error
    if not isinstance(value, dict):
        raise ConvergenceError(f"{label} must contain a JSON object")
    return value


def finite(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ConvergenceError(f"{context} is not numeric") from error
    if not math.isfinite(result):
        raise ConvergenceError(f"{context} is not finite")
    return result


def nonnegative_integer(value: str, context: str) -> int:
    numeric = finite(value, context)
    result = int(numeric)
    if numeric != result or result < 0:
        raise ConvergenceError(f"{context} is not a non-negative integer")
    return result


def read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or ())
            if not required <= fields:
                raise ConvergenceError(
                    f"{path.name} is missing columns: "
                    + ",".join(sorted(required - fields))
                )
            rows = list(reader)
    except OSError as error:
        raise ConvergenceError(f"cannot read {path}: {error}") from error
    if not rows:
        raise ConvergenceError(f"{path.name} has no rows")
    return rows


def load_controller(
    output: Path,
    expected_updates: int,
) -> dict[str, float | int | bool | None]:
    path = output / "current_source.csv"
    required = {
        "macro_particles_created",
        "represented_particles_created",
        "control_updates",
        "reverse_diagnostics_start_step",
        "reverse_demand_steps",
        "reverse_demand_step_fraction",
        "cumulative_reverse_demand_macroparticles",
        "maximum_reverse_demand_macroparticles",
        "cumulative_monitored_negative_charge",
        "cumulative_monitored_positive_charge",
        "cumulative_processed_monitored_charge",
    }
    final = read_csv(path, required)[-1]
    distribution_keys = {
        "reverse_distribution_start_step",
        "reverse_distribution_steps",
        "reverse_one_macro_steps",
        "reverse_two_macro_steps",
        "reverse_multi_macro_steps",
        "distributed_reverse_demand_macroparticles",
        "mean_reverse_demand_macroparticles",
        "rms_reverse_demand_macroparticles",
        "reverse_monitored_negative_charge",
        "reverse_monitored_positive_charge",
        "reverse_monitored_net_charge",
    }
    observed_distribution = distribution_keys & set(final)
    if observed_distribution and observed_distribution != distribution_keys:
        raise ConvergenceError(
            f"{path.name} has a partial reverse-distribution schema"
        )
    created = nonnegative_integer(
        final["macro_particles_created"],
        f"{path.name} macro particles created",
    )
    represented = finite(
        final["represented_particles_created"],
        f"{path.name} represented particles created",
    )
    updates = nonnegative_integer(
        final["control_updates"], f"{path.name} control updates"
    )
    start_step = nonnegative_integer(
        final["reverse_diagnostics_start_step"],
        f"{path.name} reverse diagnostic start step",
    )
    reverse_steps = nonnegative_integer(
        final["reverse_demand_steps"],
        f"{path.name} reverse-demand steps",
    )
    fraction = finite(
        final["reverse_demand_step_fraction"],
        f"{path.name} reverse-demand fraction",
    )
    cumulative = finite(
        final["cumulative_reverse_demand_macroparticles"],
        f"{path.name} cumulative reverse demand",
    )
    maximum = finite(
        final["maximum_reverse_demand_macroparticles"],
        f"{path.name} maximum reverse demand",
    )
    negative = finite(
        final["cumulative_monitored_negative_charge"],
        f"{path.name} monitored negative charge",
    )
    positive = finite(
        final["cumulative_monitored_positive_charge"],
        f"{path.name} monitored positive charge",
    )
    processed = finite(
        final["cumulative_processed_monitored_charge"],
        f"{path.name} processed monitored charge",
    )
    expected_fraction = reverse_steps / updates if updates else 0.0
    if (
        created <= 0
        or represented <= 0.0
        or updates != expected_updates
        or start_step != 0
        or reverse_steps > updates
        or cumulative < 0.0
        or maximum < 0.0
        or maximum > cumulative + 1e-12
        or negative > 0.0
        or positive < 0.0
        or not math.isclose(
            fraction, expected_fraction, rel_tol=1e-12, abs_tol=1e-15
        )
        or not math.isclose(
            negative + positive,
            processed,
            rel_tol=1e-12,
            abs_tol=1e-24,
        )
    ):
        raise ConvergenceError(
            f"{path.name} controller diagnostics are inconsistent "
            "or do not cover the complete run"
        )
    macro_charge = (
        represented / created * ELEMENTARY_CHARGE_C
    )
    cumulative_charge = cumulative * macro_charge
    maximum_charge = maximum * macro_charge
    monitored_absolute_charge = abs(negative) + positive
    result: dict[str, float | int | bool | None] = {
        "control_updates": updates,
        "reverse_demand_steps": reverse_steps,
        "reverse_demand_step_fraction": fraction,
        "cumulative_reverse_demand_macroparticles": cumulative,
        "maximum_reverse_demand_macroparticles": maximum,
        "macro_charge_c": macro_charge,
        "cumulative_reverse_demand_charge_c": cumulative_charge,
        "reverse_demand_charge_per_update_c":
            cumulative_charge / updates,
        "maximum_reverse_demand_charge_c": maximum_charge,
        "reverse_demand_fraction_of_absolute_monitored_charge":
            cumulative_charge / monitored_absolute_charge
            if monitored_absolute_charge > 0.0 else 0.0,
        "reverse_distribution_available":
            observed_distribution == distribution_keys,
    }
    if observed_distribution != distribution_keys:
        return result
    distribution_start = nonnegative_integer(
        final["reverse_distribution_start_step"],
        f"{path.name} reverse distribution start step",
    )
    distribution_steps = nonnegative_integer(
        final["reverse_distribution_steps"],
        f"{path.name} reverse distribution steps",
    )
    one_steps = nonnegative_integer(
        final["reverse_one_macro_steps"],
        f"{path.name} one-macro reverse steps",
    )
    two_steps = nonnegative_integer(
        final["reverse_two_macro_steps"],
        f"{path.name} two-macro reverse steps",
    )
    multi_steps = nonnegative_integer(
        final["reverse_multi_macro_steps"],
        f"{path.name} multi-macro reverse steps",
    )
    distributed = finite(
        final["distributed_reverse_demand_macroparticles"],
        f"{path.name} distributed reverse demand",
    )
    mean = finite(
        final["mean_reverse_demand_macroparticles"],
        f"{path.name} mean reverse demand",
    )
    rms = finite(
        final["rms_reverse_demand_macroparticles"],
        f"{path.name} RMS reverse demand",
    )
    reverse_negative = finite(
        final["reverse_monitored_negative_charge"],
        f"{path.name} reverse monitored negative charge",
    )
    reverse_positive = finite(
        final["reverse_monitored_positive_charge"],
        f"{path.name} reverse monitored positive charge",
    )
    reverse_net = finite(
        final["reverse_monitored_net_charge"],
        f"{path.name} reverse monitored net charge",
    )
    expected_mean = (
        distributed / distribution_steps
        if distribution_steps else 0.0
    )
    if (
        distribution_start != 0
        or distribution_steps != one_steps + two_steps + multi_steps
        or distribution_steps > reverse_steps
        or distributed < 0.0
        or mean < 0.0
        or rms < mean
        or reverse_negative > 0.0
        or reverse_positive < 0.0
        or not math.isclose(
            mean, expected_mean, rel_tol=1e-12, abs_tol=1e-15
        )
        or not math.isclose(
            reverse_negative + reverse_positive,
            reverse_net,
            rel_tol=1e-12,
            abs_tol=1e-24,
        )
    ):
        raise ConvergenceError(
            f"{path.name} reverse distribution is inconsistent "
            "or does not cover the complete run"
        )
    result.update({
        "reverse_distribution_start_step": distribution_start,
        "reverse_distribution_steps": distribution_steps,
        "reverse_one_macro_steps": one_steps,
        "reverse_two_macro_steps": two_steps,
        "reverse_multi_macro_steps": multi_steps,
        "distributed_reverse_demand_macroparticles": distributed,
        "mean_reverse_demand_macroparticles": mean,
        "rms_reverse_demand_macroparticles": rms,
        "reverse_monitored_negative_charge_c": reverse_negative,
        "reverse_monitored_positive_charge_c": reverse_positive,
        "reverse_monitored_net_charge_c": reverse_net,
    })
    return result


def atomic_json(path: Path, value: dict[str, object]) -> None:
    if path.exists():
        raise ConvergenceError(f"refusing to overwrite report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
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


def coordinate_vectors(
    rows: list[dict[str, str]],
    observables: tuple[str, ...],
    source: str,
) -> tuple[list[float], dict[str, list[float]]]:
    coordinates = [
        finite(row["coordinate"], f"{source} coordinate") for row in rows
    ]
    if coordinates != sorted(set(coordinates)):
        raise ConvergenceError(
            f"{source} coordinates must be unique and increasing"
        )
    vectors = {
        observable: [
            finite(row[observable], f"{source} {observable}")
            for row in rows
        ]
        for observable in observables
    }
    return coordinates, vectors


def deck_contract(path: Path) -> dict[str, int]:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string(
            "[global]\n" + path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, configparser.Error) as error:
        raise ConvergenceError(f"cannot read runtime deck {path}: {error}") from error
    section = parser["global"]
    try:
        return {
            "nodes": section.getint("nx"),
            "steps": section.getint("steps"),
            "start": section.getint("resolved_diagnostic_start_step"),
            "interval": section.getint("resolved_diagnostic_interval"),
            "max_mode": section.getint("resolved_max_mode"),
        }
    except (ValueError, configparser.Error) as error:
        raise ConvergenceError(
            f"runtime deck {path} has an invalid diagnostic contract"
        ) from error


def load_stage(
    output: Path,
    expected_samples: int,
    expected_nodes: int,
    expected_max_mode: int,
    expected_updates: int,
) -> dict[str, object]:
    field_path = output / "resolved_field_time_average.csv"
    field_rows = read_csv(
        field_path,
        {"samples", "profile_axis", "coordinate", *FIELD_OBSERVABLES},
    )
    if (
        len(field_rows) != expected_nodes
        or
        {row["profile_axis"] for row in field_rows} != {"x"}
        or {
            int(finite(row["samples"], "field samples"))
            for row in field_rows
        } != {expected_samples}
    ):
        raise ConvergenceError("field average has the wrong axis or samples")
    coordinates, vectors = coordinate_vectors(
        field_rows, FIELD_OBSERVABLES, field_path.name
    )

    species_path = output / "resolved_species_time_average.csv"
    species_rows = read_csv(
        species_path,
        {
            "samples", "profile_axis", "coordinate", "species",
            *SPECIES_OBSERVABLES,
        },
    )
    if (
        {row["profile_axis"] for row in species_rows} != {"x"}
        or {row["species"] for row in species_rows}
            != {"electrons", "ions"}
        or {
            int(finite(row["samples"], "species samples"))
            for row in species_rows
        } != {expected_samples}
    ):
        raise ConvergenceError(
            "species average has the wrong axis, species, or samples"
        )
    for species in ("electrons", "ions"):
        selected = [
            row for row in species_rows if row["species"] == species
        ]
        if len(selected) != len(coordinates):
            raise ConvergenceError(
                f"species average has incomplete {species} coverage"
            )
        selected_coordinates, selected_vectors = coordinate_vectors(
            selected, SPECIES_OBSERVABLES, f"{species_path.name}:{species}"
        )
        if selected_coordinates != coordinates:
            raise ConvergenceError(
                f"{species} coordinates do not match the field profile"
            )
        for observable, values in selected_vectors.items():
            vectors[f"{species}.{observable}"] = values

    mode_path = output / "resolved_modes.csv"
    mode_rows = read_csv(
        mode_path,
        {"mode", "quantity", "species", "amplitude"},
    )
    for quantity, species in MODE_OBSERVABLES:
        selected: dict[int, list[float]] = {}
        for row in mode_rows:
            if row["quantity"] != quantity or row["species"] != species:
                continue
            mode_value = finite(row["mode"], "resolved mode")
            mode = int(mode_value)
            if mode_value != mode or mode == 0:
                continue
            selected.setdefault(mode, []).append(
                finite(row["amplitude"], "resolved mode amplitude")
            )
        if not selected or sorted(selected) != list(
            range(1, expected_max_mode + 1)
        ):
            raise ConvergenceError(
                f"mode spectrum {quantity}:{species} is incomplete"
            )
        if {len(values) for values in selected.values()} != {
            expected_samples
        }:
            raise ConvergenceError(
                f"mode spectrum {quantity}:{species} has wrong sample count"
            )
        label = (
            f"modes.{species}.{quantity}"
            if species else f"modes.{quantity}"
        )
        vectors[label] = [
            statistics.fmean(selected[mode])
            for mode in sorted(selected)
        ]
    return {
        "coordinates": coordinates,
        "vectors": vectors,
        "controller": load_controller(output, expected_updates),
        "sha256": {
            "field_average": sha256(field_path),
            "species_average": sha256(species_path),
            "mode_history": sha256(mode_path),
            "current_source": sha256(output / "current_source.csv"),
        },
    }


def difference(
    candidate: list[float], baseline: list[float]
) -> dict[str, float]:
    if len(candidate) != len(baseline) or not baseline:
        raise ConvergenceError("convergence vectors have different shapes")
    delta = [
        candidate_value - baseline_value
        for candidate_value, baseline_value in zip(candidate, baseline)
    ]
    absolute_l2 = math.sqrt(statistics.fmean(
        value * value for value in delta
    ))
    absolute_linf = max(abs(value) for value in delta)
    baseline_l2 = math.sqrt(statistics.fmean(
        value * value for value in baseline
    ))
    baseline_linf = max(abs(value) for value in baseline)
    if baseline_l2 == 0.0 or baseline_linf == 0.0:
        raise ConvergenceError(
            "cannot normalize an identically zero convergence observable"
        )
    return {
        "absolute_l2": absolute_l2,
        "absolute_linf": absolute_linf,
        "relative_l2": absolute_l2 / baseline_l2,
        "relative_linf": absolute_linf / baseline_linf,
    }


def safe_ratio(candidate: float, baseline: float) -> float | None:
    if baseline > 0.0:
        return candidate / baseline
    return 0.0 if candidate == 0.0 else None


def analyze(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.convergence_manifest.resolve()
    root = manifest_path.parent
    campaign = load_json(manifest_path, "convergence manifest")
    if campaign.get("hall_convergence_version") != 2:
        raise ConvergenceError("hall_convergence_version must be 2")
    case_path = Path(str(campaign.get("case_manifest", "")))
    if sha256(case_path) != campaign.get("case_manifest_sha256"):
        raise ConvergenceError("case-manifest SHA-256 mismatch")
    runs = campaign.get("runs")
    if not isinstance(runs, list) or len(runs) != 5:
        raise ConvergenceError("convergence campaign must contain five runs")
    analyzed_axes = (
        ("population", "duration")
        if args.axis == "all" else (args.axis,)
    )
    required_stages = {"population_1"}
    if "population" in analyzed_axes:
        required_stages.update({"population_0p5", "population_2"})
    if "duration" in analyzed_axes:
        required_stages.update({"duration_0p5", "duration_2"})
    acceptance = campaign.get("acceptance")
    if not isinstance(acceptance, dict):
        raise ConvergenceError("convergence acceptance contract is missing")
    l2_tolerance = float(acceptance["relative_l2_tolerance"])
    linf_tolerance = float(acceptance["relative_linf_tolerance"])
    ratio_limit = float(
        acceptance["maximum_fine_to_coarse_change_ratio"]
    )
    reverse_charge_ratio_limit = float(
        acceptance[
            "maximum_fine_to_baseline_reverse_charge_per_update_ratio"
        ]
    )
    reverse_impulse_ratio_limit = float(
        acceptance[
            "maximum_fine_to_baseline_reverse_impulse_ratio"
        ]
    )
    reverse_macro_limit = float(
        acceptance[
            "maximum_fine_reverse_demand_macroparticles_per_update"
        ]
    )
    if any(
        not math.isfinite(value) or value <= 0.0
        for value in (
            l2_tolerance,
            linf_tolerance,
            ratio_limit,
            reverse_charge_ratio_limit,
            reverse_impulse_ratio_limit,
            reverse_macro_limit,
        )
    ):
        raise ConvergenceError(
            "convergence acceptance limits must be positive and finite"
        )

    stages: dict[str, dict[str, object]] = {}
    records: list[dict[str, object]] = []
    for run in runs:
        if not isinstance(run, dict):
            raise ConvergenceError("convergence run must be an object")
        stage = str(run.get("stage", ""))
        if not stage or stage in stages:
            raise ConvergenceError("convergence stage identity is invalid")
        if stage not in required_stages:
            continue
        deck = root / str(run.get("runtime_config", ""))
        if sha256(deck) != run.get("runtime_config_sha256"):
            raise ConvergenceError(f"{stage} runtime-config hash mismatch")
        contract = deck_contract(deck)
        if (
            contract["steps"] != int(run["steps"])
            or contract["start"] != int(run["diagnostic_start_step"])
            or contract["interval"] != int(run["diagnostic_interval"])
        ):
            raise ConvergenceError(
                f"{stage} deck and campaign cadence differ"
            )
        output = Path(str(run.get("result_dir", ""))).resolve()
        loaded = load_stage(
            output,
            int(run["diagnostic_samples"]),
            contract["nodes"],
            contract["max_mode"],
            contract["steps"],
        )
        stages[stage] = loaded
        records.append({
            "stage": stage,
            "axis": run.get("axis"),
            "factor": run.get("factor"),
            "output_dir": str(output),
            "runtime_config_sha256": run.get("runtime_config_sha256"),
            "artifacts": loaded["sha256"],
            "controller": loaded["controller"],
        })

    baseline = stages.get("population_1")
    if baseline is None:
        raise ConvergenceError("population_1 baseline is missing")
    all_axes = {
        "population": ("population_0p5", "population_2"),
        "duration": ("duration_0p5", "duration_2"),
    }
    comparisons: dict[str, object] = {}
    passed = True
    baseline_vectors = baseline["vectors"]
    assert isinstance(baseline_vectors, dict)
    for axis in analyzed_axes:
        coarse_name, fine_name = all_axes[axis]
        if (
            stages[coarse_name]["coordinates"]
                != baseline["coordinates"]
            or stages[fine_name]["coordinates"]
                != baseline["coordinates"]
        ):
            raise ConvergenceError(
                f"{axis} stage coordinates differ from the baseline"
            )
        coarse_vectors = stages[coarse_name]["vectors"]
        fine_vectors = stages[fine_name]["vectors"]
        if (
            set(coarse_vectors) != set(baseline_vectors)
            or set(fine_vectors) != set(baseline_vectors)
        ):
            raise ConvergenceError(
                f"{axis} stages expose different observables"
            )
        observable_results: dict[str, object] = {}
        for observable in sorted(baseline_vectors):
            coarse_change = difference(
                coarse_vectors[observable],
                baseline_vectors[observable],
            )
            fine_change = difference(
                fine_vectors[observable],
                baseline_vectors[observable],
            )
            denominator = coarse_change["relative_l2"]
            ratio = (
                fine_change["relative_l2"] / denominator
                if denominator > 0.0
                else (0.0 if fine_change["relative_l2"] == 0.0 else None)
            )
            observable_passed = (
                fine_change["relative_l2"] <= l2_tolerance
                and fine_change["relative_linf"] <= linf_tolerance
                and ratio is not None
                and ratio <= ratio_limit
            )
            passed = passed and observable_passed
            observable_results[observable] = {
                "coarse_to_baseline": coarse_change,
                "fine_to_baseline": fine_change,
                "fine_to_coarse_change_ratio": ratio,
                "passed": observable_passed,
            }
        comparisons[axis] = {
            "coarse_stage": coarse_name,
            "baseline_stage": "population_1",
            "fine_stage": fine_name,
            "observables": observable_results,
            "passed": all(
                value["passed"] for value in observable_results.values()
            ),
        }
    if "population" in analyzed_axes:
        controller_stages: dict[
            str, dict[str, float | int | bool | None]
        ] = {}
        for stage_name in (
            "population_0p5", "population_1", "population_2"
        ):
            controller = stages[stage_name]["controller"]
            if not isinstance(controller, dict):
                raise ConvergenceError(
                    f"{stage_name} controller diagnostics are missing"
                )
            controller_stages[stage_name] = controller
        baseline_controller = controller_stages["population_1"]
        fine_controller = controller_stages["population_2"]
        charge_ratio = safe_ratio(
            float(
                fine_controller[
                    "reverse_demand_charge_per_update_c"
                ]
            ),
            float(
                baseline_controller[
                    "reverse_demand_charge_per_update_c"
                ]
            ),
        )
        impulse_ratio = safe_ratio(
            float(
                fine_controller["maximum_reverse_demand_charge_c"]
            ),
            float(
                baseline_controller[
                    "maximum_reverse_demand_charge_c"
                ]
            ),
        )
        controller_passed = (
            charge_ratio is not None
            and charge_ratio <= reverse_charge_ratio_limit
            and impulse_ratio is not None
            and impulse_ratio <= reverse_impulse_ratio_limit
            and float(
                fine_controller[
                    "maximum_reverse_demand_macroparticles"
                ]
            ) <= reverse_macro_limit
        )
        passed = passed and controller_passed
        comparisons["controller_population"] = {
            "coarse_stage": "population_0p5",
            "baseline_stage": "population_1",
            "fine_stage": "population_2",
            "stages": controller_stages,
            "fine_to_baseline_reverse_charge_per_update_ratio":
                charge_ratio,
            "fine_to_baseline_reverse_impulse_ratio": impulse_ratio,
            "event_frequency_is_acceptance_metric": False,
            "passed": controller_passed,
            "interpretation": (
                "Acceptance uses represented charge because macro-particle "
                "weight changes with population. Reverse-event frequency is "
                "reported but is not required to decrease."
            ),
        }
    return {
        "schema_version": 2,
        "case_id": campaign.get("case_id"),
        "physics_claim": "none",
        "passed": passed,
        "analyzed_axes": list(analyzed_axes),
        "convergence_manifest": str(manifest_path),
        "convergence_manifest_sha256": sha256(manifest_path),
        "acceptance": acceptance,
        "runs": records,
        "comparisons": comparisons,
        "limitations": [
            "This is a same-seed fixed-grid population/duration study.",
            "The controller test diagnoses population scaling over this "
            "early-time window; it does not validate the cathode model.",
            "Grid, timestep, random-seed, and published-reference "
            "convergence remain independent requirements.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze completed Hall convergence outputs"
    )
    parser.add_argument("convergence_manifest", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--axis",
        choices=("all", "population", "duration"),
        default="all",
        help=(
            "analyze all stages or one completed convergence axis; "
            "the shared population_1 baseline is always required"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args)
        atomic_json(args.report, report)
    except (ConvergenceError, KeyError, TypeError, ValueError) as error:
        print(f"Hall convergence analysis error: {error}", file=sys.stderr)
        return 2
    print(
        f"Hall convergence analysis {'passed' if report['passed'] else 'failed'}: "
        f"{args.report}"
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
