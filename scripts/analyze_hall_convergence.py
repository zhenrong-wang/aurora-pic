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
        "sha256": {
            "field_average": sha256(field_path),
            "species_average": sha256(species_path),
            "mode_history": sha256(mode_path),
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


def analyze(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.convergence_manifest.resolve()
    root = manifest_path.parent
    campaign = load_json(manifest_path, "convergence manifest")
    if campaign.get("hall_convergence_version") != 1:
        raise ConvergenceError("hall_convergence_version must be 1")
    case_path = Path(str(campaign.get("case_manifest", "")))
    if sha256(case_path) != campaign.get("case_manifest_sha256"):
        raise ConvergenceError("case-manifest SHA-256 mismatch")
    runs = campaign.get("runs")
    if not isinstance(runs, list) or len(runs) != 5:
        raise ConvergenceError("convergence campaign must contain five runs")
    acceptance = campaign.get("acceptance")
    if not isinstance(acceptance, dict):
        raise ConvergenceError("convergence acceptance contract is missing")
    l2_tolerance = float(acceptance["relative_l2_tolerance"])
    linf_tolerance = float(acceptance["relative_linf_tolerance"])
    ratio_limit = float(
        acceptance["maximum_fine_to_coarse_change_ratio"]
    )

    stages: dict[str, dict[str, object]] = {}
    records: list[dict[str, object]] = []
    for run in runs:
        if not isinstance(run, dict):
            raise ConvergenceError("convergence run must be an object")
        stage = str(run.get("stage", ""))
        if not stage or stage in stages:
            raise ConvergenceError("convergence stage identity is invalid")
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
        )
        stages[stage] = loaded
        records.append({
            "stage": stage,
            "axis": run.get("axis"),
            "factor": run.get("factor"),
            "output_dir": str(output),
            "runtime_config_sha256": run.get("runtime_config_sha256"),
            "artifacts": loaded["sha256"],
        })

    baseline = stages.get("population_1")
    if baseline is None:
        raise ConvergenceError("population_1 baseline is missing")
    axes = {
        "population": ("population_0p5", "population_2"),
        "duration": ("duration_0p5", "duration_2"),
    }
    comparisons: dict[str, object] = {}
    passed = True
    baseline_vectors = baseline["vectors"]
    assert isinstance(baseline_vectors, dict)
    for axis, (coarse_name, fine_name) in axes.items():
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
    return {
        "schema_version": 1,
        "case_id": campaign.get("case_id"),
        "physics_claim": "none",
        "passed": passed,
        "convergence_manifest": str(manifest_path),
        "convergence_manifest_sha256": sha256(manifest_path),
        "acceptance": acceptance,
        "runs": records,
        "comparisons": comparisons,
        "limitations": [
            "This is a same-seed fixed-grid population/duration study.",
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
