#!/usr/bin/env python3
"""Compute the published Turner ion-density X-squared statistic."""

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


EXPECTED_NODES = {1: 129, 2: 257, 3: 513, 4: 513}
RANGES_95 = {1: (55.0, 303.0), 2: (177.0, 435.0),
             3: (405.0, 693.0), 4: (417.0, 665.0)}
RANGES_99 = {1: (48.0, 405.0), 2: (160.0, 548.0),
             3: (382.0, 798.0), 4: (392.0, 730.0)}
STEPS_PER_RF_CYCLE = {1: 400, 2: 800, 3: 1600, 4: 3200}
TOTAL_RF_CYCLES = {1: 1280, 2: 5120, 3: 5120, 4: 15360}
RF_FREQUENCY_HZ = 13.56e6
GAP_LENGTH_M = 0.067
# The publisher result files print mesh coordinates with approximately six
# significant digits. Their largest Case 1 discrepancy from the prescribed
# uniform grid is 1.91e-4 cell widths.
REFERENCE_COORDINATE_TOLERANCE_CELLS = 2.5e-4
class ComparisonError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ComparisonError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_report(path: Path, report: dict[str, object]) -> None:
    require(not path.exists(),
            f"refusing to overwrite existing report: {path}")
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


def load_columns(path: Path, columns: tuple[str, ...]) -> list[dict[str, float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, f"{path} has no CSV header")
            missing = [name for name in columns if name not in reader.fieldnames]
            require(not missing, f"{path} is missing columns {missing}")
            rows: list[dict[str, float]] = []
            for line, source in enumerate(reader, 2):
                try:
                    row = {name: float(source[name]) for name in columns}
                except (TypeError, ValueError) as error:
                    raise ComparisonError(
                        f"{path}:{line}: invalid numeric value"
                    ) from error
                require(all(math.isfinite(value) for value in row.values()),
                        f"{path}:{line}: non-finite numeric value")
                rows.append(row)
    except OSError as error:
        raise ComparisonError(f"cannot read {path}: {error}") from error
    return rows


def load_candidate(path: Path, species: str) -> list[dict[str, float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, f"{path} has no CSV header")
            prefix = "electron" if species == "electrons" else "ion"
            direct_density = f"{prefix}_density_mean_m-3"
            if "x_m" in reader.fieldnames and direct_density in reader.fieldnames:
                rows = []
                for line, source in enumerate(reader, 2):
                    try:
                        x = float(source["x_m"])
                        density = float(source[direct_density])
                    except (TypeError, ValueError) as error:
                        raise ComparisonError(
                            f"{path}:{line}: invalid numeric value"
                        ) from error
                    require(math.isfinite(x) and math.isfinite(density),
                            f"{path}:{line}: non-finite numeric value")
                    rows.append({
                        "x_m": x,
                        "density_mean_m-3": density,
                    })
                return rows

            required = {
                "species", "x_m", "number_density_mean_m-3"
            }
            missing = sorted(required - set(reader.fieldnames))
            require(not missing, f"{path} is missing columns {missing}")
            rows = []
            for line, source in enumerate(reader, 2):
                if source["species"] != species:
                    continue
                try:
                    x = float(source["x_m"])
                    density = float(source["number_density_mean_m-3"])
                except (TypeError, ValueError) as error:
                    raise ComparisonError(
                        f"{path}:{line}: invalid numeric value"
                    ) from error
                require(math.isfinite(x) and math.isfinite(density),
                        f"{path}:{line}: non-finite numeric value")
                rows.append({
                    "x_m": x,
                    "density_mean_m-3": density,
                })
            require(rows, f"{path} has no rows for species {species!r}")
            return rows
    except OSError as error:
        raise ComparisonError(f"cannot read {path}: {error}") from error


def compare(case: int, reference: Path, candidate: Path,
            audit_path: Path, metadata_path: Path,
            species: str = "ions",
            post_benchmark_window: bool = False,
            numerical_sensitivity: bool = False) -> dict[str, object]:
    require(
        not (post_benchmark_window and numerical_sensitivity),
        "post-benchmark and numerical-sensitivity modes are mutually exclusive",
    )
    require(species in ("electrons", "ions"),
            "Turner density species must be 'electrons' or 'ions'")
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ComparisonError(f"cannot read normalization audit: {error}") from error
    require(audit.get("turner_normalization_version") in (1, 2),
            "normalization audit version must be 1")
    files = audit.get("normalized_files")
    require(isinstance(files, dict) and reference.name in files,
            "reference is not present in the normalization audit")
    reference_lock = files[reference.name]
    require(isinstance(reference_lock, dict)
            and reference_lock.get("sha256") == sha256(reference),
            "reference SHA-256 differs from the normalization audit")
    expected_name = f"turner_case{case}_benchmark.csv"
    require(reference.name == expected_name,
            f"Case {case} statistical comparison requires {expected_name}")
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ComparisonError(
            f"cannot read candidate averaging metadata: {error}"
        ) from error
    published_steps_per_cycle = STEPS_PER_RF_CYCLE[case]
    if numerical_sensitivity:
        try:
            metadata_dt = float(metadata.get("dt", 0.0))
            metadata_frequency = float(metadata.get("rf_frequency", 0.0))
            steps_per_cycle_value = 1.0 / (metadata_frequency * metadata_dt)
        except (TypeError, ValueError, ZeroDivisionError) as error:
            raise ComparisonError(
                "numerical-sensitivity metadata has invalid RF timing"
            ) from error
        steps_per_cycle = int(round(steps_per_cycle_value))
        require(
            steps_per_cycle >= published_steps_per_cycle
            and steps_per_cycle % published_steps_per_cycle == 0
            and math.isclose(
                steps_per_cycle_value, steps_per_cycle,
                rel_tol=1e-12, abs_tol=1e-12,
            ),
            "numerical-sensitivity timestep must be an integer refinement "
            "of the published steps per RF cycle",
        )
    else:
        steps_per_cycle = published_steps_per_cycle
    end_step = steps_per_cycle * TOTAL_RF_CYCLES[case]
    samples = 32 * steps_per_cycle
    start_step = end_step - samples + 1
    metadata_version = metadata.get("spatial_average_version")
    require(
        metadata_version in (1, 2, 3, 4, 5, 6, 7),
        "candidate averaging metadata 'spatial_average_version' must be 1 through 7",
    )
    if metadata_version >= 6:
        require(
            metadata.get("sampling_order") == "post_collision",
            "Turner candidate sampling order must be 'post_collision'",
        )
    common_metadata = {
        "unit_system": "si",
        "interval": 1,
        "samples": samples,
        "expected_samples": samples,
        "rf_cycles": 32,
        "complete": True,
    }
    for key, expected in common_metadata.items():
        require(metadata.get(key) == expected,
                f"candidate averaging metadata {key!r} must be {expected!r}")
    if post_benchmark_window:
        diagnostic_start = metadata.get("start_step")
        diagnostic_end = metadata.get("end_step")
        require(
            isinstance(diagnostic_start, int)
            and isinstance(diagnostic_end, int)
            and diagnostic_start > end_step
            and diagnostic_end - diagnostic_start + 1 == samples
            and diagnostic_end % steps_per_cycle == 0
            and metadata.get("final_step") == diagnostic_end
            and metadata.get("reset_on_restart") is True,
            "post-benchmark candidate must be a complete reset 32-cycle "
            "window after the published duration",
        )
    else:
        expected_metadata = {
            "start_step": start_step,
            "end_step": end_step,
            "final_step": end_step,
        }
        for key, expected in expected_metadata.items():
            require(
                metadata.get(key) == expected,
                f"candidate averaging metadata {key!r} must be "
                f"{expected!r}",
            )
    require(
        math.isclose(
            float(metadata.get("rf_frequency", 0.0)),
            RF_FREQUENCY_HZ, rel_tol=1e-14, abs_tol=0.0),
        "candidate averaging metadata has the wrong RF frequency",
    )
    require(
        math.isclose(
            float(metadata.get("dt", 0.0)),
            1.0 / (RF_FREQUENCY_HZ * steps_per_cycle),
            rel_tol=1e-14, abs_tol=0.0),
        "candidate averaging metadata has the wrong timestep",
    )
    metadata_species = metadata.get("species")
    require(isinstance(metadata_species, list) and species in metadata_species,
            f"candidate averaging metadata does not contain species {species!r}")

    reference_prefix = "electron" if species == "electrons" else "ion"
    reference_mean_column = f"{reference_prefix}_density_mean_m-3"
    reference_stddev_column = (
        f"{reference_prefix}_density_population_stddev_m-3"
    )
    reference_rows = load_columns(reference, (
        "x_m", reference_mean_column, reference_stddev_column,
    ))
    candidate_rows = load_candidate(candidate, species)
    expected_nodes = EXPECTED_NODES[case]
    require(len(reference_rows) == expected_nodes,
            f"Case {case} reference must have {expected_nodes} nodes")
    candidate_nodes = len(candidate_rows)
    if numerical_sensitivity:
        require(
            candidate_nodes >= expected_nodes
            and (candidate_nodes - 1) % (expected_nodes - 1) == 0,
            f"Case {case} numerical-sensitivity candidate grid must be an "
            "integer refinement of the published node grid",
        )
        grid_refinement_ratio = (
            (candidate_nodes - 1) // (expected_nodes - 1)
        )
        comparison_rows = candidate_rows[::grid_refinement_ratio]
    else:
        require(candidate_nodes == expected_nodes,
                f"Case {case} candidate must have {expected_nodes} nodes")
        grid_refinement_ratio = 1
        comparison_rows = candidate_rows

    terms: list[float] = []
    squared_error: list[float] = []
    squared_reference: list[float] = []
    relative_errors: list[float] = []
    cell_width = GAP_LENGTH_M / (expected_nodes - 1)
    maximum_reference_coordinate_error = 0.0
    maximum_candidate_coordinate_error = 0.0
    for index, (ref, value) in enumerate(zip(reference_rows, comparison_rows)):
        x_ref = ref["x_m"]
        x_value = value["x_m"]
        x_expected = index * cell_width
        reference_error = abs(x_ref - x_expected)
        candidate_error = abs(x_value - x_expected)
        maximum_reference_coordinate_error = max(
            maximum_reference_coordinate_error, reference_error
        )
        maximum_candidate_coordinate_error = max(
            maximum_candidate_coordinate_error, candidate_error
        )
        require(
            reference_error
            <= REFERENCE_COORDINATE_TOLERANCE_CELLS * cell_width,
            f"reference coordinate is off the prescribed grid at node {index}",
        )
        require(
            candidate_error <= 1e-12 * max(1.0, abs(x_expected)),
            f"candidate coordinate is off the prescribed grid at node {index}",
        )
        mean = ref[reference_mean_column]
        sigma = ref[reference_stddev_column]
        density = value["density_mean_m-3"]
        require(mean > 0.0 and sigma > 0.0 and density >= 0.0,
                f"invalid density/statistical value at node {index}")
        delta = density - mean
        terms.append((delta / sigma) ** 2)
        squared_error.append(delta * delta)
        squared_reference.append(mean * mean)
        relative_errors.append(abs(delta) / mean)

    statistic = math.fsum(terms)
    interval_95 = RANGES_95[case]
    interval_99 = RANGES_99[case]
    statistic_report = {
        "name": (
            "Turner electron-density descriptive X-squared"
            if species == "electrons"
            else "Turner ion-density X-squared"
        ),
        "formula_variance": "population_standard_deviation_squared",
        "x_squared": statistic,
    }
    published_ion_acceptance = (
        species == "ions"
        and not post_benchmark_window
        and not numerical_sensitivity
    )
    if species == "ions":
        statistic_report.update({
            "range_95_percent": list(interval_95),
            "range_99_percent": list(interval_99),
        })
    if species == "electrons":
        statistic_report.update({
            "published_acceptance_applicable": False,
            "acceptance_reason": "published_ranges_are_for_ion_density",
        })
    elif not published_ion_acceptance:
        statistic_report.update({
            "published_acceptance_applicable": False,
            "within_published_95_percent_range":
                interval_95[0] <= statistic <= interval_95[1],
            "within_published_99_percent_range":
                interval_99[0] <= statistic <= interval_99[1],
        })
    else:
        statistic_report.update({
            "accepted_95_percent":
                interval_95[0] <= statistic <= interval_95[1],
            "accepted_99_percent":
                interval_99[0] <= statistic <= interval_99[1],
        })
    return {
        "turner_comparison_version": 1,
        "case": case,
        "reference": {
            "path": str(reference.resolve()),
            "sha256": sha256(reference),
            "normalization_audit_sha256": sha256(audit_path),
        },
        "candidate": {
            "path": str(candidate.resolve()),
            "sha256": sha256(candidate),
            "species": species,
            "averaging_metadata_path": str(metadata_path.resolve()),
            "averaging_metadata_sha256": sha256(metadata_path),
        },
        "nodes": expected_nodes,
        **({"candidate_nodes": candidate_nodes}
           if numerical_sensitivity else {}),
        "statistic": statistic_report,
        "secondary_metrics": {
            "relative_l2": math.sqrt(
                math.fsum(squared_error) / math.fsum(squared_reference)
            ),
            "maximum_pointwise_relative_error": max(relative_errors),
        },
        "coordinate_contract": {
            "mapping": "ordered_prescribed_grid_no_interpolation",
            "gap_length_m": GAP_LENGTH_M,
            "cell_width_m": cell_width,
            "reference_rounding_tolerance_cell_widths":
                REFERENCE_COORDINATE_TOLERANCE_CELLS,
            "maximum_reference_error_m":
                maximum_reference_coordinate_error,
            "maximum_candidate_error_m":
                maximum_candidate_coordinate_error,
        },
        "comparison_scope": (
            "published_baseline_electron_density_diagnostic_only"
            if species == "electrons"
            else
            "post_benchmark_density_diagnostic_only"
            if post_benchmark_window
            else "numerical_sensitivity_density_diagnostic_only"
            if numerical_sensitivity
            else "published_baseline_ion_density_statistic_only"
        ),
        **({
            "numerical_sensitivity_contract": {
                "published_steps_per_rf_cycle": published_steps_per_cycle,
                "candidate_steps_per_rf_cycle": steps_per_cycle,
                "time_refinement_ratio":
                    steps_per_cycle // published_steps_per_cycle,
                "grid_refinement_ratio": grid_refinement_ratio,
                "published_numerical_contract_changed": True,
            }
        } if numerical_sensitivity else {}),
        "averaging_contract_verified": True,
        "physics_claim": (
            "none_published_electron_density_descriptive_only"
            if species == "electrons"
            else
            "none_post_benchmark_window_outside_published_duration"
            if post_benchmark_window
            else "none_changed_published_numerical_contract"
            if numerical_sensitivity
            else "none_without_complete_solver_run_contract"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, choices=range(1, 5), required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--candidate-metadata", type=Path, required=True)
    parser.add_argument("--species", choices=("electrons", "ions"), default="ions")
    parser.add_argument("--normalization-audit", type=Path, required=True)
    parser.add_argument(
        "--post-benchmark-window", action="store_true",
        help="compare a reset 32-cycle diagnostic after the published run",
    )
    parser.add_argument(
        "--numerical-sensitivity", action="store_true",
        help="compare a complete refined diagnostic without published acceptance",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = compare(
            args.case, args.reference, args.candidate,
            args.normalization_audit, args.candidate_metadata,
            args.species, args.post_benchmark_window,
            args.numerical_sensitivity,
        )
        write_report(args.output, report)
    except (ComparisonError, OSError) as error:
        print(f"Turner comparison error: {error}", file=sys.stderr)
        return 2
    print(f"Turner Case {args.case} comparison written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
