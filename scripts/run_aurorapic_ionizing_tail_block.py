#!/usr/bin/env python3
"""Run the predeclared regional, phase-resolved ionizing-tail continuation."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, finite, integer,
    insert_global, run_process, set_global, sha256, table,
)


RULE_SHA256S = {
    "a8cded31a57af98b6c32dda816d122cda8ebf9f7daf31d906b516d1a5b12b9f2",
    "e0216347692759a4c775f3cc5b932ce5c36c62a6a0a45bee364ac0bae5380704",
    "3ef905d4e89f9503e3c5f6fa616bbe180d47fdf39f7301043540d6203c6a2e26",
}
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_ONE_BOUNDED_TAIL_BLOCK"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"


class TailBlockError(RuntimeError):
    pass


def build_deck(base: str, output: Path, checkpoint: Path,
               rule: dict[str, object]) -> str:
    execution = rule["execution_contract"]
    diagnostic = rule["diagnostic_contract"]
    regions = ",".join(
        f"{region['name']}:{region['x_min_m']}:{region['x_max_m']}"
        for region in diagnostic["regions"])
    values = {
        "steps": str(execution["end_step"]),
        "output_interval": str(execution["output_interval_steps"]),
        "output_dir": str(output),
        "spatial_average": "true",
        "spatial_average_interval": str(
            diagnostic["spatial_average_interval_steps"]),
        "spatial_average_start_step": str(execution["start_step"] + 1),
        "spatial_average_end_step": str(execution["end_step"]),
        "spatial_average_rf_frequency": "13560000",
        "spatial_average_rf_cycles": str(execution["cycles"]),
        "spatial_average_phase_bins": str(diagnostic["phase_bins"]),
        "spatial_average_reset_on_restart": "true",
        "spatial_average_sampling_order": "pre_collision",
        "phase_eedf": "true",
        "phase_eedf_species": str(diagnostic["phase_eedf_species"]),
        "phase_eedf_energy_bins": str(diagnostic["energy_bins"]),
        "phase_eedf_energy_max": str(diagnostic["energy_max_eV"]),
        "phase_eedf_regions": regions,
        "wall_impact_spectrum": "true",
        "wall_impact_reset_on_restart": "true",
        "checkpoint_interval": str(execution["end_step"]),
        "runtime_backend": "serial",
        "runtime_threads": "1",
        "restart_path": str(checkpoint),
    }
    result = base
    for key, value in values.items():
        result = set_global(result, key, value)
    result = insert_global(
        result, "phase_eedf_tail_threshold",
        str(diagnostic.get("tail_threshold_eV", 0.0)))
    surface = diagnostic.get("surface_flux")
    if surface is not None:
        additions = {
            "phase_surface_flux": "true",
            "phase_surface_flux_reset_on_restart": str(
                bool(surface["reset_on_restart"])).lower(),
            "phase_surface_flux_species": str(surface["species"]),
            "phase_surface_flux_positions": ",".join(
                str(position) for position in surface["positions_m"]),
            "phase_surface_flux_energy_bins": str(surface["energy_bins"]),
            "phase_surface_flux_energy_max": str(surface["energy_max_eV"]),
        }
        for key, value in additions.items():
            result = insert_global(result, key, value)
    return result


def analyze_surface_flux(output: Path, diagnostic: dict[str, object],
                         metadata: dict[str, object]) -> dict[str, object]:
    surface = diagnostic["surface_flux"]
    phases = int(diagnostic["phase_bins"])
    positions = [float(value) for value in surface["positions_m"]]
    bins = int(surface["energy_bins"])
    directions = list(surface["direction_order"])
    contract_ok = (
        metadata.get("phase_surface_flux_enabled") is True and
        metadata.get("phase_surface_flux_species") == surface["species"] and
        metadata.get("phase_surface_flux_energy_bins") == bins and
        metadata.get("phase_surface_flux_energy_max") ==
            float(surface["energy_max_eV"]) and
        metadata.get("phase_surface_flux_positions") == positions)

    summary = table(output / "phase_surface_flux_summary.csv")
    expected_summary_rows = phases * len(positions) * len(directions)
    summary_shape = len(summary) == expected_summary_rows
    summary_counts: dict[tuple[int, int, str], tuple[float, float]] = {}
    totals = [0 for _ in positions]
    finite_summary = True
    for row in summary:
        phase = integer(row, "phase_bin", "surface-flux summary")
        surface_id = integer(row, "surface_id", "surface-flux summary")
        direction = row.get("direction", "")
        represented = finite(
            row, "represented_crossings", "surface-flux summary")
        overflow = finite(row, "overflow_fraction", "surface-flux summary")
        macro = integer(row, "macro_crossings", "surface-flux summary")
        finite_summary = finite_summary and all(math.isfinite(float(row[key]))
            for key in ("represented_particle_flux_m-2_s-1",
                        "kinetic_energy_flux_W_m-2"))
        if (phase < 0 or phase >= phases or surface_id < 0 or
                surface_id >= len(positions) or direction not in directions or
                represented < 0.0 or overflow < 0.0 or overflow > 1.0):
            summary_shape = False
            continue
        key = (phase, surface_id, direction)
        if key in summary_counts:
            summary_shape = False
        summary_counts[key] = (represented, represented * overflow)
        totals[surface_id] += macro

    histogram_sums: dict[tuple[int, int, str], float] = {}
    histogram_rows = 0
    finite_histogram = True
    with (output / "phase_surface_flux.csv").open(
            newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            histogram_rows += 1
            phase = integer(row, "phase_bin", "surface-flux histogram")
            surface_id = integer(row, "surface_id", "surface-flux histogram")
            direction = row.get("direction", "")
            count = finite(
                row, "represented_crossings", "surface-flux histogram")
            density = finite(
                row, "probability_density", "surface-flux histogram")
            finite_histogram = finite_histogram and count >= 0.0 and density >= 0.0
            key = (phase, surface_id, direction)
            histogram_sums[key] = histogram_sums.get(key, 0.0) + count
    expected_histogram_rows = expected_summary_rows * bins
    closure = summary_shape and len(histogram_sums) == expected_summary_rows
    if closure:
        for key, (represented, overflow) in summary_counts.items():
            difference = abs(histogram_sums.get(key, -1.0) + overflow - represented)
            if difference > 1e-9 * max(1.0, represented):
                closure = False
                break
    minimum = int(surface["minimum_total_macro_crossings_per_surface"])
    return {
        "contract": contract_ok,
        "shape": summary_shape and histogram_rows == expected_histogram_rows,
        "finite": finite_summary and finite_histogram,
        "histogram_closure": closure,
        "crossing_sufficiency": all(total >= minimum for total in totals),
        "total_macro_crossings_per_surface": totals,
        "summary_rows": len(summary),
        "histogram_rows": histogram_rows,
    }


def validate_inputs(args: argparse.Namespace,
                    rule: dict[str, object]) -> tuple[int, int]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise TailBlockError("missing exact bounded-run acknowledgement")
    locked = rule["locked_initial_state"]
    paths = {
        "solver_sha256": args.executable.resolve(),
        "base_config_sha256": args.base_config.resolve(),
        "checkpoint_sha256": args.checkpoint.resolve(),
        "prior_report_sha256": args.prior_report.resolve(),
    }
    for key, path in paths.items():
        if sha256(path) != locked[key]:
            raise TailBlockError(f"locked input differs: {key}")
    execution = rule["execution_contract"]
    start = int(execution["start_step"])
    end = int(execution["end_step"])
    if (start != int(locked["step"]) or end <= start or
            end - start != int(execution["cycles"]) *
            int(execution["steps_per_cycle"])):
        raise TailBlockError("execution horizon is inconsistent")
    return start, end


def analyze_output(output: Path, rule: dict[str, object],
                   resources: dict[str, object]) -> dict[str, object]:
    execution = rule["execution_contract"]
    diagnostic = rule["diagnostic_contract"]
    start, end = int(execution["start_step"]), int(execution["end_step"])
    metadata_path = output / "spatial_average_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected_samples = int(diagnostic["spatial_average_samples"])
    expected_per_phase = int(diagnostic["samples_per_phase_bin"])
    metadata_ok = (
        metadata.get("start_step") == start + 1 and
        metadata.get("end_step") == end and
        metadata.get("interval") ==
            int(diagnostic["spatial_average_interval_steps"]) and
        metadata.get("sampling_order") == diagnostic["sampling_order"] and
        metadata.get("samples") == expected_samples and
        metadata.get("moment_samples") == expected_samples and
        metadata.get("phase_bins") == int(diagnostic["phase_bins"]) and
        metadata.get("phase_eedf_tail_threshold") == float(
            diagnostic.get("tail_threshold_eV", 0.0)) and
        set(metadata.get("phase_bin_samples", [])) == {expected_per_phase} and
        metadata.get("complete") is True)

    scalars = table(output / "scalars.csv")
    if (integer(scalars[0], "step", "initial scalar") != start or
            integer(scalars[-1], "step", "final scalar") != end):
        raise TailBlockError("scalar endpoints differ")
    maximum_particles = max(integer(row, "live_particles", "scalar")
                            for row in scalars)
    finite_scalars = all(math.isfinite(float(row[key])) for row in scalars
                         for key in ("time", "kinetic_energy", "field_energy",
                                     "total_energy"))
    maximum_field = 0.0
    field_files = sorted(output.glob("fields_*.csv"))
    if not field_files:
        raise TailBlockError("field snapshots are missing")
    for path in field_files:
        maximum_field = max(
            maximum_field,
            *(abs(finite(row, "E", path.name)) for row in table(path)))

    moments_path = output / "phase_eedf_moments.csv"
    with moments_path.open(newline="", encoding="utf-8") as stream:
        moments = list(csv.DictReader(stream))
    phase_bins = int(diagnostic["phase_bins"])
    regions = diagnostic["regions"]
    if len(moments) != phase_bins * len(regions):
        raise TailBlockError("phase EEDF moment shape differs")
    critical_regions = set(diagnostic.get("critical_regions", []))
    observation_rows = ([row for row in moments
                         if row["region"] in critical_regions]
                        if critical_regions else moments)
    minimum_observations = min(int(row["macro_observations"])
                               for row in observation_rows)
    maximum_overflow = max(float(row["overflow_fraction"])
                           for row in moments)
    finite_moments = all(
        math.isfinite(float(row[key])) for row in moments for key in (
            "represented_observations", "overflow_fraction", "mean_energy",
            "energy_standard_deviation", "mean_velocity_x", "mean_velocity_y",
            "mean_velocity_z", "drift_separated_temperature", "temperature_x",
            "temperature_y", "temperature_z", "tail_threshold",
            "tail_represented_observations", "tail_positive_x_fraction",
            "tail_negative_x_fraction",
            "tail_directional_population_imbalance", "tail_mean_velocity_x",
            "tail_longitudinal_energy_fraction"))
    region_names = sorted({row["region"] for row in moments})
    expected_names = sorted(str(region["name"]) for region in regions)
    shape_ok = region_names == expected_names

    checkpoint = output / f"checkpoint_{end}.apc"
    required = [
        "phase_eedf.csv", "phase_eedf_moments.csv",
        "spatial_phase_collision_power.csv",
        "spatial_phase_collision_rate.csv", "spatial_phase_moments.csv",
        "spatial_phase_fields.csv", "scalars.csv", "collisions.csv",
        "boundary_losses.csv", "spatial_average_metadata.json"]
    surface_result = None
    if "surface_flux" in diagnostic:
        required.extend((
            "phase_surface_flux.csv", "phase_surface_flux_summary.csv"))
    if not checkpoint.is_file() or any(not (output / name).is_file()
                                       for name in required):
        raise TailBlockError("required output is missing")
    gates = {
        "sampling_contract": metadata_ok,
        "phase_eedf_shape": shape_ok,
        "phase_eedf_observations": minimum_observations >= int(
            diagnostic["minimum_macro_observations_per_region_phase_bin"]),
        "phase_eedf_overflow": maximum_overflow <= float(
            diagnostic["maximum_overflow_fraction"]),
        "finite_diagnostics": finite_scalars and finite_moments,
        "particle_cap": maximum_particles <= int(
            execution["maximum_total_particles"]),
        "absolute_field": maximum_field <= float(
            execution["maximum_absolute_field_V_m"]),
        "resident_memory": int(resources["peak_resident_set_kib"]) <= int(
            execution["maximum_peak_resident_set_kib"]),
    }
    if "surface_flux" in diagnostic:
        surface_result = analyze_surface_flux(output, diagnostic, metadata)
        gates.update({f"surface_flux_{key}": value for key, value in
                      surface_result.items() if isinstance(value, bool)})
    return {
        "sampling": {
            "samples": metadata.get("samples"),
            "phase_bins": metadata.get("phase_bins"),
            "samples_per_phase_bin": expected_per_phase,
            "minimum_macro_observations_per_region_phase_bin":
                minimum_observations,
            "maximum_phase_region_overflow_fraction": maximum_overflow,
        },
        "safety": {
            "maximum_live_particles": maximum_particles,
            "maximum_sampled_absolute_field_V_m": maximum_field,
            "field_snapshot_files": len(field_files),
        },
        "surface_flux": surface_result,
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "final_checkpoint_sha256": sha256(checkpoint),
        "output_hashes": {name: sha256(output / name) for name in required},
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve()
    rule_sha256 = sha256(rule_path)
    if rule_sha256 not in RULE_SHA256S:
        raise TailBlockError("ionizing-tail rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    validate_inputs(args, rule)
    work = args.work_dir.resolve()
    if work.exists():
        raise TailBlockError(f"refusing to overwrite {work}")
    execution = rule["execution_contract"]
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    if available_memory < int(execution["minimum_available_memory_kib"]):
        raise TailBlockError("available memory is below the launch floor")
    if available_disk < int(execution["minimum_available_disk_kib"]):
        raise TailBlockError("available disk is below the launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck = work / "input.cfg"
    atomic_text(deck, build_deck(
        args.base_config.read_text(encoding="utf-8"), output,
        args.checkpoint.resolve(), rule))
    resources = run_process([
        str(args.executable.resolve()), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(deck),
    ], work / "stdout.txt", work / "stderr.txt",
        timeout_seconds=float(execution["timeout_seconds"]))
    result = analyze_output(output, rule, resources)
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_spatial_phase_ionizing_tail_block",
        "rule_sha256": rule_sha256,
        "inputs": {
            "solver_sha256": sha256(args.executable.resolve()),
            "base_config_sha256": sha256(args.base_config.resolve()),
            "prior_report_sha256": sha256(args.prior_report.resolve()),
            "input_checkpoint_sha256": sha256(args.checkpoint.resolve()),
            "deck_sha256": sha256(deck),
        },
        "resources": {
            **resources,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk,
        },
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    })
    atomic_json(work / "block-report.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_config", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("prior_report", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--acknowledge-cost", required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (TailBlockError, PilotError, OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps({
        "report": str(args.work_dir.resolve() / "block-report.json"),
        "all_gates_passed": result["all_gates_passed"],
        "sampling": result["sampling"],
        "resources": result["resources"],
    }, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
