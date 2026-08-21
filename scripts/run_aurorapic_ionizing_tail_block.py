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
    run_process, set_global, sha256, table,
)


RULE_SHA256 = (
    "a8cded31a57af98b6c32dda816d122cda8ebf9f7daf31d906b516d1a5b12b9f2")
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
    return result


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
    minimum_observations = min(int(row["macro_observations"])
                               for row in moments)
    maximum_overflow = max(float(row["overflow_fraction"])
                           for row in moments)
    finite_moments = all(
        math.isfinite(float(row[key])) for row in moments for key in (
            "represented_observations", "overflow_fraction", "mean_energy",
            "energy_standard_deviation", "mean_velocity_x", "mean_velocity_y",
            "mean_velocity_z", "drift_separated_temperature"))
    region_names = sorted({row["region"] for row in moments})
    expected_names = sorted(str(region["name"]) for region in regions)
    shape_ok = region_names == expected_names

    checkpoint = output / f"checkpoint_{end}.apc"
    required = (
        "phase_eedf.csv", "phase_eedf_moments.csv",
        "spatial_phase_collision_rate.csv", "spatial_phase_moments.csv",
        "spatial_phase_fields.csv", "scalars.csv", "collisions.csv",
        "boundary_losses.csv", "spatial_average_metadata.json")
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
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "final_checkpoint_sha256": sha256(checkpoint),
        "output_hashes": {name: sha256(output / name) for name in required},
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve()
    if sha256(rule_path) != RULE_SHA256:
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
        "rule_sha256": RULE_SHA256,
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
