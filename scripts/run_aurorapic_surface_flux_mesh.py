#!/usr/bin/env python3
"""Run one locked checkpoint-continuation arm of surface-flux mesh refinement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import sys

from run_aurorapic_edupic_pilot import (
    atomic_json, atomic_text, available_memory_kib, insert_global, integer,
    run_process, set_global, sha256, table,
)
from run_aurorapic_ionizing_tail_block import analyze_surface_flux


RULE_SHA256 = (
    "4bea77b968db89ca6a2e066a599d3e85b99c480de2f0cb6e56b12bdaeb891f54")
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_SURFACE_FLUX_MESH_RUN"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"


class MeshRunError(RuntimeError):
    pass


def build_deck(base: str, output: Path, checkpoint: Path,
               rule: dict[str, object], branch: str) -> str:
    fixed = rule["fixed_inputs"]
    configured = rule["branches"][branch]
    diagnostic = rule["diagnostic_contract"]
    values = {
        "nx": configured["nodes"],
        "steps": fixed["end_step"],
        "output_interval": 200,
        "output_dir": output,
        "spatial_average_interval": diagnostic[
            "spatial_average_interval_steps"],
        "spatial_average_start_step": int(fixed["start_step"]) + 1,
        "spatial_average_end_step": fixed["end_step"],
        "spatial_average_rf_cycles": fixed["measurement_cycles"],
        "spatial_average_phase_bins": diagnostic["phase_bins"],
        "checkpoint_interval": fixed["end_step"],
        "runtime_backend": "serial",
        "runtime_threads": 1,
    }
    result = base
    for key, value in values.items():
        result = set_global(result, key, str(value))
    additions = {
        "spatial_average_reset_on_restart": "true",
        "spatial_average_sampling_order": diagnostic[
            "spatial_sampling_order"],
        "phase_eedf": "true",
        "phase_eedf_species": "electrons",
        "phase_eedf_energy_bins": diagnostic["phase_eedf_energy_bins"],
        "phase_eedf_energy_max": diagnostic["phase_eedf_energy_max_eV"],
        "phase_eedf_regions": diagnostic["phase_eedf_regions"],
        "phase_surface_flux": "true",
        "phase_surface_flux_reset_on_restart": "true",
        "phase_surface_flux_species": diagnostic["surface_flux_species"],
        "phase_surface_flux_positions": ",".join(
            str(value) for value in diagnostic["surface_positions_m"]),
        "phase_surface_flux_energy_bins": diagnostic["surface_energy_bins"],
        "phase_surface_flux_energy_max": diagnostic["surface_energy_max_eV"],
        "restart_path": checkpoint,
    }
    for key, value in additions.items():
        result = insert_global(result, key, str(value))
    return result


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MeshRunError("mesh continuation cost was not acknowledged")
    rule_path = args.rule.resolve()
    if sha256(rule_path) != RULE_SHA256:
        raise MeshRunError("surface-flux mesh rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    if args.branch not in rule["branches"]:
        raise MeshRunError("unknown mesh branch")
    branch = rule["branches"][args.branch]
    for path, expected, label in (
        (args.executable.resolve(), rule["fixed_inputs"]["solver_sha256"],
         "solver"),
        (args.base_deck.resolve(), rule["fixed_inputs"]["base_deck_sha256"],
         "base deck"),
        (args.checkpoint.resolve(), branch["input_checkpoint_sha256"],
         "checkpoint"),
        (args.prior_report.resolve(), branch["prior_branch_report_sha256"],
         "prior branch report"),
    ):
        if sha256(path) != expected:
            raise MeshRunError(f"{label} SHA-256 differs")
    work = args.work_dir.resolve()
    if work.exists():
        raise MeshRunError(f"refusing to overwrite {work}")
    execution = rule["execution_contract"]
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    if available_memory < int(execution["minimum_available_memory_kib"]):
        raise MeshRunError("available memory is below the launch floor")
    if available_disk < int(execution["minimum_available_disk_kib"]):
        raise MeshRunError("available disk is below the launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck = work / "input.cfg"
    atomic_text(deck, build_deck(
        args.base_deck.read_text(encoding="utf-8"), output,
        args.checkpoint.resolve(), rule, args.branch))
    resources = run_process([
        str(args.executable.resolve()), "--allow-large-run",
        CLI_ACKNOWLEDGEMENT, str(deck)], work / "stdout.txt",
        work / "stderr.txt", float(branch["timeout_seconds"]))
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    fixed = rule["fixed_inputs"]
    diagnostic = rule["diagnostic_contract"]
    surface = analyze_surface_flux(output, {
        "phase_bins": diagnostic["phase_bins"],
        "surface_flux": {
            "species": diagnostic["surface_flux_species"],
            "positions_m": diagnostic["surface_positions_m"],
            "energy_bins": diagnostic["surface_energy_bins"],
            "energy_max_eV": diagnostic["surface_energy_max_eV"],
            "minimum_total_macro_crossings_per_surface": diagnostic[
                "minimum_total_macro_crossings_per_surface"],
            "direction_order": ["left_to_right", "right_to_left"],
        }}, metadata)
    scalars = table(output / "scalars.csv")
    fields = list(output.glob("fields_*.csv"))
    maximum_field = max(abs(float(row["E"]))
                        for path in fields for row in table(path))
    maximum_particles = max(integer(row, "live_particles", "scalar")
                            for row in scalars)
    gates = {
        "measurement_window": metadata.get("complete") is True and
            metadata.get("start_step") == int(fixed["start_step"]) + 1 and
            metadata.get("end_step") == int(fixed["end_step"]) and
            integer(scalars[0], "step", "initial scalar") ==
                int(fixed["start_step"]) and
            integer(scalars[-1], "step", "final scalar") ==
                int(fixed["end_step"]),
        "absolute_field": math.isfinite(maximum_field) and maximum_field <=
            float(execution["maximum_absolute_field_V_m"]),
        "particle_cap": maximum_particles <= int(
            execution["maximum_total_particles"]),
        "resident_memory": int(resources["peak_resident_set_kib"]) <= int(
            execution["maximum_peak_resident_set_kib"]),
        **{f"surface_flux_{key}": value for key, value in surface.items()
           if isinstance(value, bool)},
    }
    required = (
        "phase_surface_flux.csv", "phase_surface_flux_summary.csv",
        "spatial_phase_moments.csv", "spatial_phase_fields.csv",
        "spatial_phase_collision_power.csv", "spatial_average_metadata.json",
        "phase_eedf.csv", "phase_eedf_moments.csv", "scalars.csv")
    checkpoint = output / f"checkpoint_{fixed['end_step']}.apc"
    if not checkpoint.is_file() or any(not (output / name).is_file()
                                       for name in required):
        raise MeshRunError("required branch output is missing")
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "surface_flux_mesh_refinement_branch",
        "branch": args.branch,
        "rule_sha256": RULE_SHA256,
        "inputs": {
            "solver_sha256": sha256(args.executable.resolve()),
            "base_deck_sha256": sha256(args.base_deck.resolve()),
            "input_checkpoint_sha256": sha256(args.checkpoint.resolve()),
            "prior_branch_report_sha256": sha256(args.prior_report.resolve()),
            "deck_sha256": sha256(deck),
        },
        "numerics": {"nodes": branch["nodes"], "cells": branch["cells"],
                     "start_step": fixed["start_step"],
                     "end_step": fixed["end_step"]},
        "resources": {**resources,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk},
        "safety": {"maximum_live_particles": maximum_particles,
                   "maximum_sampled_absolute_field_V_m": maximum_field},
        "surface_flux": surface,
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "final_checkpoint_sha256": sha256(checkpoint),
        "output_hashes": {name: sha256(output / name) for name in required},
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    }
    atomic_json(work / "branch-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("prior_report", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--acknowledge-cost")
    try:
        report = execute(parser.parse_args())
    except (MeshRunError, OSError, ValueError, KeyError) as error:
        print(f"surface-flux mesh run rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
