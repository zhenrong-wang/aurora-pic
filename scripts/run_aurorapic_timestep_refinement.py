#!/usr/bin/env python3
"""Run one locked common-state arm of the argon timestep refinement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from run_aurorapic_edupic_pilot import (
    atomic_json, atomic_text, available_memory_kib, insert_global, run_analyzer,
    run_process, set_global, sha256, table,
)
from run_aurorapic_initialization_ab import set_species_value


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_TIMESTEP_REFINEMENT_RUN"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256 = (
    "d02b420adb3f457b5f6f471db6260fc127fda2435a83ddf7e307d22ab6640174")


class RefinementError(RuntimeError):
    pass


def common_deck(base: str, rule: dict[str, object], branch: str,
                output: Path) -> tuple[str, int, int]:
    state = rule["common_stationary_state"]
    config = rule["branches"][branch]
    steps_per_cycle = int(config["steps_per_rf_cycle"])
    equilibration_steps = int(config["equilibration_cycles"]) * steps_per_cycle
    measurement_steps = int(config["measurement_cycles"]) * steps_per_cycle
    result = base
    for species in ("electrons", "ions"):
        result = set_species_value(result, species, "particles", int(state[species]))
    for key, value in {
        "dt": config["timestep_s"],
        "steps": equilibration_steps,
        "output_interval": max(1, steps_per_cycle // 40),
        "output_dir": output,
        "spatial_average_start_step": 1,
        "spatial_average_end_step": equilibration_steps,
        "spatial_average_rf_cycles": int(config["equilibration_cycles"]),
        "spatial_average_phase_bins": 16,
        "checkpoint_interval": equilibration_steps,
        "runtime_backend": "serial",
        "runtime_threads": 1,
    }.items():
        result = set_global(result, key, str(value))
    return result, equilibration_steps, measurement_steps


def initial_deck(base: str, rule: dict[str, object], branch: str,
                 output: Path, state_path: Path) -> tuple[str, int, int]:
    result, equilibration_steps, measurement_steps = common_deck(
        base, rule, branch, output)
    state = rule["common_stationary_state"]
    result = insert_global(result, "initial_state_path", str(state_path))
    result = insert_global(
        result, "initial_state_signature", str(state["particle_state_signature"]))
    return result, equilibration_steps, measurement_steps


def measurement_deck(base: str, rule: dict[str, object], branch: str,
                     output: Path, checkpoint: Path) -> str:
    result, equilibration_steps, measurement_steps = common_deck(
        base, rule, branch, output)
    end_step = equilibration_steps + measurement_steps
    config = rule["branches"][branch]
    diagnostics = rule["fresh_measurement_contract"]
    for key, value in {
        "steps": end_step,
        "spatial_average_start_step": equilibration_steps + 1,
        "spatial_average_end_step": end_step,
        "spatial_average_rf_cycles": int(config["measurement_cycles"]),
        "checkpoint_interval": measurement_steps,
    }.items():
        result = set_global(result, key, str(value))
    for key, value in {
        "spatial_average_reset_on_restart": "true",
        "phase_eedf": "true",
        "phase_eedf_species": "electrons",
        "phase_eedf_energy_bins": diagnostics["phase_eedf_energy_bins"],
        "phase_eedf_energy_max": diagnostics["phase_eedf_energy_max_eV"],
        "phase_eedf_regions": "full_gap:0:0.025",
        "wall_impact_spectrum": "true",
        "wall_impact_reset_on_restart": "true",
        "wall_impact_energy_bins": diagnostics["wall_impact_raw_energy_bins"],
        "wall_impact_energy_max": diagnostics["wall_impact_energy_max_eV"],
        "restart_path": checkpoint,
    }.items():
        result = insert_global(result, key, str(value))
    return result


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise RefinementError("timestep refinement cost was not acknowledged")
    rule_path = args.rule.resolve()
    if sha256(rule_path) != APPROVED_RULE_SHA256:
        raise RefinementError("rule is not the approved timestep refinement")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    branch = args.branch
    if branch not in rule["branches"]:
        raise RefinementError("unknown refinement branch")
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    state_path = args.particle_state.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise RefinementError(f"refusing to overwrite {work}")
    fixed = rule["fixed_inputs"]
    state = rule["common_stationary_state"]
    for path, expected, label in (
        (executable, fixed["solver_sha256"], "solver"),
        (base_path, fixed["base_deck_sha256"], "base deck"),
        (state_path, state["particle_state_sha256"], "particle state"),
    ):
        if sha256(path) != expected:
            raise RefinementError(f"{label} SHA-256 differs")
    execution = rule["execution_contract"]
    available = available_memory_kib()
    if available < int(execution["minimum_available_memory_kib"]):
        raise RefinementError("available memory is below the launch floor")
    config = rule["branches"][branch]
    base = base_path.read_text(encoding="utf-8")
    work.mkdir(parents=True)
    equilibration = work / "equilibration"
    measurement = work / "measurement"
    equilibration.mkdir()
    measurement.mkdir()
    initial, equilibration_steps, measurement_steps = initial_deck(
        base, rule, branch, equilibration / "output", state_path)
    initial_path = equilibration / "input.cfg"
    atomic_text(initial_path, initial)
    equilibration_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(initial_path)], equilibration / "stdout.txt",
        equilibration / "stderr.txt", float(config["equilibration_timeout_seconds"]))
    checkpoint = equilibration / "output" / f"checkpoint_{equilibration_steps}.apc"
    measured = measurement_deck(
        base, rule, branch, measurement / "output", checkpoint)
    measured_path = measurement / "input.cfg"
    atomic_text(measured_path, measured)
    measurement_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(measured_path)], measurement / "stdout.txt", measurement / "stderr.txt",
        float(config["measurement_timeout_seconds"]))
    maximum_rss = max(int(equilibration_resources["peak_resident_set_kib"]),
                      int(measurement_resources["peak_resident_set_kib"]))
    if maximum_rss > int(execution["maximum_peak_resident_set_kib"]):
        raise RefinementError("branch exceeded its predeclared memory ceiling")

    output = measurement / "output"
    scripts = Path(__file__).resolve().parent
    energy = output / "energy-budget.json"
    phase = output / "phase-eedf-analysis.json"
    run_analyzer([sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
                  str(output), "--json", str(energy)], "energy analyzer")
    run_analyzer([sys.executable, str(scripts / "analyze_phase_eedf.py"),
                  str(output), "--threshold", "ionization:15.76",
                  "--max-overflow", "0.001", "--json", str(phase)],
                 "phase EEDF analyzer")
    energy_report = json.loads(energy.read_text(encoding="utf-8"))
    phase_report = json.loads(phase.read_text(encoding="utf-8"))
    fields = list(output.glob("fields_*.csv"))
    maximum_field = max(abs(float(row["E"]))
                        for path in fields for row in table(path))
    metadata = json.loads((output / "spatial_average_metadata.json").read_text())
    end_step = equilibration_steps + measurement_steps
    scalars = table(output / "scalars.csv")
    gates = {
        "energy_closure": energy_report.get("passes") is True and abs(float(
            energy_report["relative_closure_residual"])) <=
            float(execution["maximum_relative_energy_residual"]),
        "absolute_field": math.isfinite(maximum_field) and maximum_field <=
            float(execution["maximum_absolute_field_V_m"]),
        "measurement_window": metadata.get("complete") is True and
            metadata.get("start_step") == equilibration_steps + 1 and
            metadata.get("end_step") == end_step and
            int(scalars[0]["step"]) == equilibration_steps and
            int(scalars[-1]["step"]) == end_step,
        "phase_eedf": phase_report.get("passes") is True,
    }
    if not all(gates.values()):
        raise RefinementError("branch failed a solver/diagnostic gate")
    final_checkpoint = output / f"checkpoint_{end_step}.apc"
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "common_state_timestep_refinement_branch",
        "physics_claim": "paired_numerical_sensitivity_evidence_only",
        "branch": branch,
        "rule_sha256": APPROVED_RULE_SHA256,
        "inputs": {
            "solver_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "particle_state_sha256": sha256(state_path),
            "initial_config_sha256": sha256(initial_path),
            "measurement_config_sha256": sha256(measured_path),
        },
        "numerics": {
            "timestep_s": config["timestep_s"],
            "steps_per_rf_cycle": config["steps_per_rf_cycle"],
            "equilibration_steps": equilibration_steps,
            "measurement_steps": measurement_steps,
        },
        "resources": {
            "equilibration": equilibration_resources,
            "measurement": measurement_resources,
            "maximum_peak_resident_set_kib": maximum_rss,
            "available_memory_before_launch_kib": available,
        },
        "gates": gates,
        "all_gates_passed": True,
        "maximum_sampled_absolute_field_V_m": maximum_field,
        "energy_analysis_sha256": sha256(energy),
        "phase_eedf_analysis_sha256": sha256(phase),
        "equilibration_checkpoint_sha256": sha256(checkpoint),
        "final_checkpoint_sha256": sha256(final_checkpoint),
        "output_hashes": {name: sha256(output / name) for name in (
            "spatial_average.csv", "spatial_phase_fields.csv",
            "spatial_phase_moments.csv", "spatial_collision_rate.csv",
            "spatial_phase_collision_rate.csv", "phase_eedf.csv",
            "phase_eedf_moments.csv", "wall_impact_spectrum.csv",
            "wall_impact_spectrum_summary.csv")},
        "claim_boundary": (
            "One branch is not a convergence result; both locked branches "
            "must complete and pass the prospective paired comparison."),
    }
    atomic_json(work / "branch-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("particle_state", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--branch", choices=("baseline_dt", "half_dt"), required=True)
    parser.add_argument("--acknowledge-cost")
    try:
        report = execute(parser.parse_args())
    except (RefinementError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"timestep refinement rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
