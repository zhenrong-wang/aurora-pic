#!/usr/bin/env python3
"""Run one locked matched-half-step AuroraPIC threshold-traffic branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, insert_global,
    run_process, set_global, sha256,
)
from run_aurorapic_initialization_ab import set_species_value
from run_aurorapic_ionizing_tail_block import analyze_output


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_MATCHED_HALF_STEP_THRESHOLD_RUN"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256S = {
    "960ae236fd92f32d24b382236e0fe2264370ad6072fcba681a140a7757a74fbf",
}


class MatchedHalfStepError(RuntimeError):
    pass


def remove_global(text: str, key: str) -> str:
    pattern = re.compile(rf"(?m)^\s*{re.escape(key)}\s*=.*\n?")
    result, count = pattern.subn("", text, count=1)
    if count != 1:
        raise MatchedHalfStepError(
            f"base deck does not contain exactly one {key!r}")
    return result


def common_deck(base: str, rule: dict[str, object], state: dict[str, object],
                output: Path) -> str:
    execution = rule["execution_contract"]
    result = base
    result = insert_global(
        result, "collision_velocity_sampling", "leapfrog_half_step")
    for species in ("electrons", "ions"):
        result = set_species_value(
            result, species, "particles", int(state[species]))
    for key, value in {
        "output_interval": execution["output_interval_steps"],
        "output_dir": output,
        "runtime_backend": "serial",
        "runtime_threads": 1,
        "seed": state["seed"],
    }.items():
        result = set_global(result, key, str(value))
    return result


def equilibration_deck(base: str, rule: dict[str, object],
                       state: dict[str, object], output: Path,
                       state_path: Path) -> str:
    execution = rule["execution_contract"]
    end = int(execution["equilibration_steps"])
    result = common_deck(base, rule, state, output)
    for key, value in {
        "steps": end,
        "spatial_average_interval": execution["equilibration_sampling_interval_steps"],
        "spatial_average_start_step": 1,
        "spatial_average_end_step": end,
        "spatial_average_rf_cycles": execution["equilibration_cycles"],
        "spatial_average_phase_bins": execution["equilibration_phase_bins"],
        "checkpoint_interval": end,
        "initial_state_path": state_path,
        "initial_state_signature": state["particle_state_signature"],
    }.items():
        result = set_global(result, key, str(value))
    return result


def measurement_deck(base: str, rule: dict[str, object],
                     state: dict[str, object], output: Path,
                     checkpoint: Path) -> str:
    execution = rule["execution_contract"]
    diagnostic = rule["diagnostic_contract"]
    start, end = int(execution["start_step"]), int(execution["end_step"])
    regions = ",".join(
        f"{item['name']}:{item['x_min_m']}:{item['x_max_m']}"
        for item in diagnostic["regions"])
    result = common_deck(base, rule, state, output)
    result = remove_global(result, "initial_state_path")
    result = remove_global(result, "initial_state_signature")
    for key, value in {
        "steps": end,
        "spatial_average_interval": diagnostic["spatial_average_interval_steps"],
        "spatial_average_start_step": start + 1,
        "spatial_average_end_step": end,
        "spatial_average_rf_cycles": execution["measurement_cycles"],
        "spatial_average_phase_bins": diagnostic["phase_bins"],
        "checkpoint_interval": end,
    }.items():
        result = set_global(result, key, str(value))
    for key, value in {
        "spatial_average_reset_on_restart": "true",
        "spatial_average_sampling_order": diagnostic["sampling_order"],
        "phase_eedf": "true",
        "phase_eedf_species": diagnostic["phase_eedf_species"],
        "phase_eedf_energy_bins": diagnostic["energy_bins"],
        "phase_eedf_energy_max": diagnostic["energy_max_eV"],
        "phase_eedf_regions": regions,
        "phase_eedf_tail_threshold": diagnostic["tail_threshold_eV"],
        "phase_eedf_history": "true",
        "wall_impact_spectrum": "true",
        "wall_impact_reset_on_restart": "true",
        "wall_impact_energy_bins": 200,
        "wall_impact_energy_max": 500.0,
        "restart_path": checkpoint,
    }.items():
        result = insert_global(result, key, str(value).lower()
                               if isinstance(value, bool) else str(value))
    return result


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MatchedHalfStepError("matched-half-step cost was not acknowledged")
    rule_path = args.rule.resolve()
    rule_hash = sha256(rule_path)
    if rule_hash not in APPROVED_RULE_SHA256S:
        raise MatchedHalfStepError("rule is not an approved matched-half-step campaign")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    matches = [item for item in rule["locked_initial_states"]
               if item["id"] == args.initial_state_id]
    if len(matches) != 1:
        raise MatchedHalfStepError("initial-state id is not locked by the rule")
    state = matches[0]
    executable = args.executable.resolve()
    base_path = args.base_config.resolve()
    state_path = args.particle_state.resolve()
    manifest_path = args.export_manifest.resolve()
    for path, expected, label in (
        (executable, state["solver_sha256"], "solver"),
        (base_path, state["base_config_sha256"], "base config"),
        (state_path, state["particle_state_sha256"], "particle state"),
        (manifest_path, state["export_manifest_sha256"], "export manifest"),
    ):
        if sha256(path) != expected:
            raise MatchedHalfStepError(f"locked {label} differs")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("particle_state_sha256") != state["particle_state_sha256"] or
            manifest.get("particle_state_signature") !=
                state["particle_state_signature"] or
            manifest.get("source_checkpoint_sha256") !=
                state["source_checkpoint_sha256"]):
        raise MatchedHalfStepError("particle-state export provenance differs")
    execution = rule["execution_contract"]
    if (int(execution["start_step"]) != int(execution["equilibration_steps"]) or
            int(execution["end_step"]) - int(execution["start_step"]) !=
                int(execution["measurement_cycles"]) *
                int(execution["steps_per_cycle"])):
        raise MatchedHalfStepError("execution horizon is inconsistent")
    work = args.work_dir.resolve()
    if work.exists():
        raise MatchedHalfStepError(f"refusing to overwrite {work}")
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    if available_memory < int(execution["minimum_available_memory_kib"]):
        raise MatchedHalfStepError("available memory is below the launch floor")
    if available_disk < int(execution["minimum_available_disk_kib"]):
        raise MatchedHalfStepError("available disk is below the launch floor")
    work.mkdir(parents=True)
    base = base_path.read_text(encoding="utf-8")
    equilibration = work / "equilibration"
    measurement = work / "measurement"
    equilibration.mkdir()
    measurement.mkdir()
    equilibration_output = equilibration / "output"
    equilibration_input = equilibration / "input.cfg"
    atomic_text(equilibration_input, equilibration_deck(
        base, rule, state, equilibration_output, state_path))
    equilibration_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(equilibration_input)], equilibration / "stdout.txt",
        equilibration / "stderr.txt", float(execution["timeout_seconds"]))
    checkpoint = equilibration_output / (
        f"checkpoint_{int(execution['equilibration_steps'])}.apc")
    if not checkpoint.is_file():
        raise MatchedHalfStepError("equilibration checkpoint is missing")
    measurement_output = measurement / "output"
    measurement_input = measurement / "input.cfg"
    atomic_text(measurement_input, measurement_deck(
        base, rule, state, measurement_output, checkpoint))
    measurement_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(measurement_input)], measurement / "stdout.txt",
        measurement / "stderr.txt", float(execution["timeout_seconds"]))
    peak = max(int(equilibration_resources["peak_resident_set_kib"]),
               int(measurement_resources["peak_resident_set_kib"]))
    analysis_resources = {
        **measurement_resources, "peak_resident_set_kib": peak}
    result = analyze_output(measurement_output, rule, analysis_resources)
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_matched_half_step_threshold_branch",
        "rule_sha256": rule_hash,
        "inputs": {
            "initial_state_id": state["id"],
            "solver_sha256": sha256(executable),
            "base_config_sha256": sha256(base_path),
            "particle_state_sha256": sha256(state_path),
            "export_manifest_sha256": sha256(manifest_path),
            "source_checkpoint_sha256": state["source_checkpoint_sha256"],
            "equilibration_deck_sha256": sha256(equilibration_input),
            "measurement_deck_sha256": sha256(measurement_input),
        },
        "algorithm_contract": {
            "collision_velocity_sampling": "leapfrog_half_step",
            "portable_state_velocity_staggering": "time_centered",
            "half_step_rebuilt_self_consistently": True,
            "equilibration_cycles_after_rebuild":
                execution["equilibration_cycles"],
        },
        "resources": {
            "equilibration": equilibration_resources,
            "measurement": measurement_resources,
            "maximum_peak_resident_set_kib": peak,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk,
        },
        "equilibration_checkpoint_sha256": sha256(checkpoint),
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    })
    atomic_json(work / "branch-report.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_config", type=Path)
    parser.add_argument("particle_state", type=Path)
    parser.add_argument("export_manifest", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--initial-state-id", required=True)
    parser.add_argument("--acknowledge-cost", required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (MatchedHalfStepError, PilotError, OSError, ValueError,
            KeyError) as error:
        parser.error(str(error))
    print(json.dumps({
        "report": str(args.work_dir.resolve() / "branch-report.json"),
        "all_gates_passed": result["all_gates_passed"],
        "sampling": result["sampling"],
        "resources": result["resources"],
    }, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
