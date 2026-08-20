#!/usr/bin/env python3
"""Run one predeclared AuroraPIC source/loss-stationarity block."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, finite, integer,
    run_process, set_global, sha256, table,
)


RULE_SHA256 = "e2371b9ac8df48091f02ee8bacc7b04b6e81f80c23ae1146b0802bc8b1c70c42"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_ONE_BOUNDED_BALANCE_BLOCK"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"


class BalanceBlockError(RuntimeError):
    pass


def normalized_slope(values: list[float]) -> float:
    if len(values) < 2:
        raise ValueError("slope needs at least two samples")
    mean = statistics.fmean(values)
    if mean == 0.0:
        return math.inf if any(value != 0.0 for value in values) else 0.0
    center = 0.5 * (len(values) - 1)
    denominator = math.fsum((index - center) ** 2
                            for index in range(len(values)))
    slope = math.fsum((index - center) * (value - mean)
                      for index, value in enumerate(values)) / denominator
    return slope / abs(mean)


def build_deck(base: str, output: Path, checkpoint: Path,
               end_step: int, interval: int) -> str:
    values = {
        "steps": str(end_step),
        "output_interval": str(interval),
        "output_dir": str(output),
        "spatial_average": "false",
        "phase_eedf": "false",
        "wall_impact_spectrum": "false",
        "checkpoint_interval": str(end_step),
        "runtime_backend": "serial",
        "runtime_threads": "1",
        "restart_path": str(checkpoint),
    }
    result = base
    for key, value in values.items():
        result = set_global(result, key, value)
    disabled_options = (
        "spatial_average_interval", "spatial_average_start_step",
        "spatial_average_end_step", "spatial_average_rf_frequency",
        "spatial_average_rf_cycles", "spatial_average_phase_bins",
        "spatial_average_reset_on_restart", "spatial_average_sampling_order",
        "phase_eedf_species", "phase_eedf_energy_bins",
        "phase_eedf_energy_max", "phase_eedf_regions",
        "wall_impact_reset_on_restart", "wall_impact_energy_bins",
        "wall_impact_energy_max",
    )
    for key in disabled_options:
        result, count = re.subn(
            rf"(?m)^\s*{re.escape(key)}\s*=.*\n?", "", result, count=1)
        if count != 1:
            raise BalanceBlockError(
                f"base config does not contain exactly one {key!r}")
    return result


def cumulative_delta(rows: list[dict[str, str]], key: str) -> int:
    return integer(rows[-1], key, "final counter") - integer(
        rows[0], key, "initial counter")


def field_manifest(output: Path) -> dict[str, object]:
    paths = sorted(output.glob("fields_*.csv"))
    if not paths:
        raise BalanceBlockError("block produced no field snapshots")
    digest = hashlib.sha256()
    maximum = 0.0
    for path in paths:
        file_hash = sha256(path)
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\n")
        maximum = max(maximum, *(abs(finite(row, "E", path.name))
                                 for row in table(path)))
    return {"files": len(paths), "manifest_sha256": digest.hexdigest(),
            "maximum_absolute_field_V_m": maximum}


def analyze_output(output: Path, start_step: int, end_step: int,
                   steps_per_cycle: int, rule: dict[str, object],
                   resources: dict[str, object]) -> dict[str, object]:
    scalars = table(output / "scalars.csv")
    collisions = table(output / "collisions.csv")
    boundaries = table(output / "boundary_losses.csv")
    for name, rows in (("scalars", scalars), ("collisions", collisions),
                       ("boundaries", boundaries)):
        if integer(rows[0], "step", name) != start_step or integer(
                rows[-1], "step", name) != end_step:
            raise BalanceBlockError(f"{name} has the wrong block endpoints")
    scalar_by_step = {integer(row, "step", "scalars"): row for row in scalars}
    collision_by_step = {
        integer(row, "step", "collisions"): row for row in collisions}
    endpoint_steps = list(range(start_step, end_step + 1, steps_per_cycle))
    if (any(step not in scalar_by_step for step in endpoint_steps) or
            any(step not in collision_by_step for step in endpoint_steps)):
        raise BalanceBlockError("whole-cycle endpoints are incomplete")

    ionization_key = "cumulative_collisions_electron_mcc.ionization"
    created = cumulative_delta(collisions, ionization_key)
    electron_losses = sum(cumulative_delta(boundaries, key) for key in (
        "absorbed_left_count_electrons", "absorbed_right_count_electrons"))
    ion_losses = sum(cumulative_delta(boundaries, key) for key in (
        "absorbed_left_count_ions", "absorbed_right_count_ions"))
    electron_change = integer(
        scalars[-1], "live_particles_electrons", "final scalar") - integer(
            scalars[0], "live_particles_electrons", "initial scalar")
    ion_change = integer(
        scalars[-1], "live_particles_ions", "final scalar") - integer(
            scalars[0], "live_particles_ions", "initial scalar")
    if electron_change != created - electron_losses:
        raise BalanceBlockError("electron particle ledger does not close")
    if ion_change != created - ion_losses:
        raise BalanceBlockError("ion particle ledger does not close")
    if created <= 0:
        raise BalanceBlockError("block has no ionization events")

    electron_populations = [float(scalar_by_step[step][
        "live_particles_electrons"]) for step in endpoint_steps]
    ion_populations = [float(scalar_by_step[step][
        "live_particles_ions"]) for step in endpoint_steps]
    field_energies = [float(scalar_by_step[step]["field_energy"])
                      for step in endpoint_steps]
    ionization_per_cycle = [
        integer(collision_by_step[high], ionization_key, "cycle collision") -
        integer(collision_by_step[low], ionization_key, "cycle collision")
        for low, high in zip(endpoint_steps, endpoint_steps[1:])
    ]
    ionization_mean = statistics.fmean(ionization_per_cycle)
    field = field_manifest(output)
    limits = rule["prospective_block_stationarity"]
    metrics = {
        "electron_source_loss_relative_imbalance":
            (created - electron_losses) / created,
        "ion_source_loss_relative_imbalance": (created - ion_losses) / created,
        "normalized_electron_population_slope_per_cycle":
            normalized_slope(electron_populations),
        "normalized_ion_population_slope_per_cycle":
            normalized_slope(ion_populations),
        "normalized_field_energy_slope_per_cycle":
            normalized_slope(field_energies),
        "ionization_count_coefficient_of_variation":
            statistics.pstdev(ionization_per_cycle) / ionization_mean,
    }
    stationarity_gates = {
        "electron_source_loss_balance": abs(metrics[
            "electron_source_loss_relative_imbalance"]) <= float(limits[
                "maximum_absolute_electron_source_loss_relative_imbalance"]),
        "ion_source_loss_balance": abs(metrics[
            "ion_source_loss_relative_imbalance"]) <= float(limits[
                "maximum_absolute_ion_source_loss_relative_imbalance"]),
        "electron_population_slope": abs(metrics[
            "normalized_electron_population_slope_per_cycle"]) <= float(limits[
                "maximum_absolute_normalized_electron_population_slope_per_cycle"]),
        "ion_population_slope": abs(metrics[
            "normalized_ion_population_slope_per_cycle"]) <= float(limits[
                "maximum_absolute_normalized_ion_population_slope_per_cycle"]),
        "field_energy_slope": abs(metrics[
            "normalized_field_energy_slope_per_cycle"]) <= float(limits[
                "maximum_absolute_normalized_field_energy_slope_per_cycle"]),
        "ionization_variation": metrics[
            "ionization_count_coefficient_of_variation"] <= float(limits[
                "maximum_ionization_count_coefficient_of_variation"]),
    }
    execution = rule["execution_contract"]
    maximum_particles = max(int(row["live_particles"]) for row in scalars)
    hard_gates = {
        "particle_ledger_closure": True,
        "finite_scalar_history": all(
            math.isfinite(float(row[key])) for row in scalars
            for key in ("time", "kinetic_energy", "field_energy", "total_energy")),
        "particle_cap": maximum_particles <= int(
            execution["maximum_total_particles"]),
        "absolute_field": field["maximum_absolute_field_V_m"] <= float(
            execution["maximum_absolute_field_V_m"]),
        "resident_memory": int(resources["peak_resident_set_kib"]) <= int(
            execution["maximum_peak_resident_set_kib"]),
    }
    checkpoint = output / f"checkpoint_{end_step}.apc"
    if not checkpoint.is_file():
        raise BalanceBlockError("final checkpoint is missing")
    return {
        "window": {"start_step": start_step, "end_step": end_step,
                   "cycles": len(endpoint_steps) - 1},
        "particle_ledger": {
            "ionization_macro_events": created,
            "electron_wall_loss_macro_events": electron_losses,
            "ion_wall_loss_macro_events": ion_losses,
            "electron_live_particle_change": electron_change,
            "ion_live_particle_change": ion_change,
            "exact_species_closure": True,
        },
        "cycle_endpoint_populations": {
            "electrons": [int(value) for value in electron_populations],
            "ions": [int(value) for value in ion_populations],
        },
        "ionization_macro_events_per_cycle": ionization_per_cycle,
        "metrics": metrics,
        "stationarity_gates": stationarity_gates,
        "stationarity_block_passed": all(stationarity_gates.values()),
        "hard_safety_gates": hard_gates,
        "all_hard_safety_gates_passed": all(hard_gates.values()),
        "field_diagnostics": field,
        "maximum_live_particles": maximum_particles,
        "final_checkpoint_sha256": sha256(checkpoint),
        "output_sha256": {
            name: sha256(output / name) for name in (
                "scalars.csv", "collisions.csv", "boundary_losses.csv")
        },
    }


def validate_inputs(args: argparse.Namespace, rule: dict[str, object]) -> tuple[int, int]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise BalanceBlockError("missing exact bounded-run acknowledgement")
    execution = rule["execution_contract"]
    if not 1 <= args.block_index <= int(execution["maximum_blocks"]):
        raise BalanceBlockError("block index is outside the campaign")
    initial = rule["locked_initial_state"]
    if sha256(args.executable.resolve()) != initial["solver_sha256"]:
        raise BalanceBlockError("solver differs from the locked binary")
    if sha256(args.base_config.resolve()) != initial["base_config_sha256"]:
        raise BalanceBlockError("base config differs from the locked input")
    prior_hash = sha256(args.prior_report.resolve())
    checkpoint_hash = sha256(args.checkpoint.resolve())
    if args.block_index == 1:
        if (prior_hash != initial["prior_report_sha256"] or
                checkpoint_hash != initial["checkpoint_sha256"]):
            raise BalanceBlockError("first block does not continue the locked state")
    else:
        prior = json.loads(args.prior_report.read_text(encoding="utf-8"))
        if (prior.get("scope") != "aurorapic_source_loss_stationarity_block" or
                prior.get("rule_sha256") != RULE_SHA256 or
                prior.get("block_index") != args.block_index - 1 or
                prior.get("all_hard_safety_gates_passed") is not True or
                prior.get("final_checkpoint_sha256") != checkpoint_hash):
            raise BalanceBlockError("continuation chain differs from prior block")
    cycles = int(execution["cycles_per_block"])
    steps_per_cycle = int(execution["steps_per_cycle"])
    start_step = int(initial["step"]) + (args.block_index - 1) * cycles * steps_per_cycle
    return start_step, start_step + cycles * steps_per_cycle


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve()
    if sha256(rule_path) != RULE_SHA256:
        raise BalanceBlockError("stationarity rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    start_step, end_step = validate_inputs(args, rule)
    work = args.work_dir.resolve()
    if work.exists():
        raise BalanceBlockError(f"refusing to overwrite {work}")
    available = available_memory_kib()
    execution = rule["execution_contract"]
    if available < int(execution["minimum_available_memory_kib"]):
        raise BalanceBlockError("available memory is below the launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck = work / "input.cfg"
    atomic_text(deck, build_deck(
        args.base_config.read_text(encoding="utf-8"), output,
        args.checkpoint.resolve(), end_step,
        int(execution["output_interval_steps"])))
    resources = run_process([
        str(args.executable.resolve()), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(deck),
    ], work / "stdout.txt", work / "stderr.txt",
        timeout_seconds=float(execution["timeout_seconds_per_block"]))
    result = analyze_output(
        output, start_step, end_step, int(execution["steps_per_cycle"]),
        rule, resources)
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_source_loss_stationarity_block",
        "block_index": args.block_index,
        "rule_sha256": RULE_SHA256,
        "inputs": {
            "solver_sha256": sha256(args.executable.resolve()),
            "base_config_sha256": sha256(args.base_config.resolve()),
            "prior_report_sha256": sha256(args.prior_report.resolve()),
            "input_checkpoint_sha256": sha256(args.checkpoint.resolve()),
            "deck_sha256": sha256(deck),
        },
        "resources": {**resources,
                      "available_memory_before_launch_kib": available},
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
    parser.add_argument("--block-index", type=int, required=True)
    parser.add_argument("--acknowledge-cost", required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (BalanceBlockError, PilotError, OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps({
        "report": str(args.work_dir.resolve() / "block-report.json"),
        "hard_safety_passed": result["all_hard_safety_gates_passed"],
        "stationarity_passed": result["stationarity_block_passed"],
        "metrics": result["metrics"],
    }, indent=2, sort_keys=True))
    return 0 if result["all_hard_safety_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
