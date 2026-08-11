#!/usr/bin/env python3
"""Run the predeclared eduPIC-argon initialization-density A/B diagnostic."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import statistics
import sys

from run_aurorapic_edupic_pilot import (
    STEPS_PER_CYCLE, analyze_stage, atomic_json, atomic_text,
    available_memory_kib, global_integer, insert_global, run_analyzer,
    run_process, set_global, sha256,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_AN_INITIALIZATION_DIAGNOSTIC"


class StudyError(RuntimeError):
    pass


def set_species_value(text: str, species: str, key: str, value: object) -> str:
    section = re.compile(
        rf"(?ms)(^\[species\.{re.escape(species)}\]\s*$)(.*?)(?=^\[|\Z)")
    match = section.search(text)
    if match is None:
        raise StudyError(f"base deck is missing species '{species}'")
    body = match.group(2)
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    replaced, count = pattern.subn(rf"\g<1>{value}", body, count=1)
    if count != 1:
        raise StudyError(f"species '{species}' does not contain '{key}'")
    return text[:match.start(2)] + replaced + text[match.end(2):]


def branch_base(base: str, rule: dict[str, object], branch: str) -> str:
    state = rule["common_initial_state"]
    branch_config = rule["branches"][branch]
    result = base
    for species in ("electrons", "ions"):
        result = set_species_value(
            result, species, "particles", int(state[species]))
        result = set_species_value(
            result, species, "weight",
            float(branch_config[
                "electron_macro_weight"
                if species == "electrons" else
                "ion_macro_weight"]))
    return result


def initial_deck(base: str, output: Path, state_path: Path,
                 signature: int) -> str:
    result = base
    for key, value in {
        "steps": STEPS_PER_CYCLE,
        "output_interval": 100,
        "output_dir": output,
        "spatial_average_start_step": 1,
        "spatial_average_end_step": STEPS_PER_CYCLE,
        "spatial_average_rf_cycles": 1,
        "spatial_average_phase_bins": 16,
        "checkpoint_interval": STEPS_PER_CYCLE,
        "runtime_backend": "serial",
        "runtime_threads": 1,
    }.items():
        result = set_global(result, key, str(value))
    result = insert_global(result, "initial_state_path", str(state_path))
    return insert_global(result, "initial_state_signature", str(signature))


def restart_deck(base: str, cycle: int, output: Path, checkpoint: Path) -> str:
    result = base
    for key, value in {
        "steps": (cycle + 1) * STEPS_PER_CYCLE,
        "output_interval": 100,
        "output_dir": output,
        "spatial_average_start_step": cycle * STEPS_PER_CYCLE + 1,
        "spatial_average_end_step": (cycle + 1) * STEPS_PER_CYCLE,
        "spatial_average_rf_cycles": 1,
        "spatial_average_phase_bins": 16,
        "checkpoint_interval": STEPS_PER_CYCLE,
        "runtime_backend": "serial",
        "runtime_threads": 1,
    }.items():
        result = set_global(result, key, str(value))
    result = insert_global(result, "spatial_average_reset_on_restart", "true")
    return insert_global(result, "restart_path", str(checkpoint))


def normalized_slope(values: list[float]) -> float:
    coordinates = list(range(len(values)))
    x_mean = statistics.mean(coordinates)
    y_mean = statistics.mean(values)
    slope = math.fsum(
        (x - x_mean) * (y - y_mean)
        for x, y in zip(coordinates, values)
    ) / math.fsum((x - x_mean) ** 2 for x in coordinates)
    return slope / y_mean


def summarize(stages: list[dict[str, object]]) -> dict[str, object]:
    first = stages[0]["population"]
    populations = [first["initial_electrons"] + first["initial_ions"]]
    field_energies = [stages[0]["state"]["initial_field_energy_J_m-2"]]
    peak_fields = []
    ionizations = []
    for stage in stages:
        pop = stage["population"]
        populations.append(pop["final_electrons"] + pop["final_ions"])
        field_energies.append(stage["state"]["final_field_energy_J_m-2"])
        peak_fields.append(stage["maximum_sampled_absolute_electric_field_V_m"])
        ionizations.append(pop["ionization_pairs"])
    ionization_total = sum(ionizations)
    electron_losses = sum(s["population"]["electron_wall_losses"] for s in stages)
    ion_losses = sum(s["population"]["ion_wall_losses"] for s in stages)
    return {
        "normalized_total_population_slope_per_cycle": normalized_slope(populations),
        "electron_source_loss_relative_imbalance":
            1.0 - electron_losses / ionization_total,
        "ion_source_loss_relative_imbalance":
            1.0 - ion_losses / ionization_total,
        "normalized_field_energy_slope_per_cycle": normalized_slope(field_energies),
        "normalized_peak_field_slope_per_cycle": normalized_slope(peak_fields),
        "normalized_ionization_slope_per_cycle": normalized_slope(ionizations),
        "ionization_coefficient_of_variation":
            statistics.pstdev(ionizations) / statistics.mean(ionizations),
        "ionization_pairs": ionization_total,
        "electron_wall_losses": electron_losses,
        "ion_wall_losses": ion_losses,
        "initial_total_particles": populations[0],
        "final_total_particles": populations[-1],
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise StudyError("study requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    state_path = args.particle_state.resolve()
    rule_path = args.rule.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise StudyError(f"refusing to overwrite study directory: {work}")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    execution = rule["execution_contract"]
    state_contract = rule["common_initial_state"]
    if (rule.get("scope") != "predeclared_initialization_density_ab_diagnostic"
            or execution["cycles_per_branch"] != 4):
        raise StudyError("rule identity or cycle contract is invalid")
    for path, expected, label in (
        (executable, execution["solver_sha256"], "solver"),
        (base_path, execution["base_deck_sha256"], "base deck"),
        (state_path, state_contract["particle_state_sha256"], "particle state"),
    ):
        if sha256(path) != expected:
            raise StudyError(f"{label} SHA-256 does not match the rule")
    base = base_path.read_text(encoding="utf-8")
    if global_integer(base, "steps") != STEPS_PER_CYCLE:
        raise StudyError("base deck is not the one-cycle contract")
    total_cap = 2 * global_integer(base, "max_particles_per_species")
    scripts = Path(__file__).resolve().parent
    work.mkdir(parents=True)
    branches = {}
    for branch in ("control", "warm"):
        configured_base = branch_base(base, rule, branch)
        branch_dir = work / branch
        branch_dir.mkdir()
        previous = {
            "electrons": int(state_contract["electrons"]),
            "ions": int(state_contract["ions"]),
        }
        checkpoint = None
        stages = []
        for cycle in range(4):
            if available_memory_kib() < int(execution["minimum_available_memory_kib"]):
                raise StudyError("available memory is below the rule floor")
            stage_dir = branch_dir / f"cycle-{cycle + 1:04d}"
            output = stage_dir / "output"
            stage_dir.mkdir()
            deck = (initial_deck(
                        configured_base, output, state_path,
                        int(state_contract["particle_state_signature"]))
                    if cycle == 0 else
                    restart_deck(configured_base, cycle, output, checkpoint))
            deck_path = stage_dir / "input.cfg"
            atomic_text(deck_path, deck)
            resources = run_process(
                [str(executable), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(deck_path)],
                stage_dir / "stdout.txt", stage_dir / "stderr.txt",
                timeout_seconds=float(execution["maximum_wall_seconds_per_cycle"]))
            energy = output / "energy-budget.json"
            spatial = output / "spatial-collision.json"
            run_analyzer([
                sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
                str(output), "--json", str(energy)], "energy analyzer")
            run_analyzer([
                sys.executable, str(scripts / "analyze_1d_spatial_collision.py"),
                str(output), "--boundary", "dirichlet", "--json", str(spatial)],
                "spatial analyzer")
            result = analyze_stage(output, cycle, previous, total_cap)
            result["resources"] = resources
            result["deck_sha256"] = sha256(deck_path)
            result["energy_analysis_sha256"] = sha256(energy)
            result["spatial_analysis_sha256"] = sha256(spatial)
            checkpoint = output / f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc"
            result["output_checkpoint_sha256"] = sha256(checkpoint)
            atomic_json(stage_dir / "stage-report.json", result)
            if not result["passes"]:
                raise StudyError(f"{branch} cycle {cycle + 1} failed a hard gate")
            previous = {
                "electrons": result["population"]["final_electrons"],
                "ions": result["population"]["final_ions"],
            }
            stages.append(result)
        branches[branch] = {"stages": stages, "metrics": summarize(stages)}
    control = branches["control"]["metrics"]
    warm = branches["warm"]["metrics"]
    comparisons = {
        "electron_imbalance_reduction":
            control["electron_source_loss_relative_imbalance"] -
            warm["electron_source_loss_relative_imbalance"],
        "ion_imbalance_reduction":
            control["ion_source_loss_relative_imbalance"] -
            warm["ion_source_loss_relative_imbalance"],
        "population_slope_reduction":
            control["normalized_total_population_slope_per_cycle"] -
            warm["normalized_total_population_slope_per_cycle"],
    }
    supported = (
        comparisons["electron_imbalance_reduction"] >= 0.10 and
        comparisons["ion_imbalance_reduction"] >= 0.10 and
        comparisons["population_slope_reduction"] > 0.0)
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "completed_initialization_density_ab_diagnostic",
        "physics_claim": "none_initialization_diagnostic_only",
        "rule_sha256": sha256(rule_path),
        "inputs": {
            "solver_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "particle_state_sha256": sha256(state_path),
            "particle_state_signature": state_contract["particle_state_signature"],
        },
        "branches": branches,
        "comparison": comparisons,
        "warm_start_hypothesis_supported": supported,
        "equilibrium_established": False,
        "measurement_authorized": False,
        "claim_boundary": rule["claim_boundary"],
    }
    atomic_json(work / "study-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("particle_state", type=Path)
    parser.add_argument("rule", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--acknowledge-cost")
    try:
        report = execute(parser.parse_args())
    except (StudyError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"initialization A/B study rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
