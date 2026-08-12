#!/usr/bin/env python3
"""Run a predeclared, bounded warm-state transient follow-up."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from run_aurorapic_edupic_pilot import (
    STEPS_PER_CYCLE, analyze_stage, atomic_json, atomic_text,
    available_memory_kib, global_integer, run_analyzer, run_process, sha256,
)
from run_aurorapic_initialization_ab import (
    initial_deck, restart_deck, set_species_value, summarize,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_AN_INITIALIZATION_DIAGNOSTIC"


class FollowupError(RuntimeError):
    pass


def configured_base(text: str, rule: dict[str, object]) -> str:
    state = rule["initial_state"]
    result = text
    for species in ("electrons", "ions"):
        result = set_species_value(
            result, species, "particles", int(state[species]))
        result = set_species_value(
            result, species, "weight", float(state["macro_weight"]))
    return result


def evaluate(stages: list[dict[str, object]],
             rule: dict[str, object]) -> dict[str, object]:
    first = stages[0]
    population = first["population"]
    observed = {
        "first_cycle_ionization_pairs": population["ionization_pairs"],
        "first_cycle_electron_wall_losses": population["electron_wall_losses"],
        "first_cycle_total_growth_factor": population["total_growth_factor"],
        "maximum_sampled_absolute_field_V_m": max(
            stage["maximum_sampled_absolute_electric_field_V_m"]
            for stage in stages),
    }
    thresholds = rule["decision_rule"]["prospective_thresholds"]
    gates = {
        "all_numerical_hard_gates": all(stage["passes"] for stage in stages),
        "bounded_first_cycle_ionization":
            observed["first_cycle_ionization_pairs"] <=
            thresholds["maximum_first_cycle_ionization_pairs"],
        "nonvanishing_first_cycle_electron_loss":
            observed["first_cycle_electron_wall_losses"] >=
            thresholds["minimum_first_cycle_electron_wall_losses"],
        "bounded_first_cycle_population_growth":
            observed["first_cycle_total_growth_factor"] <=
            thresholds["maximum_first_cycle_total_growth_factor"],
        "bounded_transient_field":
            observed["maximum_sampled_absolute_field_V_m"] <=
            thresholds["maximum_sampled_absolute_field_V_m"],
    }
    return {"observed": observed, "gates": gates,
            "passes": all(gates.values())}


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise FollowupError(
            "follow-up requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    state_path = args.particle_state.resolve()
    rule_path = args.rule.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise FollowupError(f"refusing to overwrite study directory: {work}")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    execution = rule["execution_contract"]
    state = rule["initial_state"]
    if (rule.get("scope") !=
            "predeclared_field_consistent_initialization_followup" or
            execution["cycles"] != 4 or
            execution["steps_per_cycle"] != STEPS_PER_CYCLE):
        raise FollowupError("rule identity or cycle contract is invalid")
    for path, expected, label in (
        (executable, execution["solver_sha256"], "solver"),
        (base_path, execution["base_deck_sha256"], "base deck"),
        (state_path, state["particle_state_sha256"], "particle state"),
    ):
        if sha256(path) != expected:
            raise FollowupError(f"{label} SHA-256 does not match the rule")
    base = base_path.read_text(encoding="utf-8")
    if global_integer(base, "steps") != STEPS_PER_CYCLE:
        raise FollowupError("base deck is not the one-cycle contract")
    base = configured_base(base, rule)
    total_cap = 2 * global_integer(base, "max_particles_per_species")
    scripts = Path(__file__).resolve().parent
    work.mkdir(parents=True)
    previous = {name: int(state[name]) for name in ("electrons", "ions")}
    checkpoint = None
    stages = []
    for cycle in range(4):
        available = available_memory_kib()
        if available < int(execution["minimum_available_memory_kib"]):
            raise FollowupError("available memory is below the rule floor")
        stage_dir = work / f"cycle-{cycle + 1:04d}"
        output = stage_dir / "output"
        stage_dir.mkdir()
        deck = (initial_deck(
                    base, output, state_path,
                    int(state["particle_state_signature"]))
                if cycle == 0 else
                restart_deck(base, cycle, output, checkpoint))
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
        result.update({
            "resources": resources,
            "available_memory_before_launch_kib": available,
            "deck_sha256": sha256(deck_path),
            "energy_analysis_sha256": sha256(energy),
            "spatial_analysis_sha256": sha256(spatial),
        })
        checkpoint = output / f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc"
        result["output_checkpoint_sha256"] = sha256(checkpoint)
        atomic_json(stage_dir / "stage-report.json", result)
        if not result["passes"]:
            raise FollowupError(f"cycle {cycle + 1} failed a hard gate")
        previous = {
            "electrons": result["population"]["final_electrons"],
            "ions": result["population"]["final_ions"],
        }
        stages.append(result)
    decision = evaluate(stages, rule)
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "completed_field_consistent_initialization_followup",
        "physics_claim": "none_initialization_diagnostic_only",
        "rule_sha256": sha256(rule_path),
        "inputs": {
            "solver_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "particle_state_sha256": sha256(state_path),
            "particle_state_signature": state["particle_state_signature"],
        },
        "stages": stages,
        "metrics": summarize(stages),
        "prospective_decision": decision,
        "admissible_for_density_diagnostic": decision["passes"],
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
    except (FollowupError, OSError, ValueError, KeyError,
            json.JSONDecodeError) as error:
        print(f"initialization follow-up rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
