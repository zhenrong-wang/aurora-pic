#!/usr/bin/env python3
"""Run a predeclared control/+density warm-state initialization bracket."""

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


class BracketError(RuntimeError):
    pass


def configured_base(text: str, branch: dict[str, object]) -> str:
    result = text
    for species in ("electrons", "ions"):
        result = set_species_value(
            result, species, "particles", int(branch[species]))
        result = set_species_value(
            result, species, "weight", float(branch["macro_weight"]))
    return result


def branch_gates(stages: list[dict[str, object]],
                 rule: dict[str, object]) -> dict[str, object]:
    first_population = stages[0]["population"]
    thresholds = rule["decision_rule"]["prospective_safety_thresholds"]
    observed = {
        "first_cycle_ionization_pairs":
            first_population["ionization_pairs"],
        "first_cycle_electron_wall_losses":
            first_population["electron_wall_losses"],
        "first_cycle_total_growth_factor":
            first_population["total_growth_factor"],
        "maximum_sampled_absolute_field_V_m": max(
            stage["maximum_sampled_absolute_electric_field_V_m"]
            for stage in stages),
    }
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


def select_branch(branches: dict[str, dict[str, object]],
                  rule: dict[str, object]) -> dict[str, object]:
    decision = rule["decision_rule"]
    control_metrics = branches["control"]["metrics"]
    control_score = max(abs(control_metrics[key]) for key in (
        "electron_source_loss_relative_imbalance",
        "ion_source_loss_relative_imbalance"))
    candidates: dict[str, dict[str, object]] = {}
    for name in decision["candidate_branches"]:
        result = branches[name]
        metrics = result["metrics"]
        score = max(abs(metrics[key]) for key in (
            "electron_source_loss_relative_imbalance",
            "ion_source_loss_relative_imbalance"))
        imbalance_improvement = control_score - score
        slope_improvement = (
            abs(control_metrics["normalized_total_population_slope_per_cycle"])
            - abs(metrics["normalized_total_population_slope_per_cycle"]))
        eligible = (
            result["safety_decision"]["passes"] and
            imbalance_improvement >=
            decision["minimum_absolute_imbalance_improvement_vs_control"] and
            slope_improvement > 0.0)
        candidates[name] = {
            "maximum_absolute_source_loss_imbalance": score,
            "absolute_imbalance_improvement_vs_control":
                imbalance_improvement,
            "absolute_population_slope_improvement_vs_control":
                slope_improvement,
            "eligible": eligible,
        }
    eligible_names = [name for name, result in candidates.items()
                      if result["eligible"]]
    selected = None
    if eligible_names:
        best_score = min(candidates[name][
            "maximum_absolute_source_loss_imbalance"]
            for name in eligible_names)
        near_best = [name for name in eligible_names
                     if candidates[name][
                         "maximum_absolute_source_loss_imbalance"] <=
                     best_score + decision["parsimony_tolerance"]]
        selected = min(
            near_best,
            key=lambda name: rule["branches"][name]["added_pairs"])
    return {
        "control_maximum_absolute_source_loss_imbalance": control_score,
        "candidates": candidates,
        "selected_branch": selected,
        "density_acceleration_supported": selected is not None,
    }


def parse_states(values: list[str]) -> dict[str, Path]:
    states = {}
    for value in values:
        name, separator, path = value.partition("=")
        if not separator or not name or not path or name in states:
            raise BracketError("each --state must be a unique NAME=PATH")
        states[name] = Path(path).resolve()
    return states


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise BracketError(
            "bracket requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    rule_path = args.rule.resolve()
    work = args.work_dir.resolve()
    states = parse_states(args.state)
    if work.exists():
        raise BracketError(f"refusing to overwrite study directory: {work}")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    execution = rule["execution_contract"]
    branch_contracts = rule["branches"]
    if (rule.get("scope") != "predeclared_cycle128_density_bracket" or
            execution["cycles_per_branch"] != 4 or
            execution["steps_per_cycle"] != STEPS_PER_CYCLE or
            execution["branch_order"] != list(branch_contracts) or
            list(branch_contracts) != ["control", "plus25", "plus50"] or
            set(states) != set(branch_contracts)):
        raise BracketError("rule identity, branch order, or cycle contract invalid")
    if sha256(Path(__file__).resolve()) != rule["provenance"]["runner_sha256"]:
        raise BracketError("runner SHA-256 does not match the rule")
    for path, expected, label in (
        (executable, execution["solver_sha256"], "solver"),
        (base_path, execution["base_deck_sha256"], "base deck"),
    ):
        if sha256(path) != expected:
            raise BracketError(f"{label} SHA-256 does not match the rule")
    for name, contract in branch_contracts.items():
        if sha256(states[name]) != contract["particle_state_sha256"]:
            raise BracketError(f"{name} particle-state SHA-256 mismatch")
    base = base_path.read_text(encoding="utf-8")
    if (global_integer(base, "steps") != STEPS_PER_CYCLE or
            global_integer(base, "seed") != execution["seed"]):
        raise BracketError("base deck cycle or seed contract is invalid")
    total_cap = 2 * global_integer(base, "max_particles_per_species")
    scripts = Path(__file__).resolve().parent
    work.mkdir(parents=True)
    branches = {}
    for name, contract in branch_contracts.items():
        branch_dir = work / name
        branch_dir.mkdir()
        branch_base = configured_base(base, contract)
        previous = {species: int(contract[species])
                    for species in ("electrons", "ions")}
        checkpoint = None
        stages = []
        for cycle in range(execution["cycles_per_branch"]):
            available = available_memory_kib()
            if available < int(execution["minimum_available_memory_kib"]):
                raise BracketError("available memory is below the rule floor")
            stage_dir = branch_dir / f"cycle-{cycle + 1:04d}"
            output = stage_dir / "output"
            stage_dir.mkdir()
            deck = (initial_deck(
                        branch_base, output, states[name],
                        int(contract["particle_state_signature"]))
                    if cycle == 0 else
                    restart_deck(branch_base, cycle, output, checkpoint))
            deck_path = stage_dir / "input.cfg"
            atomic_text(deck_path, deck)
            resources = run_process(
                [str(executable), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(deck_path)],
                stage_dir / "stdout.txt", stage_dir / "stderr.txt",
                timeout_seconds=float(
                    execution["maximum_wall_seconds_per_cycle"]))
            energy = output / "energy-budget.json"
            spatial = output / "spatial-collision.json"
            run_analyzer([
                sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
                str(output), "--json", str(energy)], "energy analyzer")
            run_analyzer([
                sys.executable, str(scripts / "analyze_1d_spatial_collision.py"),
                str(output), "--boundary", "dirichlet", "--json",
                str(spatial)], "spatial analyzer")
            result = analyze_stage(output, cycle, previous, total_cap)
            result.update({
                "resources": resources,
                "available_memory_before_launch_kib": available,
                "deck_sha256": sha256(deck_path),
                "energy_analysis_sha256": sha256(energy),
                "spatial_analysis_sha256": sha256(spatial),
            })
            checkpoint = output / (
                f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc")
            result["output_checkpoint_sha256"] = sha256(checkpoint)
            atomic_json(stage_dir / "stage-report.json", result)
            if not result["passes"]:
                raise BracketError(f"{name} cycle {cycle + 1} failed hard gate")
            previous = {
                "electrons": result["population"]["final_electrons"],
                "ions": result["population"]["final_ions"],
            }
            stages.append(result)
        branches[name] = {
            "stages": stages,
            "metrics": summarize(stages),
            "safety_decision": branch_gates(stages, rule),
        }
        atomic_json(branch_dir / "branch-report.json", branches[name])
    decision = select_branch(branches, rule)
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "completed_cycle128_density_bracket",
        "physics_claim": "none_initialization_diagnostic_only",
        "rule_sha256": sha256(rule_path),
        "inputs": {
            "solver_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "particle_state_sha256": {
                name: sha256(path) for name, path in states.items()},
        },
        "branches": branches,
        "prospective_decision": decision,
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
    parser.add_argument("rule", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--state", action="append", default=[],
                        metavar="NAME=PATH")
    parser.add_argument("--acknowledge-cost")
    try:
        report = execute(parser.parse_args())
    except (BracketError, OSError, ValueError, KeyError,
            json.JSONDecodeError) as error:
        print(f"density bracket rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
