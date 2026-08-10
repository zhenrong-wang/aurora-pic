#!/usr/bin/env python3
"""Extend the AuroraPIC eduPIC startup in one gated stationarity block."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics
import sys

from run_aurorapic_edupic_pilot import (
    HARD_TIMEOUT_SECONDS,
    MAX_ABSOLUTE_FIELD_V_M,
    MAX_PARTICLE_GROWTH_FACTOR,
    MAX_RELATIVE_ENERGY_RESIDUAL,
    MAX_SPATIAL_PHASE_RESIDUAL_J_M2,
    MAX_TOTAL_PARTICLE_CAP_FRACTION,
    MIN_AVAILABLE_MEMORY_KIB,
    PilotError,
    STEPS_PER_CYCLE,
    analyze_stage,
    atomic_json,
    atomic_text,
    available_memory_kib,
    global_integer,
    integer,
    run_analyzer,
    run_process,
    sha256,
    stage_deck,
    table,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_AURORAPIC_HORIZON_BLOCK"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
HARD_MAX_CYCLES_PER_BLOCK = 4
HARD_MAX_END_CYCLE = 16
APPROVED_EXTENSION_RULE_SHA256 = (
    "14584b56a28bdd54e1865f919cfbdf87d26f932bd7aff7731e07ca8d91aa78a9"
)
MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_FIELD_ENERGY_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_PEAK_FIELD_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_IONIZATION_SLOPE_PER_CYCLE = 0.02
MAX_IONIZATION_COEFFICIENT_OF_VARIATION = 0.05


class HorizonError(RuntimeError):
    pass


def solver_command(executable: Path, deck: Path) -> list[str]:
    return [
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT, str(deck)
    ]


def authorized_end_cycle(args: argparse.Namespace) -> int:
    if args.extension_rule is None:
        return HARD_MAX_END_CYCLE
    path = args.extension_rule.resolve()
    if sha256(path) != APPROVED_EXTENSION_RULE_SHA256:
        raise HorizonError("extension rule SHA-256 is not approved")
    rule = json.loads(path.read_text(encoding="utf-8"))
    execution = rule.get("execution_contract")
    stationarity_rule = rule.get("stationarity_contract")
    baseline = rule.get("baseline")
    if (
        rule.get("schema_version") != 1
        or rule.get("case_id") != "edupic-1.0-default-argon-ccp"
        or rule.get("scope")
        != "predeclared_aurorapic_edupic_equilibration_extension"
        or not isinstance(execution, dict)
        or not isinstance(stationarity_rule, dict)
        or not isinstance(baseline, dict)
    ):
        raise HorizonError("extension rule has the wrong contract identity")
    expected_execution = {
        "first_cycle": 17,
        "maximum_cycle": 64,
        "cycles_per_block": HARD_MAX_CYCLES_PER_BLOCK,
        "maximum_blocks_per_invocation": 1,
        "serial": True,
        "maximum_wall_seconds_per_cycle": HARD_TIMEOUT_SECONDS,
        "minimum_available_memory_kib": MIN_AVAILABLE_MEMORY_KIB,
        "maximum_particle_growth_factor_per_cycle": MAX_PARTICLE_GROWTH_FACTOR,
        "maximum_total_particle_cap_fraction": MAX_TOTAL_PARTICLE_CAP_FRACTION,
        "maximum_absolute_field_V_m": MAX_ABSOLUTE_FIELD_V_M,
        "maximum_relative_energy_residual": MAX_RELATIVE_ENERGY_RESIDUAL,
        "maximum_spatial_phase_residual_J_m-2":
            MAX_SPATIAL_PHASE_RESIDUAL_J_M2,
    }
    expected_stationarity = {
        "window_cycles": 5,
        "maximum_absolute_normalized_total_population_slope_per_cycle":
            MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE,
        "maximum_absolute_normalized_field_energy_slope_per_cycle":
            MAX_NORMALIZED_FIELD_ENERGY_SLOPE_PER_CYCLE,
        "maximum_absolute_normalized_peak_field_slope_per_cycle":
            MAX_NORMALIZED_PEAK_FIELD_SLOPE_PER_CYCLE,
        "maximum_absolute_normalized_ionization_slope_per_cycle":
            MAX_NORMALIZED_IONIZATION_SLOPE_PER_CYCLE,
        "maximum_ionization_coefficient_of_variation":
            MAX_IONIZATION_COEFFICIENT_OF_VARIATION,
        "consecutive_passing_blocks_required_before_measurement": 2,
    }
    if execution != expected_execution or stationarity_rule != expected_stationarity:
        raise HorizonError("extension rule differs from the built-in safety contract")
    if args.start_cycle < int(baseline.get("cycle", -1)):
        raise HorizonError("extension cannot precede its frozen baseline")
    if args.start_cycle == int(baseline.get("cycle", -1)) and (
        args.expected_prior_report_sha256.lower()
        != baseline.get("horizon_report_sha256")
        or args.expected_input_checkpoint_sha256.lower()
        != baseline.get("checkpoint_sha256")
    ):
        raise HorizonError("first extension block differs from its frozen baseline")
    return int(execution["maximum_cycle"])


def endpoint(output: Path, cycle: int) -> dict[str, float | int]:
    scalars = table(output / "scalars.csv")
    collisions = table(output / "collisions.csv")
    expected_step = cycle * STEPS_PER_CYCLE
    if integer(scalars[-1], "step", "horizon endpoint") != expected_step:
        raise HorizonError("input endpoint does not match the declared cycle")
    ionization_key = "cumulative_collisions_electron_mcc.ionization"
    ionizations = integer(collisions[-1], ionization_key, "final collisions")
    if len(collisions) > 1:
        ionizations -= integer(collisions[0], ionization_key, "initial collisions")
    maximum_field = 0.0
    for path in output.glob("fields_*.csv"):
        for row in table(path):
            try:
                value = abs(float(row["E"]))
            except (KeyError, ValueError) as error:
                raise HorizonError("field output contains invalid E") from error
            if not math.isfinite(value):
                raise HorizonError("field output contains non-finite E")
            maximum_field = max(maximum_field, value)
    return {
        "cycle": cycle,
        "electrons": integer(
            scalars[-1], "live_particles_electrons", "horizon endpoint"
        ),
        "ions": integer(scalars[-1], "live_particles_ions", "horizon endpoint"),
        "total_particles": integer(
            scalars[-1], "live_particles", "horizon endpoint"
        ),
        "field_energy_J_m-2": float(scalars[-1]["field_energy"]),
        "maximum_sampled_absolute_field_V_m": maximum_field,
        "ionization_pairs_in_cycle": ionizations,
    }


def normalized_slope(values: list[float]) -> float:
    mean = statistics.fmean(values)
    if mean == 0.0:
        return math.inf if any(value != 0.0 for value in values) else 0.0
    x_mean = 0.5 * (len(values) - 1)
    denominator = sum((index - x_mean) ** 2 for index in range(len(values)))
    slope = sum(
        (index - x_mean) * (value - mean)
        for index, value in enumerate(values)
    ) / denominator
    return slope / abs(mean)


def stationarity(endpoints: list[dict[str, float | int]]) -> dict[str, object]:
    if len(endpoints) < 5:
        raise HorizonError("stationarity block requires an input plus four cycles")
    totals = [float(item["total_particles"]) for item in endpoints]
    fields = [float(item["field_energy_J_m-2"]) for item in endpoints]
    maxima = [
        float(item["maximum_sampled_absolute_field_V_m"])
        for item in endpoints
    ]
    ionizations = [
        float(item["ionization_pairs_in_cycle"]) for item in endpoints[1:]
    ]
    metrics = {
        "normalized_total_population_slope_per_cycle": normalized_slope(totals),
        "normalized_field_energy_slope_per_cycle": normalized_slope(fields),
        "normalized_peak_field_slope_per_cycle": normalized_slope(maxima),
        "normalized_ionization_slope_per_cycle": normalized_slope(ionizations),
        "ionization_coefficient_of_variation": (
            statistics.pstdev(ionizations) / statistics.fmean(ionizations)
            if statistics.fmean(ionizations) != 0.0 else math.inf
        ),
    }
    gates = {
        "total_population_slope": abs(
            metrics["normalized_total_population_slope_per_cycle"]
        ) <= MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE,
        "field_energy_slope": abs(
            metrics["normalized_field_energy_slope_per_cycle"]
        ) <= MAX_NORMALIZED_FIELD_ENERGY_SLOPE_PER_CYCLE,
        "peak_field_slope": abs(
            metrics["normalized_peak_field_slope_per_cycle"]
        ) <= MAX_NORMALIZED_PEAK_FIELD_SLOPE_PER_CYCLE,
        "ionization_slope": abs(
            metrics["normalized_ionization_slope_per_cycle"]
        ) <= MAX_NORMALIZED_IONIZATION_SLOPE_PER_CYCLE,
        "ionization_variation": metrics[
            "ionization_coefficient_of_variation"
        ] <= MAX_IONIZATION_COEFFICIENT_OF_VARIATION,
    }
    return {
        "window_cycles": [int(item["cycle"]) for item in endpoints],
        "thresholds": {
            "maximum_absolute_normalized_total_population_slope_per_cycle":
                MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE,
            "maximum_absolute_normalized_field_energy_slope_per_cycle":
                MAX_NORMALIZED_FIELD_ENERGY_SLOPE_PER_CYCLE,
            "maximum_absolute_normalized_peak_field_slope_per_cycle":
                MAX_NORMALIZED_PEAK_FIELD_SLOPE_PER_CYCLE,
            "maximum_absolute_normalized_ionization_slope_per_cycle":
                MAX_NORMALIZED_IONIZATION_SLOPE_PER_CYCLE,
            "maximum_ionization_coefficient_of_variation":
                MAX_IONIZATION_COEFFICIENT_OF_VARIATION,
        },
        "metrics": metrics,
        "gates": gates,
        "passed": all(gates.values()),
        "claim_boundary": (
            "An internal endpoint-trend screen is not proof of physical "
            "equilibrium or cross-code agreement."
        ),
    }


def report_end_cycle(report: dict[str, object]) -> int:
    completed = report.get("completed_through_cycle")
    if isinstance(completed, int) and report.get("all_gates_passed") is True:
        return completed
    block = report.get("block")
    if (
        isinstance(block, dict)
        and isinstance(block.get("end_cycle"), int)
        and block.get("hard_safety_gates_passed") is True
    ):
        return int(block["end_cycle"])
    raise HorizonError("prior report is not a completed safe pilot or horizon block")


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise HorizonError(
            "horizon block requires --acknowledge-cost " + ACKNOWLEDGEMENT
        )
    if args.cycles != HARD_MAX_CYCLES_PER_BLOCK:
        raise HorizonError(
            f"horizon block requires exactly {HARD_MAX_CYCLES_PER_BLOCK} cycles"
        )
    maximum_end_cycle = authorized_end_cycle(args)
    if args.start_cycle < 4 or args.start_cycle + args.cycles > maximum_end_cycle:
        raise HorizonError("requested horizon exceeds its built-in cycle bounds")
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    input_output = args.input_output_dir.resolve()
    prior_report_path = args.prior_report.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise HorizonError(f"refusing to overwrite horizon directory: {work}")
    if sha256(executable) != args.expected_binary_sha256.lower():
        raise HorizonError("solver SHA-256 does not match the locked value")
    if sha256(base_path) != args.expected_base_deck_sha256.lower():
        raise HorizonError("base deck SHA-256 does not match the locked value")
    if sha256(prior_report_path) != args.expected_prior_report_sha256.lower():
        raise HorizonError("prior report SHA-256 does not match the locked value")
    prior_report = json.loads(prior_report_path.read_text(encoding="utf-8"))
    if report_end_cycle(prior_report) != args.start_cycle:
        raise HorizonError("prior report does not end at the requested start cycle")
    prior_stationarity_streak = int(
        prior_report.get("consecutive_stationary_blocks", 0)
    )
    if args.extension_rule is not None and prior_stationarity_streak >= 2:
        raise HorizonError("equilibration is complete; start a measurement campaign")
    prior_stage = prior_report["stages"][-1]
    checkpoint = input_output / (
        f"checkpoint_{args.start_cycle * STEPS_PER_CYCLE}.apc"
    )
    if (
        sha256(checkpoint) != args.expected_input_checkpoint_sha256.lower()
        or prior_stage.get("output_checkpoint_sha256") != sha256(checkpoint)
    ):
        raise HorizonError("input checkpoint is not the locked prior output")
    base = base_path.read_text(encoding="utf-8")
    total_cap = 2 * global_integer(base, "max_particles_per_species")
    initial_endpoint = endpoint(input_output, args.start_cycle)
    previous = {
        "electrons": int(initial_endpoint["electrons"]),
        "ions": int(initial_endpoint["ions"]),
    }
    endpoints = [initial_endpoint]
    stages = []
    scripts = Path(__file__).resolve().parent
    work.mkdir(parents=True)
    for cycle in range(args.start_cycle, args.start_cycle + args.cycles):
        available = available_memory_kib()
        if available < MIN_AVAILABLE_MEMORY_KIB:
            raise HorizonError("available memory is below the hard launch floor")
        stage = work / f"cycle-{cycle + 1:04d}"
        output = stage / "output"
        stage.mkdir()
        deck_path = stage / "input.cfg"
        atomic_text(deck_path, stage_deck(base, cycle, output, checkpoint))
        resources = run_process(
            solver_command(executable, deck_path),
            stage / "stdout.txt", stage / "stderr.txt",
        )
        energy = output / "energy-budget.json"
        spatial = output / "spatial-collision.json"
        run_analyzer([
            sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
            str(output), "--json", str(energy),
        ], "energy analyzer")
        run_analyzer([
            sys.executable, str(scripts / "analyze_1d_spatial_collision.py"),
            str(output), "--boundary", "dirichlet", "--json", str(spatial),
        ], "spatial collision analyzer")
        result = analyze_stage(output, cycle, previous, total_cap)
        next_checkpoint = output / (
            f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc"
        )
        result.update({
            "resources": resources,
            "available_memory_before_launch_kib": available,
            "input_checkpoint_sha256": sha256(checkpoint),
            "output_checkpoint_sha256": sha256(next_checkpoint),
            "deck_sha256": sha256(deck_path),
            "energy_analysis_sha256": sha256(energy),
            "spatial_analysis_sha256": sha256(spatial),
        })
        atomic_json(stage / "stage-report.json", result)
        if not result["passes"]:
            raise HorizonError(f"cycle {cycle + 1} failed a hard safety gate")
        stages.append(result)
        endpoints.append(endpoint(output, cycle + 1))
        previous = {
            "electrons": result["population"]["final_electrons"],
            "ions": result["population"]["final_ions"],
        }
        checkpoint = next_checkpoint
    screen = stationarity(endpoints)
    stationary_streak = (
        prior_stationarity_streak + 1 if screen["passed"] else 0
    )
    report = {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_aurorapic_horizon_stationarity_block",
        "physics_claim": "none_internal_stationarity_screen_only",
        "inputs": {
            "executable_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "prior_report_sha256": sha256(prior_report_path),
            "input_checkpoint_sha256": args.expected_input_checkpoint_sha256.lower(),
            "extension_rule_sha256": (
                sha256(args.extension_rule.resolve())
                if args.extension_rule is not None else None
            ),
        },
        "block": {
            "start_cycle": args.start_cycle,
            "end_cycle": args.start_cycle + args.cycles,
            "cycles": args.cycles,
            "hard_safety_gates_passed": True,
        },
        "endpoints": endpoints,
        "stages": stages,
        "stationarity_screen": screen,
        "consecutive_stationary_blocks": stationary_streak,
        "comparison_measurement_eligible": (
            args.extension_rule is not None and stationary_streak >= 2
        ),
        "production_launch_authorized": False,
    }
    atomic_json(work / "horizon-report.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("input_output_dir", type=Path)
    parser.add_argument("prior_report", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--start-cycle", type=int, required=True)
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-base-deck-sha256", required=True)
    parser.add_argument("--expected-prior-report-sha256", required=True)
    parser.add_argument("--expected-input-checkpoint-sha256", required=True)
    parser.add_argument("--extension-rule", type=Path)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        report = execute(parse_args())
    except (HorizonError, PilotError, OSError, json.JSONDecodeError) as error:
        print(f"AuroraPIC eduPIC horizon rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
