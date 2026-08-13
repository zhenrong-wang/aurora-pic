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
APPROVED_PRODUCTION_RULE_SHA256 = (
    "c3b58592dbaa631aa1a3da549828247e7394d5d847e5101b1782f040aa489a8d"
)
APPROVED_STRICT_RULE_SHA256 = (
    "9034bd795245f9662ff02298b9e38ba56938a6383622aff1295823a587f21144"
)
APPROVED_POST_TREND_RULE_SHA256 = (
    "3e6e29425325e9e70557fa1a17545893fae816e49ff351bac14dc4dd82f37b27"
)
APPROVED_DENSITY_ACCELERATED_RULE_SHA256 = (
    "41ec8075a076fea6de40ffed410605a8687547d6890a37b4690bfb8402aedf98"
)
APPROVED_REBRACKETED_RULE_SHA256 = (
    "f8300ee6440dbd539bcc8e4a252ada94e471cfa892ed1c4e769c62fc7a98a878"
)
MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_FIELD_ENERGY_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_PEAK_FIELD_SLOPE_PER_CYCLE = 0.01
MAX_NORMALIZED_IONIZATION_SLOPE_PER_CYCLE = 0.02
MAX_IONIZATION_COEFFICIENT_OF_VARIATION = 0.05
PRODUCTION_WALL_IMPACT_ORIGIN_CYCLE = 64
STRICT_WALL_IMPACT_ORIGIN_CYCLE = 76
STRICT_MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE = 0.001
STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE = 0.05


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
    rule_sha256 = sha256(path)
    if rule_sha256 not in {
        APPROVED_EXTENSION_RULE_SHA256, APPROVED_PRODUCTION_RULE_SHA256,
        APPROVED_STRICT_RULE_SHA256, APPROVED_POST_TREND_RULE_SHA256,
        APPROVED_DENSITY_ACCELERATED_RULE_SHA256,
        APPROVED_REBRACKETED_RULE_SHA256,
    }:
        raise HorizonError("extension rule SHA-256 is not approved")
    rule = json.loads(path.read_text(encoding="utf-8"))
    execution = rule.get("execution_contract")
    stationarity_rule = rule.get("stationarity_contract")
    baseline = rule.get("baseline")
    scopes = {
        APPROVED_EXTENSION_RULE_SHA256:
            "predeclared_aurorapic_edupic_equilibration_extension",
        APPROVED_PRODUCTION_RULE_SHA256:
            "predeclared_aurorapic_edupic_production_equilibration_extension",
        APPROVED_STRICT_RULE_SHA256:
            "predeclared_aurorapic_edupic_strict_source_loss_equilibration",
        APPROVED_POST_TREND_RULE_SHA256:
            "predeclared_aurorapic_edupic_post_trend_strict_equilibration",
        APPROVED_DENSITY_ACCELERATED_RULE_SHA256:
            ("predeclared_aurorapic_edupic_density_accelerated_"
             "strict_equilibration"),
        APPROVED_REBRACKETED_RULE_SHA256:
            "predeclared_aurorapic_edupic_rebracketed_strict_equilibration",
    }
    expected_scope = scopes[rule_sha256]
    if (
        rule.get("schema_version") != 1
        or rule.get("case_id") != "edupic-1.0-default-argon-ccp"
        or rule.get("scope") != expected_scope
        or not isinstance(execution, dict)
        or not isinstance(stationarity_rule, dict)
        or not isinstance(baseline, dict)
    ):
        raise HorizonError("extension rule has the wrong contract identity")
    post_trend = rule_sha256 == APPROVED_POST_TREND_RULE_SHA256
    density_accelerated = (
        rule_sha256 == APPROVED_DENSITY_ACCELERATED_RULE_SHA256)
    rebracketed = rule_sha256 == APPROVED_REBRACKETED_RULE_SHA256
    dense_campaign = density_accelerated or rebracketed
    production = rule_sha256 == APPROVED_PRODUCTION_RULE_SHA256
    strict = rule_sha256 in {
        APPROVED_STRICT_RULE_SHA256, APPROVED_POST_TREND_RULE_SHA256,
        APPROVED_DENSITY_ACCELERATED_RULE_SHA256,
        APPROVED_REBRACKETED_RULE_SHA256}
    expected_execution = {
        "first_cycle": 5 if dense_campaign else (117 if post_trend else
            (81 if strict else (65 if production else 17))),
        "maximum_cycle": (28 if rebracketed else 36) if dense_campaign else (
            148 if post_trend else
            (112 if strict else (96 if production else 64))),
        "cycles_per_block": HARD_MAX_CYCLES_PER_BLOCK,
        "maximum_blocks_per_invocation": 1,
        "serial": True,
        "maximum_wall_seconds_per_cycle": (
            120 if dense_campaign else HARD_TIMEOUT_SECONDS),
        "minimum_available_memory_kib": (
            256 * 1024 if dense_campaign else MIN_AVAILABLE_MEMORY_KIB),
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
            (STRICT_MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE if strict else
             MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE),
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
    if strict:
        expected_stationarity.update({
            "minimum_population_efolding_time_cycles": 1000.0,
            "maximum_electron_source_loss_relative_imbalance":
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
            "maximum_ion_source_loss_relative_imbalance":
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
        })
    if execution != expected_execution or stationarity_rule != expected_stationarity:
        raise HorizonError("extension rule differs from the built-in safety contract")
    if production and rule.get("wall_impact_diagnostic_contract") != {
        "enabled": True,
        "energy_bins": 200,
        "energy_max_eV": 500.0,
        "accumulation_origin_cycle": PRODUCTION_WALL_IMPACT_ORIGIN_CYCLE,
        "required_electrodes": ["left", "right"],
        "required_species": ["electrons", "ions"],
        "exact_macro_count_closure": True,
        "represented_count_and_energy_closure_required": True,
        "overflow_is_reported_separately": True,
        "equilibration_spectra_are_not_comparison_measurements": True,
    }:
        raise HorizonError("production wall-impact contract differs from the built-in contract")
    if post_trend and rule.get("diagnostic_contract") != {
        "wall_impact_spectrum_continues_from_cycle": 112,
        "cycle_113_through_116_measurement_remains_diagnostic_only": True,
        "future_production_measurement_must_start_with_another_fresh_window": True,
    }:
        raise HorizonError("post-trend diagnostic contract differs")
    if dense_campaign and rule.get("diagnostic_contract") != {
        "wall_impact_spectrum_enabled": False,
        "equilibration_outputs_are_not_comparison_measurements": True,
        "future_production_measurement_must_start_with_a_fresh_window": True,
    }:
        raise HorizonError("dense-campaign diagnostic contract differs")
    if args.start_cycle < int(baseline.get("cycle", -1)):
        raise HorizonError("extension cannot precede its frozen baseline")
    if args.start_cycle == int(baseline.get("cycle", -1)) and (
        args.expected_prior_report_sha256.lower()
        != baseline.get(
            "horizon_report_sha256", baseline.get("measurement_report_sha256"))
        or args.expected_input_checkpoint_sha256.lower()
        != baseline.get("checkpoint_sha256")
    ):
        raise HorizonError("first extension block differs from its frozen baseline")
    return int(execution["maximum_cycle"])


def production_rule(path: Path | None) -> bool:
    return (
        path is not None
        and sha256(path.resolve()) in {
            APPROVED_PRODUCTION_RULE_SHA256, APPROVED_STRICT_RULE_SHA256,
            APPROVED_POST_TREND_RULE_SHA256}
    )


def strict_rule(path: Path | None) -> bool:
    return path is not None and sha256(path.resolve()) in {
        APPROVED_STRICT_RULE_SHA256, APPROVED_POST_TREND_RULE_SHA256,
        APPROVED_DENSITY_ACCELERATED_RULE_SHA256,
        APPROVED_REBRACKETED_RULE_SHA256}


def execution_limits(path: Path | None) -> tuple[int, int]:
    if (path is not None and sha256(path.resolve()) in {
            APPROVED_DENSITY_ACCELERATED_RULE_SHA256,
            APPROVED_REBRACKETED_RULE_SHA256}):
        return 120, 256 * 1024
    return HARD_TIMEOUT_SECONDS, MIN_AVAILABLE_MEMORY_KIB


def post_trend_rule(path: Path | None) -> bool:
    return (path is not None and
            sha256(path.resolve()) == APPROVED_POST_TREND_RULE_SHA256)


def wall_impact_origin_cycle(path: Path | None) -> int:
    return (112 if post_trend_rule(path) else
            STRICT_WALL_IMPACT_ORIGIN_CYCLE if strict_rule(path) else
            PRODUCTION_WALL_IMPACT_ORIGIN_CYCLE)


def wall_impact_diagnostic(output: Path, origin_step: int) -> dict[str, object]:
    summary_path = output / "wall_impact_spectrum_summary.csv"
    spectrum_path = output / "wall_impact_spectrum.csv"
    summary = table(summary_path)
    spectrum = table(spectrum_path)
    identities = {
        (row.get("species"), row.get("electrode")) for row in summary
    }
    expected = {
        ("electrons", "left"), ("electrons", "right"),
        ("ions", "left"), ("ions", "right"),
    }
    if identities != expected or len(summary) != 4 or len(spectrum) != 800:
        raise HorizonError("wall-impact diagnostic has the wrong shape or identity")
    for row in summary + spectrum:
        if integer(row, "origin_step", "wall-impact diagnostic") != origin_step:
            raise HorizonError("wall-impact diagnostic has the wrong origin")
    if any(integer(row, "count_closure", "wall-impact summary") != 1
           for row in summary):
        raise HorizonError("wall-impact macro counts do not close")
    result = {
        "origin_step": origin_step,
        "summary_sha256": sha256(summary_path),
        "spectrum_sha256": sha256(spectrum_path),
        "macro_impacts": sum(
            integer(row, "macro_impacts", "wall-impact summary")
            for row in summary
        ),
        "overflow_macro_impacts": sum(
            integer(row, "overflow_macro_impacts", "wall-impact summary")
            for row in summary
        ),
        "maximum_absolute_energy_closure_residual_J_m-2": max(
            abs(float(row["energy_closure_residual"])) for row in summary
        ),
    }
    return result


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
    result = {
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
    return result


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


def stationarity(
    endpoints: list[dict[str, float | int]],
    populations: list[dict[str, object]] | None = None,
    strict: bool = False,
) -> dict[str, object]:
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
    population_limit = (STRICT_MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE
                        if strict else MAX_NORMALIZED_POPULATION_SLOPE_PER_CYCLE)
    gates = {
        "total_population_slope": abs(
            metrics["normalized_total_population_slope_per_cycle"]
        ) <= population_limit,
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
    if strict:
        if not populations or len(populations) != 4:
            raise HorizonError("strict stationarity requires four source-loss records")
        ionization_total = sum(int(item["ionization_pairs"])
                               for item in populations)
        electron_losses = sum(int(item["electron_wall_losses"])
                              for item in populations)
        ion_losses = sum(int(item["ion_wall_losses"])
                         for item in populations)
        if ionization_total <= 0:
            raise HorizonError("strict stationarity requires positive ionization")
        metrics.update({
            "population_efolding_time_cycles": (
                1.0 / abs(metrics["normalized_total_population_slope_per_cycle"])
                if metrics["normalized_total_population_slope_per_cycle"] != 0.0
                else math.inf),
            "ionization_pairs_in_block": ionization_total,
            "electron_wall_losses_in_block": electron_losses,
            "ion_wall_losses_in_block": ion_losses,
            "electron_source_loss_relative_imbalance":
                abs(ionization_total - electron_losses) / ionization_total,
            "ion_source_loss_relative_imbalance":
                abs(ionization_total - ion_losses) / ionization_total,
        })
        gates.update({
            "electron_source_loss_balance": metrics[
                "electron_source_loss_relative_imbalance"] <=
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
            "ion_source_loss_balance": metrics[
                "ion_source_loss_relative_imbalance"] <=
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
        })
    result = {
        "window_cycles": [int(item["cycle"]) for item in endpoints],
        "thresholds": {
            "maximum_absolute_normalized_total_population_slope_per_cycle":
                population_limit,
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
    if strict:
        result["thresholds"].update({
            "minimum_population_efolding_time_cycles": 1000.0,
            "maximum_electron_source_loss_relative_imbalance":
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
            "maximum_ion_source_loss_relative_imbalance":
                STRICT_MAX_SOURCE_LOSS_RELATIVE_IMBALANCE,
        })
    return result


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
    window = report.get("window")
    if (isinstance(window, dict) and
            isinstance(window.get("end_cycle"), int) and
            report.get("all_gates_passed") is True):
        return int(window["end_cycle"])
    stages = report.get("stages")
    safety = report.get("safety_decision")
    if (isinstance(stages, list) and stages and
            isinstance(safety, dict) and safety.get("passes") is True and
            all(isinstance(stage, dict) and stage.get("passes") is True
                for stage in stages)):
        end_step = stages[-1].get("end_step")
        if isinstance(end_step, int) and end_step % STEPS_PER_CYCLE == 0:
            return end_step // STEPS_PER_CYCLE
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
    timeout_seconds, minimum_available_memory_kib = execution_limits(
        args.extension_rule)
    collect_wall_impacts = production_rule(args.extension_rule)
    strict_campaign = strict_rule(args.extension_rule)
    if (args.start_cycle < 4 or
            args.start_cycle + args.cycles > maximum_end_cycle):
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
    prior_checkpoint_sha256 = prior_report.get("final_checkpoint_sha256")
    if prior_checkpoint_sha256 is None:
        prior_checkpoint_sha256 = prior_report["stages"][-1].get(
            "output_checkpoint_sha256")
    checkpoint = input_output / (
        f"checkpoint_{args.start_cycle * STEPS_PER_CYCLE}.apc"
    )
    if (
        sha256(checkpoint) != args.expected_input_checkpoint_sha256.lower()
        or prior_checkpoint_sha256 != sha256(checkpoint)
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
        if available < minimum_available_memory_kib:
            raise HorizonError("available memory is below the hard launch floor")
        stage = work / f"cycle-{cycle + 1:04d}"
        output = stage / "output"
        stage.mkdir()
        deck_path = stage / "input.cfg"
        atomic_text(deck_path, stage_deck(
            base, cycle, output, checkpoint,
            wall_impact_spectrum=collect_wall_impacts,
        ))
        resources = run_process(
            solver_command(executable, deck_path),
            stage / "stdout.txt", stage / "stderr.txt",
            timeout_seconds=timeout_seconds,
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
        if collect_wall_impacts:
            result["wall_impact_diagnostic"] = wall_impact_diagnostic(
                output,
                wall_impact_origin_cycle(args.extension_rule) * STEPS_PER_CYCLE,
            )
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
    screen = stationarity(
        endpoints, [stage["population"] for stage in stages], strict_campaign)
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
        "wall_impact_accumulation_is_equilibration_only": collect_wall_impacts,
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
