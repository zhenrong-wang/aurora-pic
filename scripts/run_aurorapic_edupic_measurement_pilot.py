#!/usr/bin/env python3
"""Run the predeclared four-cycle AuroraPIC eduPIC measurement pilot."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from extend_aurorapic_edupic_horizon import wall_impact_diagnostic
from run_aurorapic_edupic_pilot import (
    MAX_ABSOLUTE_FIELD_V_M, MAX_RELATIVE_ENERGY_RESIDUAL,
    MIN_AVAILABLE_MEMORY_KIB, PilotError, atomic_json, atomic_text,
    available_memory_kib, insert_global, integer, run_analyzer, run_process,
    set_global, sha256, table, STEPS_PER_CYCLE,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_FRESH_AURORAPIC_MEASUREMENT_PILOT"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULES = {
    "8b662fc493df3229fc123d034cd55c602628d0dcc1b2e447e119d1a6647a9387":
        (76, 240, MIN_AVAILABLE_MEMORY_KIB),
    "f8d4a05fc5fb359cc5808fb35d08d39af3392b8d5108a2474716ee57a2c845be":
        (112, 240, MIN_AVAILABLE_MEMORY_KIB),
    "6da07382d4723853d976bcfb00d89c766b468d17e21502f8e885996628014986":
        (20, 600, 256 * 1024),
    "f9a5e33683986432f3c2050515ad6e7de02316b14dd35c14d3e6a05694e5a216":
        (24, 600, 256 * 1024),
}
MEASUREMENT_CYCLES = 4
class MeasurementPilotError(RuntimeError):
    pass


def deck(base: str, output: Path, checkpoint: Path,
         start_cycle: int = 76) -> str:
    start_step = start_cycle * STEPS_PER_CYCLE
    end_step = (start_cycle + MEASUREMENT_CYCLES) * STEPS_PER_CYCLE
    values = {
        "steps": end_step, "output_interval": 100,
        "output_dir": output, "spatial_average_start_step": start_step + 1,
        "spatial_average_end_step": end_step,
        "spatial_average_rf_cycles": MEASUREMENT_CYCLES,
        "spatial_average_phase_bins": 16,
        "checkpoint_interval": MEASUREMENT_CYCLES * STEPS_PER_CYCLE,
        "runtime_backend": "serial", "runtime_threads": 1,
    }
    result = base
    for key, value in values.items():
        result = set_global(result, key, str(value))
    additions = {
        "spatial_average_reset_on_restart": "true",
        "phase_eedf": "true", "phase_eedf_species": "electrons",
        "phase_eedf_energy_bins": "2000",
        "phase_eedf_energy_max": "500",
        "phase_eedf_regions": "full_gap:0:0.025",
        "wall_impact_spectrum": "true",
        "wall_impact_reset_on_restart": "true",
        "wall_impact_energy_bins": "200",
        "wall_impact_energy_max": "500",
        "restart_path": str(checkpoint),
    }
    for key, value in additions.items():
        result = insert_global(result, key, value)
    return result


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MeasurementPilotError("measurement pilot cost was not acknowledged")
    rule = args.rule.resolve()
    rule_hash = sha256(rule)
    if rule_hash not in APPROVED_RULES:
        raise MeasurementPilotError("measurement rule SHA-256 is not approved")
    start_cycle, timeout_seconds, minimum_available_memory_kib = (
        APPROVED_RULES[rule_hash])
    rule_data = json.loads(rule.read_text(encoding="utf-8"))
    execution = rule_data.get("execution_contract", {})
    if execution != {
        "first_cycle": start_cycle + 1,
        "last_cycle": start_cycle + MEASUREMENT_CYCLES,
        "measurement_cycles": MEASUREMENT_CYCLES,
        "serial": True,
        "maximum_wall_seconds": timeout_seconds,
        "minimum_available_memory_kib": minimum_available_memory_kib,
        "maximum_absolute_field_V_m": MAX_ABSOLUTE_FIELD_V_M,
        "maximum_relative_energy_residual": MAX_RELATIVE_ENERGY_RESIDUAL,
    }:
        raise MeasurementPilotError("measurement execution contract differs")
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    checkpoint = args.checkpoint.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise MeasurementPilotError(f"refusing to overwrite {work}")
    for path, expected, name in (
        (executable, args.expected_binary_sha256, "binary"),
        (base_path, args.expected_base_deck_sha256, "base deck"),
        (checkpoint, args.expected_checkpoint_sha256, "checkpoint"),
    ):
        if sha256(path) != expected.lower():
            raise MeasurementPilotError(f"{name} SHA-256 differs")
    available = available_memory_kib()
    if available < minimum_available_memory_kib:
        raise MeasurementPilotError("available memory is below the launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck_path = work / "input.cfg"
    atomic_text(deck_path, deck(base_path.read_text(encoding="utf-8"),
                                output, checkpoint, start_cycle))
    resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(deck_path),
    ], work / "stdout.txt", work / "stderr.txt", timeout_seconds)
    scripts = Path(__file__).resolve().parent
    energy = output / "energy-budget.json"
    spatial = output / "spatial-collision.json"
    phase = output / "phase-eedf-analysis.json"
    run_analyzer([sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
                  str(output), "--json", str(energy)], "energy analyzer")
    run_analyzer([sys.executable, str(scripts / "analyze_1d_spatial_collision.py"),
                  str(output), "--boundary", "dirichlet", "--json", str(spatial)],
                 "spatial analyzer")
    run_analyzer([sys.executable, str(scripts / "analyze_phase_eedf.py"),
                  str(output), "--threshold", "ionization:15.76",
                  "--max-overflow", "0.001", "--json", str(phase)],
                 "phase EEDF analyzer")
    scalars = table(output / "scalars.csv")
    start_step = start_cycle * STEPS_PER_CYCLE
    end_step = (start_cycle + MEASUREMENT_CYCLES) * STEPS_PER_CYCLE
    if (integer(scalars[0], "step", "measurement start") != start_step or
            integer(scalars[-1], "step", "measurement end") != end_step):
        raise MeasurementPilotError("measurement scalar window differs")
    maximum_field = max(abs(float(row["E"]))
                        for path in output.glob("fields_*.csv")
                        for row in table(path))
    energy_report = json.loads(energy.read_text(encoding="utf-8"))
    phase_report = json.loads(phase.read_text(encoding="utf-8"))
    metadata = json.loads((output / "spatial_average_metadata.json").read_text())
    gates = {
        "energy_closure": energy_report.get("passes") is True and
            abs(float(energy_report["relative_closure_residual"])) <=
            MAX_RELATIVE_ENERGY_RESIDUAL,
        "absolute_field": math.isfinite(maximum_field) and
            maximum_field <= MAX_ABSOLUTE_FIELD_V_M,
        "spatial_window_complete": metadata.get("complete") is True and
            metadata.get("start_step") == start_step + 1 and
            metadata.get("end_step") == end_step and
            metadata.get("samples") == MEASUREMENT_CYCLES * STEPS_PER_CYCLE,
        "phase_eedf": phase_report.get("passes") is True,
    }
    impact = wall_impact_diagnostic(output, start_step)
    final_checkpoint = output / f"checkpoint_{end_step}.apc"
    if not all(gates.values()):
        raise MeasurementPilotError("measurement pilot failed a diagnostic gate")
    report = {
        "schema_version": 1, "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "fresh_window_aurorapic_measurement_pilot",
        "physics_claim": "none_descriptive_pilot_only",
        "window": {"start_cycle": start_cycle,
                   "end_cycle": start_cycle + MEASUREMENT_CYCLES,
                   "measurement_cycles": MEASUREMENT_CYCLES,
                   "equilibration_statistics_excluded": True},
        "inputs": {"rule_sha256": rule_hash,
                   "binary_sha256": sha256(executable),
                   "base_deck_sha256": sha256(base_path),
                   "checkpoint_sha256": sha256(checkpoint)},
        "resources": resources, "available_memory_before_launch_kib": available,
        "gates": gates, "all_gates_passed": True,
        "maximum_sampled_absolute_field_V_m": maximum_field,
        "energy_analysis_sha256": sha256(energy),
        "spatial_analysis_sha256": sha256(spatial),
        "phase_eedf_analysis_sha256": sha256(phase),
        "wall_impact_diagnostic": impact,
        "final_checkpoint_sha256": sha256(final_checkpoint),
        "output_hashes": {
            name: sha256(output / name) for name in (
                "spatial_average.csv", "spatial_phase_fields.csv",
                "spatial_phase_moments.csv",
                "spatial_collision_rate.csv",
                "spatial_phase_collision_rate.csv", "phase_eedf.csv",
                "phase_eedf_moments.csv", "wall_impact_spectrum.csv",
                "wall_impact_spectrum_summary.csv")},
        "claim_boundary": "This four-cycle pilot checks fresh-window diagnostics; it does not establish statistical convergence or cross-code agreement.",
    }
    atomic_json(work / "measurement-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-base-deck-sha256", required=True)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--acknowledge-cost")
    try:
        report = execute(parser.parse_args())
    except (MeasurementPilotError, PilotError, OSError, ValueError,
            json.JSONDecodeError) as error:
        print(f"AuroraPIC measurement pilot rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
