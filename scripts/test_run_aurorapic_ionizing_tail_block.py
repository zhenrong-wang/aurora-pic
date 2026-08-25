#!/usr/bin/env python3
"""Focused deck-generation tests for the ionizing-tail continuation."""

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

from run_aurorapic_ionizing_tail_block import (
    ACKNOWLEDGEMENT, TailBlockError, analyze_surface_flux, build_deck,
    validate_inputs,
)


def main() -> None:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-argon-ionizing-tail-rule-20260821.json"
    ).read_text(encoding="utf-8"))
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/"
        "matched-heating-long-window-20260819/measurement/input.cfg"
    ).read_text(encoding="utf-8")
    deck = build_deck(base, Path("output"), Path("checkpoint.apc"), rule)
    required = (
        "steps = 72000", "spatial_average_interval = 2",
        "spatial_average_start_step = 56001",
        "spatial_average_phase_bins = 200",
        "phase_eedf_energy_bins = 320", "phase_eedf_energy_max = 80.0",
        "phase_eedf_tail_threshold = 0.0",
        "phase_eedf_history = false",
        "phase_eedf_regions = x000_010:0.0:0.0025",
        "runtime_backend = serial", "runtime_threads = 1",
        "restart_path = checkpoint.apc",
    )
    for text in required:
        assert text in deck, text
    assert deck.count("phase_eedf_regions =") == 1
    assert deck.count("phase_eedf_tail_threshold =") == 1
    assert deck.count("phase_eedf_history =") == 1

    surface_rule = json.loads(Path(
        "benchmarks/ccp/edupic-argon-surface-flux-rule-20260822.json"
    ).read_text(encoding="utf-8"))
    surface_deck = build_deck(
        base, Path("surface-output"), Path("checkpoint-72000.apc"),
        surface_rule)
    surface_required = (
        "steps = 88000", "spatial_average_start_step = 72001",
        "phase_surface_flux = true",
        "phase_surface_flux_reset_on_restart = true",
        "phase_surface_flux_species = electrons",
        "phase_surface_flux_positions = 0.005,0.015",
        "phase_surface_flux_energy_bins = 320",
        "phase_surface_flux_energy_max = 80.0",
    )
    for text in surface_required:
        assert text in surface_deck, text
    assert surface_deck.count("phase_surface_flux =") == 1

    diagnostic = {
        "phase_bins": 2,
        "surface_flux": {
            "species": "electrons", "positions_m": [0.5],
            "energy_bins": 2, "energy_max_eV": 2.0,
            "direction_order": ["left_to_right", "right_to_left"],
            "minimum_total_macro_crossings_per_surface": 1,
        },
    }
    metadata = {
        "phase_surface_flux_enabled": True,
        "phase_surface_flux_species": "electrons",
        "phase_surface_flux_positions": [0.5],
        "phase_surface_flux_energy_bins": 2,
        "phase_surface_flux_energy_max": 2.0,
    }
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory)
        summary_header = (
            "phase_bin,surface_id,direction,macro_crossings,"
            "represented_crossings,overflow_fraction,"
            "represented_particle_flux_m-2_s-1,"
            "kinetic_energy_flux_W_m-2\n")
        histogram_header = (
            "phase_bin,surface_id,direction,represented_crossings,"
            "probability_density\n")
        summary_rows = []
        histogram_rows = []
        for phase in range(2):
            for direction in ("left_to_right", "right_to_left"):
                summary_rows.append(
                    f"{phase},0,{direction},1,2,0,1,1\n")
                histogram_rows.extend((
                    f"{phase},0,{direction},1,0.5\n",
                    f"{phase},0,{direction},1,0.5\n",
                ))
        (output / "phase_surface_flux_summary.csv").write_text(
            summary_header + "".join(summary_rows), encoding="utf-8")
        (output / "phase_surface_flux.csv").write_text(
            histogram_header + "".join(histogram_rows), encoding="utf-8")
        result = analyze_surface_flux(output, diagnostic, metadata)
        assert all(result[key] for key in (
            "contract", "shape", "finite", "histogram_closure",
            "crossing_sufficiency")), result

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        paths = []
        for name in ("solver", "config", "checkpoint", "report"):
            path = root / name
            path.write_text(name, encoding="utf-8")
            paths.append(path)
        digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        locked = {
            "id": "state_b", "step": 10,
            "solver_sha256": digest(paths[0]),
            "base_config_sha256": digest(paths[1]),
            "checkpoint_sha256": digest(paths[2]),
            "prior_report_sha256": digest(paths[3]),
        }
        replication_rule = {
            "locked_initial_states": [{**locked, "id": "state_a"}, locked],
            "execution_contract": {
                "start_step": 10, "end_step": 18,
                "cycles": 2, "steps_per_cycle": 4,
            },
        }
        args = argparse.Namespace(
            acknowledge_cost=ACKNOWLEDGEMENT, initial_state_id="state_b",
            executable=paths[0], base_config=paths[1], checkpoint=paths[2],
            prior_report=paths[3])
        assert validate_inputs(args, replication_rule) == (10, 18)
        args.initial_state_id = "unknown"
        try:
            validate_inputs(args, replication_rule)
            raise AssertionError("unlocked initial state was accepted")
        except TailBlockError:
            pass


if __name__ == "__main__":
    main()
    print("ionizing-tail runner tests passed")
