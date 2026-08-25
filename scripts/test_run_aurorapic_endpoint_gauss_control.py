#!/usr/bin/env python3

import json
from pathlib import Path

from run_aurorapic_endpoint_gauss_control import relaxation_deck


def main() -> int:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-endpoint-gauss-control-rule-20260825.json"
    ).read_text(encoding="utf-8"))
    base = """steps = 1
output_interval = 1
output_dir = old
spatial_average_interval = 1
spatial_average_start_step = 1
spatial_average_end_step = 1
spatial_average_rf_cycles = 1
spatial_average_phase_bins = 1
phi_left_phase = 1.5707963267948966
runtime_backend = serial
runtime_threads = 1
seed = 1
checkpoint_interval = 1
initial_state_path = state.aps
initial_state_signature = 1
[species.electrons]
particles = 1
[species.ions]
particles = 1
"""
    state = rule["locked_continuation_states"][0]
    deck = relaxation_deck(
        base, rule, state, Path("output"), Path("checkpoint.apc"))
    assert "phi_left_phase = 1.5707963267948966" in deck
    assert "steps = 40000" in deck
    assert "spatial_average_start_step = 36001" in deck
    assert "restart_path = checkpoint.apc" in deck
    assert "collision_velocity_sampling = leapfrog_half_step" in deck
    assert "initial_state_path" not in deck
    print("endpoint Gauss-control runner tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
