#!/usr/bin/env python3
"""Run locked collision-free AuroraPIC common-state horizons serially."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time


class RunError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def deck(rule: dict[str, object], state: Path, horizon: int,
         output: Path) -> str:
    physics = rule["physics_contract"]
    locked = rule["locked_inputs"]
    populations = locked["source_populations"]
    return f"""config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {physics['nodes']}
length = {physics['length_m']:.17g}
dt = {physics['electron_timestep_s']:.17g}
steps = {horizon}
output_interval = {max(1, horizon)}
output_dir = {output}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_left_amplitude = {physics['electrode_voltage_amplitude_V']:.17g}
phi_left_frequency = {physics['rf_frequency_hz']:.17g}
phi_left_phase = {physics['aurorapic_initial_phase_rad']:.17g}
seed = 13507
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 1000000
initial_state_path = {state}
initial_state_signature = {locked['particle_state_signature']}
collision_velocity_sampling = leapfrog_half_step

[species.electrons]
charge = -1.6021766200000001e-19
mass = 9.1093835599999998e-31
weight = {physics['macro_weight_aurorapic_1d']}
particles = {populations['electrons']}
thermal_velocity = 0
loading = random

[species.ions]
charge = 1.6021766200000001e-19
mass = 6.6335209000000003e-26
weight = {physics['macro_weight_aurorapic_1d']}
particles = {populations['ions']}
thermal_velocity = 0
loading = random
timestep_multiplier = {physics['ion_timestep_multiplier']}
"""


def parse_resource(path: Path) -> tuple[float, int]:
    text = path.read_text(encoding="utf-8")
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if elapsed is None or rss is None:
        raise RunError("resource report is incomplete")
    parts = [float(item) for item in elapsed.group(1).split(":")]
    seconds = parts[-1] + (parts[-2] * 60 if len(parts) > 1 else 0)
    return seconds, int(rss.group(1))


def summarize_member(member: Path, horizon: int) -> dict[str, object]:
    output = member / "output"
    field = output / f"fields_{horizon}.csv"
    resources = member / "resources.txt"
    if not field.is_file() or not resources.is_file():
        raise RunError(f"horizon {horizon} has an incomplete existing directory")
    wall, rss = parse_resource(resources)
    with (output / "scalars.csv").open(newline="", encoding="utf-8") as stream:
        scalar = list(csv.DictReader(stream))[-1]
    return {
        "horizon": horizon, "field_sha256": sha256(field),
        "electron_population": int(scalar["live_particles_electrons"]),
        "ion_population": int(scalar["live_particles_ions"]),
        "wall_seconds": wall, "peak_resident_set_kib": rss,
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve(); state = args.state.resolve()
    binary = args.binary.resolve(); root = args.output_root.resolve()
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    locked = rule["locked_inputs"]; contract = rule["execution_contract"]
    if sha256(state) != locked["particle_state_sha256"]:
        raise RunError("particle-state hash differs")
    if sha256(binary) != locked["aurorapic_binary_sha256"]:
        raise RunError("AuroraPIC binary hash differs")
    root.mkdir(parents=True, exist_ok=True)
    members = []
    for horizon in rule["sampling_contract"]["matching_aurorapic_post_step_horizons"]:
        member = root / f"horizon-{horizon:04d}"
        output = member / "output"
        if member.exists():
            members.append(summarize_member(member, horizon))
            continue
        member.mkdir()
        config = member / "input.cfg"
        config.write_text(deck(rule, state, horizon, output), encoding="utf-8")
        stdout = member / "stdout.txt"; resources = member / "resources.txt"
        limit_bytes = int(contract["address_space_limit_kib"]) * 1024
        def limits() -> None:
            resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
            os.nice(10)
        started = time.monotonic()
        with stdout.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(resources),
                 str(binary), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)], stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=int(contract["timeout_seconds_each_aurorapic_horizon"]),
                preexec_fn=limits, check=False)
        if completed.returncode != (1 if horizon == 0 else 0):
            raise RunError(f"horizon {horizon} returned {completed.returncode}")
        summary = summarize_member(member, horizon)
        summary["orchestration_wall_seconds"] = time.monotonic() - started
        members.append(summary)
    return {
        "schema_version": 1,
        "scope": "aurorapic_collision_free_common_state_horizons",
        "rule_sha256": sha256(rule_path),
        "particle_state_sha256": sha256(state),
        "binary_sha256": sha256(binary),
        "members": members,
        "all_resource_gates_passed": all(
            item["peak_resident_set_kib"] <=
            int(contract["maximum_peak_resident_set_kib"])
            for item in members),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (RunError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        parser.error(str(error))
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
