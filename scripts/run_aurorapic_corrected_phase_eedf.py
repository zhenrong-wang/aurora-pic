#!/usr/bin/env python3
"""Run the locked corrected-cadence common-state phase-EEDF ensemble."""

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


def deck(rule: dict[str, object], state: Path, seed: int,
         electron_gas: Path, ion_gas: Path, output: Path) -> str:
    p = rule["physics_contract"]
    d = rule["diagnostic_contract"]
    locked = rule["locked_inputs"]
    populations = locked["source_populations"]
    regions = ",".join(d["regions"])
    return f"""config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {p['nodes']}
length = {p['length_m']:.17g}
dt = {p['electron_timestep_s']:.17g}
steps = {d['electron_pushes']}
output_interval = {d['electron_pushes']}
output_dir = {output}
spatial_average = true
spatial_average_interval = {d['sample_every_timesteps']}
spatial_average_start_step = 1
spatial_average_end_step = {d['electron_pushes']}
spatial_average_rf_frequency = {p['rf_frequency_hz']:.17g}
spatial_average_rf_cycles = {d['measurement_cycles']}
spatial_average_phase_bins = {d['phase_bins']}
spatial_average_sampling_order = pre_collision
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_left_amplitude = {p['electrode_voltage_amplitude_V']:.17g}
phi_left_frequency = {p['rf_frequency_hz']:.17g}
phi_left_phase = {p['aurorapic_initial_phase_rad']:.17g}
seed = {seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 1000000
initial_state_path = {state}
initial_state_signature = {locked['particle_state_signature']}
collision_velocity_sampling = leapfrog_half_step
subcycle_charge_deposition = pre_push_held
phase_eedf = true
phase_eedf_species = electrons
phase_eedf_energy_bins = {d['energy_bins']}
phase_eedf_energy_max = {d['energy_max_eV']:.17g}
phase_eedf_regions = {regions}
phase_eedf_tail_threshold = 15.8

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = {p['neutral_density_m3']:.17g}
neutral_temperature = 0
max_frequency = 1000000000
max_candidates_per_particle = 16
gas_data_file = {electron_gas}

[collisions.electron_mcc.channel.ionization]
secondary_species = electrons
ion_species = ions

[collisions.ion_mcc]
model = null_collision
species = ions
neutral_density = {p['neutral_density_m3']:.17g}
neutral_temperature = {p['neutral_temperature_K']:.17g}
max_frequency = 100000000
max_candidates_per_particle = 16
gas_data_file = {ion_gas}

[species.electrons]
charge = -1.6021766200000001e-19
mass = 9.1093835599999998e-31
weight = {p['macro_weight_aurorapic_1d']}
particles = {populations['electrons']}
thermal_velocity = 0
loading = random

[species.ions]
charge = 1.6021766200000001e-19
mass = 6.6335209000000003e-26
weight = {p['macro_weight_aurorapic_1d']}
particles = {populations['ions']}
thermal_velocity = 0
loading = random
timestep_multiplier = {p['ion_timestep_multiplier']}
"""


def parse_resource(path: Path) -> tuple[float, int]:
    value = path.read_text(encoding="utf-8")
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", value)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", value)
    if elapsed is None or rss is None:
        raise RunError("resource report is incomplete")
    parts = [float(item) for item in elapsed.group(1).split(":")]
    seconds = parts[-1] + (parts[-2] * 60.0 if len(parts) > 1 else 0.0)
    return seconds, int(rss.group(1))


def summarize(member: Path, seed: int, rule: dict[str, object]) -> dict[str, object]:
    output = member / "output"
    histogram = output / "phase_eedf.csv"
    moments = output / "phase_eedf_moments.csv"
    if not histogram.is_file() or not moments.is_file():
        raise RunError(f"phase-EEDF member {seed} is incomplete")
    with moments.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    critical = set(rule["diagnostic_contract"]["critical_regions"])
    selected = [row for row in rows if row["region"] in critical]
    if not selected:
        raise RunError("phase-EEDF critical moments are absent")
    wall, rss = parse_resource(member / "resources.txt")
    return {
        "seed": seed,
        "phase_eedf_sha256": sha256(histogram),
        "phase_eedf_moments_sha256": sha256(moments),
        "critical_minimum_macro_observations_per_region_phase_bin": min(
            int(row["macro_observations"]) for row in selected),
        "maximum_histogram_overflow_fraction": max(
            float(row["overflow_fraction"]) for row in rows),
        "wall_seconds": wall,
        "peak_resident_set_kib": rss,
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve(); lock_path = args.lock.resolve()
    binary = args.binary.resolve(); state = args.state.resolve()
    electron_gas = args.electron_gas.resolve(); ion_gas = args.ion_gas.resolve()
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock["rule_sha256"] != sha256(rule_path):
        raise RunError("rule hash differs from execution lock")
    if lock["aurorapic_binary_sha256"] != sha256(binary):
        raise RunError("AuroraPIC binary hash differs")
    checks = (
        (state, rule["locked_inputs"]["particle_state_sha256"]),
        (electron_gas, rule["locked_inputs"]["electron_gas_manifest_sha256"]),
        (ion_gas, rule["locked_inputs"]["ion_gas_manifest_sha256"]),
    )
    if any(sha256(path) != expected for path, expected in checks):
        raise RunError("state or gas input hash differs")
    root = args.output_root.resolve(); root.mkdir(parents=True, exist_ok=True)
    contract = rule["execution_contract"]
    members = []
    for seed in rule["diagnostic_contract"]["candidate_seeds"]:
        member = root / f"seed-{seed}"
        if member.exists():
            members.append(summarize(member, seed, rule))
            continue
        member.mkdir()
        output = member / "output"
        config = member / "input.cfg"
        config.write_text(deck(rule, state, seed, electron_gas, ion_gas, output),
                          encoding="utf-8")
        limit_bytes = int(contract["address_space_limit_kib"]) * 1024
        def limits() -> None:
            resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
            os.nice(int(contract["nice_increment"]))
        started = time.monotonic()
        with (member / "stdout.txt").open("w", encoding="utf-8") as stream:
            result = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(member / "resources.txt"),
                 str(binary), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)],
                cwd=member, stdout=stream, stderr=subprocess.STDOUT,
                timeout=int(contract["timeout_seconds_each_member"]),
                preexec_fn=limits, check=False)
        if result.returncode != 0:
            raise RunError(f"candidate seed {seed} returned {result.returncode}")
        item = summarize(member, seed, rule)
        item["orchestration_wall_seconds"] = time.monotonic() - started
        members.append(item)
    return {
        "schema_version": 1,
        "scope": "corrected_cadence_common_state_phase_eedf_execution",
        "rule_sha256": sha256(rule_path),
        "execution_lock_sha256": sha256(lock_path),
        "aurorapic_binary_sha256": sha256(binary),
        "particle_state_sha256": sha256(state),
        "members": members,
        "all_three_members_complete": len(members) == 3,
        "all_resource_gates_passed": all(
            item["peak_resident_set_kib"] <=
            int(contract["maximum_peak_resident_set_kib"]) for item in members),
        "total_member_wall_seconds": sum(item["wall_seconds"] for item in members),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--electron-gas", type=Path, required=True)
    parser.add_argument("--ion-gas", type=Path, required=True)
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
