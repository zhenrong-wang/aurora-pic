#!/usr/bin/env python3
"""Run the locked collision-enabled common-state pilot serially."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import time


class RunError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_resource(path: Path) -> tuple[float, int]:
    value = path.read_text(encoding="utf-8")
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", value)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", value)
    if elapsed is None or rss is None:
        raise RunError(f"incomplete resource report: {path}")
    parts = [float(item) for item in elapsed.group(1).split(":")]
    seconds = parts[-1] + (parts[-2] * 60.0 if len(parts) > 1 else 0.0)
    return seconds, int(rss.group(1))


def last_csv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RunError(f"empty CSV: {path}")
    return rows[-1]


def deck(rule: dict[str, object], state: Path, seed: int,
         electron_gas: Path, ion_gas: Path, output: Path) -> str:
    physics = rule["physics_contract"]
    locked = rule["locked_inputs"]
    populations = locked["initial_populations"]
    return f"""config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {physics['nodes']}
length = {physics['length_m']:.17g}
dt = {physics['electron_timestep_s']:.17g}
steps = {rule['ensemble_contract']['electron_pushes']}
output_interval = {rule['ensemble_contract']['electron_pushes']}
output_dir = {output}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_left_amplitude = {physics['electrode_voltage_amplitude_V']:.17g}
phi_left_frequency = {physics['rf_frequency_hz']:.17g}
phi_left_phase = {physics['aurorapic_initial_phase_rad']:.17g}
seed = {seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 1000000
initial_state_path = {state}
initial_state_signature = {locked['aurorapic_particle_state_signature']}
collision_velocity_sampling = leapfrog_half_step
subcycle_charge_deposition = pre_push_held

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = {physics['neutral_density_m3']:.17g}
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
neutral_density = {physics['neutral_density_m3']:.17g}
neutral_temperature = {physics['neutral_temperature_K']:.17g}
max_frequency = 100000000
max_candidates_per_particle = 16
gas_data_file = {ion_gas}

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


def summarize_native(member: Path, seed: int) -> dict[str, object]:
    metrics = last_csv(member / "edupic_collision_endpoint_metrics.csv")
    field = member / "edupic_collision_endpoint_field.csv"
    wall, rss = parse_resource(member / "resources.txt")
    return {
        "implementation": "edupic", "seed": seed,
        "endpoint": {key: int(value) for key, value in metrics.items()},
        "field_sha256": sha256(field), "wall_seconds": wall,
        "peak_resident_set_kib": rss,
    }


def summarize_aurora(member: Path, seed: int, horizon: int) -> dict[str, object]:
    output = member / "output"
    collision = last_csv(output / "collisions.csv")
    boundary = last_csv(output / "boundary_losses.csv")
    scalar = last_csv(output / "scalars.csv")
    field = output / f"fields_{horizon}.csv"
    wall, rss = parse_resource(member / "resources.txt")
    endpoint = {
        "pre_push_step": horizon + 1,
        "electrons": int(scalar["live_particles_electrons"]),
        "ions": int(scalar["live_particles_ions"]),
        "electron_elastic": int(collision[
            "cumulative_collisions_electron_mcc.elastic"]),
        "electron_excitation": int(collision[
            "cumulative_collisions_electron_mcc.excitation"]),
        "electron_ionization": int(collision[
            "cumulative_collisions_electron_mcc.ionization"]),
        "ion_isotropic": int(collision[
            "cumulative_collisions_ion_mcc.isotropic"]),
        "ion_backward": int(collision[
            "cumulative_collisions_ion_mcc.backward"]),
        "electron_absorbed_left": int(boundary[
            "absorbed_left_count_electrons"]),
        "electron_absorbed_right": int(boundary[
            "absorbed_right_count_electrons"]),
        "ion_absorbed_left": int(boundary["absorbed_left_count_ions"]),
        "ion_absorbed_right": int(boundary["absorbed_right_count_ions"]),
    }
    return {
        "implementation": "aurorapic", "seed": seed,
        "endpoint": endpoint, "field_sha256": sha256(field),
        "wall_seconds": wall, "peak_resident_set_kib": rss,
    }


def run_process(command: list[str], member: Path, contract: dict[str, object]) -> None:
    resources = member / "resources.txt"
    stdout = member / "stdout.txt"
    limit_bytes = int(contract["address_space_limit_kib"]) * 1024

    def limits() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
        os.nice(int(contract["nice_increment"]))

    with stdout.open("w", encoding="utf-8") as stream:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(resources), *command],
            cwd=member, stdout=stream, stderr=subprocess.STDOUT,
            timeout=int(contract["timeout_seconds_each_member"]),
            preexec_fn=limits, check=False)
    if result.returncode != 0:
        raise RunError(f"member {member} returned {result.returncode}")


def execute(args: argparse.Namespace) -> dict[str, object]:
    paths = {name: value.resolve() for name, value in {
        "rule": args.rule, "lock": args.lock, "edupic_binary": args.edupic_binary,
        "aurorapic_binary": args.aurorapic_binary,
        "checkpoint": args.edupic_checkpoint, "state": args.aurorapic_state,
        "electron_gas": args.electron_gas, "ion_gas": args.ion_gas}.items()}
    rule = json.loads(paths["rule"].read_text(encoding="utf-8"))
    lock = json.loads(paths["lock"].read_text(encoding="utf-8"))
    parent_report = None
    parent_members: dict[tuple[str, int], dict[str, object]] = {}
    if args.parent_report is not None:
        parent_path = args.parent_report.resolve()
        expected_parent = rule["basis"]["parent_four_period_execution_sha256"]
        if sha256(parent_path) != expected_parent:
            raise RunError("parent execution report hash differs")
        parent_report = json.loads(parent_path.read_text(encoding="utf-8"))
        parent_members = {
            (item["implementation"], int(item["seed"])): item
            for item in parent_report["members"]
        }
    expected = {
        "rule_sha256": sha256(paths["rule"]),
        "edupic_binary_sha256": sha256(paths["edupic_binary"]),
        "aurorapic_binary_sha256": sha256(paths["aurorapic_binary"]),
    }
    if any(lock.get(key) != value for key, value in expected.items()):
        raise RunError("execution lock differs from selected rule or binaries")
    locked = rule["locked_inputs"]
    checks = {
        "checkpoint": (paths["checkpoint"], locked["edupic_checkpoint_sha256"]),
        "state": (paths["state"], locked["aurorapic_particle_state_sha256"]),
        "electron_gas": (paths["electron_gas"], locked["electron_gas_manifest_sha256"]),
        "ion_gas": (paths["ion_gas"], locked["ion_gas_manifest_sha256"]),
    }
    for label, (path, wanted) in checks.items():
        if sha256(path) != wanted:
            raise RunError(f"{label} hash differs")

    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    horizon = int(rule["ensemble_contract"]["electron_pushes"])
    contract = rule["execution_contract"]
    members: list[dict[str, object]] = []
    reused_parent_members_verified = True
    for implementation in ("edupic", "aurorapic"):
        for seed in rule["ensemble_contract"]["seeds_each_implementation"]:
            member = root / implementation / f"seed-{seed}"
            if member.exists():
                summary = (summarize_native(member, seed)
                           if implementation == "edupic" else
                           summarize_aurora(member, seed, horizon))
                if parent_report is not None and (implementation, seed) in parent_members:
                    parent = parent_members[(implementation, seed)]
                    reused_parent_members_verified &= (
                        summary["endpoint"] == parent["endpoint"] and
                        summary["field_sha256"] == parent["field_sha256"])
                members.append(summary)
                continue
            member.mkdir(parents=True)
            started = time.monotonic()
            if implementation == "edupic":
                shutil.copy2(paths["checkpoint"], member / "picdata.bin")
                run_process(
                    [str(paths["edupic_binary"]),
                     str(rule["ensemble_contract"].get("rf_periods", 1)),
                     "pilot", str(seed)],
                    member, contract)
                summary = summarize_native(member, seed)
            else:
                output = member / "output"
                config = member / "input.cfg"
                config.write_text(deck(
                    rule, paths["state"], seed, paths["electron_gas"],
                    paths["ion_gas"], output), encoding="utf-8")
                run_process([
                    str(paths["aurorapic_binary"]), "--allow-large-run",
                    "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)],
                    member, contract)
                summary = summarize_aurora(member, seed, horizon)
            summary["orchestration_wall_seconds"] = time.monotonic() - started
            members.append(summary)
    maximum_rss = int(contract["maximum_peak_resident_set_kib"])
    expected_members = (2 * int(
        rule["ensemble_contract"]["members_each_implementation"]))
    complete = len(members) == expected_members
    result = {
        "schema_version": 1,
        "scope": ("collision_enabled_common_state_ensemble_execution" if
                  int(rule["ensemble_contract"].get("rf_periods", 1)) == 1 else
                  "collision_enabled_common_state_four_period_execution"),
        **expected,
        "execution_lock_sha256": sha256(paths["lock"]),
        "members": members,
        "all_members_complete": complete,
        "all_resource_gates_passed": all(
            item["peak_resident_set_kib"] <= maximum_rss for item in members),
        "total_member_wall_seconds": sum(item["wall_seconds"] for item in members),
    }
    if expected_members == 10:
        result["all_ten_members_complete"] = complete
    if parent_report is not None:
        result["parent_execution_sha256"] = sha256(args.parent_report.resolve())
        result["reused_parent_members_verified"] = reused_parent_members_verified
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--edupic-binary", type=Path, required=True)
    parser.add_argument("--aurorapic-binary", type=Path, required=True)
    parser.add_argument("--edupic-checkpoint", type=Path, required=True)
    parser.add_argument("--aurorapic-state", type=Path, required=True)
    parser.add_argument("--electron-gas", type=Path, required=True)
    parser.add_argument("--ion-gas", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--parent-report", type=Path)
    args = parser.parse_args()
    try:
        report = execute(args)
    except (RunError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        parser.error(str(error))
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
