#!/usr/bin/env python3
"""Run the preregistered Turner subcycle-policy no-op control."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time


class ControlError(RuntimeError):
    """A locked input, execution, or comparison gate failed."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deck(rule: dict[str, object], checkpoint: Path, electron_gas: Path,
         ion_gas: Path, output: Path, policy: str) -> str:
    contract = rule["paired_contract"]
    end_step = int(contract["end_step"])
    start_step = int(contract["start_step"])
    return f"""# Prospective Turner Case 1 subcycle-policy invariance control.
config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = 129
length = 0.067000000000000004
dt = 1.8436578171091445e-10
steps = {end_step}
output_interval = 400
output_dir = {output}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_right_amplitude = 450
phi_right_frequency = 13560000
seed = {contract['seed']}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 262144
restart_path = {checkpoint}
checkpoint_output = false
spatial_average = true
spatial_average_reset_on_restart = true
spatial_average_interval = 1
spatial_average_start_step = {start_step}
spatial_average_end_step = {end_step}
spatial_average_rf_frequency = 13560000
spatial_average_rf_cycles = 1
subcycle_charge_deposition = {policy}

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = 9.64e+20
neutral_temperature = 300
max_frequency = 150000000
max_candidates_per_particle = 16
gas_data_file = {electron_gas}

[collisions.electron_mcc.channel.ionization]
secondary_species = electrons
ion_species = ions

[collisions.ion_mcc]
model = null_collision
species = ions
neutral_density = 9.64e+20
neutral_temperature = 300
max_frequency = 120000000
max_candidates_per_particle = 16
gas_data_file = {ion_gas}

[species.electrons]
charge = -1.602176487e-19
mass = 9.1089999999999993e-31
weight = 261718750.00000003
particles = 65536
thermal_velocity = 674321.70332745416
loading = quiet_start

[species.ions]
charge = 1.602176487e-19
mass = 6.6700000000000003e-27
weight = 261718750.00000003
particles = 65536
thermal_velocity = 788.02423116649607
loading = quiet_start
"""


def peak_rss(path: Path) -> int:
    prefix = "Maximum resident set size (kbytes):"
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith(prefix):
            return int(line.split(":", 1)[1].strip())
    raise ControlError(f"maximum RSS absent from {path}")


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path = args.rule.resolve()
    binary = args.binary.resolve()
    checkpoint = args.checkpoint.resolve()
    case_manifest = args.case_manifest.resolve()
    electron_gas = args.electron_gas.resolve()
    ion_gas = args.ion_gas.resolve()
    root = args.output_root.resolve()
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    locked = rule["locked_inputs"]
    inputs = {
        binary: locked["solver_binary_sha256"],
        checkpoint: locked["source_checkpoint_sha256"],
        case_manifest: locked["turner_case_manifest_sha256"],
        electron_gas: locked["electron_gas_manifest_sha256"],
        ion_gas: locked["ion_gas_manifest_sha256"],
    }
    for path, expected in inputs.items():
        if sha256(path) != expected:
            raise ControlError(f"locked hash differs: {path}")
    if root.exists():
        raise ControlError(f"refusing to overwrite output root: {root}")
    root.mkdir(parents=True)
    execution = rule["execution_contract"]
    branches: dict[str, object] = {}
    policies = {"current_position": "current_position",
                "pre_push_held": "pre_push_held"}
    for name, policy in policies.items():
        member = root / name
        output = member / "output"
        member.mkdir()
        config = member / "input.cfg"
        config.write_text(deck(rule, checkpoint, electron_gas, ion_gas,
                               output, policy), encoding="utf-8")
        stdout = member / "stdout.txt"
        resources = member / "resources.txt"
        limit_bytes = int(execution["address_space_limit_kib"]) * 1024

        def limits() -> None:
            resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
            os.nice(int(execution["nice_increment"]))

        started = time.monotonic()
        with stdout.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(resources),
                 str(binary), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)], stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=int(execution["timeout_seconds_each"]),
                preexec_fn=limits, check=False)
        if completed.returncode != 0:
            raise ControlError(
                f"branch {name} returned {completed.returncode}; see {stdout}")
        branches[name] = {
            "policy": policy,
            "config_sha256": sha256(config),
            "wall_seconds": time.monotonic() - started,
            "peak_resident_set_kib": peak_rss(resources),
        }

    comparisons = []
    all_identical = True
    for filename in rule["acceptance"]["required_identical_csv_outputs"]:
        left = root / "current_position" / "output" / filename
        right = root / "pre_push_held" / "output" / filename
        if not left.is_file() or not right.is_file():
            raise ControlError(f"required comparison output absent: {filename}")
        left_hash, right_hash = sha256(left), sha256(right)
        identical = left_hash == right_hash
        all_identical = all_identical and identical
        comparisons.append({"path": filename, "byte_identical": identical,
                            "current_position_sha256": left_hash,
                            "pre_push_held_sha256": right_hash})
    rss_passed = all(
        item["peak_resident_set_kib"] <=
        int(rule["acceptance"]["maximum_peak_resident_set_kib_each"])
        for item in branches.values())
    passed = all_identical and rss_passed
    return {
        "schema_version": 1,
        "scope": rule["scope"],
        "rule_sha256": sha256(rule_path),
        "locked_input_hashes_passed": True,
        "branches": branches,
        "comparisons": comparisons,
        "all_required_physical_outputs_byte_identical": all_identical,
        "all_resource_gates_passed": rss_passed,
        "passed": passed,
        "classification": (
            "turner_subcycle_policy_invariance_established" if passed else
            "turner_subcycle_policy_invariance_failed"),
        "interpretation": (
            "The held-density correction is inactive when all Turner species "
            "advance every electron step; the existing published-duration "
            "density discrepancy remains unchanged and unresolved."),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--case-manifest", type=Path, required=True)
    parser.add_argument("--electron-gas", type=Path, required=True)
    parser.add_argument("--ion-gas", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (ControlError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        parser.error(str(error))
    args.report.write_text(json.dumps(result, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
