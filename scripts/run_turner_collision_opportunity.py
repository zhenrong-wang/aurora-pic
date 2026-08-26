#!/usr/bin/env python3
"""Run and analyze the preregistered Turner collision-opportunity ensemble."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time


class RunError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deck(rule: dict[str, object], state: Path, electron_gas: Path,
         ion_gas: Path, output: Path, seed: int, mode: str) -> str:
    p = rule["physics_contract"]
    e = rule["ensemble_contract"]
    pops = rule["locked_inputs"]["source_populations"]
    return f"""# Turner common-state collision-opportunity discriminator.
config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {p['nodes']}
length = {p['length_m']:.17g}
dt = {p['timestep_s']:.17g}
steps = {e['steps']}
output_interval = 400
output_dir = {output}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_right_amplitude = {p['voltage_amplitude_V']:.17g}
phi_right_frequency = {p['rf_frequency_hz']:.17g}
seed = {seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = 262144
initial_state_path = {state}
initial_state_signature = {rule['locked_inputs']['particle_state_signature']}
checkpoint_output = false
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 1
spatial_average_end_step = {e['steps']}
spatial_average_rf_frequency = {p['rf_frequency_hz']:.17g}
spatial_average_rf_cycles = {e['measurement_cycles']}

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = {p['neutral_density_m3']:.17g}
neutral_temperature = {p['neutral_temperature_K']:.17g}
max_frequency = {p['electron_max_frequency_s']:.17g}
opportunity_sampling = {mode}
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
max_frequency = {p['ion_max_frequency_s']:.17g}
opportunity_sampling = {mode}
max_candidates_per_particle = 16
gas_data_file = {ion_gas}

[species.electrons]
charge = -1.602176487e-19
mass = 9.1089999999999993e-31
weight = {p['macro_weight_m2']:.17g}
particles = {pops['electrons']}
thermal_velocity = 0
loading = random

[species.ions]
charge = 1.602176487e-19
mass = 6.6700000000000003e-27
weight = {p['macro_weight_m2']:.17g}
particles = {pops['ions']}
thermal_velocity = 0
loading = random
"""


def final_row(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RunError(f"empty CSV: {path}")
    result = {key: float(value) for key, value in rows[-1].items()}
    if not all(math.isfinite(value) for value in result.values()):
        raise RunError(f"non-finite CSV value: {path}")
    return result


def densities(path: Path) -> dict[str, float]:
    grouped: dict[str, list[tuple[float, float]]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            grouped.setdefault(row["species"], []).append(
                (float(row["x_m"]), float(row["number_density_mean_m-3"])))
    result = {}
    for species, points in grouped.items():
        points.sort()
        integral = sum(
            0.5 * (left[1] + right[1]) * (right[0] - left[0])
            for left, right in zip(points, points[1:]))
        if not math.isfinite(integral):
            raise RunError(f"non-finite density integral: {species}")
        result[species] = integral
    return result


def peak_rss(path: Path) -> int:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("Maximum resident set size (kbytes):"):
            return int(line.split(":", 1)[1])
    raise RunError(f"maximum RSS absent: {path}")


def summarize(member: Path) -> dict[str, object]:
    output = member / "output"
    collisions = final_row(output / "collisions.csv")
    losses = final_row(output / "boundary_losses.csv")
    power = final_row(output / "power_transfer.csv")
    scalars = final_row(output / "scalars.csv")
    return {
        "config_sha256": sha256(member / "input.cfg"),
        "peak_resident_set_kib": peak_rss(member / "resources.txt"),
        "wall_seconds": float((member / "wall-seconds.txt").read_text()),
        "candidates": int(collisions["cumulative_candidates"]),
        "null_collisions": int(collisions["cumulative_null_collisions"]),
        "channels": {
            key.removeprefix("cumulative_collisions_"): int(value)
            for key, value in collisions.items()
            if key.startswith("cumulative_collisions_")
        },
        "endpoint_populations": {
            "electrons": int(scalars["live_particles_electrons"]),
            "ions": int(scalars["live_particles_ions"]),
        },
        "integrated_density_m2": densities(output / "spatial_average.csv"),
        "wall_losses": {
            key: int(value) for key, value in losses.items()
            if key.startswith("absorbed_") and key.endswith(("_electrons", "_ions"))
            and "count" in key
        },
        "electrical_work_J_m2": {
            "electrons": power["electric_work_electrons_J_m-2"],
            "ions": power["electric_work_ions_J_m-2"],
        },
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def execute(args: argparse.Namespace) -> dict[str, object]:
    rule_path, binary, state = (args.rule.resolve(), args.binary.resolve(),
                                args.state.resolve())
    electron_gas, ion_gas = (args.electron_gas.resolve(),
                             args.ion_gas.resolve())
    root = args.output_root.resolve()
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    locked = rule["locked_inputs"]
    for path, expected in ((binary, locked["solver_binary_sha256"]),
                           (state, locked["particle_state_sha256"]),
                           (electron_gas, locked["electron_gas_manifest_sha256"]),
                           (ion_gas, locked["ion_gas_manifest_sha256"])):
        if sha256(path) != expected:
            raise RunError(f"locked hash differs: {path}")
    if root.exists() and not args.reuse_completed:
        raise RunError(f"refusing to overwrite output root: {root}")
    root.mkdir(parents=True, exist_ok=args.reuse_completed)
    execution = rule["execution_contract"]
    members: dict[str, dict[str, object]] = {}
    for mode in rule["ensemble_contract"]["modes"]:
        members[mode] = {}
        for seed in rule["ensemble_contract"]["seeds"]:
            member = root / mode / f"seed-{seed}"
            output = member / "output"
            expected_deck = deck(
                rule, state, electron_gas, ion_gas, output, seed, mode)
            if args.reuse_completed and (member / "wall-seconds.txt").is_file():
                if (member / "input.cfg").read_text(encoding="utf-8") != expected_deck:
                    raise RunError(f"completed member deck differs: {member}")
                members[mode][str(seed)] = summarize(member)
                continue
            member.mkdir(parents=True, exist_ok=args.reuse_completed)
            config = member / "input.cfg"
            config.write_text(expected_deck, encoding="utf-8")
            resources = member / "resources.txt"
            stdout = member / "stdout.txt"
            limit_bytes = int(execution["address_space_limit_kib"]) * 1024

            def limits() -> None:
                resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
                os.nice(int(execution["nice_increment"]))

            started = time.monotonic()
            with stdout.open("w", encoding="utf-8") as stream:
                completed = subprocess.run(
                    ["/usr/bin/time", "-v", "-o", str(resources),
                     str(binary), "--allow-large-run",
                     "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)],
                    stdout=stream, stderr=subprocess.STDOUT,
                    timeout=int(execution["timeout_seconds_each"]),
                    preexec_fn=limits, check=False)
            elapsed = time.monotonic() - started
            (member / "wall-seconds.txt").write_text(f"{elapsed}\n")
            if completed.returncode != 0:
                raise RunError(f"{mode} seed {seed} returned {completed.returncode}")
            members[mode][str(seed)] = summarize(member)
    modes = {}
    for mode, values in members.items():
        entries = list(values.values())
        modes[mode] = {
            "mean_candidates": mean([x["candidates"] for x in entries]),
            "mean_channels": {
                channel: mean([x["channels"][channel] for x in entries])
                for channel in entries[0]["channels"]
            },
            "mean_endpoint_populations": {
                species: mean([x["endpoint_populations"][species] for x in entries])
                for species in ("electrons", "ions")
            },
            "mean_integrated_density_m2": {
                species: mean([x["integrated_density_m2"][species] for x in entries])
                for species in ("electrons", "ions")
            },
        }
    poisson, single = modes["poisson_clock"], modes["single_bernoulli"]
    candidate_ratio = poisson["mean_candidates"] / single["mean_candidates"]
    ionization_key = "electron_mcc.ionization"
    ionization_ratio = (poisson["mean_channels"][ionization_key] /
                        single["mean_channels"][ionization_key])
    ion_density_ratio = (
        single["mean_integrated_density_m2"]["ions"] /
        poisson["mean_integrated_density_m2"]["ions"])
    p = rule["physics_contract"]
    weight = p["macro_weight_m2"]
    opportunity_means = {
        "poisson_clock": {
            "electrons": p["electron_max_frequency_s"] * p["timestep_s"],
            "ions": p["ion_max_frequency_s"] * p["timestep_s"]},
        "single_bernoulli": {
            "electrons": -math.expm1(
                -p["electron_max_frequency_s"] * p["timestep_s"]),
            "ions": -math.expm1(
                -p["ion_max_frequency_s"] * p["timestep_s"])} }
    exposure_closure = {}
    for mode in ("poisson_clock", "single_bernoulli"):
        expected = rule["ensemble_contract"]["steps"] * sum(
            modes[mode]["mean_integrated_density_m2"][species] / weight *
            opportunity_means[mode][species]
            for species in ("electrons", "ions"))
        exposure_closure[mode] = {
            "expected_candidates_from_mean_population_exposure": expected,
            "observed_to_expected_ratio": modes[mode]["mean_candidates"] / expected,
            "relative_residual": modes[mode]["mean_candidates"] / expected - 1.0,
        }
    bounds = rule["acceptance"]["poisson_to_single_total_candidate_ratio"]
    resource_pass = all(
        x["peak_resident_set_kib"] <=
        rule["acceptance"]["maximum_peak_resident_set_kib_each"]
        for mode in members.values() for x in mode.values())
    if "maximum_absolute_exposure_normalized_candidate_residual" in rule["acceptance"]:
        tolerance = rule["acceptance"][
            "maximum_absolute_exposure_normalized_candidate_residual"]
        implementation_pass = resource_pass and all(
            abs(value["relative_residual"]) <= tolerance
            for value in exposure_closure.values())
    else:
        implementation_pass = bounds[0] <= candidate_ratio <= bounds[1] and resource_pass
    direction = (ionization_ratio > 1.0 and ion_density_ratio < 1.0)
    classification = (
        "collision_opportunity_execution_failed" if not implementation_pass else
        "turner_bias_direction_supported" if direction else
        "turner_bias_direction_not_supported")
    return {
        "schema_version": 1,
        "scope": rule["scope"],
        "rule_sha256": sha256(rule_path),
        "locked_input_hashes_passed": True,
        "completed_members_reused_without_rerun": args.reuse_completed,
        "members": members,
        "ensemble_means": modes,
        "ratios": {
            "poisson_to_single_candidates": candidate_ratio,
            "poisson_to_single_ionization": ionization_ratio,
            "single_to_poisson_integrated_ion_density": ion_density_ratio,
        },
        "candidate_exposure_closure": exposure_closure,
        "implementation_passed": implementation_pass,
        "resource_gates_passed": resource_pass,
        "classification": classification,
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--electron-gas", type=Path, required=True)
    parser.add_argument("--ion-gas", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reuse-completed", action="store_true")
    args = parser.parse_args()
    try:
        result = execute(args)
    except (RunError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        parser.error(str(error))
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["implementation_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
