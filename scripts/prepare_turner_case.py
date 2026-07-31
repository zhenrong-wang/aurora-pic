#!/usr/bin/env python3
"""Audit local Turner data and generate, but never run, the Case 1 deck."""

from __future__ import annotations

import argparse
import bisect
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


class PreparationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise PreparationError(f"cannot read {path}: {error}") from error
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    if path.exists():
        raise PreparationError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def load_case(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise PreparationError(f"cannot read case manifest {path}: {error}") from error
    for section in (
        "global", "reference", "physics", "numerics", "collision_guard",
        "authorization", "reported_case1_characteristics",
    ):
        require(section in parser, f"case manifest is missing [{section}]")
    require(
        parser["global"].getint("case_manifest_version") == 1,
        "unsupported case manifest version",
    )
    return parser


def load_audit(normalized_dir: Path, case: configparser.ConfigParser) -> dict:
    path = normalized_dir / "audit.json"
    require(
        sha256(path) == case["reference"]["normalized_audit_sha256"],
        "normalization audit checksum mismatch",
    )
    try:
        audit = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PreparationError(f"cannot parse normalization audit: {error}") from error
    require(
        audit.get("turner_normalization_version") == 2,
        "unsupported Turner normalization version",
    )
    require(
        audit.get("case_id") == case["reference"]["source_case_id"],
        "normalization audit case identity mismatch",
    )
    require(
        audit.get("source_artifact", {}).get("sha256")
        == case["reference"]["source_artifact_sha256"],
        "normalization audit source identity mismatch",
    )
    files = audit.get("normalized_files")
    require(isinstance(files, dict) and files, "normalization audit has no files")
    for name, identity in files.items():
        require(
            isinstance(identity, dict) and "sha256" in identity and "bytes" in identity,
            f"incomplete normalized identity for {name}",
        )
        artifact = normalized_dir / name
        try:
            size = artifact.stat().st_size
        except OSError as error:
            raise PreparationError(f"cannot inspect {artifact}: {error}") from error
        require(size == identity["bytes"], f"normalized byte count mismatch: {name}")
        require(
            sha256(artifact) == identity["sha256"],
            f"normalized checksum mismatch: {name}",
        )
    reference_name = case["reference"]["reference_file"]
    require(reference_name in files, "Case 1 reference is absent from audit")
    return audit


def read_table(path: Path, energy_scale: float, cross_section_scale: float) -> list:
    rows = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise PreparationError(f"cannot read collision table {path}: {error}") from error
    for line_number, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        require(len(fields) == 2, f"{path}:{line_number}: expected two columns")
        try:
            energy = float(fields[0]) * energy_scale
            cross_section = float(fields[1]) * cross_section_scale
        except ValueError as error:
            raise PreparationError(f"{path}:{line_number}: invalid number") from error
        require(
            math.isfinite(energy) and math.isfinite(cross_section)
            and energy >= 0 and cross_section >= 0,
            f"{path}:{line_number}: invalid collision value",
        )
        require(not rows or energy > rows[-1][0], f"{path}: energy is not increasing")
        rows.append((energy, cross_section))
    require(len(rows) >= 2, f"{path}: collision table is too short")
    return rows


def interpolate(table: list, energy: float) -> float:
    energies = [row[0] for row in table]
    index = bisect.bisect_right(energies, energy)
    if index == 0:
        return table[0][1]
    if index == len(table):
        return table[-1][1]
    x0, y0 = table[index - 1]
    x1, y1 = table[index]
    return y0 + (y1 - y0) * (energy - x0) / (x1 - x0)


def load_channels(path: Path, normalized_dir: Path) -> tuple[float, list]:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    # ConfigParser needs a synthetic section around the leading metadata block.
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise PreparationError(f"cannot parse gas manifest {path}: {error}") from error
    neutral_mass = parser["global"].getfloat("neutral_mass")
    channels = []
    for section in parser.sections():
        if not section.startswith("collision."):
            continue
        item = parser[section]
        energy_scale = item.getfloat("energy_scale", fallback=1.0)
        cross_scale = item.getfloat("cross_section_scale", fallback=1.0)
        threshold = item.getfloat("threshold_energy", fallback=0.0)
        frame = item.get("energy_frame", fallback="projectile")
        require(
            frame in ("projectile", "center_of_mass"),
            f"{path}: unsupported energy frame {frame}",
        )
        table_path = normalized_dir / item["cross_section_file"]
        channels.append((
            read_table(table_path, energy_scale, cross_scale),
            threshold,
            frame,
        ))
    require(channels, f"{path}: no collision channels")
    return neutral_mass, channels


def sampled_thermal_bound(
    particle_mass: float,
    neutral_mass: float,
    neutral_density: float,
    neutral_temperature: float,
    boltzmann: float,
    energy_ceiling_ev: float,
    elementary_charge: float,
    sigma_limit: float,
    intervals: int,
    channels: list,
) -> tuple[float, float]:
    thermal_limit = sigma_limit * math.sqrt(
        boltzmann * neutral_temperature / neutral_mass
    )
    maximum_speed = math.sqrt(
        2.0 * energy_ceiling_ev * elementary_charge / particle_mass
    )
    peak = 0.0
    peak_speed = 0.0
    for index in range(intervals + 1):
        speed = maximum_speed * index / intervals
        maximum_relative_speed = speed + thermal_limit
        minimum_relative_speed = max(0.0, speed - thermal_limit)
        total = 0.0
        for table, threshold, frame in channels:
            energy_mass = particle_mass
            if frame == "center_of_mass":
                energy_mass = particle_mass * neutral_mass / (
                    particle_mass + neutral_mass
                )
            maximum_energy = (
                0.5 * energy_mass * maximum_relative_speed ** 2
            )
            if maximum_energy < threshold:
                continue
            minimum_energy = max(
                threshold, 0.5 * energy_mass * minimum_relative_speed ** 2
            )
            maximum_cross_section = max(
                interpolate(table, minimum_energy),
                interpolate(table, maximum_energy),
                *(
                    cross_section for energy, cross_section in table
                    if minimum_energy <= energy <= maximum_energy
                ),
            )
            total += (
                neutral_density * maximum_cross_section * maximum_relative_speed
            )
        if total > peak:
            peak = total
            peak_speed = speed
    return peak, peak_speed


def round_up_significant(value: float, digits: int) -> float:
    require(value > 0 and digits > 0, "invalid majorant rounding request")
    exponent = math.floor(math.log10(value)) - digits + 1
    scale = 10.0 ** exponent
    return math.ceil(value / scale) * scale


def prepare(args: argparse.Namespace) -> tuple[Path, Path]:
    case_path = args.case_manifest.resolve()
    normalized_dir = args.normalized_dir.resolve()
    output = args.output.resolve()
    report_path = (
        args.report.resolve()
        if args.report is not None
        else output.with_suffix(".preflight.json")
    )
    require(not output.exists(), f"refusing to overwrite existing file: {output}")
    require(
        not report_path.exists(),
        f"refusing to overwrite existing file: {report_path}",
    )
    case = load_case(case_path)
    authorization = case["authorization"]
    reported = case["reported_case1_characteristics"]
    acknowledgement = authorization["acknowledgement"]
    require(
        args.acknowledge_cost == acknowledgement,
        f"deck generation requires --acknowledge-cost {acknowledgement}",
    )
    require(
        reported.getint("total_macro_particles") > 0,
        "reported Case 1 macro-particle population must be positive",
    )
    audit = load_audit(normalized_dir, case)
    physics = case["physics"]
    numerics = case["numerics"]
    guard = case["collision_guard"]

    electron_gas = normalized_dir / "turner_he_electron.gas"
    ion_gas = normalized_dir / "turner_he_ion.gas"
    electron_neutral_mass, electron_channels = load_channels(
        electron_gas, normalized_dir
    )
    ion_neutral_mass, ion_channels = load_channels(ion_gas, normalized_dir)
    helium_mass = physics.getfloat("helium_mass_kg")
    require(
        electron_neutral_mass == helium_mass and ion_neutral_mass == helium_mass,
        "gas-manifest neutral mass differs from the benchmark contract",
    )

    elementary_charge = physics.getfloat("elementary_charge_c")
    boltzmann = physics.getfloat("boltzmann_constant_j_k")
    neutral_density = physics.getfloat("neutral_density_m3")
    neutral_temperature = physics.getfloat("neutral_temperature_k")
    electron_mass = physics.getfloat("electron_mass_kg")
    intervals = guard.getint("scan_intervals")
    sigma_limit = guard.getfloat("thermal_neutral_sigma_limit")
    electron_bound, electron_peak_speed = sampled_thermal_bound(
        electron_mass, helium_mass, neutral_density, neutral_temperature,
        boltzmann, guard.getfloat("electron_energy_ceiling_ev"),
        elementary_charge, sigma_limit, intervals, electron_channels,
    )
    ion_bound, ion_peak_speed = sampled_thermal_bound(
        helium_mass, helium_mass, neutral_density, neutral_temperature,
        boltzmann, guard.getfloat("ion_energy_ceiling_ev"),
        elementary_charge, sigma_limit, intervals, ion_channels,
    )
    safety_factor = guard.getfloat("safety_factor")
    digits = guard.getint("rounding_significant_digits")
    electron_majorant = round_up_significant(electron_bound * safety_factor, digits)
    ion_majorant = round_up_significant(ion_bound * safety_factor, digits)

    frequency = physics.getfloat("rf_frequency_hz")
    steps_per_cycle = numerics.getint("steps_per_rf_cycle")
    cycles = numerics.getint("rf_cycles")
    steps = steps_per_cycle * cycles
    require(steps == 512000, "Case 1 production duration drifted")
    dt = 1.0 / (frequency * steps_per_cycle)
    cells = numerics.getint("cells")
    nodes = numerics.getint("nodes")
    require(nodes == cells + 1, "node/cell contract is inconsistent")
    particles_per_species = (
        cells * numerics.getint("particles_per_cell_per_species")
    )
    capacity = numerics.getint("max_particles_per_species")
    require(capacity >= particles_per_species, "particle capacity is too small")
    initial_density = physics.getfloat("initial_density_m3")
    length = physics.getfloat("length_m")
    weight = initial_density * length / particles_per_species
    electron_thermal = math.sqrt(
        boltzmann * physics.getfloat("electron_temperature_k") / electron_mass
    )
    ion_thermal = math.sqrt(
        boltzmann * physics.getfloat("ion_temperature_k") / helium_mass
    )
    averaging_cycles = numerics.getint("averaging_rf_cycles")
    averaging_samples = averaging_cycles * steps_per_cycle
    averaging_start = steps - averaging_samples + 1
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else output.parent / "output"
    )
    seed_argument = getattr(args, "seed", None)
    seed = numerics.getint("seed") if seed_argument is None else seed_argument
    require(
        isinstance(seed, int) and 0 <= seed <= 4_294_967_295,
        "seed must be an unsigned 32-bit integer",
    )
    deck = f"""# Generated exact Turner et al. helium CCP benchmark Case 1 deck.
# The preparer audited local restricted inputs and did not launch this run.
config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {nodes}
length = {length:.17g}
dt = {dt:.17g}
steps = {steps}
output_interval = {numerics.getint('output_interval_steps')}
output_dir = {output_dir}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_right_amplitude = {physics.getfloat('voltage_amplitude_v'):.17g}
phi_right_frequency = {frequency:.17g}
seed = {seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = {capacity}
checkpoint_output = true
checkpoint_interval = {numerics.getint('checkpoint_interval_steps')}
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = {averaging_start}
spatial_average_end_step = {steps}
spatial_average_rf_frequency = {frequency:.17g}
spatial_average_rf_cycles = {averaging_cycles}

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = {neutral_density:.17g}
neutral_temperature = {neutral_temperature:.17g}
max_frequency = {electron_majorant:.17g}
max_candidates_per_particle = {guard.getint('max_candidates_per_particle')}
gas_data_file = {electron_gas}

[collisions.electron_mcc.channel.ionization]
secondary_species = electrons
ion_species = ions

[collisions.ion_mcc]
model = null_collision
species = ions
neutral_density = {neutral_density:.17g}
neutral_temperature = {neutral_temperature:.17g}
max_frequency = {ion_majorant:.17g}
max_candidates_per_particle = {guard.getint('max_candidates_per_particle')}
gas_data_file = {ion_gas}

[species.electrons]
charge = {-elementary_charge:.17g}
mass = {electron_mass:.17g}
weight = {weight:.17g}
particles = {particles_per_species}
thermal_velocity = {electron_thermal:.17g}
loading = quiet_start

[species.ions]
charge = {elementary_charge:.17g}
mass = {helium_mass:.17g}
weight = {weight:.17g}
particles = {particles_per_species}
thermal_velocity = {ion_thermal:.17g}
loading = quiet_start
"""
    deck_digest = hashlib.sha256(deck.encode("utf-8")).hexdigest()
    report = {
        "turner_case_preflight_version": 1,
        "case_id": case["global"]["case_id"],
        "physics_claim": authorization["physics_claim_before_complete_comparison"],
        "full_run_launched": False,
        "production_launch_authorized": False,
        "provenance": {
            "case_manifest": str(case_path),
            "case_manifest_sha256": sha256(case_path),
            "normalization_audit": str(normalized_dir / "audit.json"),
            "normalization_audit_sha256": sha256(normalized_dir / "audit.json"),
            "source_artifact_sha256": audit["source_artifact"]["sha256"],
            "generated_deck_sha256": deck_digest,
        },
        "contract": {
            "steps": steps,
            "rf_cycles": cycles,
            "steps_per_rf_cycle": steps_per_cycle,
            "timestep_s": dt,
            "nodes": nodes,
            "particles_per_species": particles_per_species,
            "macro_particle_weight": weight,
            "averaging_start_step": averaging_start,
            "averaging_end_step": steps,
            "averaging_samples": averaging_samples,
            "seed": seed,
        },
        "collision_guard": {
            "method": "sampled production-kernel thermal bound plus safety factor",
            "runtime_behavior": "fail fast if an encountered exact bound exceeds majorant",
            "scan_intervals": intervals,
            "thermal_neutral_sigma_limit": sigma_limit,
            "safety_factor": safety_factor,
            "electron": {
                "energy_ceiling_ev": guard.getfloat("electron_energy_ceiling_ev"),
                "sampled_peak_frequency_s": electron_bound,
                "sampled_peak_speed_m_s": electron_peak_speed,
                "configured_majorant_s": electron_majorant,
                "majorant_dt": electron_majorant * dt,
            },
            "ion": {
                "energy_ceiling_ev": guard.getfloat("ion_energy_ceiling_ev"),
                "sampled_peak_frequency_s": ion_bound,
                "sampled_peak_speed_m_s": ion_peak_speed,
                "configured_majorant_s": ion_majorant,
                "majorant_dt": ion_majorant * dt,
            },
        },
        "resource_floor": {
            "initial_live_particles": 2 * particles_per_species,
            "initial_particle_updates": 2 * particles_per_species * steps,
            "capacity_particles": 2 * capacity,
            "note": "wall time is intentionally not projected without a bounded runtime measurement",
        },
        "reported_case1_characteristics": {
            "source": reported["source"],
            "time_averaged_midplane_ion_density_m3": reported.getfloat(
                "time_averaged_midplane_ion_density_m3"
            ),
            "time_averaged_electron_temperature_ev": reported.getfloat(
                "time_averaged_electron_temperature_ev"
            ),
            "time_averaged_electron_power_w_m2": reported.getfloat(
                "time_averaged_electron_power_w_m2"
            ),
            "time_averaged_ion_power_w_m2": reported.getfloat(
                "time_averaged_ion_power_w_m2"
            ),
            "time_averaged_ion_current_a_m2": reported.getfloat(
                "time_averaged_ion_current_a_m2"
            ),
            "total_macro_particles": reported.getint(
                "total_macro_particles"
            ),
        },
    }
    atomic_text(output, deck)
    atomic_text(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    return output, report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit local Turner data and generate, but never run, Case 1"
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        output, report = prepare(parse_args())
    except PreparationError as error:
        print(f"Turner Case 1 preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Generated production deck without launching it: {output}")
    print(f"Wrote preflight report: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
