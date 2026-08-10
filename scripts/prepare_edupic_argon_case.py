#!/usr/bin/env python3
"""Validate a local eduPIC package and generate a bounded AuroraPIC preflight."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
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
    require(not path.exists(), f"refusing to overwrite existing file: {path}")
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
        raise PreparationError(f"cannot parse case manifest: {error}") from error
    required_sections = (
        "global", "reference", "package", "physics", "numerics",
        "compatibility", "authorization",
    )
    for section in required_sections:
        require(section in parser, f"case manifest is missing [{section}]")
    require(
        parser["global"].getint("case_manifest_version") == 1,
        "unsupported case manifest version",
    )
    require(
        parser["global"]["case_id"] == "edupic-1.0-default-argon-ccp",
        "unexpected case identity",
    )
    require(
        not parser["authorization"].getboolean(
            "production_launch_authorized"
        ),
        "this preparer must not authorize production",
    )
    return parser


def validate_package(
    package_dir: Path, case: configparser.ConfigParser
) -> dict[str, object]:
    audit_path = package_dir / "audit.json"
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PreparationError(f"cannot parse package audit: {error}") from error
    reference = case["reference"]
    require(
        audit.get("case_id") == case["global"]["case_id"],
        "package audit case identity mismatch",
    )
    source = audit.get("source", {})
    require(
        source.get("commit") == reference["commit"]
        and source.get("implementation_sha256")
        == reference["implementation_sha256"]
        and source.get("table_sha256")
        == reference["generated_cross_sections_sha256"],
        "package source identity mismatch",
    )
    artifact_names = {
        "electron_manifest": case["package"]["electron_manifest"],
        "ion_manifest": case["package"]["ion_manifest"],
        "electron_elastic": "electron_elastic.dat",
        "electron_excitation": "electron_excitation.dat",
        "electron_ionization": "electron_ionization.dat",
        "ion_isotropic": "ion_isotropic.dat",
        "ion_backward": "ion_backward.dat",
    }
    audit_artifacts = audit.get("artifacts", {})
    for key, name in artifact_names.items():
        path = package_dir / name
        expected = case["package"][f"{key}_sha256"]
        require(sha256(path) == expected, f"package checksum mismatch: {name}")
        require(
            audit_artifacts.get(name, {}).get("sha256") == expected,
            f"package audit identity mismatch: {name}",
        )
    contract = audit.get("contract", {})
    require(
        contract.get("rows") == 1_000_000
        and contract.get("energy_step_ev") == 0.001
        and contract.get("ionization_kinematics") == "opal_beaty_peterson"
        and contract.get("inelastic_transform")
        == "finite_mass_center_of_mass"
        and contract.get("ion_energy_frame") == "center_of_mass"
        and contract.get("cross_section_interpolation") == "lower_bin",
        "package collision contract mismatch",
    )
    return audit


def prepare(args: argparse.Namespace) -> tuple[Path, Path]:
    case_path = args.case_manifest.resolve()
    package_dir = args.package_dir.resolve()
    output = args.output.resolve()
    report = (
        args.report.resolve()
        if args.report is not None
        else output.with_suffix(".preflight.json")
    )
    case = load_case(case_path)
    audit = validate_package(package_dir, case)
    authorization = case["authorization"]
    maximum_steps = authorization.getint("maximum_preflight_steps")
    require(
        1 <= args.steps <= maximum_steps,
        f"preflight steps must be in [1, {maximum_steps}]",
    )
    physics = case["physics"]
    numerics = case["numerics"]
    frequency = physics.getfloat("rf_frequency_hz")
    steps_per_cycle = numerics.getint("steps_per_rf_cycle")
    ion_timestep_multiplier = numerics.getint(
        "ion_subcycling_reference"
    )
    require(
        ion_timestep_multiplier > 0
        and steps_per_cycle % ion_timestep_multiplier == 0,
        "ion subcycling must be positive and divide one RF cycle",
    )
    dt = 1.0 / (frequency * steps_per_cycle)
    neutral_density = physics.getfloat("neutral_density_m3")
    electron_majorant = numerics.getfloat("electron_max_frequency_s")
    ion_majorant = numerics.getfloat("ion_max_frequency_s")
    rate_envelope = audit["grid_rate_envelope"]
    require(
        electron_majorant > rate_envelope["electron_peak_frequency_s"]
        and ion_majorant > rate_envelope["ion_peak_frequency_s"],
        "configured majorants do not exceed imported grid envelopes",
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else output.parent / "output"
    )
    electron_manifest = package_dir / case["package"]["electron_manifest"]
    ion_manifest = package_dir / case["package"]["ion_manifest"]
    seed = numerics.getint("seed") if args.seed is None else args.seed
    require(0 <= seed <= 4_294_967_295, "seed must be unsigned 32-bit")
    phase = math.pi / 2.0
    deck = f"""# Generated bounded eduPIC 1.0 argon CCP contract preflight.
# This is not a production or validation run and carries no physics claim.
config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = {numerics.getint('nodes')}
length = {physics.getfloat('length_m'):.17g}
dt = {dt:.17g}
steps = {args.steps}
output_interval = {args.steps}
output_dir = {output_dir}
boundary = dirichlet
mode = transient
phi_left = 0
phi_right = 0
phi_left_amplitude = {physics.getfloat('voltage_amplitude_v'):.17g}
phi_left_frequency = {frequency:.17g}
phi_left_phase = {phase:.17g}
seed = {seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = {numerics.getint('max_particles_per_species')}
checkpoint_output = true
checkpoint_interval = {args.steps}

[collisions.electron_mcc]
model = null_collision
species = electrons
neutral_density = {neutral_density:.17g}
neutral_temperature = 0
max_frequency = {electron_majorant:.17g}
max_candidates_per_particle = {numerics.getint('max_candidates_per_particle')}
gas_data_file = {electron_manifest}

[collisions.electron_mcc.channel.ionization]
secondary_species = electrons
ion_species = ions

[collisions.ion_mcc]
model = null_collision
species = ions
neutral_density = {neutral_density:.17g}
neutral_temperature = {physics.getfloat('neutral_temperature_k'):.17g}
max_frequency = {ion_majorant:.17g}
max_candidates_per_particle = {numerics.getint('max_candidates_per_particle')}
gas_data_file = {ion_manifest}

[species.electrons]
charge = {-physics.getfloat('elementary_charge_c'):.17g}
mass = {physics.getfloat('electron_mass_kg'):.17g}
weight = {physics.getfloat('macro_particle_weight'):.17g}
particles = {physics.getint('initial_particles_per_species')}
thermal_velocity = 0
loading = random

[species.ions]
charge = {physics.getfloat('elementary_charge_c'):.17g}
mass = {physics.getfloat('argon_mass_kg'):.17g}
weight = {physics.getfloat('macro_particle_weight'):.17g}
particles = {physics.getint('initial_particles_per_species')}
thermal_velocity = 0
loading = random
timestep_multiplier = {ion_timestep_multiplier}
"""
    compatibility = dict(case["compatibility"])
    unresolved = [
        key for key, value in compatibility.items()
        if value.startswith("unmatched") or value.startswith("failed")
    ]
    preflight = {
        "edupic_argon_preflight_version": 1,
        "case_id": case["global"]["case_id"],
        "physics_claim": authorization["physics_claim"],
        "production_launch_authorized": False,
        "run_launched_by_preparer": False,
        "case_manifest": {
            "path": str(case_path),
            "sha256": sha256(case_path),
        },
        "package_audit": {
            "path": str(package_dir / "audit.json"),
            "sha256": sha256(package_dir / "audit.json"),
            "source_table_sha256": audit["source"]["table_sha256"],
        },
        "deck": {
            "path": str(output),
            "sha256": hashlib.sha256(deck.encode()).hexdigest(),
            "steps": args.steps,
            "rf_cycles": args.steps / steps_per_cycle,
            "timestep_s": dt,
            "seed": seed,
        },
        "compatibility": compatibility,
        "unresolved_contract_items": unresolved,
        "claim_boundary": (
            "A successful bounded execution proves package loading and runtime "
            "integration only, not cross-code agreement or physical validation."
        ),
    }
    atomic_text(output, deck)
    atomic_text(report, json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    return output, report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("package_dir", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--steps", required=True, type=int)
    parser.add_argument("--seed", type=int)
    return parser.parse_args()


def main() -> int:
    try:
        output, report = prepare(parse_args())
    except PreparationError as error:
        print(f"eduPIC argon preflight error: {error}", file=os.sys.stderr)
        return 2
    print(f"Generated bounded deck without launching it: {output}")
    print(f"Wrote preflight report: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
