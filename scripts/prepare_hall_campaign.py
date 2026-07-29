#!/usr/bin/env python3
"""Generate a resource-tiered LANDMARK Case 2 campaign deck without running it."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tempfile


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_RUN"
WORKSTATION_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_AN_OPT_IN_WORKSTATION_RUN"


class CampaignError(RuntimeError):
    pass


def load_manifest(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise CampaignError(f"cannot read {path}: {error}") from error
    for section in ("global", "reference", "magnetic_field", "pair_source",
                    "cathode_control", "diagnostics"):
        if section not in parser:
            raise CampaignError(f"manifest is missing [{section}]")
    return parser


def atomic_text(path: Path, text: str) -> None:
    if path.exists():
        raise CampaignError(f"refusing to overwrite existing file: {path}")
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


def prepare(args: argparse.Namespace) -> Path:
    manifest_path = args.case_manifest.resolve()
    manifest = load_manifest(manifest_path)
    reference = manifest["reference"]
    magnetic = manifest["magnetic_field"]
    source = manifest["pair_source"]
    cathode = manifest["cathode_control"]
    diagnostics = manifest["diagnostics"]
    tier_name = f"campaign.{args.tier}"
    if tier_name not in manifest:
        raise CampaignError(f"manifest is missing [{tier_name}]")
    tier = manifest[tier_name]
    authorization = tier["authorization"]
    required_acknowledgement = {
        "micro": None,
        "workstation": WORKSTATION_ACKNOWLEDGEMENT,
        "production": ACKNOWLEDGEMENT,
    }[args.tier]
    if (authorization == "explicit_cost_acknowledgement"
            and args.acknowledge_cost != required_acknowledgement):
        raise CampaignError(
            f"{args.tier} deck generation requires --acknowledge-cost "
            f"{required_acknowledgement}"
        )
    cells_x = tier.getint("cells_x")
    cells_y = tier.getint("cells_y")
    nodes_x = cells_x + 1
    nodes_y = cells_y
    if args.tier == "production" and (
        cells_x != reference.getint("production_cells_x")
        or cells_y != reference.getint("production_cells_y")
        or nodes_x != reference.getint("aurorapic_nodes_x")
        or nodes_y != reference.getint("aurorapic_nodes_y")
    ):
        raise CampaignError("production tier drifted from original LANDMARK Case 2")
    particles_per_cell = tier.getint("particles_per_cell_per_species")
    particles_per_species = cells_x * cells_y * particles_per_cell
    max_particles_per_species = tier.getint("max_particles_per_species")
    if max_particles_per_species < particles_per_species:
        raise CampaignError("particle capacity is below the initial population")

    source_asset = manifest_path.parent / magnetic["file"]
    try:
        digest = hashlib.sha256(source_asset.read_bytes()).hexdigest()
    except OSError as error:
        raise CampaignError(f"cannot read magnetic profile: {error}") from error
    if digest != magnetic["sha256"]:
        raise CampaignError("magnetic profile checksum mismatch")

    output = args.output.resolve()
    if output.exists():
        raise CampaignError(f"refusing to overwrite existing file: {output}")
    asset = output.parent / source_asset.name
    if asset.exists():
        try:
            asset_digest = hashlib.sha256(asset.read_bytes()).hexdigest()
        except OSError as error:
            raise CampaignError(f"cannot inspect existing asset: {error}") from error
        if asset_digest != digest:
            raise CampaignError(f"existing campaign asset has wrong checksum: {asset}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_asset, asset)

    start_step = tier.getint("diagnostic_start_step")
    steps = tier.getint("steps")
    if start_step > steps:
        raise CampaignError("diagnostic start step exceeds the production duration")
    checkpoint_interval = tier.getint("checkpoint_interval")
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else f"output_hall_landmark_{args.tier}"
    )
    deck = f"""# Generated LANDMARK Case 2 {args.tier} campaign tier.
# Purpose: {tier['purpose']}; physics claim: {tier['physics_claim']}.
# Generation authorization: {authorization}; no run has been launched.
# Selected grid: {cells_x} x {cells_y} cells.
# AuroraPIC structured grid: {nodes_x} x {nodes_y} nodes.
config_version = 1
dimension = 2
units = si
nx = {nodes_x}
ny = {nodes_y}
length_x = {reference['domain_x_m']}
length_y = {reference['domain_y_m']}
out_of_plane_depth = {source['out_of_plane_depth_m']}
dt = {reference['production_dt_s']}
steps = {tier['steps']}
mode = transient
boundary = dirichlet
boundary_x = dirichlet
boundary_y = periodic
phi_left = {reference['discharge_voltage_v']}
phi_right = 0
current_source_species = {cathode['emitted_species']}
current_source_control_mode = {cathode['control_mode']}
current_source_monitor_boundary = {cathode['monitor_boundary']}
current_source_emission_boundary = {cathode['emission_boundary']}
current_source_emission_inset = 0.001
current_source_temperature_ev = {cathode['emission_temperature_ev']}
potential_reference_axis = {cathode['potential_reference_axis']}
potential_reference_coordinate = {cathode['potential_reference_coordinate_m']}
potential_reference_target = {cathode['potential_reference_target_v']}
potential_reference_correction = {cathode['potential_reference_correction']}
output_interval = {tier['diagnostic_interval']}
output_dir = {output_dir}
resolved_diagnostics = true
resolved_diagnostic_interval = {tier['diagnostic_interval']}
resolved_diagnostic_start_step = {start_step}
resolved_profile_axis = {diagnostics['profile_axis']}
resolved_mode_axis = {diagnostics['mode_axis']}
resolved_max_mode = {tier['max_mode']}
checkpoint_output = {'true' if checkpoint_interval > 0 else 'false'}
checkpoint_interval = {checkpoint_interval}
vtk_output = false
particle_output = false
magnetic_field_profile_file = {asset.name}
magnetic_field_profile_axis = {magnetic['axis']}
seed = {args.seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = {max_particles_per_species}

[species.electrons]
charge = -1.602176634e-19
mass = 9.1093837139e-31
density = {reference['initial_density_m3']}
particles = {particles_per_species}
temperature_ev = 10
loading = random

[species.ions]
charge = 1.602176634e-19
mass = 2.1801714e-25
density = {reference['initial_density_m3']}
particles = {particles_per_species}
temperature_ev = 0.5
loading = random

[source.channel_pair_seed]
first_species = electrons
second_species = ions
peak_volumetric_pair_rate = {source['peak_volumetric_pair_rate_m3_s']}
x_min = {source['x_min_m']}
x_max = {source['x_max_m']}
y_min = {source['y_min_m']}
y_max = {source['y_max_m']}
first_temperature_ev = {source['electron_temperature_ev']}
second_temperature_ev = {source['ion_temperature_ev']}
density_profile = sinusoidal
profile_amplitude = -1
profile_mode_x = 1
"""
    atomic_text(output, deck)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate, but never launch, a LANDMARK Case 2 deck"
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--tier",
        choices=("micro", "workstation", "production"),
        default="production",
    )
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument(
        "--acknowledge-cost",
        "--acknowledge-production-cost",
        dest="acknowledge_cost",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        path = prepare(args)
    except CampaignError as error:
        print(f"Hall campaign preparation error: {error}", file=sys.stderr)
        return 2
    print(
        f"Generated {args.tier} campaign deck without launching it: {path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
