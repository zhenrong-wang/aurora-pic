#!/usr/bin/env python3
"""Generate an explicitly authorized, published-scale LANDMARK Case 2 deck."""

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
    if args.acknowledge_production_cost != ACKNOWLEDGEMENT:
        raise CampaignError(
            "production deck generation requires "
            f"--acknowledge-production-cost {ACKNOWLEDGEMENT}"
        )
    manifest_path = args.case_manifest.resolve()
    manifest = load_manifest(manifest_path)
    reference = manifest["reference"]
    magnetic = manifest["magnetic_field"]
    source = manifest["pair_source"]
    cathode = manifest["cathode_control"]
    diagnostics = manifest["diagnostics"]

    cells_x = reference.getint("production_cells_x")
    cells_y = reference.getint("production_cells_y")
    nodes_x = reference.getint("aurorapic_nodes_x")
    nodes_y = reference.getint("aurorapic_nodes_y")
    if (cells_x, cells_y, nodes_x, nodes_y) != (500, 256, 501, 256):
        raise CampaignError("manifest does not describe original LANDMARK Case 2")
    particles_per_species = cells_x * cells_y * args.particles_per_cell
    if args.max_particles_per_species < particles_per_species:
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

    start_step = args.diagnostic_start_step
    steps = reference.getint("production_steps")
    if start_step > steps:
        raise CampaignError("diagnostic start step exceeds the production duration")
    deck = f"""# Generated LANDMARK Case 2 production candidate.
# Generation was explicitly acknowledged; no run has been launched.
# Original reference grid: {cells_x} x {cells_y} cells.
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
steps = {steps}
mode = transient
boundary = dirichlet
boundary_x = dirichlet
boundary_y = periodic
phi_left = {reference['discharge_voltage_v']}
phi_right = 0
current_source_species = {cathode['emitted_species']}
current_source_monitor_boundary = {cathode['monitor_boundary']}
current_source_emission_boundary = {cathode['emission_boundary']}
current_source_emission_inset = 0.001
current_source_temperature_ev = {cathode['emission_temperature_ev']}
potential_reference_axis = {cathode['potential_reference_axis']}
potential_reference_coordinate = {cathode['potential_reference_coordinate_m']}
potential_reference_target = {cathode['potential_reference_target_v']}
output_interval = {args.diagnostic_interval}
output_dir = {args.output_dir}
resolved_diagnostics = true
resolved_diagnostic_interval = {args.diagnostic_interval}
resolved_diagnostic_start_step = {start_step}
resolved_profile_axis = {diagnostics['profile_axis']}
resolved_mode_axis = {diagnostics['mode_axis']}
resolved_max_mode = {args.max_mode}
checkpoint_output = true
checkpoint_interval = {args.checkpoint_interval}
vtk_output = false
particle_output = false
magnetic_field_profile_file = {asset.name}
magnetic_field_profile_axis = {magnetic['axis']}
seed = {args.seed}
runtime_backend = serial
runtime_threads = 1
max_particles_per_species = {args.max_particles_per_species}

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
    parser.add_argument("--output-dir", default="output_hall_landmark_case2")
    parser.add_argument("--particles-per-cell", type=int, default=75)
    parser.add_argument("--max-particles-per-species", type=int, default=80_000_000)
    parser.add_argument("--diagnostic-start-step", type=int, default=3_200_000)
    parser.add_argument("--diagnostic-interval", type=int, default=5000)
    parser.add_argument("--checkpoint-interval", type=int, default=500_000)
    parser.add_argument("--max-mode", type=int, default=128)
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--acknowledge-production-cost")
    args = parser.parse_args()
    for name in ("particles_per_cell", "max_particles_per_species",
                 "diagnostic_interval", "checkpoint_interval"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.diagnostic_start_step < 0 or args.max_mode < 0:
        parser.error("diagnostic start and max mode must be non-negative")
    if args.max_mode > 128:
        parser.error("--max-mode exceeds the 256-cell azimuthal Nyquist limit")
    return args


def main() -> int:
    try:
        path = prepare(parse_args())
    except CampaignError as error:
        print(f"Hall campaign preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Generated production candidate without launching it: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
