#!/usr/bin/env python3
"""Validate the pinned reduced LANDMARK case and source-rate derivation."""

from __future__ import annotations

import configparser
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_with_global(path: Path) -> configparser.ConfigParser:
    text = path.read_text(encoding="utf-8")
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string("[global]\n" + text)
    return parser


def number(section: configparser.SectionProxy, key: str) -> float:
    value = section.getfloat(key)
    require(math.isfinite(value), f"{section.name}.{key} must be finite")
    return value


def main() -> int:
    manifest = load_with_global(MANIFEST)
    global_section = manifest["global"]
    require(
        global_section.getint("case_manifest_version") == 1,
        "unsupported Hall case manifest version",
    )
    require(
        global_section["status"] == "reduced_integration_only",
        "Hall case must not claim production status",
    )

    reference = manifest["reference"]
    source = manifest["pair_source"]
    width = number(source, "x_max_m") - number(source, "x_min_m")
    height = number(source, "y_max_m") - number(source, "y_min_m")
    mean = number(source, "normalized_profile_mean")
    area = width * height * mean
    require(
        math.isclose(
            area, number(source, "effective_profile_area_m2"),
            rel_tol=0.0, abs_tol=1e-18,
        ),
        "Hall source effective area is inconsistent",
    )
    total_rate = (
        number(source, "peak_volumetric_pair_rate_m3_s")
        * area
        * number(source, "out_of_plane_depth_m")
    )
    require(
        math.isclose(
            total_rate,
            number(source, "derived_represented_pair_rate_s"),
            rel_tol=1e-14,
        ),
        "Hall source peak-to-integral conversion is inconsistent",
    )

    magnetic = manifest["magnetic_field"]
    magnetic_path = MANIFEST.parent / magnetic["file"]
    digest = hashlib.sha256(magnetic_path.read_bytes()).hexdigest()
    require(digest == magnetic["sha256"], "Hall magnetic profile checksum mismatch")

    runtime_path = MANIFEST.parent / global_section["runtime_config"]
    runtime = load_with_global(runtime_path)
    runtime_global = runtime["global"]
    runtime_source = runtime["source.channel_pair_seed"]
    cathode = manifest["cathode_control"]
    emitted_species = runtime[
        "species." + cathode["emitted_species"]
    ]
    require(
        math.isclose(
            runtime_global.getfloat("length_x"),
            reference.getfloat("domain_x_m"),
        )
        and math.isclose(
            runtime_global.getfloat("length_y"),
            reference.getfloat("domain_y_m"),
        )
        and math.isclose(
            runtime_global.getfloat("out_of_plane_depth"),
            number(source, "out_of_plane_depth_m"),
        ),
        "Hall runtime geometry drifted from its manifest",
    )
    require(
        math.isclose(
            runtime_global.getfloat("dt"),
            reference.getfloat("production_dt_s"),
        )
        and runtime_global["boundary_x"] == "dirichlet"
        and runtime_global["boundary_y"] == "periodic"
        and math.isclose(
            runtime_global.getfloat("phi_left"),
            reference.getfloat("discharge_voltage_v"),
        )
        and math.isclose(runtime_global.getfloat("phi_right"), 0.0),
        "Hall runtime timestep or field boundary contract drifted",
    )
    require(
        runtime_global["magnetic_field_profile_file"]
            == magnetic["file"]
        and runtime_global["magnetic_field_profile_axis"]
            == magnetic["axis"],
        "Hall runtime magnetic profile linkage drifted",
    )
    require(
        runtime_global["current_source_species"]
            == cathode["emitted_species"]
        and runtime_global["current_source_monitor_boundary"]
            == cathode["monitor_boundary"]
        and runtime_global["current_source_emission_boundary"]
            == cathode["emission_boundary"]
        and math.isclose(
            runtime_global.getfloat("length_x")
            - runtime_global.getfloat(
                "current_source_emission_inset"
            ),
            number(cathode, "emission_plane_x_m"),
            rel_tol=0.0, abs_tol=1e-15,
        )
        and math.isclose(
            runtime_global.getfloat(
                "current_source_thermal_velocity"
            ),
            number(
                cathode,
                "emission_thermal_velocity_std_m_s",
            ),
            rel_tol=1e-14,
        ),
        "Hall runtime cathode current-control linkage drifted",
    )
    expected_cathode_thermal_velocity = math.sqrt(
        number(cathode, "emission_temperature_ev")
        * abs(emitted_species.getfloat("charge"))
        / emitted_species.getfloat("mass")
    )
    require(
        math.isclose(
            expected_cathode_thermal_velocity,
            number(
                cathode,
                "emission_thermal_velocity_std_m_s",
            ),
            rel_tol=1e-14,
        ),
        "Hall cathode temperature-to-velocity conversion drifted",
    )
    require(
        runtime_global["potential_reference_axis"]
            == cathode["potential_reference_axis"]
        and math.isclose(
            runtime_global.getfloat(
                "potential_reference_coordinate"
            ),
            number(
                cathode,
                "potential_reference_coordinate_m",
            ),
        )
        and math.isclose(
            runtime_global.getfloat(
                "potential_reference_target"
            ),
            number(cathode, "potential_reference_target_v"),
        ),
        "Hall runtime potential-reference linkage drifted",
    )
    for key, manifest_key in (
        ("x_min", "x_min_m"),
        ("x_max", "x_max_m"),
        ("y_min", "y_min_m"),
        ("y_max", "y_max_m"),
        ("peak_volumetric_pair_rate", "peak_volumetric_pair_rate_m3_s"),
    ):
        require(
            math.isclose(
                runtime_source.getfloat(key), number(source, manifest_key),
                rel_tol=1e-14, abs_tol=0.0,
            ),
            f"Hall runtime source {key} drifted from its manifest",
        )
    require(
        runtime_source["density_profile"] == "sinusoidal"
        and math.isclose(
            runtime_source.getfloat("profile_amplitude"), -1.0,
        )
        and runtime_source.getint("profile_mode_x") == 1,
        "Hall runtime source profile drifted from the pinned sine-squared envelope",
    )
    reduced = manifest["reduced_contract"]
    require(
        runtime_global.getint("nx") * runtime_global.getint("ny")
        <= reduced.getint("max_cells"),
        "Hall reduced runtime exceeds its cell budget",
    )
    require(
        runtime_global.getint("steps") <= reduced.getint("max_steps")
        and runtime_global.getint("runtime_threads")
        <= reduced.getint("max_threads"),
        "Hall reduced runtime exceeds its step/thread budget",
    )
    initial_macroparticles = sum(
        runtime[section].getint("particles")
        for section in runtime.sections()
        if section.startswith("species.")
    )
    require(
        initial_macroparticles
        <= reduced.getint("max_initial_macroparticles"),
        "Hall reduced runtime exceeds its initial macro-particle budget",
    )
    require(
        reduced["physics_claim"] == "none"
        and "published_resolution" in reduced["missing_physics"]
        and "cathode_current_control"
            not in reduced["missing_physics"]
        and "potential_correction"
            not in reduced["missing_physics"],
        "Hall reduced manifest must retain its no-claim limitations",
    )
    print("Hall LANDMARK reduced-case validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
