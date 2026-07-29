#!/usr/bin/env python3
"""Validate the pinned reduced LANDMARK case and source-rate derivation."""

from __future__ import annotations

import configparser
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
ELECTRON_VOLT_J = 1.602176634e-19


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
    require(
        reference.getint("production_cells_x") == 500
        and reference.getint("production_cells_y") == 256
        and reference.getint("aurorapic_nodes_x") == 501
        and reference.getint("aurorapic_nodes_y") == 256,
        "Hall production cell/node topology drifted from LANDMARK Case 2",
    )
    require(
        reference.getint("aurorapic_nodes_x")
            == reference.getint("production_cells_x") + 1
        and reference.getint("aurorapic_nodes_y")
            == reference.getint("production_cells_y"),
        "Hall Dirichlet/periodic cell-to-node mapping is inconsistent",
    )
    provenance = manifest["provenance"]
    require(
        provenance["article_url"]
            == "https://doi.org/10.1088/1361-6595/ab46c5"
        and provenance["public_dataset_doi"] == "10.7302/5mfm-as86"
        and provenance["public_dataset_file_set_id"] == "m900nv362"
        and "512 axial by 256 azimuthal cells"
            in provenance["public_dataset_variant"]
        and provenance["public_dataset_license"] == "CC0 1.0",
        "Hall article/public-dataset provenance drifted",
    )
    expected_tiers = {
        "micro": (32, 16, 4, 200, 8_192, "none"),
        "workstation": (125, 64, 16, 5_000, 250_000, "none"),
        "production": (
            500, 256, 75, 4_000_000, 80_000_000,
            "candidate_for_reference_comparison",
        ),
    }
    previous_updates = 0
    for tier_name, expected in expected_tiers.items():
        tier = manifest[f"campaign.{tier_name}"]
        actual = (
            tier.getint("cells_x"),
            tier.getint("cells_y"),
            tier.getint("particles_per_cell_per_species"),
            tier.getint("steps"),
            tier.getint("max_particles_per_species"),
            tier["physics_claim"],
        )
        require(
            actual == expected,
            f"Hall {tier_name} campaign tier drifted",
        )
        initial_per_species = actual[0] * actual[1] * actual[2]
        require(
            actual[4] >= initial_per_species
            and tier.getint("diagnostic_start_step") <= actual[3]
            and tier.getint("max_mode") <= actual[1] // 2,
            f"Hall {tier_name} campaign bounds are inconsistent",
        )
        updates = initial_per_species * 2 * actual[3]
        require(
            updates > previous_updates,
            "Hall campaign tiers must increase monotonically in work",
        )
        previous_updates = updates
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
    initial = manifest["initial_loading"]
    diagnostics = manifest["diagnostics"]
    comparison = manifest["comparison"]
    source_registry_path = MANIFEST.parent / comparison["source_registry"]
    source_registry_digest = hashlib.sha256(
        source_registry_path.read_bytes()
    ).hexdigest()
    source_registry = load_with_global(source_registry_path)
    require(
        source_registry["global"].getint("source_registry_version") == 1
        and source_registry["global"]["case_id"]
            == global_section["case_id"]
        and source_registry_digest
            == comparison["source_registry_sha256"],
        "Hall source registry linkage drifted",
    )
    original_source = source_registry["source.original_supplement"]
    warpx_source = source_registry["source.warpx_deepblue"]
    require(
        original_source["variant"] == "original-landmark-case2-500x256"
        and original_source.getint("cells_x") == 500
        and original_source.getint("cells_y") == 256
        and original_source["doi"] == reference["doi"]
        and warpx_source["variant"] == "warpx-case2-512x256"
        and warpx_source.getint("cells_x") == 512
        and warpx_source.getint("cells_y") == 256
        and warpx_source["doi"] == provenance["public_dataset_doi"]
        and warpx_source["license"] == "CC0-1.0"
        and warpx_source["artifact_name"] == "baseline_20us.tar"
        and warpx_source["file_set_id"] == "m900nv362"
        and warpx_source["acquisition"] == "external_globus",
        "Hall reference source identities drifted",
    )
    emitted_species = runtime[
        "species." + cathode["emitted_species"]
    ]
    runtime_electrons = runtime["species.electrons"]
    runtime_ions = runtime["species.ions"]
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
                "current_source_temperature_ev"
            ),
            number(cathode, "emission_temperature_ev"),
        ),
        "Hall runtime cathode current-control linkage drifted",
    )
    expected_cathode_thermal_velocity = math.sqrt(
        number(cathode, "emission_temperature_ev")
        * ELECTRON_VOLT_J
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
    for runtime_species, prefix in (
        (runtime_electrons, "electron"),
        (runtime_ions, "ion"),
    ):
        temperature = number(
            initial, f"{prefix}_temperature_ev")
        velocity = math.sqrt(
            temperature
            * ELECTRON_VOLT_J
            / runtime_species.getfloat("mass")
        )
        require(
            math.isclose(
                runtime_species.getfloat("temperature_ev"),
                temperature,
            )
            and math.isclose(
                velocity,
                number(
                    initial,
                    f"{prefix}_thermal_velocity_std_m_s",
                ),
                rel_tol=1e-14,
            ),
            f"Hall initial {prefix} thermal contract drifted",
        )
    require(
        math.isclose(
            runtime_source.getfloat("first_temperature_ev"),
            number(source, "electron_temperature_ev"),
        )
        and math.isclose(
            runtime_source.getfloat("second_temperature_ev"),
            number(source, "ion_temperature_ev"),
        ),
        "Hall pair-source temperature contract drifted",
    )
    for runtime_species, prefix in (
        (runtime_electrons, "electron"),
        (runtime_ions, "ion"),
    ):
        velocity = math.sqrt(
            number(source, f"{prefix}_temperature_ev")
            * ELECTRON_VOLT_J
            / runtime_species.getfloat("mass")
        )
        require(
            math.isclose(
                velocity,
                number(
                    source,
                    f"{prefix}_thermal_velocity_std_m_s",
                ),
                rel_tol=1e-14,
            ),
            f"Hall pair-source {prefix} temperature conversion drifted",
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
    require(
        runtime_global.getboolean("resolved_diagnostics")
        and runtime_global.getint("resolved_diagnostic_interval") == 1
        and runtime_global.getint("resolved_diagnostic_start_step") == 0
        and runtime_global["resolved_profile_axis"]
            == diagnostics["profile_axis"]
        and runtime_global["resolved_mode_axis"]
            == diagnostics["mode_axis"]
        and runtime_global.getint("resolved_max_mode")
            == diagnostics.getint("reduced_max_mode"),
        "Hall runtime resolved-diagnostic linkage drifted",
    )
    require(
        diagnostics["field_profiles"] == "resolved_field_profiles.csv"
        and diagnostics["species_profiles"]
            == "resolved_species_profiles.csv"
        and diagnostics["mode_history"] == "resolved_modes.csv"
        and diagnostics["field_time_average"]
            == "resolved_field_time_average.csv"
        and diagnostics["species_time_average"]
            == "resolved_species_time_average.csv",
        "Hall resolved-diagnostic artifact contract drifted",
    )
    require(
        comparison.getint("reference_contract_version") == 1
        and comparison["reference_data_policy"]
            == "external_checksum_pinned"
        and comparison["comparator"] == "scripts/compare_hall.py"
        and comparison["preflight"] == "scripts/preflight_hall.py"
        and comparison["runtime_qualifier"]
            == "scripts/qualify_hall_runtime.py"
        and comparison["source_locker"] == "scripts/lock_hall_source.py"
        and comparison["reference_normalizer"]
            == "scripts/normalize_hall_reference.py"
        and comparison["ensemble_preparer"]
            == "scripts/prepare_hall_ensemble.py"
        and comparison["ensemble_aggregator"]
            == "scripts/aggregate_hall_ensemble.py",
        "Hall comparison/preflight contract drifted",
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
        and "timestep_local_cathode_control"
            in reduced["missing_physics"]
        and "affine_potential_correction"
            in reduced["missing_physics"]
        and "published_pair_thermal_loading"
            not in reduced["missing_physics"]
        and "published_initial_thermal_loading"
            not in reduced["missing_physics"],
        "Hall reduced manifest must retain its no-claim limitations",
    )
    print("Hall LANDMARK reduced-case validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
