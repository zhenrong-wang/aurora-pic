#!/usr/bin/env python3
"""Audit frozen eduPIC/AuroraPIC collision kernels without evolving plasma."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import hashlib
import json
import math
from pathlib import Path
from typing import Iterator, TextIO


PINNED_EDUPIC_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41")
EXPECTED_COLLISION_CPP_SHA256 = (
    "6a41a2ecf468847e384c4dc8915a58c76a0993c21f5f9b4adc6ccf356922e750")
EXPECTED_HASHES = {
    "electron_elastic.dat": "f660d0a7665e8edcf591162734c9173abc4c1f93c7f3c623b1ddd5ed27545f47",
    "electron_excitation.dat": "9d4ef324bc094a4ce77959dfc2e1d65dbfa8b88f9e8d06eb4a598d7073e2c96c",
    "electron_ionization.dat": "419958d75e53776ced9f8b81ff77518bf5fc5a18779d167e7231d752d1d9e7e0",
    "ion_isotropic.dat": "44da100ce415eac603bc0625d278eb0088cd5b4d707c56ff5d215d4e4327fb17",
    "ion_backward.dat": "c7cc0ce959cd73b343f7ec62ed38243f302c9db77e7eccd2bda8d3adf945f44b",
    "edupic_argon_electron.gas": "bcf8773f8f392acb256480390b8576aaa322bbd8e1bfaad2ff958d93a3665bb6",
    "edupic_argon_ion.gas": "b7ac195eae2f5a93a31d57663612f0ec464ed736f61ab28e0221864714f72731",
}
ELECTRON_MASS_KG = 9.1093835599999998e-31
ARGON_MASS_KG = 6.6335209000000003e-26
ELEMENTARY_CHARGE_C = 1.6021766200000001e-19
NEUTRAL_DENSITY_M3 = 2.0694208669001848e21
DT_E_S = 1.8436578171091445e-11
ENERGY_STEP_EV = 0.001


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_rows(stream: TextIO) -> Iterator[tuple[float, float]]:
    for line in stream:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        energy, cross_section = map(float, stripped.split())
        yield energy, cross_section


def opportunity_excess(max_frequency: float, timestep: float) -> float:
    value = max_frequency * timestep
    return (value - (1.0 - math.exp(-value))) / value


def scan_lookup_kernel(paths: list[Path], projectile_mass: float,
                       energy_mass: float, timestep: float,
                       energy_limit_ev: float = 80.0) -> dict[str, float | int]:
    with ExitStack() as stack:
        iterators = [data_rows(stack.enter_context(path.open(
            encoding="utf-8"))) for path in paths]
        previous = [next(iterator) for iterator in iterators]
        rows = 1
        compared = 0
        maximum_probability = 0.0
        maximum_probability_energy = 0.0
        maximum_relative_rate = 0.0
        maximum_relative_rate_energy = 0.0
        maximum_channel_tv = 0.0
        maximum_channel_tv_energy = 0.0
        probability_sum = 0.0
        for current in zip(*iterators, strict=True):
            rows += 1
            energies = [item[0] for item in current]
            if any(abs(energy - energies[0]) > 5e-10 for energy in energies):
                raise AuditError("cross-section channel grids differ")
            if any(abs(item[0] - previous[0][0]) > 5e-10
                   for item in previous):
                raise AuditError("cross-section channel grids differ")
            if abs(energies[0] - previous[0][0] - ENERGY_STEP_EV) > 5e-10:
                raise AuditError("cross-section energy spacing differs")
            query_ev = previous[0][0] + 0.75 * ENERGY_STEP_EV
            if query_ev > energy_limit_ev:
                break
            # At the upper quarter of a grid interval, eduPIC's nearest-bin
            # lookup selects the next row while AuroraPIC lower-bin selects
            # the previous row. This is the maximum lookup disagreement in
            # every interval.
            native_sigma = [item[1] for item in current]
            aurora_sigma = [item[1] for item in previous]
            relative_speed = math.sqrt(
                2.0 * query_ev * ELEMENTARY_CHARGE_C / energy_mass)
            native_rates = [NEUTRAL_DENSITY_M3 * value * relative_speed
                            for value in native_sigma]
            aurora_rates = [NEUTRAL_DENSITY_M3 * value * relative_speed
                            for value in aurora_sigma]
            native_total = sum(native_rates)
            aurora_total = sum(aurora_rates)
            native_probability = 1.0 - math.exp(-native_total * timestep)
            aurora_probability = 1.0 - math.exp(-aurora_total * timestep)
            probability_difference = abs(
                aurora_probability - native_probability)
            probability_sum += probability_difference
            if probability_difference > maximum_probability:
                maximum_probability = probability_difference
                maximum_probability_energy = query_ev
            if native_total > 0.0:
                relative_rate = abs(aurora_total / native_total - 1.0)
                if relative_rate > maximum_relative_rate:
                    maximum_relative_rate = relative_rate
                    maximum_relative_rate_energy = query_ev
                if aurora_total > 0.0:
                    tv = 0.5 * sum(abs(
                        a / aurora_total - n / native_total)
                        for a, n in zip(
                            aurora_rates, native_rates, strict=True))
                    if tv > maximum_channel_tv:
                        maximum_channel_tv = tv
                        maximum_channel_tv_energy = query_ev
            compared += 1
            previous = list(current)
    if compared == 0:
        raise AuditError("no cross-section intervals were compared")
    return {
        "source_rows_visited": rows,
        "upper_half_bin_queries": compared,
        "energy_limit_eV": energy_limit_ev,
        "maximum_absolute_event_probability_difference_per_step":
            maximum_probability,
        "energy_at_maximum_probability_difference_eV":
            maximum_probability_energy,
        "mean_absolute_event_probability_difference_per_step_uniform_energy":
            probability_sum / compared,
        "maximum_relative_total_rate_difference": maximum_relative_rate,
        "energy_at_maximum_relative_total_rate_difference_eV":
            maximum_relative_rate_energy,
        "maximum_channel_probability_total_variation": maximum_channel_tv,
        "energy_at_maximum_channel_probability_total_variation_eV":
            maximum_channel_tv_energy,
        "projectile_mass_kg": projectile_mass,
        "collision_energy_mass_kg": energy_mass,
        "timestep_s": timestep,
    }


def require_fragments(source: str, fragments: tuple[str, ...], label: str) -> None:
    if any(fragment not in source for fragment in fragments):
        raise AuditError(f"{label} source contract differs")


def audit(source: Path, collision_cpp: Path, package: Path) -> dict[str, object]:
    if sha256(source) != PINNED_EDUPIC_SHA256:
        raise AuditError("pinned eduPIC source hash differs")
    if sha256(collision_cpp) != EXPECTED_COLLISION_CPP_SHA256:
        raise AuditError("AuroraPIC collision implementation hash differs")
    for name, expected in EXPECTED_HASHES.items():
        if sha256(package / name) != expected:
            raise AuditError(f"locked gas artifact differs: {name}")
    native = source.read_text(encoding="utf-8")
    aurora = collision_cpp.read_text(encoding="utf-8")
    require_fragments(native, (
        "int(energy / DE_CS + 0.5)",
        "p_coll = 1 - exp(- nu * DT_E)",
        "energy = fabs(energy - E_EXC_TH * EV_TO_J)",
        "10.0 * tan(R01(MTgen) * atan(energy/EV_TO_J / 20.0))",
        "eta2 = eta + PI",
    ), "eduPIC electron collision")
    require_fragments(aurora, (
        "available_energy / (2.0 * ejected_energy_scale)",
        "primary_energy / available_energy",
        "secondary_energy / available_energy",
        "channel.config.inelastic_transform ==\n            InelasticTransformKind::FiniteMassCenterOfMass",
        "-std::log(open_unit_interval(rng)) / config_.max_frequency",
    ), "AuroraPIC collision")
    electron = scan_lookup_kernel([
        package / "electron_elastic.dat",
        package / "electron_excitation.dat",
        package / "electron_ionization.dat",
    ], ELECTRON_MASS_KG, ELECTRON_MASS_KG, DT_E_S)
    reduced_argon_mass = 0.5 * ARGON_MASS_KG
    ion = scan_lookup_kernel([
        package / "ion_isotropic.dat",
        package / "ion_backward.dat",
    ], ARGON_MASS_KG, reduced_argon_mass, 20.0 * DT_E_S)
    electron_bound = opportunity_excess(1.0e9, DT_E_S)
    ion_bound = opportunity_excess(1.0e8, 20.0 * DT_E_S)
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "frozen_collision_kernel_contract_audit",
        "provenance": {
            "pinned_edupic_source_sha256": sha256(source),
            "aurorapic_collision_cpp_sha256": sha256(collision_cpp),
            "gas_artifact_sha256": EXPECTED_HASHES,
            "analyzer_sha256": sha256(Path(__file__)),
        },
        "contract_equivalence": {
            "electron_energy_frame": "laboratory",
            "ion_energy_frame": "center_of_mass",
            "elastic_scattering": "isotropic_finite_mass_center_of_mass",
            "excitation_loss_eV": 11.5,
            "ionization_loss_eV": 15.8,
            "ionization_energy_partition": "opal_beaty_peterson_scale_10_eV",
            "ionization_pair_azimuth": "opposed_by_pi",
            "ion_isotropic_and_backward_channels": True,
        },
        "known_algorithm_differences": {
            "cross_section_lookup": {
                "edupic": "nearest 0.001 eV bin",
                "aurorapic": "lower 0.001 eV bin",
            },
            "collision_opportunities": {
                "edupic": "at most one Bernoulli collision per species step",
                "aurorapic": "Poisson null-collision clock permits repeats",
                "electron_maximum_relative_mean_excess": electron_bound,
                "ion_maximum_relative_mean_excess": ion_bound,
            },
            "random_stream": {
                "edupic": "std::mt19937",
                "aurorapic": "std::mt19937_64",
                "consequence": "trajectory-level random draws cannot be paired by seed",
            },
        },
        "upper_half_bin_worst_lookup_audit_0_to_80_eV": {
            "electron": electron,
            "ion": ion,
        },
        "assessment": {
            "source_kinematics_contracts_match": True,
            "lookup_probability_difference_below_1e_4_per_step":
                max(float(electron[
                    "maximum_absolute_event_probability_difference_per_step"]),
                    float(ion[
                    "maximum_absolute_event_probability_difference_per_step"])) < 1e-4,
            "opportunity_mean_difference_below_2_percent":
                max(electron_bound, ion_bound) < 0.02,
            "finding": (
                "The frozen collision contracts match in energy frames, losses, "
                "angular models, finite-mass transforms, and Opal partitioning. "
                "Remaining lookup and opportunity differences are quantified; "
                "a statistical product-kinematics ensemble is still required."),
        },
        "claim_boundary": (
            "This source and frozen-rate audit does not compare paired random "
            "trajectories or prove equivalence of evolved EEDFs."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("edupic_source", type=Path)
    parser.add_argument("aurorapic_collision_cpp", type=Path)
    parser.add_argument("gas_package", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        result = audit(args.edupic_source.resolve(),
                       args.aurorapic_collision_cpp.resolve(),
                       args.gas_package.resolve())
    except (AuditError, OSError, ValueError, StopIteration) as error:
        parser.error(str(error))
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
