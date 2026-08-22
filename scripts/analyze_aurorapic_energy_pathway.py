#!/usr/bin/env python3
"""Audit heating, collisional cooling, and inferred transport by RF phase."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_aurorapic_ionizing_tail import (
    BLOCK_REPORT_SHA256, LENGTH_M, NODES, PHASES, RULE_SHA256,
    region_phase_sum, selected_nodes,
)
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, read_matrix,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


IONIZING_TAIL_AUDIT_SHA256 = (
    "de14b4599a141b99692b6585cb4ad2761903c6e7d018dc5d4d5a3a80c2657d30")
RF_FREQUENCY_HZ = 13.56e6
CRITICAL_X_MIN_M = 0.005
CRITICAL_X_MAX_M = 0.015
COLLISION_POWER_SHA256 = (
    "de6235d79ce1f0dc91e3f721e4b16329902e2072b75a594920f96c75647b8570")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def periodic_derivative(values: list[float]) -> list[float]:
    if len(values) != PHASES * NODES:
        raise ValueError("energy-density field has the wrong shape")
    phase_dt = 1.0 / RF_FREQUENCY_HZ / PHASES
    return [
        (values[((phase + 1) % PHASES) * NODES + node] -
         values[((phase - 1) % PHASES) * NODES + node]) /
        (2.0 * phase_dt)
        for phase in range(PHASES) for node in range(NODES)]


def analyze(block_root: Path, block_report_path: Path,
            ionizing_tail_audit_path: Path,
            reference: Path) -> dict[str, object]:
    if sha256(block_report_path) != BLOCK_REPORT_SHA256:
        raise ValueError("ionizing-tail block report differs")
    if sha256(ionizing_tail_audit_path) != IONIZING_TAIL_AUDIT_SHA256:
        raise ValueError("ionizing-tail audit differs")
    report = json.loads(block_report_path.read_text(encoding="utf-8"))
    tail_audit = json.loads(
        ionizing_tail_audit_path.read_text(encoding="utf-8"))
    if (report.get("rule_sha256") != RULE_SHA256 or
            tail_audit.get("execution_outcome", {}).get(
                "critical_interior_scope_sufficient") is not True):
        raise ValueError("critical ionizing-tail scope is not usable")
    output = block_root / "output"
    required = ("spatial_phase_moments.csv", "spatial_phase_fields.csv")
    for name in required:
        expected = report["output_hashes"].get(name)
        if expected is None or sha256(output / name) != expected:
            raise ValueError(f"ionizing-tail output differs: {name}")
    collision_power_path = output / "spatial_phase_collision_power.csv"
    if sha256(collision_power_path) != COLLISION_POWER_SHA256:
        raise ValueError("ionizing-tail collision-power output differs")

    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    fields = rows(output / "spatial_phase_fields.csv")
    if len(moments) != PHASES * NODES or len(fields) != len(moments):
        raise ValueError("candidate phase-space shape differs")
    candidate_density = [float(row["number_density_mean_m-3"])
                         for row in moments]
    candidate_energy_density = [
        density * float(row["mean_kinetic_energy_eV"]) * ELEMENTARY_CHARGE_C
        for density, row in zip(candidate_density, moments)]
    candidate_electric_power = [
        -ELEMENTARY_CHARGE_C * density * float(moment["mean_velocity_x"]) *
        float(field["electric_field_mean_V_m"])
        for density, moment, field in zip(candidate_density, moments, fields)]

    collision_channels: dict[str, list[float]] = {}
    collision_seen: dict[str, set[int]] = {}
    with collision_power_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            channel = row["channel"]
            if not channel.startswith("electron_mcc."):
                continue
            values = collision_channels.setdefault(
                channel, [0.0] * (PHASES * NODES))
            seen = collision_seen.setdefault(channel, set())
            phase, node = int(row["phase_bin"]), int(row["node"])
            index = phase * NODES + node
            if (phase >= PHASES or node >= NODES or index in seen):
                raise ValueError("candidate collision-power rows differ")
            seen.add(index)
            values[index] = float(row["mean_power_density_W_m-3"])
    expected_channels = {
        "electron_mcc.elastic", "electron_mcc.excitation",
        "electron_mcc.ionization"}
    if set(collision_channels) != expected_channels:
        raise ValueError("candidate electron collision channels differ")
    if any(len(collision_seen[channel]) != PHASES * NODES
           for channel in expected_channels):
        raise ValueError("candidate collision-power grid is incomplete")
    candidate_collision_power = [
        math.fsum(values[index] for values in collision_channels.values())
        for index in range(PHASES * NODES)]
    candidate_storage = periodic_derivative(candidate_energy_density)
    candidate_inferred_transport_divergence = [
        electric + collision - storage
        for electric, collision, storage in zip(
            candidate_electric_power, candidate_collision_power,
            candidate_storage)]

    reference_values = {}
    for name in ("electron_density", "electron_mean_energy",
                 "electron_ohmic_power_density"):
        filename, expected = REFERENCE_FILES[name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_values[name] = flatten_phase_major(read_matrix(path))
    reference_density = reference_values["electron_density"]
    reference_energy_density = [
        density * energy * ELEMENTARY_CHARGE_C for density, energy in zip(
            reference_density, reference_values["electron_mean_energy"])]
    reference_storage = periodic_derivative(reference_energy_density)
    reference_electric_power = reference_values[
        "electron_ohmic_power_density"]

    nodes = selected_nodes(CRITICAL_X_MIN_M, CRITICAL_X_MAX_M)
    dx = LENGTH_M / (NODES - 1)

    def integrate(values: list[float], first: int, past_last: int) -> float:
        return (region_phase_sum(values, nodes, first, past_last) * dx /
                (past_last - first))

    tail_octants = tail_audit["critical_phase_octants"]

    def summarize(first: int, past_last: int,
                  ionization_ratio: float) -> dict[str, object]:
        candidate_column = integrate(candidate_density, first, past_last)
        reference_column = integrate(reference_density, first, past_last)
        candidate_heating = integrate(
            candidate_electric_power, first, past_last)
        reference_heating = integrate(
            reference_electric_power, first, past_last)
        power_ratio = ((candidate_heating / candidate_column) /
                       (reference_heating / reference_column))
        channel_power = {
            channel: integrate(values, first, past_last)
            for channel, values in sorted(collision_channels.items())}
        collision_power = math.fsum(channel_power.values())
        storage = integrate(candidate_storage, first, past_last)
        transport = integrate(
            candidate_inferred_transport_divergence, first, past_last)
        return {
            "candidate_electron_column_density_m-2": candidate_column,
            "published_edupic_electron_column_density_m-2": reference_column,
            "candidate_to_reference_density_ratio":
                candidate_column / reference_column,
            "candidate_electric_power_W_m-2": candidate_heating,
            "published_edupic_electric_power_W_m-2": reference_heating,
            "candidate_to_reference_electric_power_ratio":
                candidate_heating / reference_heating,
            "candidate_to_reference_electric_power_per_electron_ratio":
                power_ratio,
            "candidate_to_reference_effective_ionization_frequency_ratio":
                ionization_ratio,
            "ionization_ratio_divided_by_power_per_electron_ratio":
                ionization_ratio / power_ratio,
            "candidate_collision_kinetic_power_W_m-2": collision_power,
            "candidate_collision_channel_power_W_m-2": channel_power,
            "candidate_kinetic_energy_storage_rate_W_m-2": storage,
            "candidate_inferred_outward_energy_flux_divergence_W_m-2":
                transport,
            "candidate_discrete_balance_residual_W_m-2":
                candidate_heating + collision_power - storage - transport,
            "published_edupic_kinetic_energy_storage_rate_W_m-2":
                integrate(reference_storage, first, past_last),
        }

    critical_ionization = float(tail_audit[
        "critical_interior_phase_0p125_to_0p5"][
            "candidate_measured_to_published_edupic_ratio"])
    critical = summarize(25, 100, critical_ionization)
    octants = []
    for index, first in enumerate(range(0, PHASES, 25)):
        ionization_ratio = float(tail_octants[index][
            "candidate_measured_to_published_edupic_ratio"])
        value = summarize(first, first + 25, ionization_ratio)
        value.update({"lower_phase_fraction": first / PHASES,
                      "upper_phase_fraction": (first + 25) / PHASES})
        octants.append(value)
    exceptional = min(
        octants,
        key=lambda value: float(value[
            "ionization_ratio_divided_by_power_per_electron_ratio"]))
    return {
        "schema_version": 1,
        "case_id": report["case_id"],
        "scope": "descriptive_interior_electron_energy_pathway_audit",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "block_report_sha256": BLOCK_REPORT_SHA256,
            "ionizing_tail_audit_sha256": IONIZING_TAIL_AUDIT_SHA256,
            "candidate_output_sha256": {
                **{name: report["output_hashes"][name] for name in required},
                "spatial_phase_collision_power.csv": COLLISION_POWER_SHA256,
            },
            "reference_sha256": {
                REFERENCE_FILES[name][0]: REFERENCE_FILES[name][1]
                for name in ("electron_density", "electron_mean_energy",
                             "electron_ohmic_power_density")},
        },
        "comparison_contract": {
            "spatial_fraction": [0.2, 0.6],
            "critical_phase_fraction": [0.125, 0.5],
            "phase_bins": PHASES,
            "spatial_nodes": NODES,
            "candidate_collision_power_sign":
                "positive_adds_tracked_plasma_kinetic_energy",
            "transport_definition":
                "JE_plus_collision_power_minus_kinetic_storage_rate",
            "transport_approximation_boundary": (
                "The ionization channel includes newborn ion kinetic energy; "
                "the inferred term is electron-dominated but is not an exact "
                "species-separated energy flux."),
            "acceptance_thresholds_declared": False,
        },
        "critical_interior_phase_0p125_to_0p5": critical,
        "critical_phase_octants": octants,
        "most_energy_selective_octant": exceptional,
        "assessment": {
            "finding": (
                "Most phase-resolved ionization variation follows electric "
                "power per electron. Phase 0.375-0.5 has the largest extra "
                "ionization deficit relative to its power ratio, isolating "
                "an energy-selective tail, transport, or cooling effect."),
            "next_discriminator": (
                "Measure energy-resolved directional electron crossings at "
                "the 0.2 and 0.6 gap surfaces to separate tail transport from "
                "local energy-space redistribution."),
        },
        "claim_boundary": (
            "This post-hoc energy decomposition uses phase-averaged moments "
            "and an inferred transport residual. It is descriptive and does "
            "not uniquely separate transport from energy-space redistribution."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("block_root", type=Path)
    parser.add_argument("block_report", type=Path)
    parser.add_argument("ionizing_tail_audit", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.block_root.resolve(), args.block_report.resolve(),
        args.ionizing_tail_audit.resolve(), args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
