#!/usr/bin/env python3
"""Fold native-nearest and Aurora-lower collision kernels over checkpoints."""

from __future__ import annotations

import argparse
from array import array
import hashlib
import json
import math
from pathlib import Path

from audit_edupic_frozen_collision_kernel import (
    ARGON_MASS_KG, DT_E_S, ELECTRON_MASS_KG, ELEMENTARY_CHARGE_C,
    ENERGY_STEP_EV, EXPECTED_HASHES, NEUTRAL_DENSITY_M3,
)


BOLTZMANN_J_K = 1.380649e-23
NEUTRAL_TEMPERATURE_K = 350.0
LOCKED_CHECKPOINTS = {
    "microstate_51949":
        "15fab37306177b7042d4eb16592261fdbb9ed5817ee2ef05ea0304ffdb5e12b6",
    "microstate_63059":
        "f4e7d8c7a3e012903df93ef39184a4f560d99a1664f60226646d11e15e29d233",
}
CHANNELS = {
    "electrons": ("elastic", "excitation", "ionization"),
    "ions": ("isotropic", "backward"),
}
THRESHOLDS_EV = {
    "electrons": (0.0, 11.5, 15.8),
    "ions": (0.0, 0.0),
}


class FoldError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> array:
    values = array("d")
    expected_index = 0
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            energy, cross_section = map(float, stripped.split())
            if abs(energy - expected_index * ENERGY_STEP_EV) > 5e-10:
                raise FoldError(f"table energy grid differs: {path}")
            values.append(cross_section)
            expected_index += 1
    if len(values) != 1_000_000:
        raise FoldError(f"table row count differs: {path}")
    return values


def parse_checkpoint(path: Path, expected_hash: str
                     ) -> dict[str, list[tuple[float, float, float]]]:
    if sha256(path) != expected_hash:
        raise FoldError("checkpoint hash differs")
    with path.open(encoding="utf-8") as stream:
        magic = stream.readline().strip()
        if magic != "AuroraPIC-checkpoint-v21":
            raise FoldError("checkpoint format differs")
        for line in stream:
            fields = line.split()
            if fields and fields[0] == "step":
                if int(fields[1]) != 36000:
                    raise FoldError("checkpoint step differs")
                break
        else:
            raise FoldError("checkpoint step is missing")
        if stream.readline().split()[0] != "time":
            raise FoldError("checkpoint time is missing")
        species_count = stream.readline().split()
        if species_count != ["species_count", "2"]:
            raise FoldError("checkpoint species count differs")
        if stream.readline().split()[0] != "rng":
            raise FoldError("checkpoint RNG state is missing")
        result: dict[str, list[tuple[float, float, float]]] = {}
        for species_id in range(2):
            header = stream.readline().split()
            if (len(header) != 4 or header[0] != "species" or
                    int(header[1]) != species_id):
                raise FoldError("checkpoint species header differs")
            name, records = header[2], int(header[3])
            velocities: list[tuple[float, float, float]] = []
            for _ in range(records):
                fields = stream.readline().split()
                if len(fields) != 6:
                    raise FoldError("checkpoint particle record differs")
                _, _, vy, vz, vx_half = map(float, fields[:5])
                alive = int(fields[5])
                if alive:
                    velocities.append((vx_half, vy, vz))
            result[name] = velocities
    if set(result) != {"electrons", "ions"}:
        raise FoldError("checkpoint species identities differ")
    return result


def lookup_rates(energy_ev: float, relative_speed: float,
                 tables: list[array], thresholds: tuple[float, ...],
                 nearest: bool) -> list[float]:
    coordinate = energy_ev / ENERGY_STEP_EV
    index = int(coordinate + 0.5) if nearest else int(math.floor(coordinate))
    index = min(index, len(tables[0]) - 1)
    return [
        NEUTRAL_DENSITY_M3 * table[index] * relative_speed
        if nearest or energy_ev >= threshold else 0.0
        for table, threshold in zip(tables, thresholds, strict=True)]


def accumulate_event(rates: list[float], timestep: float,
                     channels: list[float]) -> tuple[float, float]:
    total = sum(rates)
    poisson_mean = total * timestep
    probability = 1.0 - math.exp(-poisson_mean)
    if total > 0.0:
        for index, rate in enumerate(rates):
            channels[index] += probability * rate / total
    return probability, poisson_mean


def fold_electrons(velocities: list[tuple[float, float, float]],
                   tables: list[array]) -> dict[str, object]:
    native_channels = [0.0] * len(tables)
    aurora_channels = [0.0] * len(tables)
    native_events = aurora_events = aurora_poisson = 0.0
    for vx, vy, vz in velocities:
        speed = math.sqrt(vx * vx + vy * vy + vz * vz)
        energy_ev = (
            0.5 * ELECTRON_MASS_KG * speed * speed /
            ELEMENTARY_CHARGE_C)
        native_rates = lookup_rates(
            energy_ev, speed, tables, THRESHOLDS_EV["electrons"], True)
        aurora_rates = lookup_rates(
            energy_ev, speed, tables, THRESHOLDS_EV["electrons"], False)
        probability, _ = accumulate_event(
            native_rates, DT_E_S, native_channels)
        native_events += probability
        probability, mean = accumulate_event(
            aurora_rates, DT_E_S, aurora_channels)
        aurora_events += probability
        aurora_poisson += mean
    return summarize(
        velocities, native_events, aurora_events, aurora_poisson,
        native_channels, aurora_channels, CHANNELS["electrons"])


def fold_ions(velocities: list[tuple[float, float, float]],
              tables: list[array]) -> dict[str, object]:
    # Three-point Gauss-Hermite rule for a standard normal: nodes 0,+/-sqrt(3)
    # with weights 2/3,1/6,1/6. The tensor product exactly integrates normal
    # moments through fifth order without stochastic sampling.
    nodes = ((0.0, 2.0 / 3.0),
             (math.sqrt(3.0), 1.0 / 6.0),
             (-math.sqrt(3.0), 1.0 / 6.0))
    neutral_stddev = math.sqrt(
        BOLTZMANN_J_K * NEUTRAL_TEMPERATURE_K / ARGON_MASS_KG)
    timestep = 20.0 * DT_E_S
    reduced_mass = 0.5 * ARGON_MASS_KG
    native_channels = [0.0] * len(tables)
    aurora_channels = [0.0] * len(tables)
    native_events = aurora_events = aurora_poisson = 0.0
    for vx, vy, vz in velocities:
        for nx, wx in nodes:
            for ny, wy in nodes:
                for nz, wz in nodes:
                    weight = wx * wy * wz
                    gx = vx - neutral_stddev * nx
                    gy = vy - neutral_stddev * ny
                    gz = vz - neutral_stddev * nz
                    speed = math.sqrt(gx * gx + gy * gy + gz * gz)
                    energy_ev = (
                        0.5 * reduced_mass * speed * speed /
                        ELEMENTARY_CHARGE_C)
                    native_rates = lookup_rates(
                        energy_ev, speed, tables,
                        THRESHOLDS_EV["ions"], True)
                    aurora_rates = lookup_rates(
                        energy_ev, speed, tables,
                        THRESHOLDS_EV["ions"], False)
                    local_native = [0.0] * len(tables)
                    local_aurora = [0.0] * len(tables)
                    probability, _ = accumulate_event(
                        native_rates, timestep, local_native)
                    native_events += weight * probability
                    probability, mean = accumulate_event(
                        aurora_rates, timestep, local_aurora)
                    aurora_events += weight * probability
                    aurora_poisson += weight * mean
                    for index in range(len(tables)):
                        native_channels[index] += weight * local_native[index]
                        aurora_channels[index] += weight * local_aurora[index]
    return summarize(
        velocities, native_events, aurora_events, aurora_poisson,
        native_channels, aurora_channels, CHANNELS["ions"],
        quadrature_nodes=27)


def summarize(velocities: list[tuple[float, float, float]],
              native_events: float, aurora_events: float,
              aurora_poisson: float, native_channels: list[float],
              aurora_channels: list[float], names: tuple[str, ...],
              quadrature_nodes: int = 1) -> dict[str, object]:
    return {
        "live_particles": len(velocities),
        "neutral_velocity_quadrature_nodes_per_particle": quadrature_nodes,
        "native_nearest_expected_first_events_per_species_step": native_events,
        "aurorapic_lower_expected_first_events_per_species_step": aurora_events,
        "aurorapic_to_native_expected_first_event_ratio":
            aurora_events / native_events,
        "aurorapic_poisson_expected_events_per_species_step": aurora_poisson,
        "poisson_mean_to_native_single_event_ratio":
            aurora_poisson / native_events,
        "channel_expected_first_events": {
            name: {
                "native_nearest": native,
                "aurorapic_lower": aurora,
                "aurorapic_to_native_ratio": aurora / native,
            }
            for name, native, aurora in zip(
                names, native_channels, aurora_channels, strict=True)
        },
    }


def execute(checkpoints: list[Path], package: Path) -> dict[str, object]:
    if len(checkpoints) != 2:
        raise FoldError("exactly two checkpoints are required")
    for name, expected in EXPECTED_HASHES.items():
        if sha256(package / name) != expected:
            raise FoldError(f"gas artifact differs: {name}")
    electron_tables = [read_table(package / name) for name in (
        "electron_elastic.dat", "electron_excitation.dat",
        "electron_ionization.dat")]
    ion_tables = [read_table(package / name) for name in (
        "ion_isotropic.dat", "ion_backward.dat")]
    states: dict[str, object] = {}
    for state_id, checkpoint in zip(
            LOCKED_CHECKPOINTS, checkpoints, strict=True):
        velocities = parse_checkpoint(
            checkpoint, LOCKED_CHECKPOINTS[state_id])
        states[state_id] = {
            "checkpoint_sha256": LOCKED_CHECKPOINTS[state_id],
            "electrons": fold_electrons(
                velocities["electrons"], electron_tables),
            "ions": fold_ions(velocities["ions"], ion_tables),
        }
    electron_ratios = [float(state["electrons"][
        "aurorapic_to_native_expected_first_event_ratio"])
        for state in states.values()]
    ion_ratios = [float(state["ions"][
        "aurorapic_to_native_expected_first_event_ratio"])
        for state in states.values()]
    poisson_ratios = [
        float(state[species]["poisson_mean_to_native_single_event_ratio"])
        for state in states.values() for species in ("electrons", "ions")]
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "actual_frozen_state_collision_kernel_fold",
        "states": states,
        "ensemble_bounds": {
            "electron_first_event_ratio_minimum": min(electron_ratios),
            "electron_first_event_ratio_maximum": max(electron_ratios),
            "ion_first_event_ratio_minimum": min(ion_ratios),
            "ion_first_event_ratio_maximum": max(ion_ratios),
        },
        "assessment": {
            "nearest_vs_lower_changes_expected_electron_events_below_0p1_percent":
                all(abs(value - 1.0) < .001 for value in electron_ratios),
            "nearest_vs_lower_changes_expected_ion_events_below_0p1_percent":
                all(abs(value - 1.0) < .001 for value in ion_ratios),
            "nearest_vs_lower_changes_all_expected_species_events_below_1_percent":
                all(abs(value - 1.0) < .01
                    for value in [*electron_ratios, *ion_ratios]),
            "lookup_plus_poisson_mean_changes_all_expected_species_events_below_1_percent":
                all(abs(value - 1.0) < .01 for value in poisson_ratios),
            "finding": (
                "In both actual frozen states, nearest-versus-lower lookup and "
                "the additional Poisson opportunity mean each change total "
                "electron and ion traffic by less than one percent. Collision-"
                "product kinematics remains a separate statistical question."),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "gas_artifact_sha256": EXPECTED_HASHES,
        },
        "claim_boundary": (
            "This deterministic fold uses frozen AuroraPIC states and first-event "
            "kernels. It does not evolve either code or pair random trajectories."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoints", nargs=2, type=Path)
    parser.add_argument("gas_package", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        result = execute(
            [path.resolve() for path in args.checkpoints],
            args.gas_package.resolve())
    except (FoldError, OSError, ValueError, ZeroDivisionError) as error:
        parser.error(str(error))
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
