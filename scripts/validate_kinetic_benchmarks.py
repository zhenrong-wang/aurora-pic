#!/usr/bin/env python3
"""Run quantitative, resource-bounded AuroraPIC kinetic benchmarks."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tempfile
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]


class BenchmarkFailure(RuntimeError):
    """Raised when a benchmark cannot be evaluated or misses its envelope."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BenchmarkFailure(message)


def bit_reverse(value: int, bits: int) -> int:
    result = 0
    for _ in range(bits):
        result = (result << 1) | (value & 1)
        value >>= 1
    return result


def perturbed_position(unit_coordinate: float, length: float, amplitude: float, wavenumber: float) -> float:
    """Invert x + amplitude*sin(k*x)/k = unit_coordinate*length."""
    target = unit_coordinate * length
    position = target
    for _ in range(8):
        residual = position + amplitude * math.sin(wavenumber * position) / wavenumber - target
        derivative = 1.0 + amplitude * math.cos(wavenumber * position)
        position -= residual / derivative
    return position % length


def write_landau_state(path: Path, particle_count: int, length: float, perturbation: float) -> None:
    require(particle_count > 0 and particle_count & (particle_count - 1) == 0,
            "Landau particle count must be a power of two")
    bits = particle_count.bit_length() - 1
    normal = statistics.NormalDist()
    wavenumber = 2.0 * math.pi / length
    with path.open("w", encoding="utf-8") as output:
        output.write(
            "AuroraPIC-particle-state-v1\n"
            "dimension 1\n"
            "units normalized\n"
            "weighting species_constant\n"
            "velocity_staggering time_centered\n"
            f"particle_count {2 * particle_count}\n"
            "records\n"
        )
        for index in range(particle_count):
            unit_position = (index + 0.5) / particle_count
            position = perturbed_position(
                unit_position, length, perturbation, wavenumber
            )
            velocity_rank = bit_reverse(index, bits)
            velocity_quantile = (velocity_rank + 0.5) / particle_count
            velocity = normal.inv_cdf(velocity_quantile)
            output.write(
                f"particle electrons {position:.17g} 0 0 {velocity:.17g} 0 0\n"
            )
        for index in range(particle_count):
            position = length * (index + 0.5) / particle_count
            output.write(
                f"particle ions {position:.17g} 0 0 0 0 0\n"
            )
        output.write("end\n")


def write_landau_config(
    path: Path,
    state_path: Path,
    output_dir: Path,
    particle_count: int,
    length: float,
) -> None:
    weight = length / particle_count
    path.write_text(
        "\n".join(
            [
                "config_version = 1",
                "units = normalized",
                "nx = 64",
                f"length = {length:.17g}",
                "dt = 0.05",
                "steps = 240",
                "output_interval = 1",
                "boundary = periodic",
                "mode = transient",
                "seed = 314159",
                f"output_dir = {output_dir.as_posix()}",
                f"initial_state_path = {state_path.as_posix()}",
                "initialization_max_relative_charge_imbalance = 1e-12",
                "initialization_max_relative_pair_imbalance = 1e-12",
                "initialization_charge_pairs = electrons:ions",
                "",
                "[species]",
                "name = electrons",
                "charge = -1",
                "mass = 1",
                f"weight = {weight:.17g}",
                f"particles = {particle_count}",
                "thermal_velocity = 0",
                "",
                "[species]",
                "name = ions",
                "charge = 1",
                "mass = 1000000",
                f"weight = {weight:.17g}",
                f"particles = {particle_count}",
                "thermal_velocity = 0",
                "",
            ]
        ),
        encoding="utf-8",
    )


def field_mode_amplitude(path: Path, wavenumber: float) -> float:
    cosine = 0.0
    sine = 0.0
    count = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["x", "rho", "phi", "E"],
                f"unexpected field header in {path}")
        for row in reader:
            position = float(row["x"])
            electric = float(row["E"])
            require(math.isfinite(position) and math.isfinite(electric),
                    f"non-finite field sample in {path}")
            cosine += electric * math.cos(wavenumber * position)
            sine += electric * math.sin(wavenumber * position)
            count += 1
    require(count > 2, f"insufficient field samples in {path}")
    return 2.0 * math.hypot(cosine, sine) / count


def linear_fit(x_values: Sequence[float], y_values: Sequence[float]) -> tuple[float, float]:
    require(len(x_values) == len(y_values) and len(x_values) >= 2,
            "linear fit requires at least two paired samples")
    mean_x = sum(x_values) / len(x_values)
    mean_y = sum(y_values) / len(y_values)
    denominator = sum((value - mean_x) ** 2 for value in x_values)
    require(denominator > 0.0, "linear fit has zero independent-variable variance")
    slope = sum(
        (x_value - mean_x) * (y_value - mean_y)
        for x_value, y_value in zip(x_values, y_values)
    ) / denominator
    intercept = mean_y - slope * mean_x
    return slope, intercept


def analyze_damped_mode(times: Sequence[float], amplitudes: Sequence[float]) -> dict[str, float | int]:
    require(len(times) == len(amplitudes) and len(times) >= 5,
            "damped-mode analysis requires at least five samples")
    candidates = [
        index
        for index in range(1, len(amplitudes) - 1)
        if amplitudes[index] >= amplitudes[index - 1]
        and amplitudes[index] > amplitudes[index + 1]
        and 1.0 <= times[index] <= 10.5
    ]
    peaks: list[int] = []
    for candidate in candidates:
        if peaks and times[candidate] - times[peaks[-1]] < 0.75:
            if amplitudes[candidate] > amplitudes[peaks[-1]]:
                peaks[-1] = candidate
            continue
        peaks.append(candidate)
    require(len(peaks) >= 4, "damped-mode analysis found fewer than four usable peaks")
    peak_times = [times[index] for index in peaks]
    peak_amplitudes = [amplitudes[index] for index in peaks]
    require(all(value > 0.0 for value in peak_amplitudes),
            "damped-mode analysis found a non-positive peak")
    damping_rate, _ = linear_fit(
        peak_times, [math.log(value) for value in peak_amplitudes]
    )
    half_periods = [
        peak_times[index + 1] - peak_times[index]
        for index in range(len(peak_times) - 1)
    ]
    require(all(value > 0.0 for value in half_periods),
            "damped-mode peak times are not strictly increasing")
    angular_frequency = math.pi / (sum(half_periods) / len(half_periods))
    return {
        "peak_count": len(peaks),
        "damping_rate": damping_rate,
        "angular_frequency": angular_frequency,
        "first_fit_peak_time": peak_times[0],
        "last_fit_peak_time": peak_times[-1],
    }


def read_scalar_energy(path: Path) -> tuple[float, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) >= 2, "scalar diagnostics require at least two samples")
    energies = [float(row["total_energy"]) for row in rows]
    require(all(math.isfinite(value) and value > 0.0 for value in energies),
            "scalar diagnostics contain invalid total energy")
    reference = energies[0]
    return reference, max(abs(value - reference) / reference for value in energies)


def run_landau(cli: Path, work: Path) -> dict[str, object]:
    particle_count = 32768
    length = 4.0 * math.pi
    perturbation = 0.01
    expected_initial_amplitude = perturbation / (2.0 * math.pi / length)
    state_path = work / "landau_1d.aps"
    config_path = work / "landau_1d.cfg"
    output_dir = work / "landau_1d_output"
    write_landau_state(
        state_path, particle_count, length, perturbation
    )
    write_landau_config(
        config_path, state_path, output_dir,
        particle_count, length
    )
    subprocess.run([str(cli), str(config_path)], check=True)

    times: list[float] = []
    amplitudes: list[float] = []
    wavenumber = 2.0 * math.pi / length
    for step in range(241):
        field_path = output_dir / f"fields_{step}.csv"
        require(field_path.is_file(), f"missing Landau field snapshot {field_path}")
        times.append(step * 0.05)
        amplitudes.append(
            field_mode_amplitude(field_path, wavenumber)
        )
    mode = analyze_damped_mode(times, amplitudes)
    initial_energy, max_energy_drift = read_scalar_energy(
        output_dir / "scalars.csv"
    )

    checks = {
        "initial_mode_amplitude": {
            "value": amplitudes[0],
            "minimum": 0.019,
            "maximum": 0.021,
        },
        "damping_rate": {
            "value": mode["damping_rate"],
            "reference": -0.1533,
            "minimum": -0.18,
            "maximum": -0.13,
        },
        "angular_frequency": {
            "value": mode["angular_frequency"],
            "reference": 1.4156,
            "minimum": 1.34,
            "maximum": 1.48,
        },
        "max_relative_total_energy_drift": {
            "value": max_energy_drift,
            "minimum": 0.0,
            "maximum": 0.001,
        },
    }
    passed = True
    for check in checks.values():
        value = float(check["value"])
        check["passed"] = (
            float(check["minimum"]) <= value <=
            float(check["maximum"])
        )
        passed = passed and bool(check["passed"])

    report: dict[str, object] = {
        "schema_version": 1,
        "benchmark": "landau_damping_1d",
        "model": "electrostatic_1D1V_Vlasov_Poisson",
        "reference": {
            "description": "standard k=0.5 Maxwellian-electron Landau damping case",
            "damping_rate": -0.1533,
            "angular_frequency": 1.4156,
            "citation": (
                "Crouseilles, Mehrenberger, and Vecil, ESAIM Proc. 32 "
                "(2011), doi:10.1051/proc/2011022; PICLas Landau tutorial"
            ),
        },
        "numerics": {
            "cells": 64,
            "particles_per_species": particle_count,
            "timestep": 0.05,
            "steps": 240,
            "length": length,
            "density_perturbation": perturbation,
            "expected_initial_mode_amplitude": expected_initial_amplitude,
        },
        "fit": mode,
        "initial_total_energy": initial_energy,
        "checks": checks,
        "passed": passed,
    }
    return report


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "cli", nargs="?", type=Path,
        default=ROOT / "build" / "aurorapic_cli",
        help="path to aurorapic_cli",
    )
    parser.add_argument(
        "--report", type=Path,
        help="optional path for the JSON report",
    )
    parser.add_argument(
        "--keep-output", action="store_true",
        help="retain generated state, config, and simulation output",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    cli = args.cli if args.cli.is_absolute() else ROOT / args.cli
    if not cli.is_file():
        print(f"kinetic benchmark failed: missing CLI {cli}", file=sys.stderr)
        return 2
    work = Path(tempfile.mkdtemp(
        prefix="aurorapic_kinetic_validation_", dir=ROOT
    ))
    keep = args.keep_output
    try:
        report = run_landau(cli, work)
        rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(rendered, encoding="utf-8")
        print(rendered, end="")
        require(
            bool(report["passed"]),
            "Landau benchmark missed one or more acceptance gates",
        )
        print("kinetic benchmark validation passed")
        return 0
    except (BenchmarkFailure, subprocess.CalledProcessError) as exc:
        print(f"kinetic benchmark failed: {exc}", file=sys.stderr)
        keep = True
        return 1
    finally:
        if keep:
            print(f"kinetic benchmark output retained at {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
