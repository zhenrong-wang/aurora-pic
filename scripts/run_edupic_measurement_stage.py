#!/usr/bin/env python3
"""Run one immutable, bounded eduPIC native-diagnostic measurement stage."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from run_edupic_stage import (
    TIMESTEPS_PER_CYCLE,
    atomic_json,
    checkpoint_state,
    child_setup,
    convergence_rows,
    sha256,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_MEASUREMENT_STAGE"
HARD_CYCLE_LIMIT = 25
HARD_INITIAL_PARTICLE_STEP_LIMIT = 25_000_000_000
HARD_TIMEOUT_SECONDS = 600

N_GRID = 400
N_XT = 200
N_EEPF = 2000
N_IFED = 200
GAP_METERS = 0.025
DE_EEPF_EV = 0.05
DE_IFED_EV = 1.0

XT_OUTPUTS = (
    "pot_xt.dat",
    "efield_xt.dat",
    "ne_xt.dat",
    "ni_xt.dat",
    "je_xt.dat",
    "ji_xt.dat",
    "powere_xt.dat",
    "poweri_xt.dat",
    "meanee_xt.dat",
    "meanei_xt.dat",
    "ioniz_xt.dat",
)

INFO_FIELDS = {
    "# of simulation cycles in this run": "measurement_cycles",
    "Electron density @ center": "electron_density_center_m3",
    "Plasma frequency @ center": "plasma_frequency_center_rad_s",
    "Debye length @ center": "debye_length_center_m",
    "Electron collision frequency": "electron_collision_frequency_s1",
    "Ion collision frequency": "ion_collision_frequency_s1",
    "Ion flux at powered electrode": "ion_flux_powered_m2_s1",
    "Ion flux at grounded electrode": "ion_flux_grounded_m2_s1",
    "Mean ion energy at powered electrode": "mean_ion_energy_powered_ev",
    "Mean ion energy at grounded electrode": "mean_ion_energy_grounded_ev",
    "Electron flux at powered electrode": "electron_flux_powered_m2_s1",
    "Electron flux at grounded electrode": "electron_flux_grounded_m2_s1",
    "Electron power density (average)": "electron_power_density_w_m3",
    "Ion power density (average)": "ion_power_density_w_m3",
    "Total power density(average)": "total_power_density_w_m3",
}


class MeasurementError(RuntimeError):
    pass


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def numeric_rows(path: Path, expected_rows: int,
                 expected_columns: int) -> tuple[list[list[float]], dict]:
    rows: list[list[float]] = []
    minimum = math.inf
    maximum = -math.inf
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise MeasurementError(f"cannot read measurement output {path}: {error}") from error
    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        if len(fields) != expected_columns:
            raise MeasurementError(
                f"{path.name} line {line_number} has {len(fields)} columns; "
                f"expected {expected_columns}")
        try:
            values = [float(field) for field in fields]
        except ValueError as error:
            raise MeasurementError(
                f"{path.name} line {line_number} has a non-numeric value") from error
        if not all(math.isfinite(value) for value in values):
            raise MeasurementError(f"{path.name} contains a non-finite value")
        rows.append(values)
        minimum = min(minimum, *values)
        maximum = max(maximum, *values)
    if len(rows) != expected_rows:
        raise MeasurementError(
            f"{path.name} has {len(rows)} rows; expected {expected_rows}")
    return rows, {
        "bytes": path.stat().st_size,
        "rows": expected_rows,
        "columns": expected_columns,
        "minimum": minimum,
        "maximum": maximum,
        "sha256": sha256(path),
    }


def require_axis(rows: list[list[float]], spacing: float, name: str) -> None:
    for index, row in enumerate(rows):
        expected = (index + 0.5) * spacing
        if not math.isclose(row[0], expected, rel_tol=1e-6, abs_tol=1e-10):
            raise MeasurementError(
                f"{name} energy axis differs at row {index + 1}: "
                f"{row[0]} != {expected}")


def measurement_outputs(stage_dir: Path) -> tuple[dict, dict]:
    outputs: dict[str, dict] = {}

    density, summary = numeric_rows(stage_dir / "density.dat", N_GRID, 3)
    for index, row in enumerate(density):
        expected_x = index * GAP_METERS / (N_GRID - 1)
        # The C reference writes this coordinate with "%8.5f" precision.
        if not math.isclose(row[0], expected_x, rel_tol=0.0, abs_tol=5.1e-6):
            raise MeasurementError(
                f"density.dat coordinate differs at row {index + 1}")
        if row[1] < 0.0 or row[2] < 0.0:
            raise MeasurementError("density.dat contains a negative density")
    summary["electron_density_minimum_m3"] = min(row[1] for row in density)
    summary["electron_density_maximum_m3"] = max(row[1] for row in density)
    summary["ion_density_minimum_m3"] = min(row[2] for row in density)
    summary["ion_density_maximum_m3"] = max(row[2] for row in density)
    outputs["density.dat"] = summary

    eepf, summary = numeric_rows(stage_dir / "eepf.dat", N_EEPF, 2)
    require_axis(eepf, DE_EEPF_EV, "eepf.dat")
    if any(row[1] < 0.0 for row in eepf):
        raise MeasurementError("eepf.dat contains a negative probability")
    eepf_normalization = sum(
        row[1] * math.sqrt(row[0]) * DE_EEPF_EV for row in eepf)
    if not math.isclose(eepf_normalization, 1.0, rel_tol=2e-5, abs_tol=2e-5):
        raise MeasurementError(
            f"eepf.dat normalization differs from one: {eepf_normalization}")
    summary["weighted_normalization"] = eepf_normalization
    outputs["eepf.dat"] = summary

    ifed, summary = numeric_rows(stage_dir / "ifed.dat", N_IFED, 3)
    require_axis(ifed, DE_IFED_EV, "ifed.dat")
    if any(row[1] < 0.0 or row[2] < 0.0 for row in ifed):
        raise MeasurementError("ifed.dat contains a negative probability")
    powered_normalization = sum(row[1] * DE_IFED_EV for row in ifed)
    grounded_normalization = sum(row[2] * DE_IFED_EV for row in ifed)
    for electrode, value in (("powered", powered_normalization),
                             ("grounded", grounded_normalization)):
        if not math.isclose(value, 1.0, rel_tol=2e-4, abs_tol=2e-4):
            raise MeasurementError(
                f"ifed.dat {electrode} normalization differs from one: {value}")
    summary["powered_normalization"] = powered_normalization
    summary["grounded_normalization"] = grounded_normalization
    outputs["ifed.dat"] = summary

    for name in XT_OUTPUTS:
        _, outputs[name] = numeric_rows(stage_dir / name, N_GRID, N_XT)

    info_path = stage_dir / "info.txt"
    try:
        info_text = info_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise MeasurementError(f"cannot read info.txt: {error}") from error
    if "CONDITION(S) VIOLATED" in info_text:
        raise MeasurementError("eduPIC reported a stability or accuracy violation")
    info: dict[str, float] = {}
    number = re.compile(r"^[ \t]*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)")
    for line in info_text.splitlines():
        if "=" not in line:
            continue
        label, value = line.split("=", 1)
        key = INFO_FIELDS.get(label.strip())
        if key is None:
            continue
        match = number.match(value)
        if match is None:
            raise MeasurementError(f"cannot parse info.txt field {label.strip()}")
        parsed = float(match.group(1))
        if not math.isfinite(parsed):
            raise MeasurementError(f"non-finite info.txt field {label.strip()}")
        info[key] = parsed
    missing = sorted(set(INFO_FIELDS.values()) - set(info))
    if missing:
        raise MeasurementError("info.txt is missing fields: " + ", ".join(missing))
    outputs["info.txt"] = {
        "bytes": info_path.stat().st_size,
        "sha256": sha256(info_path),
    }
    return outputs, info


def run(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MeasurementError(
            "measurement stage requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    if args.cycles > HARD_CYCLE_LIMIT:
        raise MeasurementError(f"cycles exceed built-in limit {HARD_CYCLE_LIMIT}")
    if args.timeout_seconds > HARD_TIMEOUT_SECONDS:
        raise MeasurementError(
            f"timeout exceeds built-in limit {HARD_TIMEOUT_SECONDS} seconds")
    if args.max_initial_particle_steps > HARD_INITIAL_PARTICLE_STEP_LIMIT:
        raise MeasurementError("requested work limit exceeds built-in hard limit")

    executable = args.executable.resolve()
    input_dir = args.input_state_dir.resolve()
    stage_dir = args.stage_dir.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise MeasurementError(f"missing executable eduPIC binary: {executable}")
    executable_sha256 = sha256(executable)
    if executable_sha256 != args.expected_binary_sha256.lower():
        raise MeasurementError("eduPIC binary SHA-256 differs from the locked value")
    if stage_dir.exists():
        raise MeasurementError(f"refusing to overwrite stage directory: {stage_dir}")

    input_checkpoint = input_dir / "picdata.bin"
    input_convergence = input_dir / "conv.dat"
    initial = checkpoint_state(input_checkpoint)
    if initial["sha256"] != args.expected_input_sha256.lower():
        raise MeasurementError("input checkpoint SHA-256 differs from the locked value")
    history = convergence_rows(input_convergence)
    if (history[-1]["cycle"] != initial["cycles"] or
            history[-1]["electrons"] != initial["electrons"] or
            history[-1]["ions"] != initial["ions"]):
        raise MeasurementError("checkpoint and convergence history state differ")
    initial_particle_steps = (
        initial["total_particles"] * args.cycles * TIMESTEPS_PER_CYCLE)
    if initial_particle_steps > args.max_initial_particle_steps:
        raise MeasurementError(
            f"stage exceeds initial-particle-step limit: {initial_particle_steps} > "
            f"{args.max_initial_particle_steps}")

    stage_dir.mkdir(parents=True)
    shutil.copy2(input_checkpoint, stage_dir / "picdata.bin")
    shutil.copy2(input_convergence, stage_dir / "conv.dat")
    environment = dict(os.environ)
    environment.update({"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"})
    command = [str(executable), str(args.cycles), "m"]
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command, cwd=stage_dir, env=environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=args.timeout_seconds,
            preexec_fn=child_setup if os.name == "posix" else None)
    except subprocess.TimeoutExpired as error:
        def timeout_text(value: str | bytes | None) -> str:
            if isinstance(value, bytes):
                return value.decode("utf-8", errors="replace")
            return value or ""
        (stage_dir / "stdout.txt").write_text(
            timeout_text(error.stdout), encoding="utf-8")
        (stage_dir / "stderr.txt").write_text(
            timeout_text(error.stderr), encoding="utf-8")
        raise MeasurementError(
            f"eduPIC measurement stage exceeded {args.timeout_seconds} second timeout") \
            from error
    wall_seconds = time.perf_counter() - started
    (stage_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (stage_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise MeasurementError(
            f"eduPIC measurement stage failed with status {completed.returncode}")

    final = checkpoint_state(stage_dir / "picdata.bin")
    final_history = convergence_rows(stage_dir / "conv.dat")
    expected_cycle = initial["cycles"] + args.cycles
    if final["cycles"] != expected_cycle or final_history[-1]["cycle"] != expected_cycle:
        raise MeasurementError("measurement stage did not reach its exact requested cycle")
    new_rows = final_history[len(history):]
    if (len(new_rows) != args.cycles or
            new_rows[0]["cycle"] != initial["cycles"] + 1 or
            new_rows[-1]["electrons"] != final["electrons"] or
            new_rows[-1]["ions"] != final["ions"]):
        raise MeasurementError("measurement convergence coverage is incomplete")

    diagnostic_outputs, info = measurement_outputs(stage_dir)
    if int(info["measurement_cycles"]) != args.cycles:
        raise MeasurementError("info.txt measurement-cycle count differs")
    report = {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_external_native_measurement_stage",
        "physics_claim": "none",
        "claim_boundary": (
            "A bounded native eduPIC measurement stage validates output and "
            "execution contracts but does not establish cross-code agreement."),
        "source_binary": {
            "path": str(executable),
            "sha256": executable_sha256,
            "locked_before_execution": True,
        },
        "resource_policy": {
            "cpu_affinity_count": 1,
            "nice_increment": 10,
            "timeout_seconds": args.timeout_seconds,
            "initial_particle_step_limit": args.max_initial_particle_steps,
            "initial_particle_steps": initial_particle_steps,
        },
        "stage": {
            "measurement_mode": True,
            "requested_cycles": args.cycles,
            "wall_seconds": wall_seconds,
            "cycles_per_wall_second": args.cycles / wall_seconds,
            "start_cycle": initial["cycles"],
            "end_cycle": final["cycles"],
        },
        "initial_state": initial,
        "final_state": final,
        "new_cycle_population": new_rows,
        "outputs": {
            "stdout_sha256": sha256(stage_dir / "stdout.txt"),
            "stderr_sha256": sha256(stage_dir / "stderr.txt"),
            "convergence_sha256": sha256(stage_dir / "conv.dat"),
            "diagnostics": diagnostic_outputs,
        },
        "reported_observables": info,
        "completed": True,
    }
    atomic_json(stage_dir / "stage-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("input_state_dir", type=Path)
    parser.add_argument("stage_dir", type=Path)
    parser.add_argument("--cycles", type=positive_integer, required=True)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--timeout-seconds", type=positive_integer, default=120)
    parser.add_argument("--max-initial-particle-steps", type=positive_integer,
                        default=1_000_000_000)
    parser.add_argument("--acknowledge-cost")
    args = parser.parse_args()
    try:
        report = run(args)
    except (MeasurementError, OSError) as error:
        print(f"eduPIC measurement stage rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
