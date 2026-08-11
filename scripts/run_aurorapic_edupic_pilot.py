#!/usr/bin/env python3
"""Run a gated AuroraPIC continuation of the pinned eduPIC argon startup."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_AURORAPIC_EDUPIC_PILOT"
STEPS_PER_CYCLE = 4000
HARD_END_CYCLE = 4
HARD_TIMEOUT_SECONDS = 60
HARD_ANALYZER_TIMEOUT_SECONDS = 45
MIN_AVAILABLE_MEMORY_KIB = 512 * 1024
MAX_PARTICLE_GROWTH_FACTOR = 4.0
MAX_TOTAL_PARTICLE_CAP_FRACTION = 0.25
MAX_ABSOLUTE_FIELD_V_M = 1.0e7
MAX_RELATIVE_ENERGY_RESIDUAL = 1.0e-10
MAX_SPATIAL_PHASE_RESIDUAL_J_M2 = 1.0e-15


class PilotError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def set_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    result, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise PilotError(f"base deck does not contain exactly one {key!r}")
    return result


def insert_global(text: str, key: str, value: str) -> str:
    if re.search(rf"(?m)^\s*{re.escape(key)}\s*=", text):
        raise PilotError(f"base deck already contains {key!r}")
    marker = text.find("\n[")
    if marker < 0:
        raise PilotError("base deck has no configuration sections")
    return text[:marker] + f"\n{key} = {value}" + text[marker:]


def global_integer(text: str, key: str) -> int:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([0-9]+)\s*$", text)
    if match is None:
        raise PilotError(f"base deck has no integer {key!r}")
    return int(match.group(1))


def stage_deck(
    base: str, cycle: int, output_dir: Path, restart: Path,
    *, wall_impact_spectrum: bool = False,
) -> str:
    start_step = cycle * STEPS_PER_CYCLE
    end_step = (cycle + 1) * STEPS_PER_CYCLE
    values = {
        "steps": end_step,
        "output_interval": 100,
        "output_dir": output_dir,
        "spatial_average_reset_on_restart": "true",
        "spatial_average_start_step": start_step + 1,
        "spatial_average_end_step": end_step,
        "spatial_average_rf_cycles": 1,
        "spatial_average_phase_bins": 16,
        "checkpoint_interval": STEPS_PER_CYCLE,
        "runtime_backend": "serial",
        "runtime_threads": 1,
    }
    result = base
    for key, value in values.items():
        if key == "spatial_average_reset_on_restart":
            result = insert_global(result, key, str(value))
        else:
            result = set_global(result, key, str(value))
    if wall_impact_spectrum:
        result = insert_global(result, "wall_impact_spectrum", "true")
        result = insert_global(result, "wall_impact_energy_bins", "200")
        result = insert_global(result, "wall_impact_energy_max", "500")
    return insert_global(result, "restart_path", str(restart))


def table(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
    except (OSError, UnicodeError, csv.Error) as error:
        raise PilotError(f"cannot read {path}: {error}") from error
    if not rows:
        raise PilotError(f"empty diagnostic table: {path}")
    return rows


def finite(row: dict[str, str], key: str, context: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as error:
        raise PilotError(f"{context} has invalid {key!r}") from error
    if not math.isfinite(value):
        raise PilotError(f"{context} has non-finite {key!r}")
    return value


def integer(row: dict[str, str], key: str, context: str) -> int:
    value = finite(row, key, context)
    result = int(value)
    if value != result:
        raise PilotError(f"{context} has non-integer {key!r}")
    return result


def available_memory_kib() -> int:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1])
    except (OSError, ValueError):
        pass
    return 2**63 - 1


def child_setup() -> None:
    try:
        os.nice(10)
    except OSError:
        pass
    if hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"):
        try:
            allowed = os.sched_getaffinity(0)
            if allowed:
                os.sched_setaffinity(0, {min(allowed)})
        except OSError:
            pass


def peak_rss_kib(pid: int) -> int:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith(("VmHWM:", "VmRSS:")):
                return int(line.split()[1])
    except (OSError, ValueError):
        pass
    return 0


def run_process(
    command: list[str], stdout: Path, stderr: Path,
    timeout_seconds: float = HARD_TIMEOUT_SECONDS,
) -> dict[str, object]:
    environment = dict(os.environ)
    environment.update({
        "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    })
    start = time.perf_counter()
    peak = 0
    with stdout.open("w", encoding="utf-8") as out, stderr.open(
        "w", encoding="utf-8"
    ) as err:
        process = subprocess.Popen(
            command, stdout=out, stderr=err, text=True, env=environment,
            preexec_fn=child_setup if os.name == "posix" else None,
        )
        while process.poll() is None:
            peak = max(peak, peak_rss_kib(process.pid))
            if time.perf_counter() - start > timeout_seconds:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                raise PilotError(
                    f"run exceeded the {timeout_seconds} second timeout"
                )
            time.sleep(0.02)
        peak = max(peak, peak_rss_kib(process.pid))
    wall = time.perf_counter() - start
    if process.returncode != 0:
        raise PilotError(f"solver failed with status {process.returncode}")
    return {"wall_seconds": wall, "peak_resident_set_kib": peak}


def run_analyzer(command: list[str], context: str) -> None:
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=HARD_ANALYZER_TIMEOUT_SECONDS,
    )
    if completed.returncode != 0:
        raise PilotError(
            f"{context} failed: {(completed.stderr or completed.stdout).strip()}"
        )


def analyze_stage(
    stage: Path, cycle: int, previous: dict[str, int], total_cap: int
) -> dict[str, object]:
    start_step = cycle * STEPS_PER_CYCLE
    end_step = (cycle + 1) * STEPS_PER_CYCLE
    scalars = table(stage / "scalars.csv")
    collisions = table(stage / "collisions.csv")
    boundaries = table(stage / "boundary_losses.csv")
    if integer(scalars[0], "step", "initial scalar") != start_step:
        raise PilotError("restart scalar history begins at the wrong step")
    if integer(scalars[-1], "step", "final scalar") != end_step:
        raise PilotError("restart scalar history ends at the wrong step")
    initial_e = integer(scalars[0], "live_particles_electrons", "initial scalar")
    initial_i = integer(scalars[0], "live_particles_ions", "initial scalar")
    if initial_e != previous["electrons"] or initial_i != previous["ions"]:
        raise PilotError("restart species populations are discontinuous")
    final_e = integer(scalars[-1], "live_particles_electrons", "final scalar")
    final_i = integer(scalars[-1], "live_particles_ions", "final scalar")
    ion_key = "cumulative_collisions_electron_mcc.ionization"
    ionizations = integer(collisions[-1], ion_key, "final collisions") - integer(
        collisions[0], ion_key, "initial collisions"
    )
    collision_totals = {}
    for key in collisions[-1]:
        if key.startswith("cumulative_") and not key.startswith(
            "cumulative_tracked_kinetic_energy_change_"
        ):
            collision_totals[key.removeprefix("cumulative_")] = (
                integer(collisions[-1], key, "final collisions")
                - integer(collisions[0], key, "initial collisions")
            )
    loss = {}
    for species in ("electrons", "ions"):
        loss[species] = sum(
            integer(boundaries[-1], f"absorbed_{side}_count_{species}", "final losses")
            - integer(boundaries[0], f"absorbed_{side}_count_{species}", "initial losses")
            for side in ("left", "right")
        )
    if initial_e + ionizations - loss["electrons"] != final_e:
        raise PilotError("electron population balance does not close")
    if initial_i + ionizations - loss["ions"] != final_i:
        raise PilotError("ion population balance does not close")

    total_initial = initial_e + initial_i
    total_final = final_e + final_i
    growth = total_final / total_initial
    cap_fraction = total_final / total_cap
    field_maximum = 0.0
    for path in stage.glob("fields_*.csv"):
        field_maximum = max(
            field_maximum,
            *(abs(finite(row, "E", "field output")) for row in table(path)),
        )
    energy = json.loads((stage / "energy-budget.json").read_text())
    spatial = json.loads((stage / "spatial-collision.json").read_text())
    relative_energy = abs(float(energy["relative_closure_residual"]))
    spatial_residual = max(
        abs(float(spatial["closure"]["maximum_spatial_global_residual_J_m-2"])),
        abs(float(spatial["closure"]["maximum_phase_spatial_residual_J_m-2"])),
    )
    gates = {
        "positive_electron_population": final_e > 0,
        "particle_growth_factor": growth <= MAX_PARTICLE_GROWTH_FACTOR,
        "particle_cap_fraction": cap_fraction <= MAX_TOTAL_PARTICLE_CAP_FRACTION,
        "absolute_electric_field": field_maximum <= MAX_ABSOLUTE_FIELD_V_M,
        "relative_energy_closure": relative_energy <= MAX_RELATIVE_ENERGY_RESIDUAL,
        "spatial_phase_closure": spatial_residual <= MAX_SPATIAL_PHASE_RESIDUAL_J_M2,
        "energy_analyzer": bool(energy.get("passes")),
        "spatial_analyzer": bool(spatial["closure"].get("passes")),
    }
    return {
        "cycle": cycle + 1,
        "start_step": start_step,
        "end_step": end_step,
        "population": {
            "initial_electrons": initial_e, "initial_ions": initial_i,
            "ionization_pairs": ionizations,
            "electron_wall_losses": loss["electrons"],
            "ion_wall_losses": loss["ions"],
            "final_electrons": final_e, "final_ions": final_i,
            "total_growth_factor": growth,
            "total_particle_cap_fraction": cap_fraction,
        },
        "collision_totals": collision_totals,
        "state": {
            "initial_kinetic_energy_J_m-2": finite(
                scalars[0], "kinetic_energy", "initial scalar"
            ),
            "final_kinetic_energy_J_m-2": finite(
                scalars[-1], "kinetic_energy", "final scalar"
            ),
            "initial_field_energy_J_m-2": finite(
                scalars[0], "field_energy", "initial scalar"
            ),
            "final_field_energy_J_m-2": finite(
                scalars[-1], "field_energy", "final scalar"
            ),
            "initial_total_energy_J_m-2": finite(
                scalars[0], "total_energy", "initial scalar"
            ),
            "final_total_energy_J_m-2": finite(
                scalars[-1], "total_energy", "final scalar"
            ),
            "final_charge_l1_C_m-2": finite(
                scalars[-1], "charge_l1", "final scalar"
            ),
        },
        "maximum_sampled_absolute_electric_field_V_m": field_maximum,
        "relative_energy_closure_residual": relative_energy,
        "maximum_spatial_phase_residual_J_m-2": spatial_residual,
        "gates": gates,
        "passes": all(gates.values()),
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise PilotError("pilot requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    if not 2 <= args.end_cycle <= HARD_END_CYCLE:
        raise PilotError(f"end cycle must be in [2, {HARD_END_CYCLE}]")
    executable = args.executable.resolve()
    base_path = args.base_deck.resolve()
    initial_dir = args.initial_output_dir.resolve()
    work = args.work_dir.resolve()
    if work.exists():
        raise PilotError(f"refusing to overwrite pilot directory: {work}")
    if sha256(executable) != args.expected_binary_sha256.lower():
        raise PilotError("solver SHA-256 does not match the locked value")
    if sha256(base_path) != args.expected_base_deck_sha256.lower():
        raise PilotError("base deck SHA-256 does not match the locked value")
    base = base_path.read_text(encoding="utf-8")
    if global_integer(base, "steps") != STEPS_PER_CYCLE:
        raise PilotError("base deck is not the exact one-cycle diagnostic deck")
    total_cap = 2 * global_integer(base, "max_particles_per_species")
    initial_scalars = table(initial_dir / "scalars.csv")
    previous = {
        "electrons": integer(initial_scalars[-1], "live_particles_electrons", "base scalar"),
        "ions": integer(initial_scalars[-1], "live_particles_ions", "base scalar"),
    }
    checkpoint = initial_dir / f"checkpoint_{STEPS_PER_CYCLE}.apc"
    if sha256(checkpoint) != args.expected_input_checkpoint_sha256.lower():
        raise PilotError("input checkpoint SHA-256 does not match the locked value")
    scripts = Path(__file__).resolve().parent
    work.mkdir(parents=True)
    stages = []
    for cycle in range(1, args.end_cycle):
        available = available_memory_kib()
        if available < MIN_AVAILABLE_MEMORY_KIB:
            raise PilotError(
                f"available memory {available} KiB is below the hard launch floor"
            )
        stage = work / f"cycle-{cycle + 1:04d}"
        output = stage / "output"
        stage.mkdir()
        deck_path = stage / "input.cfg"
        deck = stage_deck(base, cycle, output, checkpoint)
        atomic_text(deck_path, deck)
        stdout = stage / "stdout.txt"
        stderr = stage / "stderr.txt"
        resources = run_process([str(executable), str(deck_path)], stdout, stderr)
        energy_path = output / "energy-budget.json"
        spatial_path = output / "spatial-collision.json"
        run_analyzer([
            sys.executable, str(scripts / "analyze_1d_energy_budget.py"),
            str(output), "--json", str(energy_path),
        ], "energy analyzer")
        run_analyzer([
            sys.executable, str(scripts / "analyze_1d_spatial_collision.py"),
            str(output), "--boundary", "dirichlet", "--json", str(spatial_path),
        ], "spatial collision analyzer")
        result = analyze_stage(output, cycle, previous, total_cap)
        result.update({
            "resources": resources,
            "available_memory_before_launch_kib": available,
            "input_checkpoint_sha256": sha256(checkpoint),
            "output_checkpoint_sha256": sha256(
                output / f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc"
            ),
            "deck_sha256": sha256(deck_path),
            "scalars_sha256": sha256(output / "scalars.csv"),
            "collisions_sha256": sha256(output / "collisions.csv"),
            "energy_analysis_sha256": sha256(energy_path),
            "spatial_analysis_sha256": sha256(spatial_path),
        })
        atomic_json(stage / "stage-report.json", result)
        stages.append(result)
        if not result["passes"]:
            raise PilotError(f"cycle {cycle + 1} failed a continuation gate")
        previous = {
            "electrons": result["population"]["final_electrons"],
            "ions": result["population"]["final_ions"],
        }
        checkpoint = output / f"checkpoint_{(cycle + 1) * STEPS_PER_CYCLE}.apc"
    report = {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_aurorapic_multicycle_startup_pilot",
        "physics_claim": "none_equilibration_trend_only",
        "hard_gates": {
            "end_cycle": HARD_END_CYCLE,
            "timeout_seconds_per_cycle": HARD_TIMEOUT_SECONDS,
            "minimum_available_memory_kib": MIN_AVAILABLE_MEMORY_KIB,
            "maximum_particle_growth_factor_per_cycle": MAX_PARTICLE_GROWTH_FACTOR,
            "maximum_total_particle_cap_fraction": MAX_TOTAL_PARTICLE_CAP_FRACTION,
            "maximum_absolute_field_V_m": MAX_ABSOLUTE_FIELD_V_M,
            "maximum_relative_energy_residual": MAX_RELATIVE_ENERGY_RESIDUAL,
            "maximum_spatial_phase_residual_J_m-2": MAX_SPATIAL_PHASE_RESIDUAL_J_M2,
        },
        "inputs": {
            "executable_sha256": sha256(executable),
            "base_deck_sha256": sha256(base_path),
            "initial_checkpoint_sha256": args.expected_input_checkpoint_sha256.lower(),
        },
        "stages": stages,
        "completed_through_cycle": args.end_cycle,
        "all_gates_passed": True,
        "production_launch_authorized": False,
    }
    atomic_json(work / "pilot-report.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_deck", type=Path)
    parser.add_argument("initial_output_dir", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--end-cycle", type=int, default=4)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-base-deck-sha256", required=True)
    parser.add_argument("--expected-input-checkpoint-sha256", required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        report = execute(parse_args())
    except (PilotError, OSError, json.JSONDecodeError) as error:
        print(f"AuroraPIC eduPIC pilot rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
