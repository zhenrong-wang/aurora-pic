#!/usr/bin/env python3
"""Conservative regression for the Turner one-cycle startup workflow."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from run_turner_startup import StartupError, analyze


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_turner_startup.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_ONE_CYCLE_TURNER_STARTUP"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_csv(path: Path, fieldnames: list[str], data: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def scalar(step: int) -> dict:
    time = step / (400.0 * 13.56e6)
    created = step // 40
    electrons = 65536 + created - step
    ions = 65536 + created - step // 10
    return {
        "step": step,
        "time": f"{time:.17g}",
        "kinetic_energy": f"{1.0 + step * 1e-3:.17g}",
        "field_energy": f"{2.0 + step * 2e-3:.17g}",
        "total_energy": f"{3.0 + step * 3e-3:.17g}",
        "charge_l1": f"{step * 1e-9:.17g}",
        "live_particles": electrons + ions,
        "phi_left": "0",
        "phi_right": f"{450.0 * math.sin(2 * math.pi * step / 400):.17g}",
        "live_particles_electrons": electrons,
        "live_particles_ions": ions,
    }


def collision(step: int, interval_zero: bool = False) -> dict:
    candidates = step * 10
    ionizations = step // 40
    return {
        "step": step,
        "time": f"{step / (400.0 * 13.56e6):.17g}",
        "candidates": 0 if interval_zero else candidates,
        "null_collisions": 0 if interval_zero else step * 6,
        "collisions_electron_mcc.ionization":
            0 if interval_zero else ionizations,
        "cumulative_candidates": candidates,
        "cumulative_null_collisions": step * 6,
        "cumulative_collisions_electron_mcc.ionization": ionizations,
    }


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_startup_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        stage1 = work / "stage1"
        stage2 = work / "stage2"
        scalar_fields = list(scalar(0))
        write_csv(
            stage1 / "scalars.csv", scalar_fields,
            [scalar(step) for step in range(0, 201, 20)],
        )
        write_csv(
            stage2 / "scalars.csv", scalar_fields,
            [scalar(step) for step in range(200, 401, 20)],
        )
        collision_fields = list(collision(0))
        write_csv(
            stage1 / "collisions.csv", collision_fields,
            [collision(0), collision(200)],
        )
        write_csv(
            stage2 / "collisions.csv", collision_fields,
            [collision(200, interval_zero=True), collision(400)],
        )
        write_csv(
            stage2 / "fields_400.csv", ["x", "rho", "phi", "E"],
            [
                {
                    "x": index / 128,
                    "rho": 0,
                    "phi": 0,
                    "E": 10 if index in (0, 128) else 1,
                }
                for index in range(129)
            ],
        )
        report = analyze(stage1, stage2)
        balance = report["particle_balance"]
        require(
            report["startup_checks_passed"]
            and report["restart_continuity"]["integer_fields_exact"]
            and balance["ionization_pairs_created"] == 10
            and balance["inferred_electron_electrode_losses"] == 400
            and balance["inferred_ion_electrode_losses"] == 40
            and report["early_field_structure"][
                "boundary_to_bulk_rms_ratio"
            ] == 10,
            "startup analyzer produced an incomplete report",
        )

        changed = [scalar(step) for step in range(200, 401, 20)]
        changed[0]["live_particles_electrons"] -= 1
        write_csv(stage2 / "scalars.csv", scalar_fields, changed)
        try:
            analyze(stage1, stage2)
        except StartupError as error:
            require(
                "checkpoint scalar discontinuity" in str(error),
                "startup analyzer reported the wrong continuity error",
            )
        else:
            raise RuntimeError("startup analyzer accepted a restart discontinuity")

        common = [
            sys.executable, str(RUNNER), str(work / "missing-cli"),
            str(work / "missing.case"), str(work / "missing-normalized"),
            "--work-dir", str(work / "run"),
            "--report", str(work / "report.json"),
        ]
        unacknowledged = subprocess.run(
            common, cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
        )
        require(
            unacknowledged.returncode == 2
            and ACKNOWLEDGEMENT in unacknowledged.stderr,
            "startup runner bypassed its acknowledgement",
        )
        over_limit = subprocess.run(
            [
                *common, "--max-initial-updates", "1000",
                "--acknowledge-cost", ACKNOWLEDGEMENT,
            ],
            cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
        )
        require(
            over_limit.returncode == 2
            and "max-initial-updates" in over_limit.stderr,
            "startup runner bypassed its work limit",
        )
    print("Turner one-cycle startup regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
