#!/usr/bin/env python3
"""Conservative regression for the bounded Turner horizon workflow."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from extend_turner_horizon import HorizonError, analyze_cycle
from run_turner_startup import StartupError


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "extend_turner_horizon.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_HORIZON"


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
    electrons = 50000 - (step - 400) // 4
    ions = 60000 - (step - 400) // 10
    return {
        "step": step,
        "time": f"{step / (400.0 * 13.56e6):.17g}",
        "kinetic_energy": f"{1.0 + step * 1e-3:.17g}",
        "field_energy": f"{2.0 + step * 1e-3:.17g}",
        "total_energy": f"{3.0 + step * 2e-3:.17g}",
        "charge_l1": f"{step * 1e-9:.17g}",
        "live_particles": electrons + ions,
        "phi_left": "0",
        "phi_right": f"{450 * math.sin(2 * math.pi * step / 400):.17g}",
        "live_particles_electrons": electrons,
        "live_particles_ions": ions,
    }


def collision(step: int, interval_zero: bool = False) -> dict:
    return {
        "step": step,
        "time": f"{step / (400.0 * 13.56e6):.17g}",
        "candidates": 0 if interval_zero else step,
        "collisions_electron_mcc.ionization":
            0 if interval_zero else step // 20,
        "cumulative_candidates": 1000 + step,
        "cumulative_collisions_electron_mcc.ionization": 100 + step // 20,
    }


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_horizon_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        output = work / "cycle-2-output"
        scalar_fields = list(scalar(400))
        collision_fields = list(collision(400))
        scalar_rows = [scalar(step) for step in range(400, 801, 20)]
        collision_rows = [
            collision(400, interval_zero=True), collision(800)
        ]
        write_csv(output / "scalars.csv", scalar_fields, scalar_rows)
        write_csv(
            output / "collisions.csv", collision_fields, collision_rows
        )
        write_csv(
            output / "fields_800.csv", ["x", "rho", "phi", "E"],
            [
                {
                    "x": index / 128,
                    "rho": 0,
                    "phi": 0,
                    "E": 20 if index in (0, 128) else 2,
                }
                for index in range(129)
            ],
        )
        metrics, final_scalar, final_collision = analyze_cycle(
            2, scalar(400), collision(400), output
        )
        require(
            metrics["cycle"] == 2
            and metrics["population"]["ionization_pairs_created"] == 20
            and metrics["population"][
                "inferred_electron_electrode_losses"
            ] == 120
            and metrics["population"]["inferred_ion_electrode_losses"] == 60
            and metrics["phase_zero_field"][
                "boundary_to_bulk_rms_ratio"
            ] == 10
            and int(final_scalar["step"]) == 800
            and int(final_collision["cumulative_candidates"]) == 1800,
            "horizon cycle analyzer produced an incomplete report",
        )

        changed = scalar(400)
        changed["live_particles_ions"] -= 1
        try:
            analyze_cycle(2, changed, collision(400), output)
        except (HorizonError, StartupError) as error:
            require(
                "checkpoint scalar discontinuity" in str(error),
                "horizon analyzer reported the wrong continuity error",
            )
        else:
            raise RuntimeError("horizon analyzer accepted a discontinuity")

        common = [
            sys.executable, str(RUNNER), str(work / "missing-cli"),
            str(work / "missing.case"), str(work / "missing-normalized"),
            "--prior-work-dir", str(work / "prior"),
            "--prior-report", str(work / "prior.json"),
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
            "horizon runner bypassed its acknowledgement",
        )
        too_many = subprocess.run(
            [
                *common, "--additional-cycles", "4",
                "--acknowledge-cost", ACKNOWLEDGEMENT,
            ],
            cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
        )
        require(
            too_many.returncode == 2
            and "built-in limit" in too_many.stderr,
            "horizon runner bypassed its cycle limit",
        )
    print("Turner bounded horizon regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
