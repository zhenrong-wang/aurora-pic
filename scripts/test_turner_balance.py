#!/usr/bin/env python3
"""Bounded regression for the Turner source/wall balance analyzer."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_turner_balance.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_balance_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        scalars = work / "scalars.csv"
        collisions = work / "collisions.csv"
        boundaries = work / "boundary_losses.csv"
        power = work / "power_transfer.csv"
        write_csv(
            scalars,
            [
                "step", "time", "live_particles_electrons",
                "live_particles_ions",
            ],
            [[100, 1.0, 20, 30], [500, 2.0, 22, 31]],
        )
        write_csv(
            collisions,
            [
                "step", "time",
                "cumulative_collisions_electron_mcc.ionization",
            ],
            [[100, 1.0, 10], [500, 2.0, 15]],
        )
        header = ["step", "time", "counter_origin_step"]
        for species in ("electrons", "ions"):
            for side in ("left", "right"):
                header.extend([
                    f"absorbed_{side}_count_{species}",
                    f"absorbed_{side}_charge_{species}_C_m-2",
                    f"absorbed_{side}_kinetic_energy_{species}_J_m-2",
                ])
        write_csv(
            boundaries,
            header,
            [
                [100, 1.0, 100] + [0, 0.0, 0.0] * 4,
                [
                    500, 2.0, 100,
                    1, -2.0, 3.0,
                    2, -4.0, 5.0,
                    2, 4.0, 7.0,
                    2, 4.0, 9.0,
                ],
            ],
        )
        write_csv(
            power,
            [
                "step", "time", "counter_origin_step",
                "electric_work_electrons_J_m-2",
                "electric_work_ions_J_m-2",
            ],
            [
                [100, 1.0, 100, 10.0, 20.0],
                [500, 2.0, 100, 44.3, 110.6],
            ],
        )
        report = work / "report.json"
        completed = subprocess.run([
            sys.executable, str(ANALYZER),
            "--scalars", str(scalars),
            "--collisions", str(collisions),
            "--boundary-losses", str(boundaries),
            "--power-transfer", str(power),
            "--expected-steps", "400",
            "--reported-ion-current", "4",
            "--output", str(report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(completed.returncode == 0, completed.stderr)
        value = json.loads(report.read_text(encoding="utf-8"))
        require(
            value["balance_exact"]
            and value["species_balance"]["electrons"][
                "balance_residual_macro_particles"
            ] == 0
            and value["species_balance"]["ions"][
                "balance_residual_macro_particles"
            ] == 0
            and value["ion_current_context"][
                "two_electrode_mean_magnitude_A_m2"
            ] == 4.0
            and abs(value["volume_electrical_power_context"][
                "species"
            ]["electrons"]["mean_electrical_power_W_m2"] - 34.3) < 1e-12
            and abs(value["volume_electrical_power_context"][
                "species"
            ]["ions"]["mean_electrical_power_W_m2"] - 90.6) < 1e-12
            and value["physics_claim"].startswith("none_"),
            "Turner balance report is incorrect",
        )

        sensitivity_report = work / "sensitivity-report.json"
        sensitivity = subprocess.run([
            sys.executable, str(ANALYZER),
            "--scalars", str(scalars),
            "--collisions", str(collisions),
            "--boundary-losses", str(boundaries),
            "--power-transfer", str(power),
            "--expected-steps", "400",
            "--window-start-step", "100",
            "--window-end-step", "500",
            "--scope", "published_duration_numerical_sensitivity_window",
            "--reported-ion-current", "4",
            "--output", str(sensitivity_report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(sensitivity.returncode == 0, sensitivity.stderr)
        sensitivity_value = json.loads(
            sensitivity_report.read_text(encoding="utf-8")
        )
        require(
            sensitivity_value["scope"]
                == "published_duration_numerical_sensitivity_window"
            and sensitivity_value["window"]["start_step"] == 100
            and sensitivity_value["window"]["end_step"] == 500
            and sensitivity_value["physics_claim"]
                == "none_changed_published_numerical_contract",
            "Turner balance sensitivity window overstated its claim",
        )

        changed = list(csv.reader(boundaries.open(
            newline="", encoding="utf-8"
        )))
        changed[-1][2] = "99"
        with boundaries.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(changed)
        rejected = subprocess.run([
            sys.executable, str(ANALYZER),
            "--scalars", str(scalars),
            "--collisions", str(collisions),
            "--boundary-losses", str(boundaries),
            "--output", str(work / "bad.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            rejected.returncode == 2
            and "do not cover the entire" in rejected.stderr,
            "Turner balance analyzer accepted partial wall coverage",
        )
        changed[-1][2] = "100"
        with boundaries.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(changed)
        changed_power = list(csv.reader(power.open(
            newline="", encoding="utf-8"
        )))
        changed_power[-1][2] = "99"
        with power.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(changed_power)
        rejected = subprocess.run([
            sys.executable, str(ANALYZER),
            "--scalars", str(scalars),
            "--collisions", str(collisions),
            "--boundary-losses", str(boundaries),
            "--power-transfer", str(power),
            "--output", str(work / "bad-power.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            rejected.returncode == 2
            and "power counters do not cover" in rejected.stderr,
            "Turner balance analyzer accepted partial power coverage",
        )

    print("Turner source/wall balance regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
