#!/usr/bin/env python3
"""Regression for Turner pre-benchmark stationarity screening."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_turner_stationarity.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20,
    )


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_stationarity_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        startup_path = work / "startup.json"
        horizon_path = work / "horizon.json"
        startup = {
            "turner_startup_report_version": 1,
            "case_id": "turner-helium-ccp-2013-case-1",
            "diagnostics": {
                "startup_checks_passed": True,
                "particle_balance": {
                    "initial_electrons": 1000,
                    "initial_ions": 1000,
                    "final_electrons": 990,
                    "final_ions": 995,
                },
                "early_field_structure": {
                    "final_boundary_field_max_v_m": 100.0,
                },
                "energy_and_charge": {
                    "final_total_energy_j": 10.0,
                },
                "collisions": {
                    "collisions_electron_mcc.ionization": 100,
                },
            },
            "provenance": {
                "executable_sha256": "a" * 64,
                "normalization_audit_sha256": "b" * 64,
            },
        }
        write(startup_path, startup)
        cycles = []
        electrons = 990
        ions = 995
        for cycle in range(2, 6):
            electrons -= 1
            ions -= 1
            cycles.append({
                "cycle": cycle,
                "population": {
                    "electron_relative_change": -0.001,
                    "ion_relative_change": -0.001,
                    "final_electrons": electrons,
                    "final_ions": ions,
                    "ionization_pairs_created": 100,
                },
                "phase_zero_field": {
                    "boundary_field_max_v_m": 98.0 + cycle,
                },
                "energy_and_charge": {
                    "final_total_energy_j": 9.8 + 0.1 * cycle,
                },
            })
        horizon = {
            "turner_horizon_report_version": 1,
            "case_id": "turner-helium-ccp-2013-case-1",
            "cycles": cycles,
            "provenance": {
                "prior_report_sha256": sha256(startup_path),
                "executable_sha256": "a" * 64,
                "normalization_audit_sha256": "b" * 64,
            },
        }
        write(horizon_path, horizon)
        passed_path = work / "passed.json"
        passed = run([
            sys.executable, str(ANALYZER),
            str(startup_path), str(horizon_path),
            "--output", str(passed_path),
        ])
        require(
            passed.returncode == 0,
            "stationarity pass fixture failed: "
            + passed.stdout + passed.stderr,
        )
        passed_report = json.loads(passed_path.read_text(encoding="utf-8"))
        require(
            passed_report["stationarity_screen_passed"]
            and passed_report["window"]["start_cycle"] == 2
            and passed_report["window"]["end_cycle"] == 5
            and all(
                gate["passed"]
                for gate in passed_report["gates"].values()
            ),
            "stationarity pass report is incomplete",
        )

        failed_path = work / "failed.json"
        failed = run([
            sys.executable, str(ANALYZER),
            str(startup_path), str(horizon_path),
            "--max-population-change", "0.0005",
            "--output", str(failed_path),
        ])
        failed_report = json.loads(failed_path.read_text(encoding="utf-8"))
        require(
            failed.returncode == 0
            and not failed_report["stationarity_screen_passed"]
            and not failed_report["gates"]["population_change"]["passed"],
            "stationarity failure was not reported as data",
        )

        horizon["provenance"]["prior_report_sha256"] = "0" * 64
        write(horizon_path, horizon)
        rejected = run([
            sys.executable, str(ANALYZER),
            str(startup_path), str(horizon_path),
            "--output", str(work / "rejected.json"),
        ])
        require(
            rejected.returncode == 2
            and "hash-chain" in rejected.stderr,
            "stationarity analyzer accepted a broken report chain",
        )
    print("Turner stationarity-screen regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
