#!/usr/bin/env python3
"""Conservative regression for bounded Hall runtime qualification."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
QUALIFY = ROOT / "scripts" / "qualify_hall_runtime.py"
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_HALL_PROBE"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["AURORA_OPENMP_THREADS"] = "1"
    return subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError(
            "usage: test_hall_runtime_qualification.py <aurorapic_cli>"
        )
    executable = Path(sys.argv[1]).resolve()
    require(executable.is_file(), f"missing AuroraPIC CLI: {executable}")
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_qualification_test_"
    ) as temporary:
        work = Path(temporary)
        unacknowledged = run([
            sys.executable,
            str(QUALIFY),
            str(executable),
            str(CASE),
            "--tier",
            "workstation",
            "--steps",
            "1",
            "--report",
            str(work / "unacknowledged.json"),
        ])
        require(
            unacknowledged.returncode == 2
            and ACKNOWLEDGEMENT in unacknowledged.stderr,
            "workstation runtime probe bypassed its acknowledgement",
        )
        over_limit = run([
            sys.executable,
            str(QUALIFY),
            str(executable),
            str(CASE),
            "--steps",
            "4",
            "--max-initial-updates",
            "1000",
            "--report",
            str(work / "over-limit.json"),
        ])
        require(
            over_limit.returncode == 2
            and "max-initial-updates" in over_limit.stderr,
            "runtime probe bypassed its work cap",
        )
        report = work / "qualification.json"
        completed = run([
            sys.executable,
            str(QUALIFY),
            str(executable),
            str(CASE),
            "--steps",
            "4",
            "--max-initial-updates",
            "20000",
            "--timeout-seconds",
            "20",
            "--report",
            str(report),
        ])
        require(
            completed.returncode == 0,
            "bounded runtime qualification failed: "
            + completed.stdout + completed.stderr,
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        probe = result["probe"]
        production = result["projections"]["production"]
        require(
            result["physics_claim"] == "none"
            and not result["production_launch_authorized"]
            and probe["runtime_threads"] == 1
            and probe["steps"] == 4
            and probe["initial_particle_updates"] == 16384
            and probe["initial_particle_updates_per_second"] > 0
            and production["initial_particle_updates"] == 76800000000000
            and production["initial_population_only_wall_days"] > 0
            and len(result["provenance"]["executable_sha256"]) == 64,
            "runtime qualification report is incomplete",
        )
        repeated = run([
            sys.executable,
            str(QUALIFY),
            str(executable),
            str(CASE),
            "--steps",
            "4",
            "--max-initial-updates",
            "20000",
            "--report",
            str(report),
        ])
        require(
            repeated.returncode == 2
            and "overwrite" in repeated.stderr,
            "runtime qualification overwrote an existing report",
        )
    print("Hall bounded runtime qualification validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
