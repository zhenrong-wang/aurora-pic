#!/usr/bin/env python3
"""Conservative regression for bounded Turner runtime qualification."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from test_prepare_turner_case import create_fixture


ROOT = Path(__file__).resolve().parents[1]
QUALIFIER = ROOT / "scripts" / "qualify_turner_runtime.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_PROBE"


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
            "usage: test_turner_runtime_qualification.py <aurorapic_cli>"
        )
    executable = Path(sys.argv[1]).resolve()
    require(executable.is_file(), f"missing AuroraPIC CLI: {executable}")
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_runtime_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        normalized, case_path, _ = create_fixture(work)
        report = work / "qualification.json"
        common = [
            sys.executable,
            str(QUALIFIER),
            str(executable),
            str(case_path),
            str(normalized),
        ]
        base = [
            *common,
            "--steps",
            "1",
            "--report",
            str(report),
        ]
        unacknowledged = run(base)
        require(
            unacknowledged.returncode == 2
            and ACKNOWLEDGEMENT in unacknowledged.stderr
            and not report.exists(),
            "runtime probe bypassed its explicit acknowledgement",
        )
        over_limit = run([
            *common,
            "--steps",
            "4",
            "--max-initial-updates",
            "1000",
            "--report",
            str(report),
            "--acknowledge-cost",
            ACKNOWLEDGEMENT,
        ])
        require(
            over_limit.returncode == 2
            and "max-initial-updates" in over_limit.stderr,
            "runtime probe bypassed its particle-update cap",
        )
        completed = run([
            *base,
            "--acknowledge-cost",
            ACKNOWLEDGEMENT,
        ])
        require(
            completed.returncode == 0,
            "bounded runtime qualification failed: "
            + completed.stdout + completed.stderr,
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        probe = result["probe"]
        projection = result["projection"]
        require(
            result["physics_claim"] == "none"
            and not result["production_launch_authorized"]
            and probe["runtime_backend"] == "serial"
            and probe["runtime_threads"] == 1
            and probe["steps"] == 1
            and probe["initial_particles"] == 131072
            and probe["initial_particle_updates"] == 131072
            and probe["initial_particle_updates_per_second"] > 0
            and projection["initial_population_particle_updates"]
            == 67108864000
            and projection["initial_population_only_wall_seconds"] > 0
            and len(result["provenance"]["executable_sha256"]) == 64,
            "runtime qualification report is incomplete",
        )
        repeated = run([
            *base,
            "--acknowledge-cost",
            ACKNOWLEDGEMENT,
        ])
        require(
            repeated.returncode == 2 and "overwrite" in repeated.stderr,
            "runtime qualification overwrote an existing report",
        )
    print("Turner bounded runtime qualification regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
