#!/usr/bin/env python3
"""Run the resource-bounded Hall micro tier and test large-run guards."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
PREPARE = ROOT / "scripts" / "prepare_hall_campaign.py"
ANALYZE = ROOT / "scripts" / "analyze_hall_pilot.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(
    command: list[str],
    *,
    timeout: float = 30.0,
) -> subprocess.CompletedProcess[str]:
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
        timeout=timeout,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError("usage: test_hall_pilot.py <aurorapic_cli>")
    cli = Path(sys.argv[1]).resolve()
    require(cli.is_file(), f"missing AuroraPIC CLI: {cli}")
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_pilot_"
    ) as temporary:
        work = Path(temporary)
        output = work / "output"
        micro_deck = work / "micro.cfg"
        generated = run(
            [
                sys.executable,
                str(PREPARE),
                str(CASE),
                "--tier",
                "micro",
                "--output",
                str(micro_deck),
                "--output-dir",
                str(output),
            ]
        )
        require(
            generated.returncode == 0,
            f"micro deck generation failed: {generated.stderr}",
        )
        simulation = run([str(cli), str(micro_deck)])
        require(
            simulation.returncode == 0
            and "completed steps=200" in simulation.stdout,
            "bounded Hall micro pilot failed: "
            + simulation.stdout + simulation.stderr,
        )
        report = work / "analysis.json"
        analyzed = run(
            [
                sys.executable,
                str(ANALYZE),
                str(output),
                str(CASE),
                "--tier",
                "micro",
                "--report",
                str(report),
            ]
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        require(
            analyzed.returncode == 0
            and result["passed"]
            and result["physics_claim"] == "none"
            and result["metrics"]["steps"] == 200
            and result["metrics"]["resolved_samples"] == 11
            and result["metrics"]["resolved_modes"] == 9,
            f"Hall pilot analysis failed: {analyzed.stderr}",
        )

        workstation_deck = work / "workstation.cfg"
        workstation_generated = run(
            [
                sys.executable,
                str(PREPARE),
                str(CASE),
                "--tier",
                "workstation",
                "--output",
                str(workstation_deck),
                "--output-dir",
                str(work / "workstation_output"),
                "--acknowledge-cost",
                "I_UNDERSTAND_THIS_IS_AN_OPT_IN_WORKSTATION_RUN",
            ]
        )
        require(
            workstation_generated.returncode == 0,
            "workstation deck generation failed",
        )
        validated = run([str(cli), "--validate-only", str(workstation_deck)])
        require(
            validated.returncode == 0
            and "simulation not launched" in validated.stdout,
            "large deck could not be inspected safely",
        )
        blocked = run([str(cli), str(workstation_deck)])
        require(
            blocked.returncode == 1
            and "100,000,000-update CLI limit" in blocked.stderr,
            "CLI did not block an unacknowledged large Hall run",
        )

    print("Hall bounded pilot and large-run guard validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
