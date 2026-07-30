#!/usr/bin/env python3
"""Run the resource-bounded Hall micro tier and test large-run guards."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
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
            and result["metrics"]["resolved_modes"] == 9
            and result["metrics"]["reverse_diagnostics_available"]
            and result["metrics"]["reverse_diagnostics_complete"]
            and result["metrics"]["reverse_diagnostics_start_step"] == 0
            and 0.0 <= result["metrics"][
                "reverse_demand_step_fraction"
            ] <= 1.0,
            f"Hall pilot analysis failed: {analyzed.stderr}",
        )
        require(
            result["metrics"][
                "maximum_retained_negative_debt_macroparticles"
            ] == 0,
            "timestep-local Hall pilot retained negative controller debt",
        )

        saturated_output = work / "saturated_output"
        shutil.copytree(output, saturated_output)
        current_path = saturated_output / "current_source.csv"
        with current_path.open(newline="", encoding="utf-8") as stream:
            current_rows = list(csv.DictReader(stream))
            current_fields = list(current_rows[0])
        macro_weight = 5e16 * 0.025 * 0.0128 / 2048
        macro_charge = 1.602176634e-19 * macro_weight
        final_current = current_rows[-1]
        if "control_mode" in current_fields:
            current_fields.remove("control_mode")
            for row in current_rows:
                row.pop("control_mode")
        emitted_charge = float(
            final_current["cumulative_emitted_charge"]
        )
        final_current["control_macro_remainder"] = "-2"
        processed_charge = emitted_charge + 2 * macro_charge
        final_current["cumulative_processed_monitored_charge"] = str(
            processed_charge
        )
        if "cumulative_monitored_negative_charge" in final_current:
            final_current[
                "cumulative_monitored_negative_charge"
            ] = str(min(0.0, processed_charge))
            final_current[
                "cumulative_monitored_positive_charge"
            ] = str(max(0.0, processed_charge))
        final_current["charge_balance_residual"] = str(2 * macro_charge)
        if "raw_charge_balance_residual" in final_current:
            final_current["raw_charge_balance_residual"] = str(2 * macro_charge)
        if "unserved_reverse_charge" in final_current:
            final_current["unserved_reverse_charge"] = "0"
        with current_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=current_fields)
            writer.writeheader()
            writer.writerows(current_rows)
        saturated_report = work / "saturated_analysis.json"
        saturated = run(
            [
                sys.executable,
                str(ANALYZE),
                str(saturated_output),
                str(CASE),
                "--tier",
                "micro",
                "--report",
                str(saturated_report),
            ]
        )
        saturated_result = json.loads(
            saturated_report.read_text(encoding="utf-8")
        )
        require(
            saturated.returncode == 0
            and saturated_result["passed"]
            and saturated_result["metrics"][
                "controller_saturation_samples"
            ] == 1
            and saturated_result["metrics"][
                "maximum_controller_debt_macroparticles"
            ] == 2,
            "Hall pilot analyzer rejected consistent one-way actuator debt",
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
