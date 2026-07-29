#!/usr/bin/env python3
"""Bounded regression for checkpoint-chained Hall horizon preparation."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
CAMPAIGN = ROOT / "scripts" / "prepare_hall_campaign.py"
HORIZON = ROOT / "scripts" / "prepare_hall_horizon_stage.py"
WORKSTATION_ACK = "I_UNDERSTAND_THIS_IS_AN_OPT_IN_WORKSTATION_RUN"
HORIZON_ACK = (
    "I_UNDERSTAND_THIS_EXTENDS_A_HALL_RUN_FROM_A_PINNED_CHECKPOINT"
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_prior_output(path: Path, final_total: int = 260_675) -> None:
    path.mkdir(parents=True)
    final_electrons = 128_245
    final_ions = final_total - final_electrons
    (path / "scalars.csv").write_text(
        "step,time,live_particles,live_particles_electrons,"
        "live_particles_ions\n"
        "0,0,256000,128000,128000\n"
        f"5000,2.5e-8,{final_total},{final_electrons},"
        f"{final_ions}\n",
        encoding="utf-8",
    )
    (path / "checkpoint_5000.apc").write_text(
        "AuroraPIC-checkpoint-v8\n"
        "dimension 2\n"
        "units si 1 8.8541878128e-12 1\n"
        "step 5000\n"
        "time 2.5e-8\n",
        encoding="utf-8",
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError(
            "usage: test_hall_horizon_stage.py <aurorapic_cli>"
        )
    cli = Path(sys.argv[1]).resolve()
    require(cli.is_file(), "AuroraPIC CLI is missing")
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_horizon_"
    ) as temporary:
        work = Path(temporary)
        prior_output = work / "prior-output"
        prior_deck = work / "prior.cfg"
        generated = run(
            [
                sys.executable,
                str(CAMPAIGN),
                str(CASE),
                "--tier", "workstation",
                "--output", str(prior_deck),
                "--output-dir", str(prior_output),
                "--seed", "24680",
                "--acknowledge-cost", WORKSTATION_ACK,
            ]
        )
        require(
            generated.returncode == 0,
            f"could not create prior deck: {generated.stderr}",
        )
        write_prior_output(prior_output)
        guarded = run(
            [
                sys.executable,
                str(HORIZON),
                str(CASE),
                "--prior-deck", str(prior_deck),
                "--prior-output", str(prior_output),
                "--output-dir", str(work / "guarded"),
            ]
        )
        require(
            guarded.returncode == 2
            and HORIZON_ACK in guarded.stderr,
            "horizon stage bypassed its acknowledgement",
        )
        stage = work / "stage"
        prepared = run(
            [
                sys.executable,
                str(HORIZON),
                str(CASE),
                "--prior-deck", str(prior_deck),
                "--prior-output", str(prior_output),
                "--output-dir", str(stage),
                "--acknowledge-cost", HORIZON_ACK,
            ]
        )
        require(
            prepared.returncode == 0,
            f"horizon preparation failed: {prepared.stderr}",
        )
        report = json.loads(
            (stage / "horizon.json").read_text(encoding="utf-8")
        )
        require(
            report["prior_step"] == 5000
            and report["target_step"] == 20000
            and report["target_time_s"] == 1e-7
            and report["added_steps"] == 15000
            and report["diagnostic_start_step"] == 16000
            and report["diagnostic_interval"] == 400
            and report["diagnostic_samples"] == 11
            and report["runtime_threads"] == 1
            and report["runtime_backend"] == "serial"
            and report["max_particles_per_species"] == 250000
            and report["estimated_added_particle_updates_lower_bound"]
                == 260675 * 15000
            and len(report["restart_checkpoint_sha256"]) == 64
            and not report["launched"]
            and report["physics_claim"] == "none",
            "horizon stage contract is incomplete",
        )
        validated = run(
            [str(cli), "--validate-only", str(stage / "horizon.cfg")]
        )
        require(
            validated.returncode == 0,
            f"horizon deck is invalid: {validated.stderr}",
        )

        expensive_output = work / "expensive-output"
        expensive_deck = work / "expensive.cfg"
        expensive_generated = run(
            [
                sys.executable,
                str(CAMPAIGN),
                str(CASE),
                "--tier", "workstation",
                "--output", str(expensive_deck),
                "--output-dir", str(expensive_output),
                "--seed", "24680",
                "--acknowledge-cost", WORKSTATION_ACK,
            ]
        )
        require(
            expensive_generated.returncode == 0,
            "could not create budget-failure deck",
        )
        write_prior_output(expensive_output, final_total=400_000)
        over_budget = run(
            [
                sys.executable,
                str(HORIZON),
                str(CASE),
                "--prior-deck", str(expensive_deck),
                "--prior-output", str(expensive_output),
                "--output-dir", str(work / "over-budget"),
                "--acknowledge-cost", HORIZON_ACK,
            ]
        )
        require(
            over_budget.returncode == 2
            and "maximum_added_particle_updates"
                in over_budget.stderr,
            "horizon stage ignored its update budget",
        )
    print("Hall horizon stage preparation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
