#!/usr/bin/env python3
"""Conservative regression for fresh-seed Turner particle scaling."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_prepare_turner_case import create_fixture


ROOT = Path(__file__).resolve().parents[1]
PREPARER = ROOT / "scripts" / "prepare_turner_particle_scaling.py"
EXECUTABLE = ROOT / "build" / "aurorapic_cli"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LONG_TURNER_PARTICLE_SCALING_CAMPAIGN"
SEEDS = "4043517607,2304641002,3809002602"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=ROOT, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          timeout=60)


def main() -> int:
    require(EXECUTABLE.is_file(), "usage: build aurorapic_cli first")
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="aurorapic_turner_particle_scaling_",
            dir=ROOT / "tmp") as temporary:
        work = Path(temporary)
        normalized, case_path, _ = create_fixture(work)
        destination = work / "campaign"
        command = [
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--executable", str(EXECUTABLE), "--output-dir", str(destination),
            "--seeds", SEEDS, "--burn-in-rf-cycles", "4096",
            "--measurement-blocks", "64",
        ]
        rejected = run(command)
        require(rejected.returncode == 2
                and ACKNOWLEDGEMENT in rejected.stderr
                and not destination.exists(),
                "particle scaling bypassed its cost acknowledgement")
        completed = run(command + ["--acknowledge-cost", ACKNOWLEDGEMENT])
        require(completed.returncode == 0, completed.stdout + completed.stderr)
        manifest = json.loads((destination / "particle_scaling.json")
                              .read_text(encoding="utf-8"))
        require(manifest["turner_particle_scaling_preparation_version"] == 1
                and manifest["seeds"] == [4043517607, 2304641002, 3809002602]
                and len(manifest["arms"]) == 6
                and manifest["execution_policy"]
                    ["maximum_concurrent_solver_processes"] == 1
                and manifest["launched"] is False,
                "particle-scaling manifest contract is incomplete")
        arms = {(arm["seed"], arm["level"]): arm
                for arm in manifest["arms"]}
        for seed in manifest["seeds"]:
            one = arms[(seed, "particles_1x")]
            two = arms[(seed, "particles_2x")]
            require(one["burn_in_steps"] == 1638400
                    and one["measurement_blocks"] == 64
                    and two["particles_per_species"]
                        == 2 * one["particles_per_species"]
                    and abs(two["macro_weight"] * 2 - one["macro_weight"]) < 1e-6
                    and abs(two["represented_initial_particles_per_species"]
                            - one["represented_initial_particles_per_species"]) < 1.0,
                    f"particle scaling arm for seed {seed} changes density")
            one_text = (destination / one["runtime_config"]).read_text()
            two_text = (destination / two["runtime_config"]).read_text()
            require("steps = 1638400" in one_text
                    and "spatial_average_start_step = 1625601" in one_text
                    and "particles = 65536" in one_text
                    and "particles = 131072" in two_text
                    and "opportunity_sampling = single_bernoulli" in one_text
                    and "opportunity_sampling = single_bernoulli" in two_text,
                    f"particle scaling deck for seed {seed} drifted")
        overwrite = run(command + ["--acknowledge-cost", ACKNOWLEDGEMENT])
        require(overwrite.returncode == 2
                and "refusing to overwrite" in overwrite.stderr,
                "particle-scaling preparer overwrote an existing campaign")
        print("Turner particle-scaling preparation regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
