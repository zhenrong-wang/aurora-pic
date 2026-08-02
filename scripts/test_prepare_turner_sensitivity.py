#!/usr/bin/env python3
"""Conservative regression for staged Turner sensitivity preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_prepare_turner_case import create_fixture


ROOT = Path(__file__).resolve().parents[1]
ENSEMBLE_PREPARER = ROOT / "scripts" / "prepare_turner_ensemble.py"
PREPARER = ROOT / "scripts" / "prepare_turner_sensitivity.py"
EXECUTABLE = ROOT / "build" / "aurorapic_cli"
ENSEMBLE_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE"
)
ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_SENSITIVITY"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=60,
    )


def main() -> int:
    require(EXECUTABLE.is_file(), "usage: build aurorapic_cli first")
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_sensitivity_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        normalized, case_path, _ = create_fixture(work)
        ensemble = work / "ensemble"
        prepared_ensemble = run([
            sys.executable, str(ENSEMBLE_PREPARER),
            str(case_path), str(normalized), "--output-dir", str(ensemble),
            "--seeds", "13507,24680,97531",
            "--acknowledge-cost", ENSEMBLE_ACKNOWLEDGEMENT,
        ])
        require(prepared_ensemble.returncode == 0, prepared_ensemble.stderr)
        manifest = ensemble / "ensemble.json"
        destination = work / "sensitivity"
        base_command = [
            sys.executable, str(PREPARER), str(manifest),
            "--baseline-seed", "13507",
            "--baseline-density-bias-percent", "2.4834268915580937",
            "--executable", str(EXECUTABLE),
            "--output-dir", str(destination),
        ]
        rejected = run(base_command)
        require(
            rejected.returncode == 2
            and ACKNOWLEDGEMENT in rejected.stderr
            and not destination.exists(),
            "Turner sensitivity bypassed its cost acknowledgement",
        )
        completed = run(base_command + ["--acknowledge-cost", ACKNOWLEDGEMENT])
        require(completed.returncode == 0, completed.stdout + completed.stderr)
        sensitivity_path = destination / "sensitivity.json"
        sensitivity = json.loads(sensitivity_path.read_text(encoding="utf-8"))
        require(
            sensitivity["turner_sensitivity_preparation_version"] == 1
            and sensitivity["baseline"]["seed"] == 13507
            and sensitivity["execution_policy"]["concurrent_runs_authorized"] == 1
            and sensitivity["predeclared_interpretation"]
                ["material_shift_percentage_points"] == 0.75
            and sensitivity["claim_boundary"]
                ["sensitivity_result_is_turner_benchmark_pass"] is False
            and sensitivity["launched"] is False,
            "Turner sensitivity manifest contract is incomplete",
        )
        variants = {variant["id"]: variant
                    for variant in sensitivity["variants"]}
        require(
            set(variants) == {
                "particles_2x", "timestep_2x",
                "grid_2x_fixed_particles", "grid_2x_same_ppc",
            }
            and variants["particles_2x"]["stage"] == 1
            and variants["timestep_2x"]["steps"] == 1024000
            and variants["grid_2x_fixed_particles"]
                ["particles_per_species"] == 65536
            and variants["grid_2x_same_ppc"]["nodes"] == 257,
            "Turner sensitivity variants drifted",
        )
        for identifier, variant in variants.items():
            deck = destination / variant["runtime_config"]
            text = deck.read_text(encoding="utf-8")
            require(
                sha256(deck) == variant["runtime_config_sha256"]
                and f"seed = {variant['seed']}" in text
                and variant["result_dir"] in text
                and variant["published_acceptance_applicable"] is False,
                f"Turner sensitivity variant {identifier} is inconsistent",
            )

        overwrite = run(base_command + ["--acknowledge-cost", ACKNOWLEDGEMENT])
        require(
            overwrite.returncode == 2 and "refusing to overwrite" in overwrite.stderr,
            "Turner sensitivity preparer overwrote an existing campaign",
        )

    print("Turner sensitivity preparation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
