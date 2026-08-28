#!/usr/bin/env python3
"""Bounded regression for corrected Turner density-window preparation."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
PREPARER = ROOT / "scripts" / "prepare_turner_density_window.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_POST_BENCHMARK_DIAGNOSTIC"


BASE = """config_version = 1
units = si
dimension = 1
velocity_dimensions = 3
nx = 129
length = 0.067
dt = 1.8436578171091445e-10
steps = 512000
output_interval = 400
output_dir = old-output
mode = steady_state
runtime_backend = serial
runtime_threads = 1
checkpoint_output = true
checkpoint_interval = 12800
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 499201
spatial_average_end_step = 512000
spatial_average_rf_frequency = 13560000
spatial_average_rf_cycles = 32
periodic_convergence = true
periodic_convergence_reset_on_restart = false
periodic_convergence_rf_frequency = 13560000
periodic_convergence_cycles_per_block = 32
periodic_convergence_minimum_blocks = 16
periodic_convergence_minimum_effective_blocks = 8

[collisions.electron_mcc]
model = null_collision
opportunity_sampling = single_bernoulli

[collisions.ion_mcc]
model = null_collision
opportunity_sampling = single_bernoulli
"""


def run(root: Path, base_text: str, suffix: str) -> subprocess.CompletedProcess[str]:
    base = root / f"base-{suffix}.cfg"
    checkpoint = root / "checkpoint_512000.apc"
    base.write_text(base_text, encoding="utf-8")
    checkpoint.write_bytes(b"checkpoint fixture")
    return subprocess.run([
        "python3", str(PREPARER), str(base), str(checkpoint),
        "--source-step", "512000",
        "--output-dir", str(root / f"output-{suffix}"),
        "--output-config", str(root / f"window-{suffix}.cfg"),
        "--report", str(root / f"window-{suffix}.json"),
        "--acknowledge-cost", ACKNOWLEDGEMENT,
    ], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_window_", dir=ROOT / "tmp"
    ) as temporary:
        root = Path(temporary)
        completed = run(root, BASE, "valid")
        assert completed.returncode == 0, completed.stderr
        deck = (root / "window-valid.cfg").read_text(encoding="utf-8")
        report = json.loads((root / "window-valid.json").read_text())
        required = (
            "mode = transient\n",
            "steps = 524800\n",
            f"restart_path = {root / 'checkpoint_512000.apc'}\n",
            "spatial_average_reset_on_restart = true\n",
            "spatial_average_start_step = 512001\n",
            "spatial_average_end_step = 524800\n",
            "spatial_average_interval = 1\n",
            "spatial_average_rf_cycles = 32\n",
            "periodic_convergence = true\n",
            "periodic_convergence_reset_on_restart = false\n",
            "periodic_convergence_rf_frequency = 13560000\n",
            "periodic_convergence_cycles_per_block = 32\n",
        )
        assert all(value in deck for value in required)
        assert deck.count("opportunity_sampling = single_bernoulli") == 2
        assert report["window"] == {
            "averaging_samples": 12800,
            "end_step": 524800,
            "reset_on_restart": True,
            "rf_cycles": 32,
            "start_step": 512001,
            "steps": 12800,
        }
        assert report["execution"]["launched"] is False
        assert report["execution"]["mode"] == "transient"
        assert report["execution"]["periodic_convergence"] == (
            "preserved_for_checkpoint_compatibility"
        )

        rejected = run(
            root, BASE.replace(
                "opportunity_sampling = single_bernoulli",
                "opportunity_sampling = poisson_clock", 1,
            ), "poisson",
        )
        assert rejected.returncode == 2
        assert "single_bernoulli" in rejected.stderr

    print("Turner density-window preparation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
