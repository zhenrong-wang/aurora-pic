#!/usr/bin/env python3
"""Synthetic regression for guarded 1D3V checkpoint particle export."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
EXPORTER = ROOT / "scripts" / "export_checkpoint_particle_state.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        checkpoint = work / "input.apc"
        state = work / "output.aps"
        manifest = work / "manifest.json"
        checkpoint.write_text(
            "AuroraPIC-checkpoint-v15\n"
            "dimension 1\n"
            "units si 1 8.854e-12\n"
            "velocity_dimensions 3\n"
            "step 40\n"
            "time 1.25\n"
            "species_count 2\n"
            "rng 1 2 3\n"
            "species 0 electrons 3\n"
            "0.1 1 2 3 0.5 1\n"
            "0.2 4 5 6 3.5 0\n"
            "0.3 7 8 9 6.5 1\n"
            "species 1 ions 1\n"
            "0.4 -1 -2 -3 -0.5 1\n",
            encoding="utf-8")
        command = [
            "python3", str(EXPORTER), str(checkpoint), str(state),
            "--manifest", str(manifest),
            "--expected-checkpoint-sha256", sha256(checkpoint),
            "--expected-step", "40",
            "--expected-species", "electrons=2",
            "--expected-species", "ions=1",
        ]
        completed = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout)
        report = json.loads(manifest.read_text(encoding="utf-8"))
        contents = state.read_text(encoding="utf-8")
        if not (
            report["particle_count"] == 3
            and report["species"]["electrons"]["discarded_inactive_records"] == 1
            and report["particle_state_signature"] > 0
            and "AuroraPIC-particle-state-v2\n" in contents
            and "velocity_dimensions 3\n" in contents
            and "particle electrons 0.10000000000000001 0 0 1 2 3\n" in contents
            and " 4 5 6 " not in contents
            and "particle ions 0.40000000000000002 0 0 -1 -2 -3\n" in contents
        ):
            raise RuntimeError("checkpoint export changed records or provenance")
        repeated = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if repeated.returncode == 0 or "refusing to overwrite" not in repeated.stderr:
            raise RuntimeError("checkpoint exporter overwrote existing evidence")
        wrong_hash = command.copy()
        wrong_hash[wrong_hash.index("--expected-checkpoint-sha256") + 1] = "0" * 64
        state.unlink()
        manifest.unlink()
        rejected = subprocess.run(
            wrong_hash, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if rejected.returncode == 0 or "SHA-256" not in rejected.stderr:
            raise RuntimeError("checkpoint exporter accepted the wrong source hash")
    print("checkpoint particle-state export regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
