#!/usr/bin/env python3
"""Synthetic regression for exact-charge neutral-pair augmentation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "augment_neutral_particle_pairs.py"


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        source = work / "source.aps"
        output = work / "output.aps"
        manifest = work / "manifest.json"
        source.write_text(
            "AuroraPIC-particle-state-v2\ndimension 1\n"
            "velocity_dimensions 3\nunits si\n"
            "weighting species_constant\nvelocity_staggering time_centered\n"
            "particle_count 5\nrecords\n"
            "particle electrons 0.1 0 0 1 2 3\n"
            "particle electrons 0.6 0 0 4 5 6\n"
            "particle ions 0.12 0 0 -1 -2 -3\n"
            "particle ions 0.62 0 0 -4 -5 -6\n"
            "particle ions 0.8 0 0 -7 -8 -9\nend\n",
            encoding="utf-8")
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        command = [
            "python3", str(SCRIPT), str(source), str(output),
            "--manifest", str(manifest), "--expected-source-sha256", digest,
            "--expected-electrons", "2", "--expected-ions", "3",
            "--length", "1", "--bins", "2", "--nodes", "5",
            "--macro-weight", "10", "--added-pairs", "3",
        ]
        completed = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout)
        report = json.loads(manifest.read_text(encoding="utf-8"))
        contents = output.read_text(encoding="utf-8")
        if not (
            report["output_counts"] == {"electrons": 5, "ions": 6}
            and report["pairing_pool_size"] == 2
            and report["charge_preservation"][
                "node_charge_relative_l1_error"] < 1e-14
            and contents.count("particle electrons 0.11 0 0 1 2 3") == 1
            and contents.count("particle ions 0.11 0 0 -1 -2 -3") == 1
        ):
            raise RuntimeError("neutral-pair augmentation contract changed")
        repeated = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if repeated.returncode == 0 or "refusing to overwrite" not in repeated.stderr:
            raise RuntimeError("neutral-pair augmentation overwrote evidence")
    print("neutral-pair augmentation regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
