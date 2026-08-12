#!/usr/bin/env python3
"""Synthetic regression for quasi-neutral APS warm-state preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_quasineutral_particle_state.py"


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        source = work / "source.aps"
        output = work / "output.aps"
        manifest = work / "manifest.json"
        source.write_text(
            "AuroraPIC-particle-state-v2\n"
            "dimension 1\nvelocity_dimensions 3\nunits si\n"
            "weighting species_constant\nvelocity_staggering time_centered\n"
            "particle_count 8\nrecords\n"
            "particle electrons 0.1 0 0 1 2 3\n"
            "particle electrons 0.2 0 0 4 5 6\n"
            "particle electrons 0.6 0 0 7 8 9\n"
            "particle ions 0.12 0 0 -1 -2 -3\n"
            "particle ions 0.22 0 0 -4 -5 -6\n"
            "particle ions 0.3 0 0 -7 -8 -9\n"
            "particle ions 0.62 0 0 -10 -11 -12\n"
            "particle ions 0.7 0 0 -13 -14 -15\nend\n",
            encoding="utf-8")
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        command = [
            "python3", str(SCRIPT), str(source), str(output),
            "--manifest", str(manifest),
            "--expected-source-sha256", digest,
            "--expected-electrons", "3", "--expected-ions", "5",
            "--length", "1", "--bins", "2", "--nodes", "5",
            "--source-weight", "10", "--weight-factor", "2",
        ]
        completed = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout)
        report = json.loads(manifest.read_text(encoding="utf-8"))
        contents = output.read_text(encoding="utf-8")
        if not (
            report["decomposition"]["paired_bulk_per_species"] == 3
            and report["decomposition"]["source_residual_particles"] == 2
            and report["decomposition"]["retained_residual_particles"] == 1
            and report["output_counts"] == {"electrons": 3, "ions": 4}
            and report["transformed_macro_weight"] == 20
            and "particle electrons 0.11 0 0 1 2 3" in contents
            and "particle ions 0.11 0 0 -1 -2 -3" in contents
        ):
            raise RuntimeError("quasi-neutral transform changed its contract")
        repeated = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if repeated.returncode == 0 or "refusing to overwrite" not in repeated.stderr:
            raise RuntimeError("quasi-neutral transform overwrote evidence")
    print("quasi-neutral particle-state regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
