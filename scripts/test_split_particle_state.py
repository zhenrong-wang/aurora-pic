#!/usr/bin/env python3
"""Synthetic regression for exact reciprocal-weight particle splitting."""

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "split_particle_state.py"


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="aurorapic_split_") as temporary:
        root = Path(temporary)
        source = root / "source.aps"
        output = root / "split.aps"
        manifest = root / "manifest.json"
        source.write_text(
            "AuroraPIC-particle-state-v2\ndimension 1\nvelocity_dimensions 3\n"
            "units si\nweighting species_constant\n"
            "velocity_staggering time_centered\nparticle_count 2\nrecords\n"
            "particle electrons 0.25 0 0 1 2 3\n"
            "particle ions 0.75 0 0 -1 -2 -3\nend\n", encoding="utf-8")
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        command = [
            "python3", str(SCRIPT), str(source), str(output),
            "--manifest", str(manifest), "--expected-source-sha256", digest,
            "--expected-electrons", "1", "--expected-ions", "1",
            "--factor", "2", "--length", "1", "--nodes", "5",
            "--source-macro-weight", "10",
        ]
        completed = subprocess.run(command, text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError(completed.stderr or completed.stdout)
        report = json.loads(manifest.read_text(encoding="utf-8"))
        assert report["output_counts"] == {"electrons": 2, "ions": 2}
        assert report["child_macro_weight"] == 5
        assert report["charge_preservation"]["node_charge_relative_l1_error"] == 0
        assert output.read_text(encoding="utf-8").count(
            "particle electrons 0.25 0 0 1 2 3") == 2


if __name__ == "__main__":
    main()
