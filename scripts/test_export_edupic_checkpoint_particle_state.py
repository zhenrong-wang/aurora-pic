#!/usr/bin/env python3
"""Focused tests for eduPIC checkpoint to half-step APS v3 export."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import struct
import tempfile

from export_checkpoint_particle_state import state_signature
from export_edupic_checkpoint_particle_state import ExportError, execute


def checkpoint_bytes() -> bytes:
    values = [0.25, 7.0, 2.0]
    values += [0.1, 0.9, 1.0, -2.0, 3.0, -4.0, 5.0, -6.0]
    values += [1.0, 0.4, 7.0, 8.0, 9.0]
    return struct.pack("<" + "d" * len(values), *values)


def main() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        checkpoint = root / "picdata.bin"
        output = root / "state.aps"
        manifest = root / "state.json"
        checkpoint.write_bytes(checkpoint_bytes())
        digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
        args = argparse.Namespace(
            checkpoint=checkpoint, output=output, manifest=manifest,
            expected_checkpoint_sha256=digest, expected_cycles=7,
            expected_time_s=0.25, expected_electrons=2, expected_ions=1)
        result = execute(args)
        expected = state_signature("si", {
            "electrons": [(0.1, 1.0, 3.0, 5.0),
                          (0.9, -2.0, -4.0, -6.0)],
            "ions": [(0.4, 7.0, 8.0, 9.0)]},
            version=3, velocity_staggering="leapfrog_half_step")
        assert result["particle_state_signature"] == expected
        assert "velocity_staggering leapfrog_half_step" in output.read_text()
        bad = argparse.Namespace(**vars(args))
        bad.output = root / "bad.aps"; bad.manifest = root / "bad.json"
        bad.expected_electrons = 3
        try:
            execute(bad)
        except ExportError:
            pass
        else:
            raise AssertionError("mismatched eduPIC population was accepted")
    print("eduPIC checkpoint particle-state exporter tests passed")


if __name__ == "__main__":
    main()
