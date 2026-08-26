#!/usr/bin/env python3
"""Focused tests for cycle phase-snapshot reductions."""

import tempfile
from pathlib import Path

from analyze_edupic_phase_snapshot_onset import (
    profile_metrics, relative_range, snapshot_set_sha256,
)


def main() -> None:
    metrics = profile_metrics([2.0, 2.0], [-3.0, -3.0], 1.0, 0.2, 0.8)
    assert abs(metrics["regional_mean_squared_field_V2_m2"] - 4.0) < 1e-12
    assert abs(metrics["regional_charge_density_rms_C_m3"] - 3.0) < 1e-12
    assert abs(relative_range([9.0, 11.0]) - 0.2) < 1e-12
    with tempfile.TemporaryDirectory() as directory:
        first, second = Path(directory) / "fields_400.csv", Path(directory) / "fields_800.csv"
        first.write_bytes(b"a"); second.write_bytes(b"b")
        assert snapshot_set_sha256([second, first]) == snapshot_set_sha256([first, second])
    print("phase snapshot onset analyzer tests passed")


if __name__ == "__main__":
    main()
