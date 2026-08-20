#!/usr/bin/env python3
"""Focused tests for the transport-balance audit."""

import csv
import math
from pathlib import Path
import tempfile

from analyze_edupic_transport_balance import (
    candidate_balance, matrix_phase_space, parse_flux, phase_resolution_metrics,
)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_candidate_ledger() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        write_csv(root / "collisions.csv", [
            "step", "cumulative_collisions_electron_mcc.ionization"], [
                {"step": 10, "cumulative_collisions_electron_mcc.ionization": 7},
                {"step": 20, "cumulative_collisions_electron_mcc.ionization": 12},
            ])
        boundary_fields = ["step", "absorbed_left_count_electrons",
                           "absorbed_right_count_electrons",
                           "absorbed_left_count_ions", "absorbed_right_count_ions"]
        write_csv(root / "boundary_losses.csv", boundary_fields, [
            dict(zip(boundary_fields, [10, 2, 3, 4, 5])),
            dict(zip(boundary_fields, [20, 4, 5, 6, 7])),
        ])
        scalar_fields = ["step", "time", "live_particles_electrons",
                         "live_particles_ions"]
        write_csv(root / "scalars.csv", scalar_fields, [
            dict(zip(scalar_fields, [10, 1.0, 100, 100])),
            dict(zip(scalar_fields, [20, 3.0, 101, 101])),
        ])
        result = candidate_balance(root)
        assert result["ionization_macro_events"] == 5
        assert result["electron_wall_loss_macro_events"] == 4
        assert result["ion_wall_loss_macro_events"] == 4
        assert result["exact_macro_particle_ledger_closure"] is True


def test_phase_coarsening() -> None:
    # A phase-constant field is unchanged by conservative phase coarsening.
    candidate = [float(node) for _phase in range(200) for node in range(400)]
    reference = [[float(node)] * 200 for node in range(400)]
    result = phase_resolution_metrics(candidate, reference)
    assert result["fine_200_phase_bins"]["relative_l2"] == 0.0
    assert result["coarse_grained_16_phase_bins"]["relative_l2"] == 0.0
    # Exercise the phase-major to space-major mapping explicitly.
    assert matrix_phase_space(candidate)[27][103] == 27.0


def test_info_flux_parser() -> None:
    text = """
Electron flux at powered electrode = 2.1e18 [m^-2 s^-1]
Electron flux at grounded electrode = 2.3e18 [m^-2 s^-1]
"""
    assert math.isclose(parse_flux(text, "Electron"), 4.4e18)


if __name__ == "__main__":
    test_candidate_ledger()
    test_phase_coarsening()
    test_info_flux_parser()
    print("transport-balance analyzer tests passed")
