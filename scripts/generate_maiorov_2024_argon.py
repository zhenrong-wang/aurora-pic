#!/usr/bin/env python3
"""Generate the argon collision tables published by Maiorov et al. (2024).

The source paper gives analytic fits in Angstrom squared as functions of
electron energy.  This script evaluates those equations on a deterministic
logarithmic grid and emits an AuroraPIC gas package.  It does not copy an
external cross-section dataset.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


IONIZATION_EV = 15.759
EXCITATION_EV = 11.50
ANGSTROM2_TO_M2 = 1.0e-20


def elastic_angstrom2(energy_ev: float) -> float:
    x = energy_ev / IONIZATION_EV
    terms = (
        (0.02, 24.0, 1.03, 2.83, 1.0),
        (7.64, -65.4, 1961.0, 1.37, 0.455),
    )
    value = sum(
        (sigma + alpha * x**delta) / (1.0 + beta * x**gamma)
        for sigma, alpha, beta, gamma, delta in terms
    )
    if value < 0.0 or not math.isfinite(value):
        raise RuntimeError(
            f"invalid elastic fit at {energy_ev:.17g} eV: {value}"
        )
    return value


def threshold_fit_angstrom2(
    energy_ev: float,
    threshold_ev: float,
    alpha: float,
    beta: float,
    gamma: float,
    delta: float,
) -> float:
    if energy_ev <= threshold_ev:
        return 0.0
    excess = energy_ev / threshold_ev - 1.0
    return alpha * excess**delta / (1.0 + beta * excess) ** gamma


def energies(max_energy_ev: float, points: int) -> list[float]:
    lower = 1.0e-5
    logarithmic = [
        lower * (max_energy_ev / lower) ** (index / (points - 1))
        for index in range(points)
    ]
    return sorted(
        set([0.0, EXCITATION_EV, IONIZATION_EV, max_energy_ev] + logarithmic)
    )


def write_table(
    path: Path,
    grid: list[float],
    evaluator,
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Generated from Maiorov et al. (2024) analytic fit\n")
        stream.write("# energy_eV cross_section_m2\n")
        for energy_ev in grid:
            cross_section_m2 = evaluator(energy_ev) * ANGSTROM2_TO_M2
            stream.write(f"{energy_ev:.17g} {cross_section_m2:.17g}\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--max-energy-ev", type=float, default=400.0)
    parser.add_argument("--points", type=int, default=6001)
    args = parser.parse_args()
    if not math.isfinite(args.max_energy_ev) or args.max_energy_ev <= 20.0:
        parser.error("--max-energy-ev must be finite and greater than 20")
    if args.points < 100:
        parser.error("--points must be at least 100")

    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    grid = energies(args.max_energy_ev, args.points)
    elastic = output / "electron_elastic.dat"
    excitation = output / "electron_excitation.dat"
    ionization = output / "electron_ionization.dat"
    write_table(elastic, grid, elastic_angstrom2)
    write_table(
        excitation,
        grid,
        lambda energy: threshold_fit_angstrom2(
            energy, EXCITATION_EV, 0.802, 0.229, 1.55, 0.702
        ),
    )
    write_table(
        ionization,
        grid,
        lambda energy: threshold_fit_angstrom2(
            energy, IONIZATION_EV, 3.19, 0.326, 1.92, 1.08
        ),
    )

    version = hashlib.sha256(
        (sha256(elastic) + sha256(excitation) + sha256(ionization)).encode()
    ).hexdigest()[:40]
    manifest = output / "maiorov_2024_argon.gas"
    manifest.write_text(
        f"""gas_data_version = 2
units = si
gas = argon
neutral_mass = 6.6335209e-26
dataset_id = maiorov.2024.argon.analytic-fits
dataset_version = {version}
data_provenance = Deterministic evaluation of published analytic fits
citation = Maiorov et al., Plasma Physics Reports 50, 1029-1041 (2024), doi:10.1134/S1063780X24601263
retrieved = 2026-08-24
license = Generated numerical evaluation of equations and parameters cited above

[collision.elastic]
type = elastic
cross_section_file = electron_elastic.dat
energy_scale = 1.602176634e-19
cross_section_interpolation = linear
angular_model = isotropic

[collision.excitation]
type = excitation
cross_section_file = electron_excitation.dat
energy_scale = 1.602176634e-19
cross_section_interpolation = linear
threshold_energy = {EXCITATION_EV * 1.602176634e-19:.17g}
inelastic_transform = finite_mass_center_of_mass

[collision.ionization]
type = ionization
cross_section_file = electron_ionization.dat
energy_scale = 1.602176634e-19
cross_section_interpolation = linear
threshold_energy = {IONIZATION_EV * 1.602176634e-19:.17g}
ionization_kinematics = equal_energy_isotropic
inelastic_transform = finite_mass_center_of_mass
""",
        encoding="utf-8",
        newline="\n",
    )
    print(f"manifest={manifest}")
    print(f"manifest_sha256={sha256(manifest)}")
    print(f"dataset_version={version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
