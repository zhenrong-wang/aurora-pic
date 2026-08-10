#!/usr/bin/env python3
"""Conservative regression for the pinned eduPIC table importer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile

from import_edupic_cross_sections import ImportFailure, convert


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expect_failure(action, text: str) -> None:
    try:
        action()
    except ImportFailure as error:
        require(text in str(error), f"unexpected failure: {error}")
        return
    raise RuntimeError("expected importer failure")


def fixture(path: Path, *, bad_grid: bool = False) -> None:
    energies = [0.0, 0.001, 0.002, 0.003]
    if bad_grid:
        energies[2] = 0.0025
    rows = [
        f"{energy:.4f} {1 + index} 0 0 {2 + index} {3 + index}\n"
        for index, energy in enumerate(energies)
    ]
    path.write_text("".join(rows), encoding="utf-8")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="aurorapic-edupic-import-") as root_text:
        root = Path(root_text)
        source = root / "cross_sections.dat"
        fixture(source)
        output = root / "argon"
        audit = convert(
            source,
            output,
            source_sha256=digest(source),
            retrieved="2026-08-10",
            expected_rows=4,
        )
        require(audit["contract"]["rows"] == 4, "row count was not audited")
        electron_manifest = (output / "edupic_argon_electron.gas").read_text()
        ion_manifest = (output / "edupic_argon_ion.gas").read_text()
        require(
            "ionization_kinematics = opal_beaty_peterson" in electron_manifest
            and "ionization_ejected_energy_scale" in electron_manifest
            and electron_manifest.count(
                "cross_section_interpolation = lower_bin"
            ) == 3,
            "electron manifest lost ionization kinematics",
        )
        require(
            electron_manifest.count(
                "inelastic_transform = finite_mass_center_of_mass"
            ) == 2,
            "electron manifest lost finite-mass inelastic transforms",
        )
        require(
            ion_manifest.count("energy_frame = center_of_mass") == 2
            and ion_manifest.count(
                "cross_section_interpolation = lower_bin"
            ) == 2
            and "angular_model = backward" in ion_manifest,
            "ion manifest lost center-of-mass scattering contracts",
        )
        persisted = json.loads((output / "audit.json").read_text())
        require(
            len(persisted["artifacts"]) == 7
            and all(
                len(identity["sha256"]) == 64
                for identity in persisted["artifacts"].values()
            ),
            "output identities are incomplete",
        )
        expect_failure(
            lambda: convert(
                source,
                root / "bad-hash",
                source_sha256="0" * 64,
                retrieved="2026-08-10",
                expected_rows=4,
            ),
            "SHA-256 differs",
        )
        malformed = root / "bad-grid.dat"
        fixture(malformed, bad_grid=True)
        expect_failure(
            lambda: convert(
                malformed,
                root / "bad-grid-output",
                source_sha256=digest(malformed),
                retrieved="2026-08-10",
                expected_rows=4,
            ),
            "energy grid mismatch",
        )
        expect_failure(
            lambda: convert(
                source,
                output,
                source_sha256=digest(source),
                retrieved="2026-08-10",
                expected_rows=4,
            ),
            "refusing to overwrite",
        )
    print("eduPIC cross-section importer regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
