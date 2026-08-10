#!/usr/bin/env python3
"""Synthetic regression for the bounded eduPIC argon case preparer."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

from prepare_edupic_argon_case import PreparationError, prepare


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="aurorapic-edupic-case-") as text:
        root = Path(text)
        package = root / "package"
        package.mkdir()
        artifacts = {
            "edupic_argon_electron.gas": "electron manifest\n",
            "edupic_argon_ion.gas": "ion manifest\n",
            "electron_elastic.dat": "0 1\n1 1\n",
            "electron_excitation.dat": "0 0\n1 0\n",
            "electron_ionization.dat": "0 0\n1 0\n",
            "ion_isotropic.dat": "0 1\n1 1\n",
            "ion_backward.dat": "0 1\n1 1\n",
        }
        identities = {}
        for name, content in artifacts.items():
            path = package / name
            path.write_text(content, encoding="utf-8")
            identities[name] = {
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        source_table_sha = "a" * 64
        implementation_sha = "b" * 64
        audit = {
            "case_id": "edupic-1.0-default-argon-ccp",
            "source": {
                "commit": "test-commit",
                "implementation_sha256": implementation_sha,
                "table_sha256": source_table_sha,
            },
            "contract": {
                "rows": 1_000_000,
                "energy_step_ev": 0.001,
                "ionization_kinematics": "opal_beaty_peterson",
                "ion_energy_frame": "center_of_mass",
            },
            "grid_rate_envelope": {
                "electron_peak_frequency_s": 5.0e8,
                "ion_peak_frequency_s": 5.0e7,
            },
            "artifacts": identities,
        }
        (package / "audit.json").write_text(
            json.dumps(audit), encoding="utf-8"
        )
        case = root / "case.case"
        case.write_text(
            f"""case_manifest_version = 1
case_id = edupic-1.0-default-argon-ccp
[reference]
commit = test-commit
implementation_sha256 = {implementation_sha}
generated_cross_sections_sha256 = {source_table_sha}
[package]
electron_manifest = edupic_argon_electron.gas
electron_manifest_sha256 = {identities['edupic_argon_electron.gas']['sha256']}
ion_manifest = edupic_argon_ion.gas
ion_manifest_sha256 = {identities['edupic_argon_ion.gas']['sha256']}
electron_elastic_sha256 = {identities['electron_elastic.dat']['sha256']}
electron_excitation_sha256 = {identities['electron_excitation.dat']['sha256']}
electron_ionization_sha256 = {identities['electron_ionization.dat']['sha256']}
ion_isotropic_sha256 = {identities['ion_isotropic.dat']['sha256']}
ion_backward_sha256 = {identities['ion_backward.dat']['sha256']}
[physics]
length_m = 0.025
neutral_temperature_k = 350
neutral_density_m3 = 2e21
rf_frequency_hz = 13.56e6
voltage_amplitude_v = 250
electron_mass_kg = 9.1e-31
argon_mass_kg = 6.6e-26
elementary_charge_c = 1.6e-19
macro_particle_weight = 7e8
initial_particles_per_species = 1000
[numerics]
nodes = 400
steps_per_rf_cycle = 4000
max_particles_per_species = 1000000
max_candidates_per_particle = 16
electron_max_frequency_s = 1e9
ion_max_frequency_s = 1e8
seed = 7
[compatibility]
geometry_and_drive = matched
ion_timestep = unmatched_test
[authorization]
maximum_preflight_steps = 4000
production_launch_authorized = false
physics_claim = none_contract_preflight_only
""",
            encoding="utf-8",
        )
        output = root / "preflight.cfg"
        report = root / "preflight.json"
        arguments = argparse.Namespace(
            case_manifest=case,
            package_dir=package,
            output=output,
            report=report,
            output_dir=root / "run",
            steps=2,
            seed=None,
        )
        prepare(arguments)
        deck = output.read_text(encoding="utf-8")
        result = json.loads(report.read_text(encoding="utf-8"))
        require(
            "phi_left_phase = 1.5707963267948966" in deck
            and "neutral_temperature = 0" in deck
            and "neutral_temperature = 350" in deck,
            "generated deck lost the matched drive or neutral contracts",
        )
        require(
            result["deck"]["steps"] == 2
            and result["unresolved_contract_items"] == ["ion_timestep"]
            and not result["production_launch_authorized"],
            "preflight report changed its authorization boundary",
        )
        arguments.output = root / "too-long.cfg"
        arguments.report = root / "too-long.json"
        arguments.steps = 4001
        try:
            prepare(arguments)
        except PreparationError as error:
            require("preflight steps" in str(error), "unexpected rejection")
        else:
            raise RuntimeError("preparer authorized more than one RF cycle")
    print("eduPIC argon case preparer regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
