#!/usr/bin/env python3
"""Conservative regression tests for the local LXCat/BOLSIG+ importer."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
IMPORTER = ROOT / "scripts" / "import_lxcat.py"
FIXTURE = ROOT / "tests" / "fixtures" / "synthetic_lxcat.txt"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def base_command(output_dir: Path, source: Path = FIXTURE) -> list[str]:
    return [
        sys.executable,
        str(IMPORTER),
        str(source),
        "--output-dir",
        str(output_dir),
        "--gas",
        "Ar",
        "--neutral-mass",
        "6.6335209e-26",
        "--dataset-id",
        "aurorapic.synthetic.lxcat",
        "--dataset-version",
        "test-1",
        "--provenance",
        "AuroraPIC synthetic parser fixture",
        "--citation",
        "AuroraPIC test fixture",
        "--retrieved",
        "2026-07-28",
        "--license",
        "Synthetic test data",
        "--neutral-density",
        "1.0e20",
    ]


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_lxcat_test_", dir=project_tmp
    ) as temporary:
        root = Path(temporary)
        output_dir = root / "argon"
        subprocess.run(base_command(output_dir), check=True)

        manifest = (output_dir / "Ar.gas").read_text(encoding="utf-8")
        require(
            "gas_data_version = 2\n" in manifest
            and "units = si\n" in manifest
            and manifest.count("[collision.") == 3
            and manifest.count("angular_model = isotropic") == 1
            and "energy_scale = 1.6021766339999999e-19" in manifest
            and "threshold_energy = 3.2043532679999998e-19" in manifest
            and "threshold_energy = 8.0108831699999997e-19" in manifest,
            "generated gas manifest is incomplete",
        )
        audit = json.loads(
            (output_dir / "audit.json").read_text(encoding="utf-8")
        )
        require(
            audit["process_count"] == 3
            and audit["units"]["source_energy"] == "eV"
            and audit["units"]["manifest"] == "si"
            and audit["processes"][0]["angular_model"] == "isotropic"
            and audit["rate_envelope"][
                "recommended_max_frequency_s"
            ]
            > 0.0,
            "generated gas audit is incomplete",
        )
        for channel in (
            "elastic_001.dat",
            "excitation_001.dat",
            "ionization_001.dat",
        ):
            require(
                (output_dir / channel).is_file(),
                f"missing converted table {channel}",
            )

        attachment_source = root / "with_attachment.txt"
        attachment_source.write_text(
            FIXTURE.read_text(encoding="utf-8")
            + "ATTACHMENT\n"
            + "Ar -> Ar-\n"
            + "COMMENT: synthetic attachment validation data\n"
            + "--------------------\n"
            + "0.0 0.0\n"
            + "1.0 2.0e-21\n"
            + "20.0 0.0\n"
            + "--------------------\n",
            encoding="utf-8",
        )
        attachment_output = root / "argon_attachment"
        subprocess.run(
            base_command(attachment_output, attachment_source),
            check=True,
        )
        attachment_manifest = (
            attachment_output / "Ar.gas"
        ).read_text(encoding="utf-8")
        require(
            attachment_manifest.count("[collision.") == 4
            and "type = attachment\n" in attachment_manifest
            and (attachment_output / "attachment_001.dat").is_file(),
            "attachment channel was not converted",
        )

        existing = subprocess.run(
            base_command(output_dir),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        require(
            existing.returncode != 0
            and "already exists" in existing.stderr,
            "importer overwrote an existing package",
        )

        incomplete_source = root / "incomplete.txt"
        incomplete_source.write_text(
            "\n".join(
                FIXTURE.read_text(encoding="utf-8").splitlines()[:9]
            )
            + "\n--------------------\n",
            encoding="utf-8",
        )
        incomplete = subprocess.run(
            base_command(root / "incomplete", incomplete_source),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        require(
            incomplete.returncode != 0
            and "incomplete electron collision set" in incomplete.stderr,
            "importer accepted an incomplete set without opt-in",
        )

        mass_mismatch_command = base_command(root / "mass_mismatch")
        neutral_mass_index = (
            mass_mismatch_command.index("--neutral-mass") + 1
        )
        mass_mismatch_command[neutral_mass_index] = "1.0e-26"
        mass_mismatch = subprocess.run(
            mass_mismatch_command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        require(
            mass_mismatch.returncode != 0
            and "elastic mass ratio disagrees" in mass_mismatch.stderr,
            "importer accepted inconsistent elastic mass metadata",
        )

        unsafe_metadata_command = base_command(root / "unsafe_metadata")
        citation_index = unsafe_metadata_command.index("--citation") + 1
        unsafe_metadata_command[citation_index] = "Citation # truncated"
        unsafe_metadata = subprocess.run(
            unsafe_metadata_command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        require(
            unsafe_metadata.returncode != 0
            and "manifest comment delimiter" in unsafe_metadata.stderr,
            "importer accepted metadata that would be truncated",
        )

    print("LXCat importer validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
