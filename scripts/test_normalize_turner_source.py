#!/usr/bin/env python3
"""Bounded end-to-end regression for Turner supplement normalization."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
NORMALIZER = ROOT / "scripts" / "normalize_turner_source.py"
NAMES = (
    "turner_benchmark_he_electron_table.dat",
    "turner_benchmark_he_ion_table.dat",
    "turner_benchmark_refined_results.dat",
    "turner_benchmark_results.dat",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def electron_table() -> bytes:
    blocks = [
        ("ELASTIC", "He", "1.36e-4", 0.0),
        ("EXCITATION", "He -> triplet", "19.82", 19.82),
        ("EXCITATION", "He -> singlet", "20.61", 20.61),
        ("IONIZATION", "He -> He+", "24.587", 24.587),
    ]
    lines = ["Generated on 28 Nov 2012.", "Biagi-v7.1"]
    for kind, target, parameter, threshold in blocks:
        lines.extend([
            kind, target, parameter, "-----",
            f"{threshold:.17g} 0",
            f"{threshold + 100.0:.17g} 1e-20",
            "-----",
        ])
    return ("\n".join(lines) + "\n").encode()


def ion_table() -> bytes:
    lines = [
        "# Centre of mass energy",
        "# Isotropic scattering",
        "# Backward scattering",
    ]
    for index in range(101):
        energy = 10000.0 * index / 100.0
        lines.append(f"{energy:.17g} 1 2")
    return ("\n".join(lines) + "\n").encode()


def result_table() -> bytes:
    lines = ["# synthetic Turner-format normalization fixture"]
    for case in range(1, 5):
        lines.extend([
            f"# Case {case}",
            "0 1 1 1 1 1 1",
            "0.067 2 1 1 2 1 1",
        ])
    return ("\n".join(lines) + "\n").encode()


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_normalize_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        members = {
            NAMES[0]: electron_table(),
            NAMES[1]: ion_table(),
            NAMES[2]: result_table(),
            NAMES[3]: result_table(),
        }
        artifact = work / "fixture.zip"
        with zipfile.ZipFile(artifact, "w", zipfile.ZIP_DEFLATED) as archive:
            for name, data in members.items():
                archive.writestr(name, data)
        kinds = {
            NAMES[0]: (
                "lxcat-electron-table",
                "process_counts = ELASTIC:1,EXCITATION:2,IONIZATION:1",
            ),
            NAMES[1]: ("turner-ion-table", "rows = 101"),
            NAMES[2]: (
                "turner-results", "case_rows = 1:2,2:2,3:2,4:2"
            ),
            NAMES[3]: (
                "turner-results", "case_rows = 1:2,2:2,3:2,4:2"
            ),
        }
        sections = []
        for name, data in members.items():
            kind, semantic = kinds[name]
            sections.append(
                f"[member.{name}]\n"
                f"bytes = {len(data)}\n"
                f"sha256 = {sha256(data)}\n"
                f"kind = {kind}\n{semantic}\n"
            )
        registry = work / "fixture.sources"
        registry.write_text(
            "source_registry_version = 1\n"
            "case_id = synthetic-turner-normalization\n"
            "doi = 10.0/synthetic\n"
            "landing_url = https://example.invalid\n\n"
            "[source.publisher_supplement]\n"
            "artifact_name = fixture.zip\n"
            f"artifact_bytes = {artifact.stat().st_size}\n"
            f"sha256 = {sha256(artifact.read_bytes())}\n"
            "acquisition = synthetic\n"
            "acquired = 2026-07-30\n"
            "license = synthetic\n"
            "redistribution = test-only\n"
            f"members = {', '.join(NAMES)}\n\n"
            + "\n".join(sections),
            encoding="utf-8",
        )

        output = work / "normalized"
        completed = subprocess.run(
            [sys.executable, str(NORMALIZER), str(registry), str(artifact),
             "--output-dir", str(output)],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        require(completed.returncode == 0,
                f"Turner normalizer rejected valid input: {completed.stderr}")
        electron_manifest = (
            output / "turner_he_electron.gas"
        ).read_text(encoding="utf-8")
        ion_manifest = (
            output / "turner_he_ion.gas"
        ).read_text(encoding="utf-8")
        audit = json.loads((output / "audit.json").read_text(encoding="utf-8"))
        require(
            "threshold_energy = 3.1755140885879999e-18" in electron_manifest
            and "threshold_energy = 3.9392716900158001e-18"
            in electron_manifest
            and electron_manifest.count("[collision.") == 4,
            "electron threshold conversion or channel set is incomplete",
        )
        require(
            ion_manifest.count("energy_frame = center_of_mass") == 2
            and "angular_model = backward" in ion_manifest
            and "cross_section_scale = 9.9999999999999995e-21"
            in ion_manifest,
            "ion unit/frame/scattering normalization is incomplete",
        )
        require(
            audit["transformations"]["interpolation"].startswith("linear")
            and audit["case_rows"]["benchmark"]["1"] == 2
            and len(audit["normalized_files"]) == 16,
            "normalization audit is incomplete",
        )

        refused = subprocess.run(
            [sys.executable, str(NORMALIZER), str(registry), str(artifact),
             "--output-dir", str(output)],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        require(refused.returncode == 2 and "refusing to overwrite" in refused.stderr,
                "Turner normalizer overwrote an existing package")

    print("Turner supplement normalization passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
