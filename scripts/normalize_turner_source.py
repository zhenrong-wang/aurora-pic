#!/usr/bin/env python3
"""Normalize a verified local Turner CCP supplement for AuroraPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import zipfile

from import_lxcat import EV_TO_J, parse_lxcat
from verify_turner_source import VerificationError, verify


HELIUM_MASS_KG = 6.67e-27
EXPECTED_THRESHOLDS_EV = (0.0, 19.82, 20.61, 24.587)
ELECTRON_MEMBER = "turner_benchmark_he_electron_table.dat"
ION_MEMBER = "turner_benchmark_he_ion_table.dat"
RESULT_MEMBERS = (
    "turner_benchmark_results.dat",
    "turner_benchmark_refined_results.dat",
)


class NormalizationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise NormalizationError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path: Path, points: list[tuple[float, float]],
                columns: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# Locally normalized Turner CCP publisher supplement\n")
        stream.write(f"# {columns}\n")
        for first, second in points:
            stream.write(f"{first:.17g} {second:.17g}\n")


def write_manifest(path: Path, channels: list[dict[str, object]],
                   dataset_id: str, acquired: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("gas_data_version = 2\n")
        stream.write("units = si\n")
        stream.write("gas = He\n")
        stream.write(f"neutral_mass = {HELIUM_MASS_KG:.17g}\n")
        stream.write(f"dataset_id = {dataset_id}\n")
        stream.write("dataset_version = Turner-2013-supplement\n")
        stream.write(
            "data_provenance = Turner et al. 2013 publisher supplement\n"
        )
        stream.write(
            "citation = M. M. Turner et al., Physics of Plasmas 20, "
            "013507 (2013), doi:10.1063/1.4775084\n"
        )
        stream.write(f"retrieved = {acquired}\n")
        stream.write(
            "license = Publisher supplement terms - redistribution not approved\n"
        )
        for channel in channels:
            stream.write(f"\n[collision.{channel['name']}]\n")
            stream.write(f"type = {channel['type']}\n")
            stream.write(
                f"cross_section_file = {channel['cross_section_file']}\n"
            )
            stream.write(f"energy_scale = {EV_TO_J:.17g}\n")
            if channel.get("cross_section_scale") is not None:
                stream.write(
                    "cross_section_scale = "
                    f"{channel['cross_section_scale']:.17g}\n"
                )
            if channel.get("threshold_energy") is not None:
                stream.write(
                    "threshold_energy = "
                    f"{channel['threshold_energy']:.17g}\n"
                )
            stream.write(
                f"angular_model = {channel.get('angular_model', 'isotropic')}\n"
            )
            if channel.get("energy_frame") is not None:
                stream.write(f"energy_frame = {channel['energy_frame']}\n")


def parse_ion(text: str) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    for line in text.splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        try:
            energy, isotropic, backward = (float(item) for item in value.split())
        except ValueError as error:
            raise NormalizationError(f"invalid Turner ion row: {value!r}") from error
        rows.append((energy, isotropic, backward))
    require(len(rows) == 101, "Turner ion table must contain 101 rows")
    return rows


def parse_results(text: str) -> dict[int, list[list[float]]]:
    cases: dict[int, list[list[float]]] = {}
    current: int | None = None
    for line in text.splitlines():
        value = line.strip()
        if value.startswith("# Case "):
            try:
                current = int(value.removeprefix("# Case "))
            except ValueError as error:
                raise NormalizationError(f"invalid result case: {value!r}") from error
            cases[current] = []
        elif current is not None and value and not value.startswith("#"):
            try:
                row = [float(item) for item in value.split()]
            except ValueError as error:
                raise NormalizationError(
                    f"invalid Turner Case {current} row: {value!r}"
                ) from error
            require(len(row) == 7, "Turner result rows must have seven columns")
            cases[current].append(row)
    require(set(cases) == {1, 2, 3, 4},
            "Turner results must contain Cases 1--4")
    return cases


def write_reference(path: Path, rows: list[list[float]]) -> None:
    headings = (
        "x_m,electron_density_mean_m-3,"
        "electron_density_mean_stddev_m-3,"
        "electron_density_population_stddev_m-3,"
        "ion_density_mean_m-3,ion_density_mean_stddev_m-3,"
        "ion_density_population_stddev_m-3\n"
    )
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(headings)
        for row in rows:
            stream.write(",".join(f"{value:.17g}" for value in row) + "\n")


def normalized_files(stage: Path) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for path in sorted(stage.iterdir()):
        if path.is_file() and path.name != "audit.json":
            result[path.name] = {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
    return result


def write_normalized(
    stage: Path,
    report: dict[str, object],
    members: dict[str, bytes],
    acquired: str,
) -> None:
    source_path = stage / ".turner-electron-source.dat"
    source_path.write_bytes(members[ELECTRON_MEMBER])
    try:
        electron = parse_lxcat(source_path)
    finally:
        source_path.unlink(missing_ok=True)
    electron = [process for process in electron if process.target == "He"]
    require(
        [process.source_type for process in electron]
        == ["ELASTIC", "EXCITATION", "EXCITATION", "IONIZATION"],
        "Turner electron process order or identity differs",
    )
    require(
        all(
            abs(process.threshold_ev - expected) < 1e-12
            for process, expected in zip(electron, EXPECTED_THRESHOLDS_EV)
        ),
        "Turner electron thresholds differ from the published contract",
    )
    electron_names = (
        "electron_elastic",
        "electron_excitation_triplet",
        "electron_excitation_singlet",
        "electron_ionization",
    )
    electron_channels: list[dict[str, object]] = []
    for name, process in zip(electron_names, electron):
        table_name = f"{name}.dat"
        write_table(
            stage / table_name, list(process.points),
            "energy_eV cross_section_m2",
        )
        channel: dict[str, object] = {
            "name": name.removeprefix("electron_"),
            "type": process.source_type.lower(),
            "cross_section_file": table_name,
            "angular_model": "isotropic",
        }
        if process.threshold_ev > 0.0:
            channel["threshold_energy"] = process.threshold_ev * EV_TO_J
        electron_channels.append(channel)
    write_manifest(
        stage / "turner_he_electron.gas", electron_channels,
        "turner.helium.ccp.electron", acquired,
    )

    ion = parse_ion(members[ION_MEMBER].decode("utf-8"))
    write_table(
        stage / "ion_isotropic.dat",
        [(row[0], row[1]) for row in ion],
        "center_of_mass_energy_eV cross_section_1e-20_m2",
    )
    write_table(
        stage / "ion_backward.dat",
        [(row[0], row[2]) for row in ion],
        "center_of_mass_energy_eV cross_section_1e-20_m2",
    )
    ion_channels = [
        {
            "name": "isotropic",
            "type": "elastic",
            "cross_section_file": "ion_isotropic.dat",
            "cross_section_scale": 1.0e-20,
            "angular_model": "isotropic",
            "energy_frame": "center_of_mass",
        },
        {
            "name": "backward",
            "type": "elastic",
            "cross_section_file": "ion_backward.dat",
            "cross_section_scale": 1.0e-20,
            "angular_model": "backward",
            "energy_frame": "center_of_mass",
        },
    ]
    write_manifest(
        stage / "turner_he_ion.gas", ion_channels,
        "turner.helium.ccp.ion", acquired,
    )

    case_rows: dict[str, dict[str, int]] = {}
    for member in RESULT_MEMBERS:
        results = parse_results(members[member].decode("utf-8"))
        label = "refined" if "refined" in member else "benchmark"
        for case, rows in results.items():
            write_reference(stage / f"turner_case{case}_{label}.csv", rows)
        case_rows[label] = {
            str(case): len(rows) for case, rows in results.items()
        }

    audit = {
        "turner_normalization_version": 1,
        "case_id": report["case_id"],
        "doi": report["doi"],
        "source_artifact": {
            "name": report["artifact_name"],
            "bytes": report["artifact_bytes"],
            "sha256": report["artifact_sha256"],
        },
        "source_members": report["members"],
        "license": report["license"],
        "redistribution": report["redistribution"],
        "acquired": acquired,
        "constants": {
            "electron_volt_J": EV_TO_J,
            "helium_mass_kg": HELIUM_MASS_KG,
        },
        "transformations": {
            "electron_energy": "numeric eV retained; manifest scales eV to J",
            "electron_cross_section": "numeric m2 retained",
            "inelastic_threshold": "eV multiplied by exact elementary charge",
            "ion_energy": "numeric center-of-mass eV retained; manifest scales to J",
            "ion_cross_section": "numeric 1e-20 m2 retained; manifest scales by 1e-20",
            "results": "seven numeric columns converted losslessly to CSV",
            "interpolation": "linear with upper endpoint clamp; no resampling",
        },
        "case_rows": case_rows,
        "normalized_files": normalized_files(stage),
    }
    (stage / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def normalize(registry: Path, artifact: Path, output: Path) -> None:
    require(not output.exists(), f"refusing to overwrite output: {output}")
    report = verify(registry, artifact)
    acquired = str(report["acquired"])
    with zipfile.ZipFile(artifact) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        write_normalized(stage, report, members, acquired)
        os.replace(stage, output)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("registry", type=Path)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        normalize(args.registry, args.artifact, args.output_dir)
    except (NormalizationError, VerificationError, OSError,
            UnicodeError, zipfile.BadZipFile) as error:
        print(f"Turner normalization error: {error}", file=sys.stderr)
        return 2
    print(f"normalized Turner supplement locally to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
