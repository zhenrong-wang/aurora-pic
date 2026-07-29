#!/usr/bin/env python3
"""Digitize the seven-code envelope in Charoy et al. (2019), Figure 6.

This is a reproducible screening aid, not a substitute for the authors'
native tables.  It writes derived data outside the source tree by default and
refuses unrecognized PDF/vector input.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


PDF_SHA256 = "3bd8999a50217450f4ce233ebd468a652b4f94250f2dca00a1903256c7598b96"
CASE_ID = "landmark-axial-azimuthal-2019"
COLORS = (
    "rgb(0%, 0%, 100%)",
    "rgb(0%, 50%, 0%)",
    "rgb(100%, 0%, 0%)",
    "rgb(75%, 75%, 0%)",
    "rgb(75%, 0%, 75%)",
    "rgb(0%, 75%, 75%)",
    "rgb(0%, 0%, 0%)",
)
NUMBER = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
)
X_ZERO = 1582.185420
X_HALF_CM = 2450.839854
CALIBRATIONS = (
    ("electric_x_v_m", 15390.431058, 961.841158, 1.0e4),
    ("ion_density_m3", 8044.376237, 781.390538 / 0.5, 1.0e17),
    ("electron_temperature_ev", 1178.960919, 1041.947743 / 10.0, 1.0),
)


class DigitizeError(RuntimeError):
    pass


def digest(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise DigitizeError(f"cannot hash {path}: {error}") from error


def case_id(path: Path) -> str:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string(
            "[global]\n" + path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, configparser.Error) as error:
        raise DigitizeError(f"cannot read case manifest {path}: {error}") from error
    return parser["global"].get("case_id", "").strip()


def render_page(pdf: Path, destination: Path) -> None:
    executable = shutil.which("pdftocairo")
    if executable is None:
        raise DigitizeError("pdftocairo is required to render the vector page")
    completed = subprocess.run(
        [
            executable, "-f", "20", "-l", "20", "-svg",
            str(pdf), str(destination),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise DigitizeError(
            f"pdftocairo failed: {completed.stderr.strip()}"
        )


def figure_paths(svg: Path) -> list[list[tuple[float, float]]]:
    try:
        root = ET.parse(svg).getroot()
    except (OSError, ET.ParseError) as error:
        raise DigitizeError(f"cannot parse rendered SVG: {error}") from error
    selected: list[tuple[str, list[tuple[float, float]]]] = []
    for element in root.iter():
        data = element.attrib.get("d", "")
        if (
            not element.tag.endswith("path")
            or element.attrib.get("stroke-width") != "30"
            or len(data) <= 1000
        ):
            continue
        values = [float(value) for value in NUMBER.findall(data)]
        if len(values) % 2:
            raise DigitizeError("Figure 6 curve has an odd coordinate count")
        selected.append(
            (
                element.attrib.get("stroke", ""),
                list(zip(values[0::2], values[1::2])),
            )
        )
    if len(selected) != 21:
        raise DigitizeError(
            f"expected 21 Figure 6 curves, found {len(selected)}"
        )
    for panel in range(3):
        observed = tuple(
            color for color, _ in selected[panel * 7:(panel + 1) * 7]
        )
        if observed != COLORS:
            raise DigitizeError(
                f"Figure 6 panel {panel + 1} color signature changed"
            )
    return [points for _, points in selected]


def interpolate(points: list[tuple[float, float]], x: float) -> float:
    ordered = sorted(points)
    if x < ordered[0][0] or x > ordered[-1][0]:
        raise DigitizeError(f"sample coordinate {x} lies outside a curve")
    for index in range(1, len(ordered)):
        x1, y1 = ordered[index - 1]
        x2, y2 = ordered[index]
        if x1 <= x <= x2:
            if x2 == x1:
                continue
            fraction = (x - x1) / (x2 - x1)
            return y1 + fraction * (y2 - y1)
    raise DigitizeError(f"could not interpolate sample coordinate {x}")


def csv_text(curves: list[list[tuple[float, float]]]) -> str:
    import io

    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    header = ["coordinate_m"]
    for name, _, _, _ in CALIBRATIONS:
        header.extend((name, name.replace("_v_m", "_uncertainty_v_m")
                      if name == "electric_x_v_m"
                      else name.replace("_m3", "_uncertainty_m3")
                      if name == "ion_density_m3"
                      else "electron_temperature_uncertainty_ev"))
    writer.writerow(header)
    # The workstation and production diagnostics both contain this 0.2 mm grid.
    for sample in range(1, 125):
        coordinate_m = sample * 0.0002
        x_raw = X_ZERO + (coordinate_m / 0.005) * X_HALF_CM
        row: list[float] = [coordinate_m]
        for panel, (_, y_zero, raw_per_unit, scale) in enumerate(CALIBRATIONS):
            values = [
                (interpolate(curves[panel * 7 + code], x_raw) - y_zero)
                / raw_per_unit * scale
                for code in range(7)
            ]
            center = 0.5 * (min(values) + max(values))
            uncertainty = 0.5 * (max(values) - min(values))
            row.extend((center, uncertainty))
        writer.writerow(row)
    return stream.getvalue()


def write_output(args: argparse.Namespace) -> Path:
    pdf = args.pdf.resolve()
    case = args.case_manifest.resolve()
    output = args.output_dir.resolve()
    if digest(pdf) != PDF_SHA256:
        raise DigitizeError("accepted-manuscript PDF SHA-256 mismatch")
    if case_id(case) != CASE_ID:
        raise DigitizeError(f"case manifest is not {CASE_ID}")
    if output.exists():
        raise DigitizeError(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    ))
    try:
        svg = temporary / "figure6.svg"
        render_page(pdf, svg)
        curves = figure_paths(svg)
        profile = temporary / "profiles.csv"
        profile.write_text(csv_text(curves), encoding="utf-8")
        manifest = temporary / "reference.hall-reference"
        manifest.write_text(
            f"""[reference]
hall_reference_version = 2
comparison_scope = digitized_profile_screening
case_id = {CASE_ID}
case_variant = charoy-2019-figure6-seven-code-envelope
case_manifest_sha256 = {digest(case)}
profile_data_file = profiles.csv
profile_data_sha256 = {digest(profile)}
profile_axis = x
coordinate_column = coordinate_m
coordinate_absolute_tolerance = 1e-12
reference_start_time_s = 1.6e-5
reference_end_time_s = 2.0e-5
provenance = Digitized seven-code envelope from accepted-manuscript Figure 6, PDF SHA-256 {PDF_SHA256}
citation = T. Charoy et al., Plasma Sources Science and Technology 28 (2019) 105010, doi:10.1088/1361-6595/ab46c5
retrieved = 2026-07-29
license = Accepted manuscript available under CC BY-NC-ND 3.0 after embargo

[profile.axial_field]
simulation_source = field
simulation_column = electric_x
reference_column = electric_x_v_m
reference_uncertainty_column = electric_x_uncertainty_v_m
relative_tolerance = 0.05
uncertainty_multiplier = 1

[profile.ion_density]
simulation_source = species
simulation_species = ions
simulation_column = number_density
reference_column = ion_density_m3
reference_uncertainty_column = ion_density_uncertainty_m3
relative_tolerance = 0.05
uncertainty_multiplier = 1

[profile.electron_temperature]
simulation_source = species
simulation_species = electrons
simulation_column = temperature_ev
reference_column = electron_temperature_ev
reference_uncertainty_column = electron_temperature_uncertainty_ev
relative_tolerance = 0.05
uncertainty_multiplier = 1
""",
            encoding="utf-8",
        )
        (temporary / "digitization.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "comparison_scope": "digitized_profile_screening",
                    "physics_claim": "none",
                    "source_pdf": str(pdf),
                    "source_pdf_sha256": PDF_SHA256,
                    "source_page": 20,
                    "figure": 6,
                    "curve_count": len(curves),
                    "code_curve_count_per_observable": 7,
                    "profile_rows": 124,
                    "envelope_method": "midpoint_and_half_range",
                    "sampling": "linear interpolation on 0.2 mm centers",
                    "profile_sha256": digest(profile),
                    "reference_manifest_sha256": digest(manifest),
                    "limitations": [
                        "digitized from a publication figure",
                        "published native tables are preferred for validation",
                        "inter-code spread is not measurement uncertainty",
                    ],
                },
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        svg.unlink()
        os.replace(temporary, output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a pinned, profile-only Figure 6 screening reference"
    )
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--case-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    try:
        output = write_output(parse_args())
    except (DigitizeError, ValueError) as error:
        print(f"Hall Figure 6 digitization error: {error}", file=sys.stderr)
        return 2
    print(f"Created Figure 6 screening reference: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
