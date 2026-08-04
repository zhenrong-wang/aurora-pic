#!/usr/bin/env python3
"""Export AuroraPIC phase/region EEDFs to a solver-neutral interchange."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path} has no rows")
    return rows


def number(value: str, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {name}: {value}")
    return result


def write_rows(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--code", default="AuroraPIC")
    parser.add_argument("--code-version", required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--species", default="electrons")
    parser.add_argument("--provenance")
    args = parser.parse_args()
    source_histogram = read_rows(args.input_dir / "phase_eedf.csv")
    source_moments = read_rows(args.input_dir / "phase_eedf_moments.csv")
    if "energy_eV" not in source_histogram[0]:
        raise ValueError("interchange export requires SI phase EEDF energy_eV")

    grouped: defaultdict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    for row in source_histogram:
        grouped[(int(row["phase_bin"]), row["region"])].append(row)
    moments = {(int(row["phase_bin"]), row["region"]): row
               for row in source_moments}
    if set(grouped) != set(moments):
        raise ValueError("histogram and moment phase/region keys differ")

    distributions: list[dict] = []
    moment_output: list[dict] = []
    max_residual = 0.0
    for key in sorted(grouped):
        bins = sorted(grouped[key], key=lambda row: int(row["energy_bin"]))
        centers = [number(row["energy_eV"], "energy center") for row in bins]
        if len(centers) < 2:
            raise ValueError("at least two energy bins are required")
        widths = [centers[index + 1] - centers[index]
                  for index in range(len(centers) - 1)]
        width = widths[0]
        if width <= 0.0 or any(abs(item - width) > 1e-10 * width for item in widths):
            raise ValueError("AuroraPIC exporter requires uniform energy bins")
        moment = moments[key]
        represented = number(moment["represented_observations"],
                             "represented observations")
        overflow = number(moment["overflow_fraction"], "overflow fraction")
        if represented < 0.0 or not 0.0 <= overflow <= 1.0:
            raise ValueError("invalid observation count or overflow fraction")
        mass_sum = 0.0
        for index, row in enumerate(bins):
            count = number(row["represented_count"], "represented count")
            mass = count / represented if represented > 0.0 else 0.0
            mass_sum += mass
            distributions.append({
                "phase_bin": key[0], "phase_fraction": row["phase_fraction"],
                "region": key[1], "x_min_m": row["x_min"],
                "x_max_m": row["x_max"], "energy_bin": index,
                "energy_lower_eV": centers[index] - 0.5 * width,
                "energy_upper_eV": centers[index] + 0.5 * width,
                "probability_mass": mass, "represented_count": count,
            })
        max_residual = max(max_residual, abs(mass_sum + overflow - 1.0)
                           if represented > 0.0 else abs(mass_sum))
        moment_output.append({
            "phase_bin": key[0], "phase_fraction": moment["phase_fraction"],
            "region": key[1], "x_min_m": moment["x_min"],
            "x_max_m": moment["x_max"],
            "represented_observations": represented,
            "overflow_fraction": overflow,
            "mean_energy_eV": moment["mean_energy"],
            "energy_standard_deviation_eV": moment["energy_standard_deviation"],
            "drift_separated_temperature_eV": moment["drift_separated_temperature"],
        })
    if max_residual > 1e-10:
        raise ValueError(f"histogram normalization residual {max_residual} exceeds 1e-10")
    phases = sorted({int(row["phase_bin"]) for row in moment_output})
    if phases != list(range(len(phases))):
        raise ValueError("phase bins must be contiguous from zero")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "distributions.csv", (
        "phase_bin", "phase_fraction", "region", "x_min_m", "x_max_m",
        "energy_bin", "energy_lower_eV", "energy_upper_eV",
        "probability_mass", "represented_count"), distributions)
    write_rows(args.output_dir / "moments.csv", (
        "phase_bin", "phase_fraction", "region", "x_min_m", "x_max_m",
        "represented_observations", "overflow_fraction", "mean_energy_eV",
        "energy_standard_deviation_eV", "drift_separated_temperature_eV"),
        moment_output)
    manifest = {
        "schema": "aurorapic.phase-eedf-interchange", "schema_version": 1,
        "case_id": args.case_id,
        "code": {"name": args.code, "version": args.code_version},
        "species": args.species,
        "phase": {"bins": len(phases), "coordinate": "cycle_fraction",
                  "sampling_point": "bin_center", "periodic": True},
        "energy": {"unit": "eV", "histogram_value": "probability_mass",
                   "overflow_recorded_separately": True},
        "space": {"coordinate": "x", "unit": "m",
                  "regions": "closed-lower-open-upper"},
        "files": {"distributions": "distributions.csv", "moments": "moments.csv"},
        "maximum_normalization_residual": max_residual,
    }
    if args.provenance:
        manifest["provenance"] = args.provenance
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir),
                      "phase_region_pairs": len(moment_output),
                      "maximum_normalization_residual": max_residual},
                     indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
