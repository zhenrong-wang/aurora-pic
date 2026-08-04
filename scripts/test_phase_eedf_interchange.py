#!/usr/bin/env python3
"""Bounded synthetic regression for phase-EEDF export and comparison."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="aurorapic_eedf_", dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        source = work / "source"
        source.mkdir()
        (source / "phase_eedf.csv").write_text(
            "phase_bin,phase_fraction,region_id,region,x_min,x_max,energy_bin,energy_eV,represented_count,probability_density\n"
            "0,0.25,0,bulk,0,1,0,0.5,2,0.4\n"
            "0,0.25,0,bulk,0,1,1,1.5,3,0.6\n"
            "1,0.75,0,bulk,0,1,0,0.5,1,0.2\n"
            "1,0.75,0,bulk,0,1,1,1.5,4,0.8\n", encoding="utf-8")
        (source / "phase_eedf_moments.csv").write_text(
            "phase_bin,phase_fraction,region_id,region,x_min,x_max,macro_observations,represented_observations,overflow_fraction,mean_energy,energy_standard_deviation,mean_velocity_x,mean_velocity_y,mean_velocity_z,drift_separated_temperature\n"
            "0,0.25,0,bulk,0,1,5,5,0,1.1,0.49,0,0,0,0.7333333333333333\n"
            "1,0.75,0,bulk,0,1,5,5,0,1.3,0.4,0,0,0,0.8666666666666667\n",
            encoding="utf-8")
        exported = work / "exported"
        export = subprocess.run([
            sys.executable, str(ROOT / "scripts/export_phase_eedf.py"),
            str(source), str(exported), "--code-version", "test",
            "--case-id", "synthetic"], text=True, capture_output=True)
        require(export.returncode == 0, export.stderr)
        identical = subprocess.run([
            sys.executable, str(ROOT / "scripts/compare_phase_eedf.py"),
            str(exported), str(exported), "--tail-eV", "1.0",
            "--max-tv", "0", "--max-mean-energy-relative", "0",
            "--max-temperature-relative", "0"], text=True, capture_output=True)
        require(identical.returncode == 0, identical.stderr)
        report = json.loads(identical.stdout)
        require(report["acceptance"]["passes"] is True and
                report["maximum"]["total_variation"] == 0.0,
                "identical interchange comparison was not exact")

        rebinned = work / "rebinned"
        rebinned.mkdir()
        for name in ("manifest.json", "moments.csv"):
            (rebinned / name).write_bytes((exported / name).read_bytes())
        with (exported / "distributions.csv").open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            original_rows = list(reader)
            fields = reader.fieldnames
        split_rows = []
        for row in original_rows:
            low = float(row["energy_lower_eV"])
            high = float(row["energy_upper_eV"])
            mass = float(row["probability_mass"])
            midpoint = 0.5 * (low + high)
            for part, (part_low, part_high) in enumerate(((low, midpoint),
                                                          (midpoint, high))):
                item = dict(row)
                item["energy_bin"] = str(2 * int(row["energy_bin"]) + part)
                item["energy_lower_eV"] = str(part_low)
                item["energy_upper_eV"] = str(part_high)
                item["probability_mass"] = str(0.5 * mass)
                item["represented_count"] = str(
                    0.5 * float(row["represented_count"]))
                split_rows.append(item)
        with (rebinned / "distributions.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(split_rows)
        rebin_compare = subprocess.run([
            sys.executable, str(ROOT / "scripts/compare_phase_eedf.py"),
            str(exported), str(rebinned), "--max-tv", "0",
            "--max-mean-energy-relative", "0",
            "--max-temperature-relative", "0"], text=True, capture_output=True)
        require(rebin_compare.returncode == 0, rebin_compare.stderr)
        require(json.loads(rebin_compare.stdout)["maximum"]["total_variation"] == 0.0,
                "conservative energy rebin changed the distribution")

        changed = work / "changed"
        changed.mkdir()
        for name in ("manifest.json", "moments.csv", "distributions.csv"):
            (changed / name).write_bytes((exported / name).read_bytes())
        with (changed / "distributions.csv").open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
            fields = reader.fieldnames
        rows[0]["probability_mass"] = "0.3"
        rows[1]["probability_mass"] = "0.7"
        with (changed / "distributions.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        failed = subprocess.run([
            sys.executable, str(ROOT / "scripts/compare_phase_eedf.py"),
            str(exported), str(changed), "--max-tv", "0.05",
            "--max-mean-energy-relative", "0.01",
            "--max-temperature-relative", "0.01"], text=True, capture_output=True)
        require(failed.returncode == 1, "perturbed histogram passed acceptance")
        changed_report = json.loads(failed.stdout)
        require(abs(changed_report["maximum"]["total_variation"] - 0.1) < 1e-14,
                "unexpected perturbed histogram TV")
    print("phase EEDF interchange regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
