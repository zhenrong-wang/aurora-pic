#!/usr/bin/env python3
"""Convert a user-supplied LXCat/BOLSIG+ file into an AuroraPIC gas package.

The converter runs locally and does not download or redistribute source data.
It supports the collision processes currently implemented by AuroraPIC:
ELASTIC, EXCITATION, IONIZATION, and ATTACHMENT.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Iterable, Sequence


ELECTRON_MASS_KG = 9.1093837139e-31
EV_TO_J = 1.602176634e-19
KNOWN_TYPES = {
    "ELASTIC",
    "EFFECTIVE",
    "MOMENTUM",
    "IONIZATION",
    "ATTACHMENT",
    "EXCITATION",
    "ROTATION",
}
SUPPORTED_TYPES = {
    "ELASTIC",
    "IONIZATION",
    "ATTACHMENT",
    "EXCITATION",
}
SEPARATOR = re.compile(r"^-{5,}\s*$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9_.:-]+$")


class ImportFailure(RuntimeError):
    """Raised when source data cannot be converted without ambiguity."""


@dataclass(frozen=True)
class Process:
    source_type: str
    target: str
    parameter_line: str
    threshold_ev: float
    points: tuple[tuple[float, float], ...]


def parse_number(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise ImportFailure(f"{context}: invalid number {value!r}") from exc
    if not math.isfinite(result):
        raise ImportFailure(f"{context}: non-finite number {value!r}")
    return result


def validate_manifest_value(name: str, value: str) -> None:
    if not value.strip():
        raise ImportFailure(f"{name} must not be empty")
    if value != value.strip():
        raise ImportFailure(f"{name} must not have leading or trailing whitespace")
    if any(character in value for character in ("\n", "\r", "#", ";")):
        raise ImportFailure(
            f"{name} contains a newline or manifest comment delimiter"
        )


def base_target(value: str) -> str:
    return value.split("<->", 1)[0].split("->", 1)[0].strip()


def parse_lxcat(path: Path) -> list[Process]:
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ImportFailure(f"cannot read source file {path}: {exc}") from exc

    result: list[Process] = []
    index = 0
    while index < len(lines):
        source_type = lines[index].strip().upper()
        if source_type not in KNOWN_TYPES:
            index += 1
            continue
        block_line = index + 1
        if index + 1 >= len(lines) or not lines[index + 1].strip():
            raise ImportFailure(
                f"{path}:{block_line}: {source_type} is missing its target"
            )
        target = base_target(lines[index + 1])
        if not target:
            raise ImportFailure(
                f"{path}:{block_line + 1}: {source_type} has an empty target"
            )

        parameter_index = index + 2
        parameter_line = ""
        if source_type != "ATTACHMENT":
            if parameter_index >= len(lines):
                raise ImportFailure(
                    f"{path}:{block_line}: {source_type} is missing parameters"
                )
            parameter_line = lines[parameter_index].strip()
            if not parameter_line:
                raise ImportFailure(
                    f"{path}:{parameter_index + 1}: empty process parameters"
                )
            search_index = parameter_index + 1
        else:
            search_index = parameter_index

        while search_index < len(lines) and not SEPARATOR.fullmatch(
            lines[search_index].strip()
        ):
            if lines[search_index].strip().upper() in KNOWN_TYPES:
                raise ImportFailure(
                    f"{path}:{block_line}: {source_type} has no data delimiter"
                )
            search_index += 1
        if search_index >= len(lines):
            raise ImportFailure(
                f"{path}:{block_line}: {source_type} has no data delimiter"
            )

        table_start = search_index + 1
        table_end = table_start
        while table_end < len(lines) and not SEPARATOR.fullmatch(
            lines[table_end].strip()
        ):
            table_end += 1
        if table_end >= len(lines):
            raise ImportFailure(
                f"{path}:{block_line}: {source_type} table is not terminated"
            )

        points: list[tuple[float, float]] = []
        for line_index in range(table_start, table_end):
            row = lines[line_index].strip()
            if not row:
                continue
            columns = row.split()
            if len(columns) != 2:
                raise ImportFailure(
                    f"{path}:{line_index + 1}: expected energy and cross section"
                )
            energy = parse_number(
                columns[0], f"{path}:{line_index + 1} energy"
            )
            cross_section = parse_number(
                columns[1], f"{path}:{line_index + 1} cross section"
            )
            if energy < 0.0 or cross_section < 0.0:
                raise ImportFailure(
                    f"{path}:{line_index + 1}: table values must be non-negative"
                )
            if points and energy <= points[-1][0]:
                raise ImportFailure(
                    f"{path}:{line_index + 1}: energies must increase strictly"
                )
            points.append((energy, cross_section))
        if len(points) < 2:
            raise ImportFailure(
                f"{path}:{block_line}: {source_type} requires at least two points"
            )

        threshold_ev = 0.0
        if source_type in {"EXCITATION", "IONIZATION"}:
            threshold_ev = parse_number(
                parameter_line.split()[0],
                f"{path}:{parameter_index + 1} threshold",
            )
            if threshold_ev <= 0.0:
                raise ImportFailure(
                    f"{path}:{parameter_index + 1}: inelastic threshold must be positive"
                )
            if any(
                energy < threshold_ev and cross_section > 0.0
                for energy, cross_section in points
            ):
                raise ImportFailure(
                    f"{path}:{block_line}: {source_type} has nonzero data below threshold"
                )

        result.append(
            Process(
                source_type=source_type,
                target=target,
                parameter_line=parameter_line,
                threshold_ev=threshold_ev,
                points=tuple(points),
            )
        )
        index = table_end + 1

    if not result:
        raise ImportFailure(f"{path}: no LXCat/BOLSIG+ collision blocks found")
    return result


def cross_section_at(process: Process, energy_ev: float) -> float:
    if energy_ev < process.threshold_ev:
        return 0.0
    points = process.points
    if energy_ev <= points[0][0]:
        return points[0][1]
    if energy_ev >= points[-1][0]:
        return points[-1][1]
    low = 0
    high = len(points) - 1
    while high - low > 1:
        middle = (low + high) // 2
        if points[middle][0] <= energy_ev:
            low = middle
        else:
            high = middle
    e0, s0 = points[low]
    e1, s1 = points[high]
    return s0 + (energy_ev - e0) * (s1 - s0) / (e1 - e0)


def maximum_rate_coefficient(
    processes: Sequence[Process],
    projectile_mass_kg: float,
    maximum_energy_ev: float,
) -> tuple[float, float]:
    breakpoints = {0.0, maximum_energy_ev}
    for process in processes:
        breakpoints.add(min(maximum_energy_ev, process.threshold_ev))
        breakpoints.update(
            energy
            for energy, _ in process.points
            if 0.0 <= energy <= maximum_energy_ev
        )
    ordered = sorted(breakpoints)
    candidates = set(ordered)
    for first, second in zip(ordered, ordered[1:]):
        if not first < second:
            continue
        sigma_first = sum(
            cross_section_at(process, first) for process in processes
        )
        sigma_second = sum(
            cross_section_at(process, second) for process in processes
        )
        slope = (sigma_second - sigma_first) / (second - first)
        intercept = sigma_first - slope * first
        if slope != 0.0:
            critical = -intercept / (3.0 * slope)
            if first < critical < second:
                candidates.add(critical)

    best_energy = 0.0
    best_value = 0.0
    for energy_ev in candidates:
        total_cross_section = sum(
            cross_section_at(process, energy_ev)
            for process in processes
        )
        speed = math.sqrt(
            2.0 * energy_ev * EV_TO_J / projectile_mass_kg
        )
        value = total_cross_section * speed
        if value > best_value:
            best_energy = energy_ev
            best_value = value
    return best_energy, best_value


def channel_names(processes: Sequence[Process]) -> list[str]:
    totals: dict[str, int] = {}
    names: list[str] = []
    for process in processes:
        base = process.source_type.lower()
        totals[base] = totals.get(base, 0) + 1
        names.append(f"{base}_{totals[base]:03d}")
    return names


def write_package(
    output_dir: Path,
    args: argparse.Namespace,
    processes: Sequence[Process],
) -> None:
    if output_dir.exists():
        raise ImportFailure(
            f"output directory already exists: {output_dir}"
        )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.",
            dir=output_dir.parent,
        )
    )
    try:
        names = channel_names(processes)
        for name, process in zip(names, processes):
            table_path = temporary / f"{name}.dat"
            with table_path.open("w", encoding="utf-8") as handle:
                handle.write("# User-supplied LXCat/BOLSIG+ conversion\n")
                handle.write("# energy_eV  cross_section_m2\n")
                for energy, cross_section in process.points:
                    handle.write(
                        f"{energy:.17g} {cross_section:.17g}\n"
                    )

        manifest_path = temporary / f"{args.gas}.gas"
        with manifest_path.open("w", encoding="utf-8") as handle:
            handle.write("gas_data_version = 2\n")
            handle.write("units = si\n")
            handle.write(f"gas = {args.gas}\n")
            handle.write(f"neutral_mass = {args.neutral_mass:.17g}\n")
            handle.write(f"dataset_id = {args.dataset_id}\n")
            handle.write(f"dataset_version = {args.dataset_version}\n")
            handle.write(f"data_provenance = {args.provenance}\n")
            handle.write(f"citation = {args.citation}\n")
            handle.write(f"retrieved = {args.retrieved}\n")
            handle.write(f"license = {args.license}\n")
            for name, process in zip(names, processes):
                handle.write(f"\n[collision.{name}]\n")
                handle.write(
                    f"type = {process.source_type.lower()}\n"
                )
                handle.write(f"cross_section_file = {name}.dat\n")
                handle.write(
                    f"energy_scale = {EV_TO_J:.17g}\n"
                )
                if process.source_type == "ELASTIC":
                    handle.write("angular_model = isotropic\n")
                if process.threshold_ev > 0.0:
                    handle.write(
                        f"threshold_energy = {process.threshold_ev:.17g}\n"
                    )

        maximum_energy = (
            args.max_energy_ev
            if args.max_energy_ev is not None
            else max(process.points[-1][0] for process in processes)
        )
        peak_energy, rate_coefficient = maximum_rate_coefficient(
            processes, args.projectile_mass, maximum_energy
        )
        audit = {
            "format": 1,
            "source_file": str(args.source.resolve()),
            "source_sha256": args.source_sha256,
            "gas": args.gas,
            "dataset_id": args.dataset_id,
            "dataset_version": args.dataset_version,
            "retrieved": args.retrieved,
            "units": {
                "source_energy": "eV",
                "source_cross_section": "m^2",
                "manifest": "si",
            },
            "process_count": len(processes),
            "elastic_mass_ratio_checks": args.mass_ratio_checks,
            "processes": [
                {
                    "channel": name,
                    "type": process.source_type.lower(),
                    "angular_model": "isotropic",
                    "target": process.target,
                    "threshold_ev": process.threshold_ev,
                    "points": len(process.points),
                    "minimum_energy_ev": process.points[0][0],
                    "maximum_energy_ev": process.points[-1][0],
                    "maximum_cross_section_m2": max(
                        value for _, value in process.points
                    ),
                }
                for name, process in zip(names, processes)
            ],
            "rate_envelope": {
                "maximum_energy_ev": maximum_energy,
                "peak_energy_ev": peak_energy,
                "maximum_sigma_v_m3_s": rate_coefficient,
                "safety_factor": args.safety_factor,
                "neutral_density_m3": args.neutral_density,
                "recommended_max_frequency_s": (
                    None
                    if args.neutral_density is None
                    else rate_coefficient
                    * args.neutral_density
                    * args.safety_factor
                ),
            },
        }
        (temporary / "audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.rename(output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gas", required=True)
    parser.add_argument("--neutral-mass", required=True, type=float)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--dataset-version", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--citation", required=True)
    parser.add_argument("--retrieved", default=date.today().isoformat())
    parser.add_argument("--license", required=True)
    parser.add_argument(
        "--projectile-mass",
        type=float,
        default=ELECTRON_MASS_KG,
        help="Projectile mass in kg (default: electron mass)",
    )
    parser.add_argument("--neutral-density", type=float)
    parser.add_argument("--safety-factor", type=float, default=1.2)
    parser.add_argument("--max-energy-ev", type=float)
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Allow a set missing elastic, excitation, or ionization",
    )
    parser.add_argument(
        "--skip-unsupported",
        action="store_true",
        help="Explicitly omit unsupported source process types",
    )
    parser.add_argument(
        "--ignore-mass-ratio",
        action="store_true",
        help="Accept elastic electron/neutral mass-ratio disagreement",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    try:
        if not IDENTIFIER.fullmatch(args.gas):
            raise ImportFailure(
                "--gas may contain only letters, digits, '_', '-', '.', ':'"
            )
        if not IDENTIFIER.fullmatch(args.dataset_id):
            raise ImportFailure(
                "--dataset-id contains unsupported characters"
            )
        if args.neutral_mass <= 0.0 or not math.isfinite(args.neutral_mass):
            raise ImportFailure("--neutral-mass must be positive and finite")
        if args.projectile_mass <= 0.0 or not math.isfinite(
            args.projectile_mass
        ):
            raise ImportFailure(
                "--projectile-mass must be positive and finite"
            )
        if args.neutral_density is not None and (
            args.neutral_density <= 0.0
            or not math.isfinite(args.neutral_density)
        ):
            raise ImportFailure(
                "--neutral-density must be positive and finite"
            )
        if args.safety_factor < 1.0 or not math.isfinite(args.safety_factor):
            raise ImportFailure("--safety-factor must be finite and at least 1")
        if args.max_energy_ev is not None and (
            args.max_energy_ev <= 0.0
            or not math.isfinite(args.max_energy_ev)
        ):
            raise ImportFailure("--max-energy-ev must be positive and finite")
        try:
            date.fromisoformat(args.retrieved)
        except ValueError as exc:
            raise ImportFailure(
                "--retrieved must be a valid YYYY-MM-DD date"
            ) from exc
        for option, value in (
            ("--dataset-version", args.dataset_version),
            ("--provenance", args.provenance),
            ("--citation", args.citation),
            ("--license", args.license),
        ):
            validate_manifest_value(option, value)

        import hashlib

        args.source_sha256 = hashlib.sha256(
            args.source.read_bytes()
        ).hexdigest()
        all_processes = parse_lxcat(args.source)
        target_processes = [
            process
            for process in all_processes
            if process.target == args.gas
        ]
        if not target_processes:
            available = sorted({process.target for process in all_processes})
            raise ImportFailure(
                f"no processes target {args.gas!r}; available: {available}"
            )
        unsupported = sorted(
            {
                process.source_type
                for process in target_processes
                if process.source_type not in SUPPORTED_TYPES
            }
        )
        if unsupported and not args.skip_unsupported:
            raise ImportFailure(
                "unsupported process types would make the conversion "
                f"incomplete: {unsupported}; use --skip-unsupported explicitly"
            )
        processes = [
            process
            for process in target_processes
            if process.source_type in SUPPORTED_TYPES
        ]
        present = {process.source_type for process in processes}
        required = {"ELASTIC", "EXCITATION", "IONIZATION"}
        missing = sorted(required - present)
        if missing and not args.allow_partial:
            raise ImportFailure(
                f"incomplete electron collision set; missing {missing}; "
                "use --allow-partial explicitly"
            )
        expected_mass_ratio = (
            args.projectile_mass / args.neutral_mass
        )
        args.mass_ratio_checks = []
        for process in processes:
            if process.source_type != "ELASTIC":
                continue
            source_mass_ratio = parse_number(
                process.parameter_line.split()[0],
                f"{args.source} elastic mass ratio",
            )
            if source_mass_ratio <= 0.0:
                raise ImportFailure(
                    "elastic electron/neutral mass ratio must be positive"
                )
            relative_error = abs(
                source_mass_ratio - expected_mass_ratio
            ) / expected_mass_ratio
            args.mass_ratio_checks.append(
                {
                    "source": source_mass_ratio,
                    "expected": expected_mass_ratio,
                    "relative_error": relative_error,
                }
            )
            if relative_error > 0.05 and not args.ignore_mass_ratio:
                raise ImportFailure(
                    "elastic mass ratio disagrees with projectile and "
                    "neutral masses by more than 5%; use "
                    "--ignore-mass-ratio explicitly"
                )
        write_package(args.output_dir, args, processes)
        print(
            f"wrote {len(processes)} channels to {args.output_dir}",
            flush=True,
        )
        return 0
    except (ImportFailure, OSError) as exc:
        print(f"gas import failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
