#!/usr/bin/env python3
"""Prepare a bounded Hall population/duration convergence campaign."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

from prepare_hall_campaign import CampaignError, load_manifest, prepare


ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_AN_OPT_IN_HALL_CONVERGENCE_PLAN"
)


class ConvergencePreparationError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def factors(section: object, key: str) -> list[float]:
    try:
        values = [
            float(item.strip())
            for item in section[key].split(",")
        ]
    except (KeyError, ValueError) as error:
        raise ConvergencePreparationError(
            f"{key} must contain comma-separated factors"
        ) from error
    if (
        len(values) != 3
        or values != sorted(values)
        or not math.isclose(values[1], 1.0)
        or any(not math.isfinite(value) or value <= 0.0 for value in values)
    ):
        raise ConvergencePreparationError(
            f"{key} must be three increasing positive factors centered on 1"
        )
    return values


def factor_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def replace_global(text: str, key: str, value: object) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=).*$")
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ConvergencePreparationError(
            f"generated deck must contain one global {key}"
        )
    return pattern.sub(rf"\g<1> {value}", text, count=1)


def replace_section(
    text: str, section: str, key: str, value: object
) -> str:
    section_pattern = re.compile(
        rf"(?ms)^\[{re.escape(section)}\]\n.*?(?=^\[|\Z)"
    )
    sections = list(section_pattern.finditer(text))
    if len(sections) != 1:
        raise ConvergencePreparationError(
            f"generated deck must contain one [{section}] section"
        )
    match = sections[0]
    block = match.group(0)
    key_pattern = re.compile(
        rf"(?m)^({re.escape(key)}\s*=).*$"
    )
    if len(list(key_pattern.finditer(block))) != 1:
        raise ConvergencePreparationError(
            f"generated deck must contain one [{section}] {key}"
        )
    replacement = key_pattern.sub(
        rf"\g<1> {value}", block, count=1
    )
    return text[:match.start()] + replacement + text[match.end():]


def stage_definition(
    axis: str,
    factor: float,
    base_ppc: int,
    base_steps: int,
) -> tuple[str, int, int]:
    if axis == "population":
        ppc_value = base_ppc * factor
        if not ppc_value.is_integer():
            raise ConvergencePreparationError(
                "population factor does not produce integral particles/cell"
            )
        return (
            f"population_{factor_label(factor)}",
            int(ppc_value),
            base_steps,
        )
    return (
        f"duration_{factor_label(factor)}",
        base_ppc,
        int(round(base_steps * factor)),
    )


def prepare_convergence(args: argparse.Namespace) -> Path:
    destination = args.output_dir.resolve()
    if destination.exists():
        raise ConvergencePreparationError(
            f"refusing to overwrite convergence directory: {destination}"
        )
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise ConvergencePreparationError(
            "convergence generation requires --acknowledge-cost "
            + ACKNOWLEDGEMENT
        )
    case_path = args.case_manifest.resolve()
    case = load_manifest(case_path)
    section_name = "convergence.workstation"
    if section_name not in case:
        raise ConvergencePreparationError(
            f"case manifest is missing [{section_name}]"
        )
    contract = case[section_name]
    if contract.get("convergence_contract_version") != "2":
        raise ConvergencePreparationError(
            "convergence_contract_version must be 2"
        )
    if contract.get("base_tier") != "workstation":
        raise ConvergencePreparationError(
            "only the workstation base tier is supported"
        )
    population_factors = factors(contract, "population_factors")
    duration_factors = factors(contract, "duration_factors")
    base = case["campaign.workstation"]
    cells_x = base.getint("cells_x")
    cells_y = base.getint("cells_y")
    base_ppc = base.getint("particles_per_cell_per_species")
    base_steps = base.getint("steps")
    base_capacity = base.getint("max_particles_per_species")
    averaging_fraction = contract.getfloat("averaging_fraction")
    samples = contract.getint("diagnostic_samples")
    if not 0.0 < averaging_fraction < 1.0 or samples < 2:
        raise ConvergencePreparationError(
            "averaging_fraction and diagnostic_samples are invalid"
        )

    definitions: list[tuple[str, str, float, int, int]] = []
    for factor in population_factors:
        name, ppc, steps = stage_definition(
            "population", factor, base_ppc, base_steps
        )
        definitions.append((name, "population", factor, ppc, steps))
    for factor in duration_factors:
        if math.isclose(factor, 1.0):
            continue
        name, ppc, steps = stage_definition(
            "duration", factor, base_ppc, base_steps
        )
        definitions.append((name, "duration", factor, ppc, steps))
    if len(definitions) > contract.getint("maximum_stages"):
        raise ConvergencePreparationError(
            "convergence plan exceeds maximum_stages"
        )

    aggregate_updates = sum(
        2 * cells_x * cells_y * ppc * steps
        for _, _, _, ppc, steps in definitions
    )
    maximum_updates = contract.getint(
        "maximum_aggregate_initial_particle_updates"
    )
    if aggregate_updates > maximum_updates:
        raise ConvergencePreparationError(
            "convergence plan exceeds its aggregate update limit"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp",
        dir=destination.parent,
    ))
    try:
        runs: list[dict[str, object]] = []
        for name, axis, factor, ppc, steps in definitions:
            deck = temporary / f"{name}.cfg"
            result_dir = destination / "results" / name
            prepare(argparse.Namespace(
                case_manifest=case_path,
                tier="workstation",
                output=deck,
                output_dir=str(result_dir),
                seed=args.seed,
                acknowledge_cost=(
                    "I_UNDERSTAND_THIS_IS_AN_OPT_IN_WORKSTATION_RUN"
                ),
            ))
            start = int(round(steps * (1.0 - averaging_fraction)))
            window = steps - start
            if window % (samples - 1) != 0:
                raise ConvergencePreparationError(
                    f"stage {name} cannot produce the requested sample cadence"
                )
            interval = window // (samples - 1)
            capacity = int(round(base_capacity * ppc / base_ppc))
            particles = cells_x * cells_y * ppc
            text = deck.read_text(encoding="utf-8")
            for key, value in (
                ("steps", steps),
                ("output_interval", interval),
                ("resolved_diagnostic_interval", interval),
                ("resolved_diagnostic_start_step", start),
                ("checkpoint_interval", max(1, steps // 2)),
                ("max_particles_per_species", capacity),
            ):
                text = replace_global(text, key, value)
            for species in ("species.electrons", "species.ions"):
                text = replace_section(text, species, "particles", particles)
            deck.write_text(text, encoding="utf-8")
            runs.append({
                "stage": name,
                "axis": axis,
                "factor": factor,
                "particles_per_cell_per_species": ppc,
                "initial_particles_per_species": particles,
                "steps": steps,
                "diagnostic_start_step": start,
                "diagnostic_interval": interval,
                "diagnostic_samples": samples,
                "max_particles_per_species": capacity,
                "initial_particle_updates": 2 * particles * steps,
                "runtime_config": deck.name,
                "runtime_config_sha256": sha256(deck),
                "result_dir": str(result_dir),
                "launched": False,
            })
        manifest = {
            "hall_convergence_version": 2,
            "case_id": case["global"]["case_id"],
            "tier": "workstation",
            "physics_claim": contract["physics_claim"],
            "case_manifest": str(case_path),
            "case_manifest_sha256": sha256(case_path),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "seed": args.seed,
            "fixed_cells_x": cells_x,
            "fixed_cells_y": cells_y,
            "aggregate_initial_particle_updates": aggregate_updates,
            "maximum_aggregate_initial_particle_updates": maximum_updates,
            "acceptance": {
                "relative_l2_tolerance":
                    contract.getfloat("relative_l2_tolerance"),
                "relative_linf_tolerance":
                    contract.getfloat("relative_linf_tolerance"),
                "maximum_fine_to_coarse_change_ratio":
                    contract.getfloat(
                        "maximum_fine_to_coarse_change_ratio"
                    ),
                "maximum_fine_to_baseline_reverse_charge_per_update_ratio":
                    contract.getfloat(
                        "maximum_fine_to_baseline_"
                        "reverse_charge_per_update_ratio"
                    ),
                "maximum_fine_to_baseline_reverse_impulse_ratio":
                    contract.getfloat(
                        "maximum_fine_to_baseline_reverse_impulse_ratio"
                    ),
                "maximum_fine_reverse_demand_macroparticles_per_update":
                    contract.getfloat(
                        "maximum_fine_reverse_demand_macroparticles_per_update"
                    ),
            },
            "launched": False,
            "runs": runs,
            "warnings": [
                "Campaign preparation did not launch a simulation.",
                "All five runs require the independent CLI large-run "
                "acknowledgement.",
                "Single-seed convergence is necessary but not sufficient "
                "for a physics claim.",
                "Controller acceptance uses represented charge, not raw "
                "macro-particle counts, because macro weight changes with "
                "population.",
            ],
        }
        (temporary / "convergence.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except (CampaignError, OSError, ValueError) as error:
        shutil.rmtree(temporary, ignore_errors=True)
        raise ConvergencePreparationError(str(error)) from error
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare, but never launch, the Hall workstation convergence plan"
        )
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        output = prepare_convergence(args)
    except ConvergencePreparationError as error:
        print(f"Hall convergence preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Hall convergence plan without launching it: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
