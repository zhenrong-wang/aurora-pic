#!/usr/bin/env python3
"""Prepare a staged, checksum-locked Turner Case 1 sensitivity matrix."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_SENSITIVITY"
)


class SensitivityPreparationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SensitivityPreparationError(message)


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise SensitivityPreparationError(f"cannot read {path}: {error}") from error


def load_json(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SensitivityPreparationError(
            f"cannot read {description}: {error}"
        ) from error
    require(isinstance(value, dict), f"{description} must be a JSON object")
    return value


def replace_value(text: str, section: str | None, key: str, value: object) -> str:
    lines = text.splitlines()
    active: str | None = None
    matches: list[int] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            active = stripped[1:-1]
            continue
        if active == section and stripped.startswith(f"{key} ="):
            matches.append(index)
    require(len(matches) == 1,
            f"expected one {key!r} value in section {section!r}")
    lines[matches[0]] = f"{key} = {value}"
    return "\n".join(lines) + "\n"


def prepare(args: argparse.Namespace) -> Path:
    require(
        args.acknowledge_cost == ACKNOWLEDGEMENT,
        f"sensitivity preparation requires --acknowledge-cost {ACKNOWLEDGEMENT}",
    )
    manifest_path = args.ensemble_manifest.resolve()
    ensemble = load_json(manifest_path, "Turner ensemble manifest")
    require(
        ensemble.get("turner_ensemble_preparation_version") == 1
        and ensemble.get("case_id") == "turner-helium-ccp-2013-case-1",
        "unsupported Turner ensemble manifest",
    )
    runs = ensemble.get("runs")
    require(isinstance(runs, list), "Turner ensemble has no prepared runs")
    matches = [run for run in runs
               if isinstance(run, dict) and run.get("seed") == args.seed]
    require(len(matches) == 1, f"baseline seed {args.seed} is not unique")
    baseline_run = matches[0]
    ensemble_dir = manifest_path.parent
    baseline_config = ensemble_dir / str(baseline_run["runtime_config"])
    baseline_preflight = ensemble_dir / str(baseline_run["preflight_report"])
    require(
        sha256(baseline_config) == baseline_run.get("runtime_config_sha256")
        and sha256(baseline_preflight) == baseline_run.get("preflight_report_sha256"),
        "baseline prepared artifacts fail their ensemble checksums",
    )
    preflight = load_json(baseline_preflight, "baseline preflight")
    contract = preflight.get("contract")
    require(
        isinstance(contract, dict)
        and contract.get("seed") == args.seed
        and contract.get("steps") == 512000
        and contract.get("steps_per_rf_cycle") == 400
        and contract.get("rf_cycles") == 1280
        and contract.get("nodes") == 129
        and contract.get("particles_per_species") == 65536,
        "baseline is not the exact published Turner Case 1 deck",
    )
    executable = args.executable.resolve()
    require(executable.is_file(), f"solver executable does not exist: {executable}")

    destination = args.output_dir.resolve()
    require(not destination.exists(),
            f"refusing to overwrite sensitivity directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    ))
    baseline_text = baseline_config.read_text(encoding="utf-8")
    variants = [
        {
            "id": "particles_2x",
            "stage": 1,
            "changes": [
                (None, "max_particles_per_species", 524288),
                ("species.electrons", "particles", 131072),
                ("species.electrons", "weight", "130859375.00000001"),
                ("species.ions", "particles", 131072),
                ("species.ions", "weight", "130859375.00000001"),
            ],
            "steps": 512000,
            "particles_per_species": 131072,
            "nodes": 129,
            "hypothesis": "tests macro-particle noise and discrete-particle heating",
        },
        {
            "id": "timestep_2x",
            "stage": 1,
            "changes": [
                (None, "dt", "9.2182890855457226e-11"),
                (None, "steps", 1024000),
                (None, "output_interval", 800),
                (None, "checkpoint_interval", 25600),
                (None, "spatial_average_start_step", 998401),
                (None, "spatial_average_end_step", 1024000),
            ],
            "steps": 1024000,
            "particles_per_species": 65536,
            "nodes": 129,
            "hypothesis": "tests particle push, RF integration, and MCC time splitting",
        },
        {
            "id": "grid_2x_fixed_particles",
            "stage": 2,
            "changes": [(None, "nx", 257)],
            "steps": 512000,
            "particles_per_species": 65536,
            "nodes": 257,
            "hypothesis": "tests field-grid resolution at fixed total particle count",
        },
        {
            "id": "grid_2x_same_ppc",
            "stage": 2,
            "changes": [
                (None, "nx", 257),
                (None, "max_particles_per_species", 524288),
                ("species.electrons", "particles", 131072),
                ("species.electrons", "weight", "130859375.00000001"),
                ("species.ions", "particles", 131072),
                ("species.ions", "weight", "130859375.00000001"),
            ],
            "steps": 512000,
            "particles_per_species": 131072,
            "nodes": 257,
            "hypothesis": "tests joint grid and particle refinement at fixed particles per cell",
        },
    ]
    records: list[dict[str, object]] = []
    try:
        for variant in variants:
            identifier = str(variant["id"])
            config = baseline_text
            result_dir = destination / "results" / identifier
            config = replace_value(config, None, "output_dir", result_dir)
            for section, key, value in variant["changes"]:
                config = replace_value(config, section, key, value)
            variant_dir = temporary / identifier
            variant_dir.mkdir()
            deck = variant_dir / "turner_case1.cfg"
            deck.write_text(
                "# Diagnostic sensitivity variant; not the published numerical contract.\n"
                + config,
                encoding="utf-8",
            )
            validation = subprocess.run(
                [str(executable), "--validate-only", str(deck)],
                cwd=executable.parent.parent,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=60,
            )
            require(
                validation.returncode == 0,
                f"variant {identifier} failed solver validation: "
                + validation.stdout + validation.stderr,
            )
            initial_live = 2 * int(variant["particles_per_species"])
            records.append({
                "id": identifier,
                "stage": variant["stage"],
                "runtime_config": f"{identifier}/turner_case1.cfg",
                "runtime_config_sha256": sha256(deck),
                "result_dir": str(result_dir),
                "seed": args.seed,
                "steps": variant["steps"],
                "rf_cycles": 1280,
                "averaging_rf_cycles": 32,
                "nodes": variant["nodes"],
                "particles_per_species": variant["particles_per_species"],
                "initial_particle_updates": initial_live * int(variant["steps"]),
                "hypothesis": variant["hypothesis"],
                "published_acceptance_applicable": False,
                "launched": False,
                "completed": False,
            })

        output = {
            "turner_sensitivity_preparation_version": 1,
            "case_id": ensemble["case_id"],
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "baseline": {
                "ensemble_manifest": str(manifest_path),
                "ensemble_manifest_sha256": sha256(manifest_path),
                "seed": args.seed,
                "runtime_config": str(baseline_config),
                "runtime_config_sha256": sha256(baseline_config),
                "observed_integrated_density_bias_percent":
                    args.baseline_density_bias_percent,
            },
            "solver": {
                "executable": str(executable),
                "executable_sha256": sha256(executable),
            },
            "variants": records,
            "execution_policy": {
                "concurrent_runs_authorized": 1,
                "stage_1_first": ["particles_2x", "timestep_2x"],
                "stage_2_only_if_needed": [
                    "grid_2x_fixed_particles", "grid_2x_same_ppc"
                ],
                "host_policy": "serial_backend_low_cpu_and_idle_io_priority",
            },
            "predeclared_interpretation": {
                "primary_observable": "integrated_ion_density_bias_percent",
                "material_shift_percentage_points": 0.75,
                "practical_equivalence_percentage_points": 0.50,
                "ambiguous_interval_percentage_points": [0.50, 0.75],
                "material_shift_rationale":
                    "approximately three standard errors of the completed "
                    "three-seed baseline mean, rounded conservatively",
                "published_x_squared_role": "reported_but_not_used_as_variant_pass_fail",
            },
            "claim_boundary": {
                "physics_claim": "none_from_preparation_or_single_variant",
                "published_contract_changed": True,
                "sensitivity_result_is_turner_benchmark_pass": False,
            },
            "launched": False,
            "warnings": [
                "No simulation was launched.",
                "Run at most one variant at a time on a shared workstation.",
                "Do not tune or select variants by their published X-squared result.",
            ],
        }
        (temporary / "sensitivity.json").write_text(
            json.dumps(output, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination / "sensitivity.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ensemble_manifest", type=Path)
    parser.add_argument("--baseline-seed", dest="seed", type=int, default=13507)
    parser.add_argument("--baseline-density-bias-percent", type=float, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        output = prepare(parse_args())
    except (SensitivityPreparationError, OSError, ValueError,
            subprocess.SubprocessError) as error:
        print(f"Turner sensitivity preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Turner sensitivity matrix without launching it: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
