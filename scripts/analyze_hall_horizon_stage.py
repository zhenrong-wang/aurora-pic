#!/usr/bin/env python3
"""Audit one Hall horizon stage and its published-profile trend."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


ELEMENTARY_CHARGE_C = 1.602176634e-19


class HorizonAnalysisError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise HorizonAnalysisError(f"cannot hash {path}: {error}") from error


def load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise HorizonAnalysisError(f"cannot read {label}: {error}") from error
    if not isinstance(value, dict):
        raise HorizonAnalysisError(f"{label} must contain an object")
    return value


def rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or ())
            if not required <= fields:
                raise HorizonAnalysisError(
                    f"{path.name} is missing {sorted(required - fields)}"
                )
            result = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise HorizonAnalysisError(f"cannot read {path}: {error}") from error
    if not result:
        raise HorizonAnalysisError(f"{path.name} has no rows")
    return result


def finite(value: object, context: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise HorizonAnalysisError(f"{context} is not numeric") from error
    if not math.isfinite(result):
        raise HorizonAnalysisError(f"{context} is not finite")
    return result


def integer(value: object, context: str) -> int:
    number = finite(value, context)
    result = int(number)
    if number != result or result < 0:
        raise HorizonAnalysisError(
            f"{context} is not a non-negative integer"
        )
    return result


def checkpoint_identity(path: Path) -> tuple[int, float]:
    try:
        with path.open(encoding="utf-8") as stream:
            magic = stream.readline().strip()
            dimension = stream.readline().split()
            units = stream.readline().split()
            step = stream.readline().split()
            time = stream.readline().split()
    except (OSError, UnicodeError) as error:
        raise HorizonAnalysisError(
            f"cannot inspect checkpoint {path}: {error}"
        ) from error
    if (
        not magic.startswith("AuroraPIC-checkpoint-v")
        or dimension != ["dimension", "2"]
        or not units
        or units[0] != "units"
        or len(step) != 2
        or step[0] != "step"
        or len(time) != 2
        or time[0] != "time"
    ):
        raise HorizonAnalysisError("final checkpoint header is invalid")
    return integer(step[1], "checkpoint step"), finite(
        time[1], "checkpoint time"
    )


def validate_hash(path: Path, expected: object, label: str) -> str:
    actual = sha256(path)
    if actual != expected:
        raise HorizonAnalysisError(f"{label} SHA-256 mismatch")
    return actual


def comparison_summary(
    report: dict[str, object],
    expected_output: Path,
    label: str,
) -> dict[str, dict[str, float | int]]:
    if (
        report.get("schema_version") != 1
        or report.get("comparison_scope")
            != "digitized_profile_screening"
        or report.get("physics_claim") != "none"
    ):
        raise HorizonAnalysisError(
            f"{label} is not a profile-screening comparison"
        )
    simulation = report.get("simulation")
    if not isinstance(simulation, dict):
        raise HorizonAnalysisError(f"{label} simulation identity is missing")
    if Path(str(simulation.get("output_dir", ""))).resolve() != expected_output:
        raise HorizonAnalysisError(f"{label} output directory mismatch")
    summaries = report.get("profile_summary")
    if not isinstance(summaries, list) or not summaries:
        raise HorizonAnalysisError(f"{label} profile summary is missing")
    result: dict[str, dict[str, float | int]] = {}
    for item in summaries:
        if not isinstance(item, dict) or not isinstance(item.get("name"), str):
            raise HorizonAnalysisError(f"{label} profile summary is invalid")
        name = str(item["name"])
        if name in result:
            raise HorizonAnalysisError(f"{label} duplicates {name}")
        result[name] = {
            "relative_l2": finite(
                item.get("relative_l2"), f"{label} {name} relative_l2"
            ),
            "relative_linf": finite(
                item.get("relative_linf"), f"{label} {name} relative_linf"
            ),
            "passed_points": integer(
                item.get("passed_points"), f"{label} {name} passed_points"
            ),
            "points": integer(
                item.get("points"), f"{label} {name} points"
            ),
        }
    return result


def controller_metrics(path: Path) -> dict[str, float | int]:
    data = rows(
        path,
        {
            "macro_particles_created", "represented_particles_created",
            "unserved_reverse_charge", "charge_balance_residual",
        },
    )
    final = data[-1]
    created = integer(
        final["macro_particles_created"],
        f"{path.name} macro_particles_created",
    )
    represented = finite(
        final["represented_particles_created"],
        f"{path.name} represented_particles_created",
    )
    if created <= 0 or represented <= 0.0:
        raise HorizonAnalysisError(
            f"{path.name} cannot establish macro charge"
        )
    macro_charge = represented / created * ELEMENTARY_CHARGE_C
    maximum_unserved = max(
        abs(finite(row["unserved_reverse_charge"], path.name))
        for row in data
    )
    maximum_residual = max(
        abs(finite(row["charge_balance_residual"], path.name))
        for row in data
    )
    return {
        "samples": len(data),
        "macro_particles_created": created,
        "macro_charge_c": macro_charge,
        "maximum_unserved_reverse_charge_c": maximum_unserved,
        "maximum_unserved_reverse_macroparticles":
            maximum_unserved / macro_charge,
        "unserved_fraction_of_cumulative_emission":
            maximum_unserved / macro_charge / created,
        "maximum_charge_balance_residual_c": maximum_residual,
    }


def atomic_json(path: Path, value: dict[str, object]) -> None:
    if path.exists():
        raise HorizonAnalysisError(f"refusing to overwrite report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def analyze(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.horizon_manifest.resolve()
    manifest = load_json(manifest_path, "horizon manifest")
    if (
        manifest.get("hall_horizon_stage_version") != 1
        or manifest.get("physics_claim") != "none"
    ):
        raise HorizonAnalysisError("unsupported horizon-stage manifest")
    case = Path(str(manifest.get("case_manifest", ""))).resolve()
    prior_deck = Path(str(manifest.get("prior_deck", ""))).resolve()
    prior_output = Path(str(manifest.get("prior_output", ""))).resolve()
    restart = Path(str(manifest.get("restart_checkpoint", ""))).resolve()
    runtime = Path(str(manifest.get("runtime_config", ""))).resolve()
    output = Path(str(manifest.get("result_dir", ""))).resolve()
    validate_hash(case, manifest.get("case_manifest_sha256"), "case manifest")
    validate_hash(prior_deck, manifest.get("prior_deck_sha256"), "prior deck")
    validate_hash(
        prior_output / "scalars.csv",
        manifest.get("prior_scalars_sha256"),
        "prior scalars",
    )
    validate_hash(
        restart,
        manifest.get("restart_checkpoint_sha256"),
        "restart checkpoint",
    )
    runtime_hash = validate_hash(
        runtime, manifest.get("runtime_config_sha256"), "runtime config"
    )

    scalar_required = {
        "step", "time", "total_energy", "charge_l1", "live_particles",
        "live_particles_electrons", "live_particles_ions",
    }
    current_scalars = rows(output / "scalars.csv", scalar_required)
    prior_scalars = rows(prior_output / "scalars.csv", scalar_required)
    prior_step = integer(manifest.get("prior_step"), "manifest prior_step")
    target_step = integer(manifest.get("target_step"), "manifest target_step")
    target_time = finite(
        manifest.get("target_time_s"), "manifest target_time_s"
    )
    first = current_scalars[0]
    final = current_scalars[-1]
    if (
        integer(first["step"], "first scalar step") != prior_step
        or integer(final["step"], "final scalar step") != target_step
        or not math.isclose(
            finite(final["time"], "final scalar time"),
            target_time,
            rel_tol=1e-10,
        )
        or any(
            first[key] != prior_scalars[-1][key]
            for key in scalar_required
        )
    ):
        raise HorizonAnalysisError(
            "scalar restart/final identity does not match the stage"
        )
    steps = [integer(row["step"], "scalar step") for row in current_scalars]
    times = [finite(row["time"], "scalar time") for row in current_scalars]
    if steps != sorted(set(steps)) or times != sorted(set(times)):
        raise HorizonAnalysisError("scalar steps/times are not increasing")
    for row in current_scalars:
        finite(row["total_energy"], "scalar total_energy")
        finite(row["charge_l1"], "scalar charge_l1")
    capacity = integer(
        manifest.get("max_particles_per_species"),
        "manifest particle capacity",
    )
    maximum_by_species = {
        species: max(
            integer(row[f"live_particles_{species}"], species)
            for row in current_scalars
        )
        for species in ("electrons", "ions")
    }
    if any(value > capacity for value in maximum_by_species.values()):
        raise HorizonAnalysisError("particle population exceeded capacity")

    final_checkpoint = output / f"checkpoint_{target_step}.apc"
    final_checkpoint_step, final_checkpoint_time = checkpoint_identity(
        final_checkpoint
    )
    if (
        final_checkpoint_step != target_step
        or not math.isclose(
            final_checkpoint_time, target_time, rel_tol=1e-10
        )
    ):
        raise HorizonAnalysisError("final checkpoint identity mismatch")

    current_report_path = args.current_comparison.resolve()
    prior_report_path = args.prior_comparison.resolve()
    current_report = load_json(current_report_path, "current comparison")
    prior_report = load_json(prior_report_path, "prior comparison")
    simulation = current_report.get("simulation")
    if not isinstance(simulation, dict):
        raise HorizonAnalysisError(
            "current comparison simulation identity is missing"
        )
    if (
        simulation.get("runtime_config_sha256") != runtime_hash
        or simulation.get("case_manifest_sha256") != sha256(case)
    ):
        raise HorizonAnalysisError(
            "current comparison is not bound to the horizon stage"
        )
    current_summary = comparison_summary(
        current_report, output, "current comparison"
    )
    prior_summary = comparison_summary(
        prior_report, prior_output, "prior comparison"
    )
    if set(current_summary) != set(prior_summary):
        raise HorizonAnalysisError("comparison observables changed")
    trend: dict[str, object] = {}
    trend_passed = True
    for name in sorted(current_summary):
        before = prior_summary[name]
        after = current_summary[name]
        l2_ratio = (
            float(after["relative_l2"]) / float(before["relative_l2"])
            if float(before["relative_l2"]) > 0.0 else None
        )
        linf_ratio = (
            float(after["relative_linf"]) / float(before["relative_linf"])
            if float(before["relative_linf"]) > 0.0 else None
        )
        improved = (
            l2_ratio is not None
            and linf_ratio is not None
            and l2_ratio <= 1.0
            and linf_ratio <= 1.0
        )
        trend_passed = trend_passed and improved
        trend[name] = {
            "prior": before,
            "current": after,
            "relative_l2_change_ratio": l2_ratio,
            "relative_linf_change_ratio": linf_ratio,
            "improved": improved,
        }

    prior_controller = controller_metrics(
        prior_output / "current_source.csv"
    )
    current_controller = controller_metrics(
        output / "current_source.csv"
    )
    controller_review = (
        float(current_controller[
            "maximum_unserved_reverse_macroparticles"
        ])
        > float(prior_controller[
            "maximum_unserved_reverse_macroparticles"
        ]) + 1e-9
    )
    potential_rows = rows(
        output / "potential_reference.csv",
        {"target", "corrected_line_mean"},
    )
    maximum_potential_error = max(
        abs(
            finite(row["corrected_line_mean"], "corrected potential")
            - finite(row["target"], "target potential")
        )
        for row in potential_rows
    )
    ready = trend_passed and not controller_review
    return {
        "schema_version": 1,
        "case_id": manifest.get("case_id"),
        "physics_claim": "none",
        "numerical_integrity_passed": True,
        "profile_trend_passed": trend_passed,
        "controller_review_required": controller_review,
        "ready_for_next_stage": ready,
        "horizon_manifest": str(manifest_path),
        "horizon_manifest_sha256": sha256(manifest_path),
        "current_comparison": str(current_report_path),
        "current_comparison_sha256": sha256(current_report_path),
        "prior_comparison": str(prior_report_path),
        "prior_comparison_sha256": sha256(prior_report_path),
        "runtime_config_sha256": runtime_hash,
        "final_checkpoint": str(final_checkpoint),
        "final_checkpoint_sha256": sha256(final_checkpoint),
        "metrics": {
            "prior_step": prior_step,
            "target_step": target_step,
            "target_time_s": target_time,
            "scalar_samples": len(current_scalars),
            "final_live_particles": integer(
                final["live_particles"], "final live_particles"
            ),
            "maximum_live_particles_by_species": maximum_by_species,
            "max_particles_per_species": capacity,
            "maximum_potential_reference_error_v":
                maximum_potential_error,
        },
        "profile_trend": trend,
        "prior_controller": prior_controller,
        "current_controller": current_controller,
        "decision": (
            "Hold the next horizon extension for cathode-controller review."
            if controller_review
            else "The stage is eligible for the next guarded extension."
        ),
        "limitations": [
            "Profile data are digitized publication curves, not native tables.",
            "The 80-100 ns window does not match the published 16-20 us window.",
            "Improvement with time is necessary but does not establish validation.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit a Hall horizon stage and published-profile trend"
    )
    parser.add_argument("horizon_manifest", type=Path)
    parser.add_argument("--current-comparison", type=Path, required=True)
    parser.add_argument("--prior-comparison", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args)
        atomic_json(args.report, report)
    except (
        HorizonAnalysisError, KeyError, TypeError, ValueError
    ) as error:
        print(f"Hall horizon analysis error: {error}", file=sys.stderr)
        return 2
    if not report["ready_for_next_stage"]:
        print(
            "Hall horizon stage is not ready for extension",
            file=sys.stderr,
        )
        return 1
    print(f"Hall horizon stage is ready for extension: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
