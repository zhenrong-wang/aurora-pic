#!/usr/bin/env python3
"""Run one preregistered AuroraPIC near-threshold field-work member."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, run_process,
    sha256,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PROMOTION_BAND_WORK_RUN"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256 = (
    "14742db5d59e61d950f625a711a6e59ff6af33d59b7554f430a5b6067ca9b019"
)


class PromotionBandRunError(RuntimeError):
    pass


def set_value(deck: str, key: str, value: object) -> str:
    prefix = f"{key} ="
    lines = deck.splitlines()
    matches = [index for index, line in enumerate(lines)
               if line.startswith(prefix)]
    if len(matches) != 1:
        raise PromotionBandRunError(f"expected one deck key: {key}")
    lines[matches[0]] = f"{key} = {value}"
    return "\n".join(lines) + "\n"


def make_deck(prior: str, output: Path, checkpoint: Path,
              rule: dict[str, object]) -> str:
    execution = rule["execution_contract"]
    diagnostic = rule["diagnostic_contract"]
    deck = prior
    for key, value in (
        ("steps", execution["aurorapic_end_step"]),
        ("output_dir", output),
        ("spatial_average_start_step",
         int(execution["aurorapic_start_step"]) + 1),
        ("spatial_average_end_step", execution["aurorapic_end_step"]),
        ("spatial_average_rf_cycles", execution["measurement_cycles"]),
        ("checkpoint_interval", execution["aurorapic_end_step"]),
        ("restart_path", checkpoint),
    ):
        deck = set_value(deck, key, value)
    if "phase_eedf_promotion_band_min" in deck:
        raise PromotionBandRunError("prior deck already contains promotion band")
    tail = f"phase_eedf_tail_threshold = {diagnostic['promotion_band_eV'][1]}\n"
    if deck.count(tail) != 1:
        raise PromotionBandRunError("tail-threshold deck contract differs")
    deck = deck.replace(
        tail,
        tail + "phase_eedf_promotion_band_min = " +
        str(diagnostic["promotion_band_eV"][0]) + "\n")
    return deck


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def analyze(output: Path, rule: dict[str, object],
            resources: dict[str, object]) -> dict[str, object]:
    diagnostic = rule["diagnostic_contract"]
    execution = rule["execution_contract"]
    metadata_path = output / "spatial_average_metadata.json"
    crossings_path = output / "phase_eedf_threshold_crossings.csv"
    scalars_path = output / "scalars.csv"
    checkpoint_path = output / f"checkpoint_{execution['aurorapic_end_step']}.apc"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    rows = table(crossings_path)
    selected = [
        row for row in rows
        if row["region"] in diagnostic["critical_regions"]
        and float(diagnostic["critical_phase_fraction"][0]) <=
            float(row["phase_fraction"]) <
            float(diagnostic["critical_phase_fraction"][1])
    ]
    integer_fields = (
        "field_push_macro_observations",
        "field_push_promotion_band_observations",
        "field_push_promotion_band_promotions",
    )
    counts = {name: sum(int(row[name]) for row in selected)
              for name in integer_fields}
    suffix = "eV"
    work_fields = {
        "signed_macro_work_sum_eV":
            f"field_push_promotion_band_signed_macro_work_sum_{suffix}",
        "positive_macro_work_sum_eV":
            f"field_push_promotion_band_positive_macro_work_sum_{suffix}",
        "negative_macro_work_sum_eV":
            f"field_push_promotion_band_negative_macro_work_sum_{suffix}",
    }
    work = {name: sum(float(row[column]) for row in selected)
            for name, column in work_fields.items()}
    closure = work["signed_macro_work_sum_eV"] - (
        work["positive_macro_work_sum_eV"] -
        work["negative_macro_work_sum_eV"])
    scalar_rows = table(scalars_path)
    maximum_particles = max(
        int(row["live_particles_electrons"]) +
        int(row["live_particles_ions"]) for row in scalar_rows)
    observations = counts["field_push_promotion_band_observations"]
    promotions = counts["field_push_promotion_band_promotions"]
    total = counts["field_push_macro_observations"]
    gates = {
        "metadata_contract":
            metadata.get("sampling_order") == diagnostic["sampling_order"] and
            metadata.get("phase_bins") == diagnostic["phase_bins"] and
            metadata.get("phase_eedf_tail_threshold") ==
                diagnostic["promotion_band_eV"][1] and
            metadata.get("phase_eedf_promotion_band_min") ==
                diagnostic["promotion_band_eV"][0],
        "table_shape": len(rows) == 200 * 7 and len(selected) == 25,
        "band_observation_population": observations >=
            diagnostic["minimum_band_observations_per_member"],
        "band_promotion_population": promotions >=
            diagnostic["minimum_band_promotions_per_member"],
        "work_closure": abs(closure) <=
            diagnostic["work_closure_relative_tolerance"] * max(
                1.0,
                work["positive_macro_work_sum_eV"] +
                work["negative_macro_work_sum_eV"]),
        "particle_cap": maximum_particles <=
            execution["maximum_total_particles"],
        "resident_memory": resources["peak_resident_set_kib"] <=
            execution["maximum_peak_resident_set_kib"],
        "final_checkpoint_present": checkpoint_path.is_file(),
    }
    return {
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "critical_scope": {
            "field_push_observations": total,
            "promotion_band_observations": observations,
            "promotion_band_promotions": promotions,
            **work,
            "work_closure_residual_eV": closure,
            "band_supply_fraction": observations / total,
            "band_promotion_probability": promotions / observations,
            "mean_signed_work_eV":
                work["signed_macro_work_sum_eV"] / observations,
            "mean_positive_work_eV":
                work["positive_macro_work_sum_eV"] / observations,
            "mean_negative_work_eV":
                work["negative_macro_work_sum_eV"] / observations,
        },
        "safety": {"maximum_total_live_particles": maximum_particles},
        "final_checkpoint_sha256": sha256(checkpoint_path),
        "output_sha256": {
            path.name: sha256(path)
            for path in (metadata_path, crossings_path, scalars_path)
        },
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise PromotionBandRunError("promotion-band run was not acknowledged")
    rule_path = args.rule.resolve()
    if sha256(rule_path) != APPROVED_RULE_SHA256:
        raise PromotionBandRunError("rule is not approved")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    states = [state for state in rule["locked_inputs"]
              ["aurorapic_continuation_states"]
              if state["id"] == args.state_id]
    if len(states) != 1:
        raise PromotionBandRunError("state id is not locked")
    state = states[0]
    executable = args.executable.resolve()
    checkpoint = args.checkpoint.resolve()
    prior_deck = args.prior_deck.resolve()
    prior_report = args.prior_report.resolve()
    for path, expected, label in (
        (executable, rule["locked_inputs"]["aurorapic_solver_sha256"], "solver"),
        (checkpoint, state["checkpoint_sha256"], "checkpoint"),
        (prior_deck, state["prior_deck_sha256"], "prior deck"),
        (prior_report, state["prior_report_sha256"], "prior report"),
    ):
        if sha256(path) != expected:
            raise PromotionBandRunError(f"locked {label} differs")
    previous = json.loads(prior_report.read_text(encoding="utf-8"))
    if previous.get("all_gates_passed") is not True:
        raise PromotionBandRunError("prior member did not pass")
    work = args.work_dir.resolve()
    if work.exists():
        raise PromotionBandRunError(f"refusing to overwrite {work}")
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    execution = rule["execution_contract"]
    if available_memory < execution["minimum_available_memory_kib"]:
        raise PromotionBandRunError("available memory is below launch floor")
    if available_disk < execution["minimum_available_disk_kib"]:
        raise PromotionBandRunError("available disk is below launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck_path = work / "input.cfg"
    atomic_text(deck_path, make_deck(
        prior_deck.read_text(encoding="utf-8"), output, checkpoint, rule))
    resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(deck_path)], work / "stdout.txt", work / "stderr.txt",
        float(execution["timeout_seconds_per_member"]))
    result = analyze(output, rule, resources)
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_promotion_band_work_member",
        "state_id": args.state_id,
        "rule_sha256": APPROVED_RULE_SHA256,
        "inputs": {
            "solver_sha256": sha256(executable),
            "checkpoint_sha256": sha256(checkpoint),
            "prior_deck_sha256": sha256(prior_deck),
            "prior_report_sha256": sha256(prior_report),
            "deck_sha256": sha256(deck_path),
        },
        "resources": {
            **resources,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk,
        },
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    })
    atomic_json(work / "member-report.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("executable", type=Path)
    parser.add_argument("prior_deck", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("prior_report", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--state-id", required=True)
    parser.add_argument("--acknowledge-cost", required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (PromotionBandRunError, PilotError, OSError, ValueError,
            KeyError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
