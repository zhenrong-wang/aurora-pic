#!/usr/bin/env python3
"""Run one locked AuroraPIC near-threshold mover-decomposition member."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, run_process,
    sha256,
)
from run_aurorapic_promotion_band_work import make_deck, table


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_MOVER_DECOMPOSITION_RUN"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256 = (
    "e9ff35fec7486e109e661bf26eb8ec54449b9cdebad02ce12cbb9a84b03854ae"
)


class MoverDecompositionRunError(RuntimeError):
    pass


def summarize(selected: list[dict[str, str]]) -> dict[str, float | int]:
    count_fields = (
        "field_push_macro_observations",
        "field_push_promotion_band_observations",
        "field_push_promotion_band_promotions",
    )
    values: dict[str, float | int] = {
        field: sum(int(row[field]) for row in selected)
        for field in count_fields
    }
    columns = {
        "signed_macro_work_sum_eV":
            "field_push_promotion_band_signed_macro_work_sum_eV",
        "positive_macro_work_sum_eV":
            "field_push_promotion_band_positive_macro_work_sum_eV",
        "negative_macro_work_sum_eV":
            "field_push_promotion_band_negative_macro_work_sum_eV",
        "origin_energy_sum_eV":
            "field_push_promotion_band_origin_macro_energy_sum_eV",
        "origin_longitudinal_energy_sum_eV":
            "field_push_promotion_band_origin_longitudinal_macro_energy_sum_eV",
        "linear_work_sum_eV":
            "field_push_promotion_band_linear_macro_work_sum_eV",
        "positive_linear_work_sum_eV":
            "field_push_promotion_band_positive_linear_macro_work_sum_eV",
        "negative_linear_work_sum_eV":
            "field_push_promotion_band_negative_linear_macro_work_sum_eV",
        "quadratic_work_sum_eV":
            "field_push_promotion_band_quadratic_macro_work_sum_eV",
    }
    values.update({
        name: sum(float(row[column]) for row in selected)
        for name, column in columns.items()
    })
    observations = int(values["field_push_promotion_band_observations"])
    origin = float(values["origin_energy_sum_eV"])
    values.update({
        "band_supply_fraction": observations /
            int(values["field_push_macro_observations"]),
        "band_promotion_probability":
            int(values["field_push_promotion_band_promotions"]) / observations,
        "mean_origin_energy_eV": origin / observations,
        "origin_longitudinal_energy_fraction":
            float(values["origin_longitudinal_energy_sum_eV"]) / origin,
        "mean_positive_linear_work_eV":
            float(values["positive_linear_work_sum_eV"]) / observations,
        "mean_quadratic_work_eV":
            float(values["quadratic_work_sum_eV"]) / observations,
        "mean_positive_work_eV":
            float(values["positive_macro_work_sum_eV"]) / observations,
        "linear_work_closure_residual_eV":
            float(values["linear_work_sum_eV"]) -
            (float(values["positive_linear_work_sum_eV"]) -
             float(values["negative_linear_work_sum_eV"])),
        "total_work_decomposition_closure_residual_eV":
            float(values["signed_macro_work_sum_eV"]) -
            (float(values["linear_work_sum_eV"]) +
             float(values["quadratic_work_sum_eV"])),
        "signed_work_closure_residual_eV":
            float(values["signed_macro_work_sum_eV"]) -
            (float(values["positive_macro_work_sum_eV"]) -
             float(values["negative_macro_work_sum_eV"])),
    })
    return values


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
        if row["region"] == diagnostic["critical_region"]
        and float(diagnostic["critical_phase_fraction"][0]) <=
            float(row["phase_fraction"]) <
            float(diagnostic["critical_phase_fraction"][1])
    ]
    critical = summarize(selected)
    scalar_rows = table(scalars_path)
    maximum_particles = max(
        int(row["live_particles_electrons"]) +
        int(row["live_particles_ions"]) for row in scalar_rows)
    tolerance = float(diagnostic["work_closure_relative_tolerance"])
    total_scale = max(
        1.0, float(critical["positive_macro_work_sum_eV"]) +
        float(critical["negative_macro_work_sum_eV"]))
    linear_scale = max(
        1.0, float(critical["positive_linear_work_sum_eV"]) +
        float(critical["negative_linear_work_sum_eV"]))
    gates = {
        "metadata_contract":
            metadata.get("sampling_order") == diagnostic["sampling_order"] and
            metadata.get("phase_bins") == diagnostic["phase_bins"] and
            metadata.get("phase_eedf_tail_threshold") ==
                diagnostic["promotion_band_eV"][1] and
            metadata.get("phase_eedf_promotion_band_min") ==
                diagnostic["promotion_band_eV"][0],
        "table_shape": len(rows) == 200 * 7 and len(selected) == 25,
        "band_observation_population":
            int(critical["field_push_promotion_band_observations"]) >=
            int(diagnostic["minimum_band_observations_per_member"]),
        "band_promotion_population":
            int(critical["field_push_promotion_band_promotions"]) >=
            int(diagnostic["minimum_band_promotions_per_member"]),
        "origin_energy_partition":
            0.0 <= float(critical[
                "origin_longitudinal_energy_fraction"]) <= 1.0,
        "signed_work_closure": abs(float(critical[
            "signed_work_closure_residual_eV"])) <= tolerance * total_scale,
        "linear_work_closure": abs(float(critical[
            "linear_work_closure_residual_eV"])) <= tolerance * linear_scale,
        "total_work_decomposition_closure": abs(float(critical[
            "total_work_decomposition_closure_residual_eV"])) <=
                tolerance * total_scale,
        "particle_cap": maximum_particles <=
            int(execution["maximum_total_particles"]),
        "resident_memory": resources["peak_resident_set_kib"] <=
            int(execution["maximum_peak_resident_set_kib"]),
        "final_checkpoint_present": checkpoint_path.is_file(),
    }
    return {
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "critical_scope": critical,
        "safety": {"maximum_total_live_particles": maximum_particles},
        "final_checkpoint_sha256": sha256(checkpoint_path),
        "output_sha256": {
            path.name: sha256(path)
            for path in (metadata_path, crossings_path, scalars_path)
        },
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MoverDecompositionRunError(
            "mover-decomposition run was not acknowledged")
    rule_path = args.rule.resolve()
    if sha256(rule_path) != APPROVED_RULE_SHA256:
        raise MoverDecompositionRunError("rule is not approved")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    states = [state for state in rule["locked_inputs"]
              ["aurorapic_continuation_states"]
              if state["id"] == args.state_id]
    if len(states) != 1:
        raise MoverDecompositionRunError("state id is not locked")
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
            raise MoverDecompositionRunError(f"locked {label} differs")
    if json.loads(prior_report.read_text(encoding="utf-8")).get(
            "all_gates_passed") is not True:
        raise MoverDecompositionRunError("prior member did not pass")
    work = args.work_dir.resolve()
    if work.exists():
        raise MoverDecompositionRunError(f"refusing to overwrite {work}")
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    execution = rule["execution_contract"]
    if available_memory < execution["minimum_available_memory_kib"]:
        raise MoverDecompositionRunError(
            "available memory is below launch floor")
    if available_disk < execution["minimum_available_disk_kib"]:
        raise MoverDecompositionRunError(
            "available disk is below launch floor")
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
        "scope": "aurorapic_near_threshold_mover_decomposition_member",
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
    except (MoverDecompositionRunError, PilotError, OSError, ValueError,
            KeyError, ZeroDivisionError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
