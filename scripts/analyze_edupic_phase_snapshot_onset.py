#!/usr/bin/env python3
"""Evaluate the preregistered four-cycle phase-snapshot onset comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

from analyze_edupic_grid_field_sampling import spatial_mean_square


METRICS = ("regional_mean_squared_field_V2_m2",
           "regional_charge_density_rms_C_m3")
PASSIVE_FILES = (
    "stdout.txt", "edupic_phase_eedf.csv", "edupic_phase_eedf_moments.csv",
    "edupic_phase_eedf_history.csv",
    "edupic_phase_eedf_threshold_crossings.csv",
    "edupic_phase_eedf_field_push_thresholds.csv",
    "edupic_phase_eedf_promotion_band_work.csv",
    "edupic_phase_eedf_mover_decomposition.csv",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot_set_sha256(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: int(item.stem.split("_")[-1])):
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def profile_metrics(electric_field: list[float], charge_density: list[float],
                    length: float, lower: float, upper: float) -> dict[str, float]:
    e2 = spatial_mean_square(electric_field, length, lower, upper)
    rho2 = spatial_mean_square(charge_density, length, lower, upper)
    return {METRICS[0]: e2, METRICS[1]: math.sqrt(rho2)}


def relative_range(values: list[float]) -> float:
    mean = math.fsum(values) / len(values)
    return (max(values) - min(values)) / max(abs(mean), 1e-300)


def read_native(path: Path, cycles: int, phases: int, nodes: int,
                length: float, lower: float, upper: float) -> dict[tuple[int, int], dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != cycles * phases * nodes:
        raise ValueError(f"native snapshot ledger has {len(rows)} rows")
    result = {}
    for cycle in range(1, cycles + 1):
        for phase in range(1, phases + 1):
            block = rows[((cycle - 1) * phases + phase - 1) * nodes:
                         ((cycle - 1) * phases + phase) * nodes]
            fields, charges = [], []
            for node, row in enumerate(block):
                if int(row["measurement_cycle"]) != cycle or int(row["phase_index"]) != phase:
                    raise ValueError("native cycle/phase ordering differs")
                if abs(float(row["phase_fraction"]) - phase / phases) > 1e-12:
                    raise ValueError("native phase fraction differs")
                if int(row["node"]) != node or abs(float(row["x_m"]) - node * length / (nodes - 1)) > 1e-12:
                    raise ValueError("native coordinate ordering differs")
                values = [float(row[key]) for key in
                          ("charge_density_C_m3", "potential_V", "electric_field_V_m")]
                if not all(math.isfinite(value) for value in values):
                    raise ValueError("native snapshot contains a non-finite value")
                charges.append(values[0])
                fields.append(values[2])
            result[cycle, phase] = profile_metrics(fields, charges, length, lower, upper)
    return result


def read_aurora(root: Path, cycles: int, phases: int, nodes: int,
                length: float, lower: float, upper: float) -> tuple[dict[tuple[int, int], dict[str, float]], list[Path]]:
    paths = sorted((path for path in root.glob("fields_*.csv")
                    if 36400 <= int(path.stem.split("_")[-1]) <= 52000),
                   key=lambda item: int(item.stem.split("_")[-1]))
    if len(paths) != cycles * phases:
        raise ValueError(f"AuroraPIC snapshot set has {len(paths)} files")
    result = {}
    for index, path in enumerate(paths):
        expected_step = 36400 + index * 400
        if int(path.stem.split("_")[-1]) != expected_step:
            raise ValueError("AuroraPIC snapshot step ordering differs")
        cycle, phase0 = divmod(index, phases)
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != nodes:
            raise ValueError("AuroraPIC snapshot node count differs")
        fields, charges = [], []
        for node, row in enumerate(rows):
            x = float(row["x"])
            values = [float(row[key]) for key in ("rho", "phi", "E")]
            if abs(x - node * length / (nodes - 1)) > 1e-12:
                raise ValueError("AuroraPIC coordinate ordering differs")
            if not all(math.isfinite(value) for value in values):
                raise ValueError("AuroraPIC snapshot contains a non-finite value")
            charges.append(values[0])
            fields.append(values[2])
        result[cycle + 1, phase0 + 1] = profile_metrics(fields, charges, length, lower, upper)
    return result, paths


def parse_resources(path: Path) -> tuple[float, int]:
    text = path.read_text(encoding="utf-8")
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if elapsed is None or rss is None:
        raise ValueError("resource report is incomplete")
    parts = [float(value) for value in elapsed.group(1).split(":")]
    seconds = parts[-1] + (parts[-2] * 60 if len(parts) > 1 else 0) + (parts[-3] * 3600 if len(parts) > 2 else 0)
    return seconds, int(rss.group(1))


def analyze(rule_path: Path, native_root: Path, native_reference_root: Path,
            aurora_roots: list[Path], instrumenter: Path, source: Path,
            instrumented_source: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    contract = rule["execution_contract"]; reduction = rule["reduction_contract"]
    locked = rule["locked_inputs"]
    cycles, phases, nodes = (int(contract[key]) for key in
                             ("native_cycles", "snapshots_per_cycle", "nodes"))
    length = float(reduction["length_m"])
    lower, upper = [length * float(value) for value in reduction["critical_x_over_L"]]
    phase_indices = [round(float(value) * phases) for value in reduction["phase_neighborhood"]]
    primary = round(float(reduction["primary_phase_fraction"]) * phases)

    native_members = []
    for seed in locked["native_seeds"]:
        member_root = native_root / f"seed-{seed}"
        ledger = member_root / "edupic_phase_snapshots.csv"
        metrics = read_native(ledger, cycles, phases, nodes, length, lower, upper)
        checkpoint_hash = sha256(member_root / "picdata.bin")
        passive = {name: sha256(member_root / name) == sha256(
            native_reference_root / f"seed-{seed}" / name) for name in PASSIVE_FILES}
        seconds, rss = parse_resources(member_root / "stderr.txt")
        native_members.append({"seed": seed, "ledger_sha256": sha256(ledger),
                               "final_checkpoint_sha256": checkpoint_hash,
                               "passive_diagnostics": passive, "wall_seconds": seconds,
                               "peak_resident_set_kib": rss, "metrics": metrics})

    aurora_members = []
    for root, expected in zip(aurora_roots, locked["aurorapic_members"], strict=True):
        metrics, paths = read_aurora(root, cycles, phases, nodes, length, lower, upper)
        aurora_members.append({"id": expected["id"],
                               "snapshot_set_sha256": snapshot_set_sha256(paths),
                               "metrics": metrics})

    def summarize(members: list[dict[str, object]]) -> tuple[dict[str, object], bool]:
        summary = {}; repeatable = True
        for cycle in range(1, cycles + 1):
            for phase in phase_indices:
                values = {metric: [member["metrics"][cycle, phase][metric] for member in members]
                          for metric in METRICS}
                ranges = {metric: relative_range(items) for metric, items in values.items()}
                means = {metric: math.fsum(items) / len(items) for metric, items in values.items()}
                summary[f"cycle_{cycle}_phase_{phase / phases:.1f}"] = {
                    "means": means, "relative_ranges": ranges}
        return summary, repeatable

    native_summary, _ = summarize(native_members); aurora_summary, _ = summarize(aurora_members)
    comparisons = {}
    for key in native_summary:
        comparisons[key] = {
            "ensemble_mean": {metric: aurora_summary[key]["means"][metric] /
                              native_summary[key]["means"][metric] for metric in METRICS},
            "aurorapic_members": {
                str(member["id"]): {metric: member["metrics"][
                    int(key.split("_")[1]), round(float(key.split("_")[3]) * phases)][metric] /
                    native_summary[key]["means"][metric] for metric in METRICS}
                for member in aurora_members},
        }

    integrity = {
        "rule_basis_hash_matches": sha256(Path("benchmarks/ccp/edupic-net-charge-sheath-result-20260826.json")) == rule["basis"]["net_charge_sheath_result_sha256"],
        "source_hash_matches": sha256(source) == rule["instrumentation"]["source_sha256"],
        "instrumenter_hash_matches": sha256(instrumenter) == rule["instrumentation"]["instrumenter_sha256"],
        "instrumented_source_hash_matches": sha256(instrumented_source) == rule["instrumentation"]["instrumented_source_sha256"],
        "native_checkpoint_hashes_match": all(member["final_checkpoint_sha256"] == locked["expected_passive_native_checkpoint_sha256"][str(member["seed"])] for member in native_members),
        "native_prior_diagnostics_byte_identical": all(all(member["passive_diagnostics"].values()) for member in native_members),
        "native_binary_hashes_match": all(sha256(native_root / f"seed-{seed}" / "edupic") == rule["instrumentation"]["instrumented_binary_sha256"] for seed in locked["native_seeds"]),
        "aurora_snapshot_hashes_match": all(member["snapshot_set_sha256"] == expected["snapshot_set_sha256"] for member, expected in zip(aurora_members, locked["aurorapic_members"], strict=True)),
        "resource_limits_pass": all(member["peak_resident_set_kib"] <= int(contract["maximum_peak_resident_set_kib"]) and member["wall_seconds"] <= int(contract["timeout_seconds_per_member"]) for member in native_members),
    }
    integrity_passed = all(integrity.values())
    aurora_repeat = all(value <= .15 for cell in aurora_summary.values() for value in cell["relative_ranges"].values())
    native_repeat = all(value <= .20 for cell in native_summary.values() for value in cell["relative_ranges"].values())
    primary_field_repeat = {
        "aurorapic": all(aurora_summary[f"cycle_{cycle}_phase_{primary / phases:.1f}"]
                         ["relative_ranges"][METRICS[0]] <= .15
                         for cycle in range(1, cycles + 1)),
        "native": all(native_summary[f"cycle_{cycle}_phase_{primary / phases:.1f}"]
                      ["relative_ranges"][METRICS[0]] <= .20
                      for cycle in range(1, cycles + 1)),
    }
    repeatability_details = {}
    for label, summary, limit in (("aurorapic", aurora_summary, .15),
                                  ("native", native_summary, .20)):
        failures = [{"cell": cell, "metric": metric, "relative_range": value}
                    for cell, item in summary.items()
                    for metric, value in item["relative_ranges"].items()
                    if value > limit]
        repeatability_details[label] = {
            "maximum_relative_range": max(
                value for item in summary.values()
                for value in item["relative_ranges"].values()),
            "failures": failures,
        }
    stability = {"aurorapic": {}, "native": {}}
    for label, members in (("aurorapic", aurora_members), ("native", native_members)):
        for member in members:
            ident = str(member.get("id", member.get("seed")))
            values = [member["metrics"][cycle, primary][METRICS[0]] for cycle in range(1, cycles + 1)]
            stability[label][ident] = relative_range(values)
    stability_passed = all(value <= .20 for group in stability.values() for value in group.values())
    member_ids = [str(member["id"]) for member in aurora_members]
    primary_ratios = {
        member_id: [comparisons[f"cycle_{cycle}_phase_{primary / phases:.1f}"]
                    ["aurorapic_members"][member_id][METRICS[0]]
                    for cycle in range(1, cycles + 1)]
        for member_id in member_ids}
    four_cycle_ratios = {
        member_id: {phase / phases: math.fsum(
            comparisons[f"cycle_{cycle}_phase_{phase / phases:.1f}"]
            ["aurorapic_members"][member_id][METRICS[0]]
            for cycle in range(1, cycles + 1)) / cycles
                    for phase in phase_indices}
        for member_id in member_ids}
    gates_passed = integrity_passed and aurora_repeat and native_repeat and stability_passed
    cycle_one = all(values[0] <= .90 for values in primary_ratios.values())
    persists = all(all(value <= .90 for value in values)
                   for values in primary_ratios.values())
    emerges = (all(.95 <= values[0] <= 1.05 for values in primary_ratios.values())
               and any(all(values[cycle] <= .90 for values in primary_ratios.values())
                       for cycle in range(1, cycles)))
    neighborhood = all(all(value <= .90 for value in values.values())
                       for values in four_cycle_ratios.values())
    if not gates_passed: outcome = "inconclusive_failed_joint_gate"
    elif persists: outcome = "primary_phase_field_deficit_present_cycle_one_and_persists"
    elif emerges: outcome = "primary_phase_field_deficit_emerges_within_window"
    elif cycle_one: outcome = "primary_phase_field_deficit_present_cycle_one_but_not_persistent"
    else: outcome = "mixed_or_intermediate_snapshot_onset_result"
    return {
        "schema_version": 1, "case_id": rule["case_id"],
        "scope": "four_cycle_phase_snapshot_onset_result",
        "inputs": {"rule_sha256": sha256(rule_path), "instrumenter_sha256": sha256(instrumenter)},
        "integrity": integrity, "integrity_gate_passed": integrity_passed,
        "repeatability": {"aurorapic_passed": aurora_repeat, "native_passed": native_repeat,
                          "aurorapic_limit": .15, "native_limit": .20,
                          "post_hoc_primary_phase_field_only": {
                              **primary_field_repeat,
                              "both_passed": all(primary_field_repeat.values()),
                              "not_a_substitute_for_locked_joint_gate": True,
                          },
                          "details": repeatability_details},
        "cycle_stability": {"relative_ranges": stability, "limit": .20, "passed": stability_passed},
        "native_members": [{key: value for key, value in member.items() if key != "metrics"} for member in native_members],
        "aurorapic_members": [{key: value for key, value in member.items() if key != "metrics"} for member in aurora_members],
        "native_summary": native_summary, "aurorapic_summary": aurora_summary,
        "aurorapic_to_native_ratios": comparisons,
        "primary_phase_field_energy_ratios_by_cycle": primary_ratios,
        "four_cycle_field_energy_ratios_by_phase": four_cycle_ratios,
        "decision": {"all_joint_gates_passed": gates_passed,
                     "deficit_present_in_cycle_one": cycle_one,
                     "deficit_persists_all_cycles": persists,
                     "deficit_emerges_within_window": emerges,
                     "phase_neighborhood_deficit": neighborhood,
                     "formal_outcome": outcome},
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--native-reference-root", type=Path, required=True)
    parser.add_argument("--aurora-root", type=Path, action="append", required=True)
    parser.add_argument("--instrumenter", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--instrumented-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.rule, args.native_root, args.native_reference_root,
                     args.aurora_root, args.instrumenter, args.source,
                     args.instrumented_source)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()
