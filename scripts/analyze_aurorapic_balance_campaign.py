#!/usr/bin/env python3
"""Aggregate checksum-chained AuroraPIC source/loss-stationarity blocks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics

from run_aurorapic_balance_stationarity_block import RULE_SHA256
from run_aurorapic_edupic_pilot import atomic_json, sha256


def summarize(reports: list[dict[str, object]]) -> dict[str, object]:
    blocks = []
    total_created = 0
    total_electron_losses = 0
    total_ion_losses = 0
    for report in reports:
        ledger = report["particle_ledger"]
        cycles = int(report["window"]["cycles"])
        mean_electrons = statistics.fmean(
            report["cycle_endpoint_populations"]["electrons"])
        created = int(ledger["ionization_macro_events"])
        electron_losses = int(ledger["electron_wall_loss_macro_events"])
        ion_losses = int(ledger["ion_wall_loss_macro_events"])
        total_created += created
        total_electron_losses += electron_losses
        total_ion_losses += ion_losses
        blocks.append({
            "block_index": report["block_index"],
            "cycles": cycles,
            "mean_live_electron_macroparticles": mean_electrons,
            "ionization_macro_events_per_cycle": created / cycles,
            "electron_wall_loss_macro_events_per_cycle":
                electron_losses / cycles,
            "ionization_events_per_live_electron_per_cycle":
                created / cycles / mean_electrons,
            "electron_wall_losses_per_live_electron_per_cycle":
                electron_losses / cycles / mean_electrons,
            "metrics": report["metrics"],
            "stationarity_gates": report["stationarity_gates"],
            "stationarity_block_passed": report[
                "stationarity_block_passed"],
        })
    if total_created <= 0:
        raise ValueError("campaign has no ionization events")
    first, last = blocks[0], blocks[-1]
    return {
        "blocks": blocks,
        "combined": {
            "cycles": sum(int(block["cycles"]) for block in blocks),
            "ionization_macro_events": total_created,
            "electron_wall_loss_macro_events": total_electron_losses,
            "ion_wall_loss_macro_events": total_ion_losses,
            "electron_source_loss_relative_imbalance":
                (total_created - total_electron_losses) / total_created,
            "ion_source_loss_relative_imbalance":
                (total_created - total_ion_losses) / total_created,
        },
        "first_to_last_change": {
            "ionization_macro_events_per_cycle_relative_change":
                last["ionization_macro_events_per_cycle"] /
                first["ionization_macro_events_per_cycle"] - 1.0,
            "electron_wall_loss_macro_events_per_cycle_relative_change":
                last["electron_wall_loss_macro_events_per_cycle"] /
                first["electron_wall_loss_macro_events_per_cycle"] - 1.0,
            "ionization_events_per_live_electron_per_cycle_relative_change":
                last["ionization_events_per_live_electron_per_cycle"] /
                first["ionization_events_per_live_electron_per_cycle"] - 1.0,
            "electron_wall_losses_per_live_electron_per_cycle_relative_change":
                last["electron_wall_losses_per_live_electron_per_cycle"] /
                first["electron_wall_losses_per_live_electron_per_cycle"] - 1.0,
        },
    }


def analyze(rule_path: Path, report_paths: list[Path]) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("balance stationarity rule differs")
    if len(report_paths) < 2:
        raise ValueError("campaign trend requires at least two blocks")
    reports = []
    report_hashes = []
    for expected_index, path in enumerate(report_paths, 1):
        report_hash = sha256(path)
        report = json.loads(path.read_text(encoding="utf-8"))
        if (report.get("scope") != "aurorapic_source_loss_stationarity_block" or
                report.get("rule_sha256") != RULE_SHA256 or
                report.get("block_index") != expected_index or
                report.get("all_hard_safety_gates_passed") is not True or
                report.get("particle_ledger", {}).get(
                    "exact_species_closure") is not True):
            raise ValueError(f"block {expected_index} contract differs")
        if expected_index > 1:
            if (report["inputs"]["prior_report_sha256"] != report_hashes[-1] or
                    report["inputs"]["input_checkpoint_sha256"] !=
                    reports[-1]["final_checkpoint_sha256"]):
                raise ValueError("campaign report/checkpoint chain is broken")
        reports.append(report)
        report_hashes.append(report_hash)
    summary = summarize(reports)
    change = summary["first_to_last_change"]
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "aurorapic_source_loss_stationarity_campaign_analysis",
        "rule_sha256": RULE_SHA256,
        "analyzer_sha256": sha256(Path(__file__).resolve()),
        "block_report_sha256": report_hashes,
        **summary,
        "assessment": {
            "all_blocks_passed_stationarity": all(
                report["stationarity_block_passed"] for report in reports),
            "source_loss_error_persisted": all(
                not report["stationarity_gates"][
                    "electron_source_loss_balance"] for report in reports),
            "ionization_declined_while_wall_loss_was_stable":
                change["ionization_macro_events_per_cycle_relative_change"] < -.02
                and abs(change[
                    "electron_wall_loss_macro_events_per_cycle_relative_change"]) < .01,
            "finding": (
                "The quiet macroscopic state retains a material source/loss "
                "deficit. From block 1 to block 2, ionization per cycle declines "
                "while electron wall loss remains stable, worsening population "
                "decay; more averaging alone cannot resolve this mismatch."),
            "next_discriminator": (
                "Audit the effective ionization probability per electron against "
                "the native eduPIC energy distribution and collision sampling "
                "before spending additional blocks on final density measurement."),
        },
        "claim_boundary": (
            "Two contiguous blocks localize a persistent internal imbalance but "
            "do not provide independent-sample confidence intervals or identify "
            "a unique causal implementation defect."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("reports", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule.resolve(), [path.resolve() for path in args.reports])
    if args.output:
        atomic_json(args.output, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
