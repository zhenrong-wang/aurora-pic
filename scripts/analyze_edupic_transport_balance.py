#!/usr/bin/env python3
"""Audit eduPIC/AuroraPIC phase resolution and particle source/loss balance."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import statistics

from analyze_aurorapic_heating_seed_ensemble import vector_mean
from analyze_aurorapic_microstate_ensemble import load_member
from compare_edupic_phase_space import (
    REFERENCE_FILES, flatten_phase_major, metrics, read_matrix, resample_matrix,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


NODES = 400
PHASES = 200
COARSE_PHASES = 16
MACRO_WEIGHT = 7.0e8
CAMPAIGN_SHA256 = (
    "6d690d2c786043935d6677139acbefd5cedb03fc2d58dc131d18e825afcf4693")
LONG_WINDOW_REPORT_SHA256 = (
    "1d87b9ebfbe3668513c844fc5314868e1bae625df259c7f659981f9ba03a9b4a")
CONTINUOUS_REPORT_HASHES = {
    "cycles_20_24": "cdcda4685741833aad2806e4c7a44cab5565c0815b1bb3ba382a3ef3ce7722f3",
    "cycles_24_28": "169d8e5cb0748d65e209dacdc21f9b278177fce3c882fbab033a0eae7017f95f",
    "cycles_28_32": "8d96f51fa9b65e3cf848f5e74c2baa2bd2dfa4948c7819e01f61486365469eb9",
}


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def counter_delta(rows: list[dict[str, str]], key: str) -> int:
    if len(rows) < 2:
        raise ValueError("counter diagnostic needs at least two rows")
    return int(rows[-1][key]) - int(rows[0][key])


def candidate_balance(output: Path) -> dict[str, object]:
    collisions = csv_rows(output / "collisions.csv")
    boundary = csv_rows(output / "boundary_losses.csv")
    scalars = csv_rows(output / "scalars.csv")
    if not (collisions[0]["step"] == boundary[0]["step"] == scalars[0]["step"] and
            collisions[-1]["step"] == boundary[-1]["step"] == scalars[-1]["step"]):
        raise ValueError("candidate ledger windows differ")
    duration = float(scalars[-1]["time"]) - float(scalars[0]["time"])
    if duration <= 0.0:
        raise ValueError("candidate ledger duration is not positive")
    created = counter_delta(
        collisions, "cumulative_collisions_electron_mcc.ionization")
    electron_loss = sum(counter_delta(boundary, key) for key in (
        "absorbed_left_count_electrons", "absorbed_right_count_electrons"))
    ion_loss = sum(counter_delta(boundary, key) for key in (
        "absorbed_left_count_ions", "absorbed_right_count_ions"))
    electron_change = int(scalars[-1]["live_particles_electrons"]) - int(
        scalars[0]["live_particles_electrons"])
    ion_change = int(scalars[-1]["live_particles_ions"]) - int(
        scalars[0]["live_particles_ions"])
    if electron_change != created - electron_loss:
        raise ValueError("electron particle ledger does not close")
    if ion_change != created - ion_loss:
        raise ValueError("ion particle ledger does not close")
    source_flux = created * MACRO_WEIGHT / duration
    electron_flux = electron_loss * MACRO_WEIGHT / duration
    ion_flux = ion_loss * MACRO_WEIGHT / duration
    return {
        "first_step": int(scalars[0]["step"]),
        "last_step": int(scalars[-1]["step"]),
        "duration_s": duration,
        "ionization_macro_events": created,
        "electron_wall_loss_macro_events": electron_loss,
        "ion_wall_loss_macro_events": ion_loss,
        "electron_live_particle_change": electron_change,
        "ion_live_particle_change": ion_change,
        "ionization_source_flux_m-2_s-1": source_flux,
        "electron_wall_flux_m-2_s-1": electron_flux,
        "ion_wall_flux_m-2_s-1": ion_flux,
        "electron_source_minus_loss_fraction":
            (source_flux - electron_flux) / source_flux,
        "ion_source_minus_loss_fraction": (source_flux - ion_flux) / source_flux,
        "input_sha256": {
            name: sha256(output / name) for name in (
                "collisions.csv", "boundary_losses.csv", "scalars.csv")
        },
        "exact_macro_particle_ledger_closure": True,
    }


def parse_flux(info: str, species: str) -> float:
    pattern = re.compile(
        rf"{species} flux at (?:powered|grounded) electrode\s*=\s*"
        rf"([0-9.+\-Ee]+)", re.IGNORECASE)
    values = [float(value) for value in pattern.findall(info)]
    if len(values) != 2:
        raise ValueError(f"expected two {species} wall fluxes")
    return math.fsum(values)


def matrix_phase_space(values: list[float]) -> list[list[float]]:
    if len(values) != PHASES * NODES:
        raise ValueError("candidate phase-space vector has wrong shape")
    return [[values[phase * NODES + node] for phase in range(PHASES)]
            for node in range(NODES)]


def phase_resolution_metrics(candidate: list[float], reference: list[list[float]]) -> dict[str, object]:
    fine = metrics(candidate, flatten_phase_major(reference))
    coarse_candidate = resample_matrix(matrix_phase_space(candidate), COARSE_PHASES)
    coarse_reference = resample_matrix(reference, COARSE_PHASES)
    coarse = metrics(flatten_phase_major(coarse_candidate),
                     flatten_phase_major(coarse_reference))
    fine_error = float(fine["relative_l2"])
    coarse_error = float(coarse["relative_l2"])
    ratio = (fine_error / coarse_error if coarse_error > 0.0 else
             (1.0 if fine_error == 0.0 else math.inf))
    return {
        "fine_200_phase_bins": fine,
        "coarse_grained_16_phase_bins": coarse,
        "fine_to_coarse_relative_l2_ratio": ratio,
    }


def reference_balance(campaign_root: Path) -> dict[str, object]:
    report_path = campaign_root / "campaign-report.json"
    if sha256(report_path) != CAMPAIGN_SHA256:
        raise ValueError("native eduPIC campaign report differs")
    campaign = json.loads(report_path.read_text(encoding="utf-8"))
    sources, electrons, ions = [], [], []
    stage_hashes = []
    for stage in campaign["stages"]:
        root = campaign_root / f"stage-{stage['start_cycle']:06d}-{stage['end_cycle']:06d}"
        stage_report_path = root / "stage-report.json"
        if sha256(stage_report_path) != stage["stage_report_sha256"]:
            raise ValueError("native eduPIC stage report differs")
        stage_report = json.loads(stage_report_path.read_text(encoding="utf-8"))
        diagnostics = stage_report["outputs"]["diagnostics"]
        for name in ("ioniz_xt.dat", "info.txt"):
            if sha256(root / name) != diagnostics[name]["sha256"]:
                raise ValueError(f"native eduPIC {name} differs")
        ionization = read_matrix(root / "ioniz_xt.dat")
        # eduPIC's 400 diagnostics are endpoint-inclusive nodal samples.
        dx = 0.025 / (NODES - 1)
        source = math.fsum(
            (0.5 if node in (0, NODES - 1) else 1.0) *
            math.fsum(ionization[node]) / PHASES
            for node in range(NODES)) * dx
        info = (root / "info.txt").read_text(encoding="utf-8")
        sources.append(source)
        electrons.append(parse_flux(info, "Electron"))
        ions.append(parse_flux(info, "Ion"))
        stage_hashes.append(stage["stage_report_sha256"])

    def summary(values: list[float]) -> dict[str, float]:
        return {"mean": statistics.fmean(values),
                "sample_standard_deviation": statistics.stdev(values),
                "minimum": min(values), "maximum": max(values)}

    electron_fraction = [(s - w) / s for s, w in zip(sources, electrons)]
    ion_fraction = [(s - w) / s for s, w in zip(sources, ions)]
    return {
        "campaign_report_sha256": CAMPAIGN_SHA256,
        "blocks": len(sources),
        "cycles_per_block": 16,
        "stage_report_sha256": stage_hashes,
        "ionization_source_flux_m-2_s-1": summary(sources),
        "electron_wall_flux_m-2_s-1": summary(electrons),
        "ion_wall_flux_m-2_s-1": summary(ions),
        "electron_source_minus_loss_fraction": summary(electron_fraction),
        "ion_source_minus_loss_fraction": summary(ion_fraction),
        "statistical_boundary": (
            "Contiguous blocks are summarized descriptively; they are not "
            "treated as independent samples."),
    }


def analyze(rule_path: Path, candidate_roots: dict[str, Path],
            continuous_roots: dict[str, Path], long_window_root: Path,
            reference_raw: Path,
            campaign_root: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    members = [load_member("locked_source_microstate",
                           candidate_roots["locked_source_microstate"], rule, True)]
    members.extend(load_member(name, candidate_roots[name], rule)
                   for name in ("microstate_51949", "microstate_63059"))
    ensemble = {name: vector_mean([member[name] for member in members])
                for name in ("density", "field", "current", "power")}
    phase_audit = {}
    for key, reference_key in (
            ("density", "electron_density"), ("field", "electric_field"),
            ("current", "electron_current_density"),
            ("power", "electron_ohmic_power_density")):
        filename, expected = REFERENCE_FILES[reference_key]
        path = reference_raw / filename
        if sha256(path) != expected:
            raise ValueError(f"locked reference differs: {filename}")
        phase_audit[key] = phase_resolution_metrics(
            ensemble[key], read_matrix(path))

    candidate_ledgers = {
        name: candidate_balance(root / "measurement" / "output")
        for name, root in candidate_roots.items()
    }
    continuous_ledgers = {}
    for name, root in continuous_roots.items():
        report_path = root / "measurement-report.json"
        if sha256(report_path) != CONTINUOUS_REPORT_HASHES[name]:
            raise ValueError(f"{name} measurement report differs")
        continuous_ledgers[name] = candidate_balance(root / "output")

    long_report_path = long_window_root / "branch-report.json"
    if sha256(long_report_path) != LONG_WINDOW_REPORT_SHA256:
        raise ValueError("long-window branch report differs")
    long_window_ledger = candidate_balance(
        long_window_root / "measurement" / "output")

    reference = reference_balance(campaign_root)
    candidate_e_imbalances = [float(value["electron_source_minus_loss_fraction"])
                              for value in candidate_ledgers.values()]
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "descriptive_phase_resolution_and_particle_balance_audit",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "microstate_rule_sha256": sha256(rule_path),
            "candidate_branch_report_sha256": {
                name: sha256(root / "branch-report.json")
                for name, root in candidate_roots.items()
            },
            "continuous_measurement_report_sha256": CONTINUOUS_REPORT_HASHES,
            "long_window_branch_report_sha256": LONG_WINDOW_REPORT_SHA256,
            "native_campaign_report_sha256": CAMPAIGN_SHA256,
        },
        "phase_resolution": phase_audit,
        "aurorapic_matched_microstate_ledgers": candidate_ledgers,
        "aurorapic_continuous_ledgers": continuous_ledgers,
        "aurorapic_long_window_ledger": long_window_ledger,
        "edupic_native_production_balance": reference,
        "assessment": {
            "fine_current_error_is_resolution_sensitive":
                phase_audit["current"]["fine_to_coarse_relative_l2_ratio"] > 2.0,
            "coarse_current_relative_l2": phase_audit["current"][
                "coarse_grained_16_phase_bins"]["relative_l2"],
            "coarse_field_relative_l2": phase_audit["field"][
                "coarse_grained_16_phase_bins"]["relative_l2"],
            "coarse_density_relative_l2": phase_audit["density"][
                "coarse_grained_16_phase_bins"]["relative_l2"],
            "mean_matched_aurorapic_electron_source_minus_loss_fraction":
                statistics.fmean(candidate_e_imbalances),
            "mean_native_edupic_electron_source_minus_loss_fraction":
                reference["electron_source_minus_loss_fraction"]["mean"],
            "long_window_aurorapic_electron_source_minus_loss_fraction":
                long_window_ledger["electron_source_minus_loss_fraction"],
            "finding": (
                "Conservative RF-phase coarse graining recovers strong current "
                "and field waveform agreement. The remaining material discrepancy "
                "is electron density and the associated source/loss balance."),
        },
        "claim_boundary": (
            "This post-hoc diagnostic localizes discrepancies in checksum-bound "
            "simulation outputs. It is not an experimental validation or a "
            "predeclared statistical acceptance test."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("microstate_51949", type=Path)
    parser.add_argument("microstate_63059", type=Path)
    parser.add_argument("continuous_cycle24", type=Path)
    parser.add_argument("continuous_cycle28", type=Path)
    parser.add_argument("continuous_cycle32", type=Path)
    parser.add_argument("long_window", type=Path)
    parser.add_argument("reference_raw", type=Path)
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.rule.resolve(), {
            "locked_source_microstate": args.baseline.resolve(),
            "microstate_51949": args.microstate_51949.resolve(),
            "microstate_63059": args.microstate_63059.resolve(),
        }, {
            "cycles_20_24": args.continuous_cycle24.resolve(),
            "cycles_24_28": args.continuous_cycle28.resolve(),
            "cycles_28_32": args.continuous_cycle32.resolve(),
        }, args.long_window.resolve(), args.reference_raw.resolve(),
        args.campaign_root.resolve())
    if args.output:
        atomic_json(args.output, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
