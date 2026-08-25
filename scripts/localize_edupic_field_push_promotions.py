#!/usr/bin/env python3
"""Localize the pooled matched-half-step field-push promotion deficit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


class PromotionLocalizationError(RuntimeError):
    pass


SCOPES = {
    "upstream_x010_to_x020_phase_0p125_to_0p5": ((1,), 0.125, 0.5),
    "x020_to_x040_phase_0p125_to_0p5": ((2,), 0.125, 0.5),
    "x040_to_x060_phase_0p125_to_0p5": ((3,), 0.125, 0.5),
    "critical_x020_to_x060_phase_0p125_to_0p25": ((2, 3), 0.125, 0.25),
    "critical_x020_to_x060_phase_0p25_to_0p375": ((2, 3), 0.25, 0.375),
    "critical_x020_to_x060_phase_0p375_to_0p5": ((2, 3), 0.375, 0.5),
    "critical_x020_to_x060_phase_0p125_to_0p5": ((2, 3), 0.125, 0.5),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {
        "phase_bin", "phase_fraction", "region_id",
        "field_push_macro_observations", "field_push_promotions",
        "field_push_demotions",
    }
    if not rows or not required.issubset(rows[0]):
        raise PromotionLocalizationError(f"invalid field-push table: {path}")
    cells = {(int(row["phase_bin"]), int(row["region_id"])) for row in rows}
    if cells != {(phase, region) for phase in range(200) for region in range(7)}:
        raise PromotionLocalizationError(f"phase-region shape differs: {path}")
    return rows


def relative_range(values: list[float]) -> float:
    mean = sum(values) / len(values)
    if mean == 0.0:
        return 0.0 if max(values) == min(values) else math.inf
    return (max(values) - min(values)) / mean


def aggregate(rows: list[dict[str, str]], regions: tuple[int, ...],
              phase_min: float, phase_max: float) -> dict[str, float | int]:
    selected = [
        row for row in rows
        if int(row["region_id"]) in regions
        and phase_min <= float(row["phase_fraction"]) < phase_max
    ]
    observations = sum(
        int(row["field_push_macro_observations"]) for row in selected)
    promotions = sum(int(row["field_push_promotions"]) for row in selected)
    demotions = sum(int(row["field_push_demotions"]) for row in selected)
    if observations <= 0 or promotions > observations or demotions > observations:
        raise PromotionLocalizationError("invalid aggregate counts")
    return {
        "rows": len(selected),
        "field_push_macro_observations": observations,
        "field_push_promotions": promotions,
        "field_push_demotions": demotions,
        "field_push_promotions_per_million_pushes":
            1.0e6 * promotions / observations,
        "field_push_demotions_per_million_pushes":
            1.0e6 * demotions / observations,
    }


def localize(result_path: Path, first_paths: list[Path],
             second_paths: list[Path], native_paths: list[Path]) -> dict[str, object]:
    if len(first_paths) != 2 or len(second_paths) != 2 or len(native_paths) != 3:
        raise PromotionLocalizationError(
            "two first blocks, two second blocks, and three native members required")
    prior = json.loads(result_path.read_text(encoding="utf-8"))
    if prior.get("all_provenance_population_and_repeatability_gates_passed") is not True:
        raise PromotionLocalizationError("pooled source result did not pass gates")
    provenance = prior["provenance"]
    expected = {
        "first": provenance["first_crossings_sha256"],
        "second": provenance["second_crossings_sha256"],
        "native": [
            "6576575c2315664f49e2834671e4d69e7197da2a0d9a2ff68cf82a5bb6a9ab00",
            "1763e851fba7156c9106aa3901ffe953d4f2cc612e26f63a0456dc95fedf8e3e",
            "21e770ed671c486cca963cc41603f1c442670a7bca422d9527efabf2d9e1235c",
        ],
    }
    paths = {"first": first_paths, "second": second_paths, "native": native_paths}
    actual = {name: [sha256(path) for path in values]
              for name, values in paths.items()}
    if actual != expected:
        raise PromotionLocalizationError("one or more locked table hashes differ")

    first = [read_rows(path) for path in first_paths]
    second = [read_rows(path) for path in second_paths]
    native = [read_rows(path) for path in native_paths]
    candidate = [first[index] + second[index] for index in range(2)]
    comparisons: dict[str, object] = {}
    for name, (regions, phase_min, phase_max) in SCOPES.items():
        candidate_values = [
            aggregate(rows, regions, phase_min, phase_max) for rows in candidate]
        native_values = [
            aggregate(rows, regions, phase_min, phase_max) for rows in native]
        promotion = "field_push_promotions_per_million_pushes"
        demotion = "field_push_demotions_per_million_pushes"
        native_mean = {
            promotion: sum(float(value[promotion]) for value in native_values) / 3,
            demotion: sum(float(value[demotion]) for value in native_values) / 3,
        }
        comparisons[name] = {
            "aurorapic_pooled_microstates": candidate_values,
            "native_edupic_members": native_values,
            "native_edupic_ensemble_mean": native_mean,
            "aurorapic_to_native_ratios": [
                {promotion: float(value[promotion]) / native_mean[promotion],
                 demotion: float(value[demotion]) / native_mean[demotion]}
                for value in candidate_values
            ],
            "relative_ranges": {
                "aurorapic_promotions": relative_range(
                    [float(value[promotion]) for value in candidate_values]),
                "aurorapic_demotions": relative_range(
                    [float(value[demotion]) for value in candidate_values]),
                "native_promotions": relative_range(
                    [float(value[promotion]) for value in native_values]),
                "native_demotions": relative_range(
                    [float(value[demotion]) for value in native_values]),
            },
        }

    critical = comparisons[
        "critical_x020_to_x060_phase_0p125_to_0p5"]
    prior_ratios = prior["critical_phase_0p125_to_0p5"][
        "pooled_aurorapic_to_native_ratios"]
    reproduced = all(
        abs(current["field_push_promotions_per_million_pushes"] -
            locked["field_push_promotions_per_million_pushes"]) < 1.0e-12
        and abs(current["field_push_demotions_per_million_pushes"] -
                locked["field_push_demotions_per_million_pushes"]) < 1.0e-12
        for current, locked in zip(
            critical["aurorapic_to_native_ratios"], prior_ratios, strict=True))
    if not reproduced:
        raise PromotionLocalizationError("critical aggregate does not reproduce")

    middle = comparisons[
        "critical_x020_to_x060_phase_0p25_to_0p375"]
    upstream = comparisons[
        "upstream_x010_to_x020_phase_0p125_to_0p5"]
    interior = comparisons[
        "x020_to_x040_phase_0p125_to_0p5"]
    middle_ratios = [value["field_push_promotions_per_million_pushes"]
                     for value in middle["aurorapic_to_native_ratios"]]
    upstream_ratios = [value["field_push_promotions_per_million_pushes"]
                       for value in upstream["aurorapic_to_native_ratios"]]
    interior_ratios = [value["field_push_promotions_per_million_pushes"]
                       for value in interior["aurorapic_to_native_ratios"]]
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "post_hoc_phase_space_localization_of_pooled_matched_half_step_promotions",
        "all_integrity_checks_passed": True,
        "comparisons": comparisons,
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "pooled_result_sha256": sha256(result_path),
            "input_sha256": actual,
        },
        "assessment": {
            "critical_aggregate_exactly_reproduced": reproduced,
            "observation": (
                "The matched-half-step promotion deficit is present in all "
                "three critical RF-phase octants and is strongest at phase "
                "0.25--0.375, where the two AuroraPIC/native ratios are "
                f"{middle_ratios[0]:.4f} and {middle_ratios[1]:.4f}. Over "
                "phase 0.125--0.5 the upstream x/L=0.1--0.2 ratios are "
                f"{upstream_ratios[0]:.4f} and {upstream_ratios[1]:.4f}, "
                "then fall in x/L=0.2--0.4 to "
                f"{interior_ratios[0]:.4f} and {interior_ratios[1]:.4f}. "
                "The field-push traffic discrepancy therefore spans the "
                "critical phase window and amplifies immediately beyond the "
                "upstream region rather than arising in one isolated phase."),
        },
        "claim_boundary": (
            "This post hoc localization reuses checksum-locked ledgers from a "
            "prospectively accepted aggregate comparison. It identifies where "
            "the existing field-push discrepancy is largest, but its spatial "
            "and phase subdivisions are descriptive rather than preregistered "
            "acceptance gates and do not establish causal onset or experimental "
            "validity."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pooled_result", type=Path)
    parser.add_argument("first", nargs=2, type=Path)
    parser.add_argument("second", nargs=2, type=Path)
    parser.add_argument("native", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        result = localize(
            args.pooled_result, args.first, args.second, args.native)
    except (PromotionLocalizationError, OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
