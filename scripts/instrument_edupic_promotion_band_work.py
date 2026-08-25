#!/usr/bin/env python3
"""Add passive near-threshold field-work attribution to pinned eduPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_field_push_thresholds import instrument as instrument_push
from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_surface_flux import replace_once


PROMOTION_BAND_MIN_EV = 11.5
TAIL_THRESHOLD_EV = 15.8


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_push(source)
    source = replace_once(
        source,
        "Ullong field_push_demotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};\n",
        """Ullong field_push_demotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong promotion_band_observations[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong promotion_band_promotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
double promotion_band_signed_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double promotion_band_positive_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double promotion_band_negative_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
""",
        "promotion-band storage",
    )
    helpers = r'''
void record_promotion_band_work(double x, int phase,
                                double energy_before_eV,
                                double energy_after_eV){
    if (energy_before_eV < 11.5 || energy_before_eV >= 15.8) return;
    const double work_eV = energy_after_eV - energy_before_eV;
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        promotion_band_observations[phase][region]++;
        if (energy_after_eV >= 15.8)
            promotion_band_promotions[phase][region]++;
        promotion_band_signed_work[phase][region] += work_eV;
        promotion_band_positive_work[phase][region] += max(0.0, work_eV);
        promotion_band_negative_work[phase][region] += max(0.0, -work_eV);
    }
}

void save_promotion_band_work(void){
    FILE *output = fopen("edupic_phase_eedf_promotion_band_work.csv", "w");
    fprintf(output, "phase_bin,phase_fraction,region_id,region,x_min,x_max,field_push_promotion_band_observations,field_push_promotion_band_promotions,field_push_promotion_band_promotion_fraction,field_push_promotion_band_signed_macro_work_sum_eV,field_push_promotion_band_positive_macro_work_sum_eV,field_push_promotion_band_negative_macro_work_sum_eV,field_push_promotion_band_mean_signed_work_eV\n");
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const double observations = static_cast<double>(
                promotion_band_observations[phase][region]);
            fprintf(output,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%llu,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                phase, phase_fraction, region,
                phase_eedf_region_name[region],
                phase_eedf_region_min[region],
                phase_eedf_region_max[region],
                promotion_band_observations[phase][region],
                promotion_band_promotions[phase][region],
                observations > 0.0 ? static_cast<double>(
                    promotion_band_promotions[phase][region]) /
                    observations : 0.0,
                promotion_band_signed_work[phase][region],
                promotion_band_positive_work[phase][region],
                promotion_band_negative_work[phase][region],
                observations > 0.0 ?
                    promotion_band_signed_work[phase][region] /
                    observations : 0.0);
        }
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void save_field_push_thresholds(void){\n",
        helpers + "void save_field_push_thresholds(void){\n",
        "promotion-band helpers",
    )
    source = replace_once(
        source,
        "    bool     out, energetic_before, energetic_after, field_push_energetic_before;\n",
        "    bool     out, energetic_before, energetic_after, field_push_energetic_before;\n"
        "    double   field_push_energy_before, field_push_energy_after;\n",
        "promotion-band locals",
    )
    source = replace_once(
        source,
        """            field_push_energetic_before = 0.5 * E_MASS *
                (vx_e[k]*vx_e[k] + vy_e[k]*vy_e[k] + vz_e[k]*vz_e[k]) /
                EV_TO_J >= 15.8;
            vx_e[k] -= e_x * FACTOR_E;
            x_e[k]  += vx_e[k] * DT_E;
            if (measurement_mode && x_e[k] >= 0.0 && x_e[k] <= L) {
                const bool field_push_energetic_after = 0.5 * E_MASS *
                    (vx_e[k]*vx_e[k] + vy_e[k]*vy_e[k] +
                     vz_e[k]*vz_e[k]) / EV_TO_J >= 15.8;
                record_field_push_threshold(
                    x_e[k], t_index, field_push_energetic_before,
                    field_push_energetic_after);
            }
""",
        """            field_push_energy_before = 0.5 * E_MASS *
                (vx_e[k]*vx_e[k] + vy_e[k]*vy_e[k] + vz_e[k]*vz_e[k]) /
                EV_TO_J;
            field_push_energetic_before = field_push_energy_before >= 15.8;
            vx_e[k] -= e_x * FACTOR_E;
            x_e[k]  += vx_e[k] * DT_E;
            if (measurement_mode && x_e[k] >= 0.0 && x_e[k] <= L) {
                field_push_energy_after = 0.5 * E_MASS *
                    (vx_e[k]*vx_e[k] + vy_e[k]*vy_e[k] +
                     vz_e[k]*vz_e[k]) / EV_TO_J;
                const bool field_push_energetic_after =
                    field_push_energy_after >= 15.8;
                record_field_push_threshold(
                    x_e[k], t_index, field_push_energetic_before,
                    field_push_energetic_after);
                record_promotion_band_work(
                    x_e[k], t_index, field_push_energy_before,
                    field_push_energy_after);
            }
""",
        "promotion-band mover hook",
    )
    source = replace_once(
        source,
        "        save_field_push_thresholds();\n",
        "        save_field_push_thresholds();\n"
        "        save_promotion_band_work();\n",
        "promotion-band output call",
    )
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if sha256(args.source) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("eduPIC source differs from pinned C implementation")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        instrument(args.source.read_text(encoding="utf-8")),
        encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 1,
        "transform": "passive_edupic_c_promotion_band_field_work_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "promotion_band_eV": [PROMOTION_BAND_MIN_EV, TAIL_THRESHOLD_EV],
        "before_state": "post_collision_leapfrog_velocity_entering_mover",
        "after_state": "pre_collision_leapfrog_velocity_leaving_mover",
        "spatial_attribution": "post_drift_position_for_surviving_particles",
        "work_definition": "after_total_kinetic_energy_minus_before_total_kinetic_energy_eV",
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
    }
    if args.audit:
        args.audit.write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
