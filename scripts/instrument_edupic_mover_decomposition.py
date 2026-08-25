#!/usr/bin/env python3
"""Add passive near-threshold mover decomposition to pinned eduPIC."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_promotion_band_work import instrument as instrument_work
from instrument_edupic_promotion_band_work import sha256
from instrument_edupic_surface_flux import replace_once


def instrument(source: str) -> str:
    source = instrument_work(source)
    source = replace_once(
        source,
        "double promotion_band_negative_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};\n",
        """double promotion_band_negative_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_origin_energy[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_origin_longitudinal_energy[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_linear_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_positive_linear_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_negative_linear_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double mover_band_quadratic_work[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
""",
        "mover-decomposition storage",
    )
    helpers = r'''
void record_mover_decomposition(double x, int phase,
                                double origin_energy_eV,
                                double origin_vx,
                                double delta_vx){
    if (origin_energy_eV < 11.5 || origin_energy_eV >= 15.8) return;
    const double origin_longitudinal_energy_eV =
        0.5 * E_MASS * origin_vx * origin_vx / EV_TO_J;
    const double linear_work_eV =
        E_MASS * origin_vx * delta_vx / EV_TO_J;
    const double quadratic_work_eV =
        0.5 * E_MASS * delta_vx * delta_vx / EV_TO_J;
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        mover_band_origin_energy[phase][region] += origin_energy_eV;
        mover_band_origin_longitudinal_energy[phase][region] +=
            origin_longitudinal_energy_eV;
        mover_band_linear_work[phase][region] += linear_work_eV;
        mover_band_positive_linear_work[phase][region] +=
            max(0.0, linear_work_eV);
        mover_band_negative_linear_work[phase][region] +=
            max(0.0, -linear_work_eV);
        mover_band_quadratic_work[phase][region] += quadratic_work_eV;
    }
}

void save_mover_decomposition(void){
    FILE *output = fopen("edupic_phase_eedf_mover_decomposition.csv", "w");
    fprintf(output, "phase_bin,phase_fraction,region_id,region,x_min,x_max,field_push_promotion_band_observations,field_push_promotion_band_origin_macro_energy_sum_eV,field_push_promotion_band_origin_longitudinal_macro_energy_sum_eV,field_push_promotion_band_linear_macro_work_sum_eV,field_push_promotion_band_positive_linear_macro_work_sum_eV,field_push_promotion_band_negative_linear_macro_work_sum_eV,field_push_promotion_band_quadratic_macro_work_sum_eV,field_push_promotion_band_mean_origin_energy_eV,field_push_promotion_band_origin_longitudinal_energy_fraction,field_push_promotion_band_mean_linear_work_eV,field_push_promotion_band_mean_positive_linear_work_eV,field_push_promotion_band_mean_quadratic_work_eV\n");
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const double observations = static_cast<double>(
                promotion_band_observations[phase][region]);
            const double origin_energy =
                mover_band_origin_energy[phase][region];
            fprintf(output,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                phase, phase_fraction, region,
                phase_eedf_region_name[region],
                phase_eedf_region_min[region],
                phase_eedf_region_max[region],
                promotion_band_observations[phase][region],
                origin_energy,
                mover_band_origin_longitudinal_energy[phase][region],
                mover_band_linear_work[phase][region],
                mover_band_positive_linear_work[phase][region],
                mover_band_negative_linear_work[phase][region],
                mover_band_quadratic_work[phase][region],
                observations > 0.0 ? origin_energy / observations : 0.0,
                origin_energy > 0.0 ?
                    mover_band_origin_longitudinal_energy[phase][region] /
                    origin_energy : 0.0,
                observations > 0.0 ?
                    mover_band_linear_work[phase][region] / observations : 0.0,
                observations > 0.0 ?
                    mover_band_positive_linear_work[phase][region] /
                    observations : 0.0,
                observations > 0.0 ?
                    mover_band_quadratic_work[phase][region] /
                    observations : 0.0);
        }
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void save_promotion_band_work(void){\n",
        helpers + "void save_promotion_band_work(void){\n",
        "mover-decomposition helpers",
    )
    source = replace_once(
        source,
        "    double   field_push_energy_before, field_push_energy_after;\n",
        "    double   field_push_energy_before, field_push_energy_after;\n"
        "    double   field_push_vx_before;\n",
        "mover-decomposition local",
    )
    source = replace_once(
        source,
        """            field_push_energetic_before = field_push_energy_before >= 15.8;
            vx_e[k] -= e_x * FACTOR_E;
""",
        """            field_push_energetic_before = field_push_energy_before >= 15.8;
            field_push_vx_before = vx_e[k];
            vx_e[k] -= e_x * FACTOR_E;
""",
        "mover-decomposition origin velocity",
    )
    source = replace_once(
        source,
        """                record_promotion_band_work(
                    x_e[k], t_index, field_push_energy_before,
                    field_push_energy_after);
""",
        """                record_promotion_band_work(
                    x_e[k], t_index, field_push_energy_before,
                    field_push_energy_after);
                record_mover_decomposition(
                    x_e[k], t_index, field_push_energy_before,
                    field_push_vx_before, vx_e[k] - field_push_vx_before);
""",
        "mover-decomposition hook",
    )
    source = replace_once(
        source,
        "        save_promotion_band_work();\n",
        "        save_promotion_band_work();\n"
        "        save_mover_decomposition();\n",
        "mover-decomposition output call",
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
        "transform": "passive_edupic_c_near_threshold_mover_decomposition_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "promotion_band_eV": [11.5, 15.8],
        "identity": "delta_K = m*v_x*delta_v_x + 0.5*m*delta_v_x^2",
        "quadratic_term_interpretation": "particle_sampled_field_strength_proportional_to_E_squared",
        "linear_term_interpretation": "longitudinal_velocity_field_alignment",
        "spatial_attribution": "post_drift_position_for_surviving_particles",
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
