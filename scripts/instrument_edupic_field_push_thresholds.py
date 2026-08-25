#!/usr/bin/env python3
"""Add passive field-push threshold attribution to pinned native eduPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_surface_flux import replace_once
from instrument_edupic_threshold_crossings import instrument as instrument_thresholds


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_thresholds(source)
    source = replace_once(
        source,
        "Ullong threshold_subthreshold_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};\n",
        """Ullong threshold_subthreshold_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong field_push_observations[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong field_push_promotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong field_push_demotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
""",
        "field-push storage",
    )
    helpers = r'''
void record_field_push_threshold(double x, int phase,
                                 bool energetic_before,
                                 bool energetic_after){
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        field_push_observations[phase][region]++;
        if (!energetic_before && energetic_after)
            field_push_promotions[phase][region]++;
        else if (energetic_before && !energetic_after)
            field_push_demotions[phase][region]++;
    }
}

void save_field_push_thresholds(void){
    FILE *output = fopen("edupic_phase_eedf_field_push_thresholds.csv", "w");
    fprintf(output, "phase_bin,phase_fraction,region_id,region,x_min,x_max,field_push_macro_observations,field_push_promotions,field_push_demotions,field_push_promotions_per_million_pushes,field_push_demotions_per_million_pushes\n");
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const double observations = static_cast<double>(
                field_push_observations[phase][region]);
            fprintf(output,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%llu,%llu,%.17g,%.17g\n",
                phase, phase_fraction, region,
                phase_eedf_region_name[region],
                phase_eedf_region_min[region],
                phase_eedf_region_max[region],
                field_push_observations[phase][region],
                field_push_promotions[phase][region],
                field_push_demotions[phase][region],
                observations > 0.0 ? 1.0e6 *
                    static_cast<double>(field_push_promotions[phase][region]) /
                    observations : 0.0,
                observations > 0.0 ? 1.0e6 *
                    static_cast<double>(field_push_demotions[phase][region]) /
                    observations : 0.0);
        }
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void save_threshold_crossings(void){\n",
        helpers + "void save_threshold_crossings(void){\n",
        "field-push helpers",
    )
    source = replace_once(
        source,
        "    bool     out, energetic_before, energetic_after;\n",
        "    bool     out, energetic_before, energetic_after, field_push_energetic_before;\n",
        "field-push local",
    )
    source = replace_once(
        source,
        """            vx_e[k] -= e_x * FACTOR_E;
            x_e[k]  += vx_e[k] * DT_E;
""",
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
        "electron field-push hook",
    )
    source = replace_once(
        source,
        "        save_threshold_crossings();\n",
        "        save_threshold_crossings();\n        save_field_push_thresholds();\n",
        "field-push output call",
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
    args.output.write_text(instrument(args.source.read_text(encoding="utf-8")),
                           encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 1,
        "transform": "passive_edupic_c_field_push_threshold_ledger_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "tail_threshold_eV": 15.8,
        "before_state": "post_collision_leapfrog_velocity_entering_mover",
        "after_state": "pre_collision_leapfrog_velocity_leaving_mover",
        "spatial_attribution": "post_drift_position_for_surviving_particles",
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
    }
    if args.audit:
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
