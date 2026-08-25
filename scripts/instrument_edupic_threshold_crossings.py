#!/usr/bin/env python3
"""Add a passive unconditional energetic-threshold ledger to native eduPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_phase_history import instrument as instrument_history
from instrument_edupic_surface_flux import replace_once


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_history(source)
    source = replace_once(
        source,
        "Ullong phase_history_born[N_XT][PHASE_EEDF_REGIONS] = {{0}};\n",
        """Ullong phase_history_born[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_electron_time[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_energetic_time[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_interstep_promotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_interstep_demotions[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_collision_promotions[N_XT][PHASE_EEDF_REGIONS][3] = {{{0}}};
Ullong threshold_collision_demotions[N_XT][PHASE_EEDF_REGIONS][3] = {{{0}}};
Ullong threshold_energetic_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong threshold_subthreshold_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};
void record_threshold_birth(double x, int phase, double energy_eV);
""",
        "threshold-crossing storage",
    )
    source = replace_once(
        source,
        "void update_electron_history(int particle_id, double energy_eV){\n",
        """void update_electron_history(int particle_id, double x,
                             double energy_eV, int phase){
""",
        "threshold update signature",
    )
    source = replace_once(
        source,
        """    ElectronHistory &history = electron_history[particle_id];
    history.age_steps++;
    const bool energetic = energy_eV >= 15.8;
""",
        """    ElectronHistory &history = electron_history[particle_id];
    const bool energetic = energy_eV >= 15.8;
    const bool has_previous = history.age_steps > 0;
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        threshold_electron_time[phase][region]++;
        if (energetic) threshold_energetic_time[phase][region]++;
        if (has_previous && !history.energetic_previous_step && energetic)
            threshold_interstep_promotions[phase][region]++;
        else if (has_previous && history.energetic_previous_step && !energetic)
            threshold_interstep_demotions[phase][region]++;
    }
    history.age_steps++;
""",
        "unconditional threshold update",
    )
    helpers = r'''
void record_threshold_collision(double x, int phase, int collision_type,
                                bool energetic_before,
                                bool energetic_after){
    if (energetic_before == energetic_after) return;
    int process = collision_type == E_ELA ? 0 :
                  (collision_type == E_EXC ? 1 : 2);
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        if (energetic_after)
            threshold_collision_promotions[phase][region][process]++;
        else threshold_collision_demotions[phase][region][process]++;
    }
}

void record_threshold_birth(double x, int phase, double energy_eV){
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        if (energy_eV >= 15.8) threshold_energetic_births[phase][region]++;
        else threshold_subthreshold_births[phase][region]++;
    }
}

'''
    source = replace_once(
        source,
        "void measure_phase_eedf(double x, double vx, double vy, double vz,\n",
        helpers +
        "void measure_phase_eedf(double x, double vx, double vy, double vz,\n",
        "threshold event helpers",
    )
    output = r'''
void save_threshold_crossings(void){
    FILE *output = fopen("edupic_phase_eedf_threshold_crossings.csv", "w");
    fprintf(output, "phase_bin,phase_fraction,region_id,region,x_min,x_max,electron_time_macro_observations,energetic_time_macro_observations,energetic_fraction,interstep_promotions,interstep_demotions,interstep_promotions_per_million_electron_steps,interstep_demotions_per_million_electron_steps,elastic_collision_promotions,elastic_collision_demotions,excitation_collision_promotions,excitation_collision_demotions,ionization_collision_promotions,ionization_collision_demotions,charge_exchange_collision_promotions,charge_exchange_collision_demotions,attachment_collision_promotions,attachment_collision_demotions,bgk_collision_promotions,bgk_collision_demotions,energetic_births,subthreshold_births\n");
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const double observations = static_cast<double>(
                threshold_electron_time[phase][region]);
            fprintf(output,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%llu,%.17g,%llu,%llu,%.17g,%.17g,%llu,%llu,%llu,%llu,%llu,%llu,0,0,0,0,0,0,%llu,%llu\n",
                phase, phase_fraction, region,
                phase_eedf_region_name[region],
                phase_eedf_region_min[region],
                phase_eedf_region_max[region],
                threshold_electron_time[phase][region],
                threshold_energetic_time[phase][region],
                observations > 0.0 ?
                    static_cast<double>(threshold_energetic_time[phase][region]) /
                    observations : 0.0,
                threshold_interstep_promotions[phase][region],
                threshold_interstep_demotions[phase][region],
                observations > 0.0 ? 1.0e6 *
                    static_cast<double>(threshold_interstep_promotions[phase][region]) /
                    observations : 0.0,
                observations > 0.0 ? 1.0e6 *
                    static_cast<double>(threshold_interstep_demotions[phase][region]) /
                    observations : 0.0,
                threshold_collision_promotions[phase][region][0],
                threshold_collision_demotions[phase][region][0],
                threshold_collision_promotions[phase][region][1],
                threshold_collision_demotions[phase][region][1],
                threshold_collision_promotions[phase][region][2],
                threshold_collision_demotions[phase][region][2],
                threshold_energetic_births[phase][region],
                threshold_subthreshold_births[phase][region]);
        }
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void save_phase_history(void){\n",
        output + "void save_phase_history(void){\n",
        "threshold output",
    )
    source = replace_once(
        source,
        "update_electron_history(k, energy);",
        "update_electron_history(k, x_e[k], energy, t_index);",
        "threshold update call",
    )
    source = replace_once(
        source,
        "int collision_electron (double xe, double *vxe, double *vye, double *vze, int eindex){",
        "int collision_electron (double xe, double *vxe, double *vye, double *vze, int eindex, int phase){",
        "collision phase argument",
    )
    source = replace_once(
        source,
        """        vz_e[N_e] = wz + F2 * gz;
        N_e++;
""",
        """        vz_e[N_e] = wz + F2 * gz;
        record_threshold_birth(xe, phase, 0.5 * E_MASS *
            (vx_e[N_e]*vx_e[N_e] + vy_e[N_e]*vy_e[N_e] +
             vz_e[N_e]*vz_e[N_e]) / EV_TO_J);
        N_e++;
""",
        "secondary threshold birth",
    )
    source = replace_once(
        source,
        "    bool     out;\n",
        "    bool     out, energetic_before, energetic_after;\n",
        "threshold locals",
    )
    source = replace_once(
        source,
        """                collision_type = collision_electron(
                    x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index);
""",
        """                energetic_before = energy >= 15.8;
                collision_type = collision_electron(
                    x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index,
                    t_index);
                energetic_after = 0.5 * E_MASS *
                    (vx_e[k]*vx_e[k] + vy_e[k]*vy_e[k] + vz_e[k]*vz_e[k]) /
                    EV_TO_J >= 15.8;
                record_threshold_collision(
                    x_e[k], t_index, collision_type,
                    energetic_before, energetic_after);
""",
        "collision threshold call",
    )
    source = replace_once(
        source,
        "        save_phase_history();\n",
        "        save_phase_history();\n        save_threshold_crossings();\n",
        "threshold output call",
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
        "transform": "passive_edupic_c_threshold_crossing_ledger_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "history_update_every_timesteps": 1,
        "tail_threshold_eV": 15.8,
        "interstep_definition":
            "consecutive synchronized pre-collision energetic states",
        "collision_definition": "accepted collision before/after state",
        "unconditional_denominator": "live electron timesteps",
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
