#!/usr/bin/env python3
"""Add passive energetic-electron histories to pinned native eduPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_phase_anisotropy import instrument as instrument_anisotropy
from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_surface_flux import replace_once


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_anisotropy(source)
    source = replace_once(
        source,
        "double phase_eedf_tail_vt2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};\n",
        """double phase_eedf_tail_vt2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};

struct ElectronHistory {
    unsigned int age_steps;
    unsigned int energetic_steps;
    unsigned int consecutive_energetic_steps;
    unsigned int tail_entries;
    unsigned int elastic_collisions;
    unsigned int excitation_collisions;
    unsigned int ionization_collisions;
    bool born_during_window;
    bool energetic_previous_step;
    ElectronHistory() : age_steps(0), energetic_steps(0),
        consecutive_energetic_steps(0), tail_entries(0),
        elastic_collisions(0), excitation_collisions(0),
        ionization_collisions(0), born_during_window(false),
        energetic_previous_step(false) {}
};
ElectronHistory electron_history[MAX_N_P];
double phase_history_age_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_energetic_steps_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_duty_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_consecutive_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_entries_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_elastic_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_excitation_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_history_ionization_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
Ullong phase_history_born[N_XT][PHASE_EEDF_REGIONS] = {{0}};
""",
        "electron-history storage",
    )
    source = replace_once(
        source,
        "void measure_phase_eedf(double x, double vx, double vy, double vz, double energy_eV, int phase){",
        """void update_electron_history(int particle_id, double energy_eV){
    ElectronHistory &history = electron_history[particle_id];
    history.age_steps++;
    const bool energetic = energy_eV >= 15.8;
    if (energetic) {
        history.energetic_steps++;
        history.consecutive_energetic_steps++;
        if (!history.energetic_previous_step) history.tail_entries++;
    } else {
        history.consecutive_energetic_steps = 0;
    }
    history.energetic_previous_step = energetic;
}

void measure_phase_eedf(double x, double vx, double vy, double vz,
                        double energy_eV, int phase, int particle_id){""",
        "history update and measurement signature",
    )
    source = replace_once(
        source,
        """            phase_eedf_tail_vt2_sum[phase][region] += vy * vy + vz * vz;
        }
""",
        """            phase_eedf_tail_vt2_sum[phase][region] += vy * vy + vz * vz;
            const ElectronHistory &history = electron_history[particle_id];
            phase_history_age_sum[phase][region] += history.age_steps;
            phase_history_energetic_steps_sum[phase][region] += history.energetic_steps;
            phase_history_duty_sum[phase][region] += history.age_steps > 0
                ? static_cast<double>(history.energetic_steps) /
                  static_cast<double>(history.age_steps) : 0.0;
            phase_history_consecutive_sum[phase][region] +=
                history.consecutive_energetic_steps;
            phase_history_entries_sum[phase][region] += history.tail_entries;
            phase_history_elastic_sum[phase][region] += history.elastic_collisions;
            phase_history_excitation_sum[phase][region] += history.excitation_collisions;
            phase_history_ionization_sum[phase][region] += history.ionization_collisions;
            if (history.born_during_window) phase_history_born[phase][region]++;
        }
""",
        "history tail accumulation",
    )
    history_output = r'''
void save_phase_history(void){
    FILE *output = fopen("edupic_phase_eedf_history.csv", "w");
    fprintf(output, "phase_bin,phase_fraction,region_id,region,x_min,x_max,tail_macro_observations,tail_mean_age_steps,tail_mean_energetic_steps,tail_mean_energetic_duty_fraction,tail_mean_consecutive_energetic_steps,tail_mean_entries,tail_mean_elastic_collisions,tail_mean_excitation_collisions,tail_mean_ionization_collisions,tail_mean_charge_exchange_collisions,tail_mean_bgk_collisions,tail_born_during_window_fraction\n");
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const double tail = static_cast<double>(
                phase_eedf_tail_macro[phase][region]);
            fprintf(output,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,0,0,%.17g\n",
                phase, phase_fraction, region,
                phase_eedf_region_name[region],
                phase_eedf_region_min[region],
                phase_eedf_region_max[region],
                phase_eedf_tail_macro[phase][region],
                tail > 0.0 ? phase_history_age_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_energetic_steps_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_duty_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_consecutive_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_entries_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_elastic_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_excitation_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? phase_history_ionization_sum[phase][region] / tail : 0.0,
                tail > 0.0 ? static_cast<double>(phase_history_born[phase][region]) / tail : 0.0);
        }
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "//---------------------------------------------------------------------//\n// simulation of one radiofrequency cycle",
        history_output +
        "//---------------------------------------------------------------------//\n// simulation of one radiofrequency cycle",
        "history output insertion",
    )
    source = replace_once(
        source,
        "if ((t % 2) == 0) measure_phase_eedf(x_e[k], mean_v, vy_e[k], vz_e[k], energy, t_index);",
        """update_electron_history(k, energy);
                if ((t % 2) == 0) measure_phase_eedf(
                    x_e[k], mean_v, vy_e[k], vz_e[k], energy, t_index, k);""",
        "history update call",
    )
    source = replace_once(
        source,
        """                vz_e[k] = vz_e[N_e-1];
                N_e--;
""",
        """                vz_e[k] = vz_e[N_e-1];
                electron_history[k] = electron_history[N_e-1];
                N_e--;
""",
        "history compaction",
    )
    source = replace_once(
        source,
        "void collision_electron (double xe, double *vxe, double *vye, double *vze, int eindex){",
        "int collision_electron (double xe, double *vxe, double *vye, double *vze, int eindex){",
        "collision process return type",
    )
    source = replace_once(
        source,
        "    double t0,t1,t2,rnd;\n",
        "    double t0,t1,t2,rnd;\n    int collision_type;\n",
        "collision process declaration",
    )
    source = replace_once(
        source,
        "if (rnd < (t0/t2)){                              // elastic scattering\n",
        "if (rnd < (t0/t2)){                              // elastic scattering\n        collision_type = E_ELA;\n",
        "elastic collision tag",
    )
    source = replace_once(
        source,
        "} else if (rnd < (t1/t2)){                       // excitation\n",
        "} else if (rnd < (t1/t2)){                       // excitation\n        collision_type = E_EXC;\n",
        "excitation collision tag",
    )
    source = replace_once(
        source,
        "} else {                                         // ionization\n",
        "} else {                                         // ionization\n        collision_type = E_ION;\n",
        "ionization collision tag",
    )
    source = replace_once(
        source,
        "        x_e[N_e]  = xe;                              // add new electron\n",
        """        electron_history[N_e] = ElectronHistory();
        electron_history[N_e].born_during_window = true;
        x_e[N_e]  = xe;                              // add new electron
""",
        "secondary history initialization",
    )
    source = replace_once(
        source,
        """    (*vze) = wz + F2 * gz;
}

//----------------------------------------------------------------------//
// Ar+ / Ar collision""",
        """    (*vze) = wz + F2 * gz;
    return collision_type;
}

//----------------------------------------------------------------------//
// Ar+ / Ar collision""",
        "collision process return",
    )
    source = replace_once(
        source,
        "    int      k, t, p, energy_index;\n",
        "    int      k, t, p, energy_index, collision_type;\n",
        "collision process local",
    )
    source = replace_once(
        source,
        """                collision_electron(x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index);
                N_e_coll++;
""",
        """                collision_type = collision_electron(
                    x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index);
                if (collision_type == E_ELA) electron_history[k].elastic_collisions++;
                else if (collision_type == E_EXC) electron_history[k].excitation_collisions++;
                else if (collision_type == E_ION) electron_history[k].ionization_collisions++;
                N_e_coll++;
""",
        "collision history call",
    )
    source = replace_once(
        source,
        "        save_phase_eedf();\n",
        "        save_phase_eedf();\n        save_phase_history();\n",
        "history output call",
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
        "transform": "passive_edupic_c_regional_phase_history_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "sample_every_timesteps": 2,
        "history_update_every_timesteps": 1,
        "tail_threshold_eV": 15.8,
        "history_origin": "measurement_continuation_start",
        "left_censoring_reported": True,
        "particle_compaction_tags_copied": True,
        "ionization_secondary_tags_reset": True,
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
