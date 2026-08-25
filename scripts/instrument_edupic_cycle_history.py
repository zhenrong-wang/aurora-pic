#!/usr/bin/env python3
"""Add passive cycle-resolved state history to pinned native eduPIC."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_surface_flux import replace_once
from instrument_edupic_threshold_crossings import (
    instrument as instrument_threshold_crossings,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_threshold_crossings(source)
    source = replace_once(
        source,
        "Ullong threshold_subthreshold_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};\n",
        """Ullong threshold_subthreshold_births[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong cycle_electron_elastic_collisions = 0;
Ullong cycle_electron_excitation_collisions = 0;
Ullong cycle_electron_ionization_collisions = 0;
""",
        "cycle collision counters",
    )
    source = replace_once(
        source,
        """                else if (collision_type == E_ION) electron_history[k].ionization_collisions++;
                N_e_coll++;
""",
        """                else if (collision_type == E_ION) electron_history[k].ionization_collisions++;
                if (collision_type == E_ELA) cycle_electron_elastic_collisions++;
                else if (collision_type == E_EXC) cycle_electron_excitation_collisions++;
                else if (collision_type == E_ION) cycle_electron_ionization_collisions++;
                N_e_coll++;
""",
        "cycle collision accounting",
    )
    diagnostic = r'''
void save_cycle_history(void){
    double electron_kinetic_J = 0.0;
    double ion_kinetic_J = 0.0;
    double energetic_kinetic_J = 0.0;
    Ullong energetic_electrons = 0;
    for (int particle=0; particle<N_e; ++particle) {
        const double v2 = vx_e[particle]*vx_e[particle] +
            vy_e[particle]*vy_e[particle] + vz_e[particle]*vz_e[particle];
        const double kinetic = 0.5 * E_MASS * v2;
        electron_kinetic_J += kinetic;
        if (kinetic / EV_TO_J >= E_ION_TH) {
            energetic_electrons++;
            energetic_kinetic_J += kinetic;
        }
    }
    for (int particle=0; particle<N_i; ++particle) {
        const double v2 = vx_i[particle]*vx_i[particle] +
            vy_i[particle]*vy_i[particle] + vz_i[particle]*vz_i[particle];
        ion_kinetic_J += 0.5 * AR_MASS * v2;
    }
    double field_energy_J = 0.0;
    double charge_l1_C = 0.0;
    double maximum_absolute_field_V_m = 0.0;
    for (int node=0; node<N_G; ++node) {
        const double quadrature =
            (node == 0 || node == N_G-1) ? 0.5 : 1.0;
        field_energy_J += quadrature * 0.5 * EPSILON0 *
            efield[node] * efield[node] * ELECTRODE_AREA * DX;
        charge_l1_C += quadrature * fabs(
            E_CHARGE * (i_density[node] - e_density[node])) *
            ELECTRODE_AREA * DX;
        maximum_absolute_field_V_m = max(
            maximum_absolute_field_V_m, fabs(efield[node]));
    }
    const bool first_cycle = cycle == cycles_done + 1;
    FILE *output = fopen(
        "edupic_cycle_history.csv", first_cycle ? "w" : "a");
    if (first_cycle) {
        fprintf(output,
            "cycle,time_s,electrons,ions,energetic_electrons,energetic_fraction,electron_kinetic_J,ion_kinetic_J,energetic_electron_kinetic_J,field_energy_J,charge_l1_C,maximum_absolute_field_V_m,electron_collisions,electron_elastic_collisions,electron_excitation_collisions,electron_ionization_collisions,ion_collisions,electron_absorbed_powered,electron_absorbed_grounded,ion_absorbed_powered,ion_absorbed_grounded,potential_powered_V,potential_grounded_V\n");
    }
    fprintf(output,
        "%d,%.17g,%d,%d,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%.17g,%.17g\n",
        cycle, Time, N_e, N_i, energetic_electrons,
        N_e > 0 ? static_cast<double>(energetic_electrons) /
            static_cast<double>(N_e) : 0.0,
        electron_kinetic_J, ion_kinetic_J, energetic_kinetic_J,
        field_energy_J, charge_l1_C, maximum_absolute_field_V_m,
        N_e_coll, cycle_electron_elastic_collisions,
        cycle_electron_excitation_collisions,
        cycle_electron_ionization_collisions, N_i_coll,
        N_e_abs_pow, N_e_abs_gnd, N_i_abs_pow, N_i_abs_gnd,
        pot[0], pot[N_G-1]);
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void save_threshold_crossings(void){\n",
        diagnostic + "void save_threshold_crossings(void){\n",
        "cycle history function",
    )
    source = replace_once(
        source,
        '    fprintf(datafile,"%8d  %8d  %8d\\n",cycle,N_e,N_i);\n',
        """    fprintf(datafile,"%8d  %8d  %8d\\n",cycle,N_e,N_i);
    if (measurement_mode) save_cycle_history();
""",
        "cycle history call",
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
        "transform": "passive_edupic_cycle_state_history_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "cycle_sample_phase": "after_boundary_and_collision_at_rf_cycle_end",
        "energetic_threshold_eV": 15.8,
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
        "composes_threshold_crossing_diagnostic": True,
    }
    if args.audit:
        args.audit.write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
