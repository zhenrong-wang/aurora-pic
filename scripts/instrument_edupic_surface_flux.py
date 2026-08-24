#!/usr/bin/env python3
"""Instrument pinned eduPIC C implementation with passive electron flux."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_SOURCE_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41"
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def replace_once(source: str, old: str, new: str, label: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one source anchor, found {count}")
    return source.replace(old, new, 1)


def instrument(source: str) -> str:
    source = replace_once(
        source,
        """std::random_device rd{}; 
std::mt19937 MTgen(rd());
""",
        """unsigned int diagnostic_rng_seed(){
    const char* value = std::getenv("EDUPIC_DIAGNOSTIC_SEED");
    if (value == nullptr) return 12345u;
    return static_cast<unsigned int>(std::strtoul(value, nullptr, 10));
}
std::mt19937 MTgen(diagnostic_rng_seed());
""",
        "deterministic continuation seed",
    )
    source = replace_once(
        source,
        "Ullong   N_i_coll                   = 0;                     // counter for ion collisions\n",
        """Ullong   N_i_coll                   = 0;                     // counter for ion collisions

// Passive phase-resolved internal electron surface-flux diagnostic.
const int SURFACE_FLUX_SURFACES = 2;
const int SURFACE_FLUX_DIRECTIONS = 2;
const int SURFACE_FLUX_ENERGY_BINS = 320;
const double SURFACE_FLUX_ENERGY_MAX_EV = 80.0;
const double surface_flux_positions[SURFACE_FLUX_SURFACES] = {0.2 * L, 0.6 * L};
Ullong surface_flux_macro_crossings[N_XT][SURFACE_FLUX_SURFACES][SURFACE_FLUX_DIRECTIONS] = {{{0}}};
Ullong surface_flux_overflow_macro_crossings[N_XT][SURFACE_FLUX_SURFACES][SURFACE_FLUX_DIRECTIONS] = {{{0}}};
double surface_flux_kinetic_energy_joule[N_XT][SURFACE_FLUX_SURFACES][SURFACE_FLUX_DIRECTIONS] = {{{0.0}}};
double surface_flux_histogram_macro[N_XT][SURFACE_FLUX_SURFACES][SURFACE_FLUX_DIRECTIONS][SURFACE_FLUX_ENERGY_BINS] = {{{{0.0}}}};
""",
        "surface accumulator insertion",
    )
    diagnostic = r'''
//-------------------------------------------------------------------------//
// passive phase-resolved internal electron surface-flux diagnostic       //
//-------------------------------------------------------------------------//

void measure_electron_surface_crossings(double old_x, double new_x,
                                        double vx, double vy, double vz,
                                        int phase){
    const double energy_joule = 0.5 * E_MASS * (vx*vx + vy*vy + vz*vz);
    const double energy_eV = energy_joule / EV_TO_J;
    for (int surface=0; surface<SURFACE_FLUX_SURFACES; ++surface) {
        int direction = SURFACE_FLUX_DIRECTIONS;
        if (old_x < surface_flux_positions[surface] &&
            new_x >= surface_flux_positions[surface]) direction = 0;
        else if (old_x > surface_flux_positions[surface] &&
                 new_x <= surface_flux_positions[surface]) direction = 1;
        if (direction == SURFACE_FLUX_DIRECTIONS) continue;
        surface_flux_macro_crossings[phase][surface][direction]++;
        surface_flux_kinetic_energy_joule[phase][surface][direction] += energy_joule;
        if (energy_eV >= SURFACE_FLUX_ENERGY_MAX_EV) {
            surface_flux_overflow_macro_crossings[phase][surface][direction]++;
        } else {
            const int bin = min(
                static_cast<int>(energy_eV /
                    (SURFACE_FLUX_ENERGY_MAX_EV /
                     static_cast<double>(SURFACE_FLUX_ENERGY_BINS))),
                SURFACE_FLUX_ENERGY_BINS - 1);
            surface_flux_histogram_macro[phase][surface][direction][bin] += 1.0;
        }
    }
}

void save_surface_flux(void){
    FILE *histogram = fopen("edupic_phase_surface_flux.csv", "w");
    FILE *summary = fopen("edupic_phase_surface_flux_summary.csv", "w");
    fprintf(histogram, "phase_bin,phase_fraction,surface_id,position_m,direction,energy_bin,energy_eV,represented_crossings_m-2,probability_density_eV-1\n");
    fprintf(summary, "phase_bin,phase_fraction,surface_id,position_m,direction,macro_crossings,overflow_macro_crossings,represented_crossings_m-2,overflow_fraction,represented_particle_flux_m-2_s-1,kinetic_energy_flux_W_m-2\n");
    const double bin_width = SURFACE_FLUX_ENERGY_MAX_EV /
        static_cast<double>(SURFACE_FLUX_ENERGY_BINS);
    const double represented_per_macro = WEIGHT / ELECTRODE_AREA;
    const double phase_duration = static_cast<double>(no_of_cycles) *
        PERIOD / static_cast<double>(N_XT);
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int surface=0; surface<SURFACE_FLUX_SURFACES; ++surface) {
            for (int direction=0; direction<SURFACE_FLUX_DIRECTIONS; ++direction) {
                const double represented = static_cast<double>(
                    surface_flux_macro_crossings[phase][surface][direction]) *
                    represented_per_macro;
                const char* direction_name =
                    direction == 0 ? "left_to_right" : "right_to_left";
                for (int bin=0; bin<SURFACE_FLUX_ENERGY_BINS; ++bin) {
                    const double represented_bin =
                        surface_flux_histogram_macro[phase][surface][direction][bin] *
                        represented_per_macro;
                    fprintf(histogram,
                        "%d,%.17g,%d,%.17g,%s,%d,%.17g,%.17g,%.17g\n",
                        phase, phase_fraction, surface,
                        surface_flux_positions[surface], direction_name, bin,
                        (static_cast<double>(bin) + 0.5) * bin_width,
                        represented_bin, represented > 0.0
                            ? represented_bin / represented / bin_width : 0.0);
                }
                const Ullong macro =
                    surface_flux_macro_crossings[phase][surface][direction];
                const Ullong overflow =
                    surface_flux_overflow_macro_crossings[phase][surface][direction];
                fprintf(summary,
                    "%d,%.17g,%d,%.17g,%s,%llu,%llu,%.17g,%.17g,%.17g,%.17g\n",
                    phase, phase_fraction, surface,
                    surface_flux_positions[surface], direction_name, macro,
                    overflow, represented, macro > 0
                        ? static_cast<double>(overflow) / static_cast<double>(macro)
                        : 0.0, represented / phase_duration,
                    surface_flux_kinetic_energy_joule[phase][surface][direction] *
                        represented_per_macro / phase_duration);
            }
        }
    }
    fclose(histogram);
    fclose(summary);
}

'''
    cycle_anchor = """//---------------------------------------------------------------------//
// simulation of one radiofrequency cycle                              //
//---------------------------------------------------------------------//

void do_one_cycle (void){
"""
    source = replace_once(source, cycle_anchor, diagnostic + cycle_anchor,
                          "surface function insertion")
    source = replace_once(
        source,
        "    double   g, g_sqr, gx, gy, gz, vx_a, vy_a, vz_a, e_x, energy, nu, p_coll, v_sqr, velocity;\n",
        "    double   g, g_sqr, gx, gy, gz, vx_a, vy_a, vz_a, e_x, energy, nu, p_coll, v_sqr, velocity, old_x;\n",
        "old position declaration",
    )
    source = replace_once(
        source,
        """            vx_e[k] -= e_x * FACTOR_E;
            x_e[k]  += vx_e[k] * DT_E;
""",
        """            old_x = x_e[k];
            vx_e[k] -= e_x * FACTOR_E;
            x_e[k]  += vx_e[k] * DT_E;
            if (measurement_mode) {
                measure_electron_surface_crossings(
                    old_x, x_e[k], vx_e[k], vy_e[k], vz_e[k], t_index);
            }
""",
        "surface crossing call",
    )
    source = replace_once(
        source,
        """    if (measurement_mode) {
        check_and_save_info();
    }
""",
        """    if (measurement_mode) {
        check_and_save_info();
        save_surface_flux();
    }
""",
        "surface output call",
    )
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if not args.source.is_file():
        parser.error(f"source does not exist: {args.source}")
    source_hash = sha256(args.source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise RuntimeError(
            f"source hash {source_hash} differs from pinned {EXPECTED_SOURCE_SHA256}"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(instrument(args.source.read_text(encoding="utf-8")),
                           encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 2,
        "transform": "passive_edupic_c_phase_surface_flux_v2",
        "source_sha256": source_hash,
        "instrumented_source_sha256": sha256(args.output),
        "surfaces_gap_fraction": [0.2, 0.6],
        "phase_bins": 200,
        "energy_bins": 320,
        "energy_max_eV": 80.0,
        "crossing_velocity": "post_kick_leapfrog_drift_velocity",
        "rng_initialization": "EDUPIC_DIAGNOSTIC_SEED environment variable; default 12345",
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
        "checkpoint_format": "native C implementation double-scalar header",
    }
    if args.audit:
        args.audit.parent.mkdir(parents=True, exist_ok=True)
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
