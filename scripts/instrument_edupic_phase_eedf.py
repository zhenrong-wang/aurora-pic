#!/usr/bin/env python3
"""Instrument pinned eduPIC C source with a passive regional phase EEDF."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_surface_flux import replace_once


EXPECTED_SOURCE_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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

// Passive regional phase-EEDF diagnostic, sampled every second step.
const int PHASE_EEDF_REGIONS = 7;
const int PHASE_EEDF_ENERGY_BINS = 320;
const double PHASE_EEDF_ENERGY_MAX_EV = 80.0;
const double phase_eedf_region_min[PHASE_EEDF_REGIONS] =
    {0.0, 0.1*L, 0.2*L, 0.4*L, 0.6*L, 0.8*L, 0.9*L};
const double phase_eedf_region_max[PHASE_EEDF_REGIONS] =
    {0.1*L, 0.2*L, 0.4*L, 0.6*L, 0.8*L, 0.9*L, L};
const char* phase_eedf_region_name[PHASE_EEDF_REGIONS] =
    {"x000_010", "x010_020", "x020_040", "x040_060", "x060_080", "x080_090", "x090_100"};
Ullong phase_eedf_macro[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong phase_eedf_overflow[N_XT][PHASE_EEDF_REGIONS] = {{0}};
double phase_eedf_energy_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_energy_squared_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_histogram[N_XT][PHASE_EEDF_REGIONS][PHASE_EEDF_ENERGY_BINS] = {{{0.0}}};
""",
        "phase EEDF accumulator insertion",
    )
    diagnostic = r'''
//-------------------------------------------------------------------------//
// passive regional phase-EEDF diagnostic                                 //
//-------------------------------------------------------------------------//

void measure_phase_eedf(double x, double energy_eV, int phase){
    for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
        if (x < phase_eedf_region_min[region] ||
            x > phase_eedf_region_max[region]) continue;
        phase_eedf_macro[phase][region]++;
        phase_eedf_energy_sum[phase][region] += energy_eV;
        phase_eedf_energy_squared_sum[phase][region] += energy_eV * energy_eV;
        if (energy_eV >= PHASE_EEDF_ENERGY_MAX_EV) {
            phase_eedf_overflow[phase][region]++;
        } else {
            const int bin = min(
                static_cast<int>(energy_eV / PHASE_EEDF_ENERGY_MAX_EV *
                                 PHASE_EEDF_ENERGY_BINS),
                PHASE_EEDF_ENERGY_BINS - 1);
            phase_eedf_histogram[phase][region][bin] += 1.0;
        }
    }
}

void save_phase_eedf(void){
    FILE *histogram = fopen("edupic_phase_eedf.csv", "w");
    FILE *moments = fopen("edupic_phase_eedf_moments.csv", "w");
    fprintf(histogram, "phase_bin,phase_fraction,region_id,region,x_min,x_max,energy_bin,energy_eV,macro_count,represented_count,probability_density_eV-1\n");
    fprintf(moments, "phase_bin,phase_fraction,region_id,region,x_min,x_max,macro_observations,represented_observations_m-2,overflow_fraction,mean_energy_eV,energy_standard_deviation_eV\n");
    const double bin_width = PHASE_EEDF_ENERGY_MAX_EV /
        static_cast<double>(PHASE_EEDF_ENERGY_BINS);
    const double represented_per_macro = WEIGHT / ELECTRODE_AREA;
    for (int phase=0; phase<N_XT; ++phase) {
        const double phase_fraction =
            (static_cast<double>(phase) + 0.5) / static_cast<double>(N_XT);
        for (int region=0; region<PHASE_EEDF_REGIONS; ++region) {
            const Ullong macro = phase_eedf_macro[phase][region];
            const double represented = static_cast<double>(macro) *
                represented_per_macro;
            for (int bin=0; bin<PHASE_EEDF_ENERGY_BINS; ++bin) {
                const double macro_bin = phase_eedf_histogram[phase][region][bin];
                fprintf(histogram,
                    "%d,%.17g,%d,%s,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g\n",
                    phase, phase_fraction, region, phase_eedf_region_name[region],
                    phase_eedf_region_min[region], phase_eedf_region_max[region],
                    bin, (static_cast<double>(bin) + 0.5) * bin_width,
                    macro_bin, macro_bin * represented_per_macro,
                    macro > 0 ? macro_bin / static_cast<double>(macro) / bin_width : 0.0);
            }
            const double mean = macro > 0
                ? phase_eedf_energy_sum[phase][region] /
                  static_cast<double>(macro) : 0.0;
            const double variance = macro > 0
                ? max(0.0, phase_eedf_energy_squared_sum[phase][region] /
                      static_cast<double>(macro) - mean * mean) : 0.0;
            fprintf(moments,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g\n",
                phase, phase_fraction, region, phase_eedf_region_name[region],
                phase_eedf_region_min[region], phase_eedf_region_max[region],
                macro, represented, macro > 0
                    ? static_cast<double>(phase_eedf_overflow[phase][region]) /
                      static_cast<double>(macro) : 0.0,
                mean, sqrt(variance));
        }
    }
    fclose(histogram);
    fclose(moments);
}

'''
    cycle_anchor = """//---------------------------------------------------------------------//
// simulation of one radiofrequency cycle                              //
//---------------------------------------------------------------------//

void do_one_cycle (void){
"""
    source = replace_once(source, cycle_anchor, diagnostic + cycle_anchor,
                          "phase EEDF function insertion")
    source = replace_once(
        source,
        """                energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
                meanee_xt[p][t_index]   += c1 * energy;
""",
        """                energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
                if ((t % 2) == 0) measure_phase_eedf(x_e[k], energy, t_index);
                meanee_xt[p][t_index]   += c1 * energy;
""",
        "pre-collision phase EEDF call",
    )
    source = replace_once(
        source,
        """    if (measurement_mode) {
        check_and_save_info();
    }
""",
        """    if (measurement_mode) {
        check_and_save_info();
        save_phase_eedf();
    }
""",
        "phase EEDF output call",
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
        "transform": "passive_edupic_c_regional_phase_eedf_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "regions": ["x000_010", "x010_020", "x020_040", "x040_060",
                    "x060_080", "x080_090", "x090_100"],
        "sample_every_timesteps": 2,
        "sampling_order": "pre_collision_synchronized_velocity_at_old_position",
        "energy_bins": 320,
        "energy_max_eV": 80.0,
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
