#!/usr/bin/env python3
"""Add passive velocity-anisotropy moments to pinned eduPIC phase EEDFs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from instrument_edupic_phase_eedf import (
    EXPECTED_SOURCE_SHA256,
    instrument as instrument_phase_eedf,
)
from instrument_edupic_surface_flux import replace_once


TAIL_THRESHOLD_EV = 15.8


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def instrument(source: str) -> str:
    source = instrument_phase_eedf(source)
    source = replace_once(
        source,
        "double phase_eedf_histogram[N_XT][PHASE_EEDF_REGIONS][PHASE_EEDF_ENERGY_BINS] = {{{0.0}}};\n",
        """double phase_eedf_histogram[N_XT][PHASE_EEDF_REGIONS][PHASE_EEDF_ENERGY_BINS] = {{{0.0}}};
double phase_eedf_vx_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_vy_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_vz_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_vx2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_vy2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_vz2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
Ullong phase_eedf_tail_macro[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong phase_eedf_tail_positive_x[N_XT][PHASE_EEDF_REGIONS] = {{0}};
Ullong phase_eedf_tail_negative_x[N_XT][PHASE_EEDF_REGIONS] = {{0}};
double phase_eedf_tail_vx_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_tail_vx2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
double phase_eedf_tail_vt2_sum[N_XT][PHASE_EEDF_REGIONS] = {{0.0}};
""",
        "velocity-space accumulator insertion",
    )
    source = replace_once(
        source,
        "void measure_phase_eedf(double x, double energy_eV, int phase){",
        "void measure_phase_eedf(double x, double vx, double vy, double vz, double energy_eV, int phase){",
        "velocity-space diagnostic signature",
    )
    source = replace_once(
        source,
        """        phase_eedf_energy_squared_sum[phase][region] += energy_eV * energy_eV;
        if (energy_eV >= PHASE_EEDF_ENERGY_MAX_EV) {
""",
        """        phase_eedf_energy_squared_sum[phase][region] += energy_eV * energy_eV;
        phase_eedf_vx_sum[phase][region] += vx;
        phase_eedf_vy_sum[phase][region] += vy;
        phase_eedf_vz_sum[phase][region] += vz;
        phase_eedf_vx2_sum[phase][region] += vx * vx;
        phase_eedf_vy2_sum[phase][region] += vy * vy;
        phase_eedf_vz2_sum[phase][region] += vz * vz;
        if (energy_eV >= 15.8) {
            phase_eedf_tail_macro[phase][region]++;
            if (vx >= 0.0) phase_eedf_tail_positive_x[phase][region]++;
            else phase_eedf_tail_negative_x[phase][region]++;
            phase_eedf_tail_vx_sum[phase][region] += vx;
            phase_eedf_tail_vx2_sum[phase][region] += vx * vx;
            phase_eedf_tail_vt2_sum[phase][region] += vy * vy + vz * vz;
        }
        if (energy_eV >= PHASE_EEDF_ENERGY_MAX_EV) {
""",
        "velocity-space accumulation",
    )
    source = replace_once(
        source,
        '    fprintf(moments, "phase_bin,phase_fraction,region_id,region,x_min,x_max,macro_observations,represented_observations_m-2,overflow_fraction,mean_energy_eV,energy_standard_deviation_eV\\n");',
        '    fprintf(moments, "phase_bin,phase_fraction,region_id,region,x_min,x_max,macro_observations,represented_observations_m-2,overflow_fraction,mean_energy_eV,energy_standard_deviation_eV,mean_velocity_x,mean_velocity_y,mean_velocity_z,drift_separated_temperature,temperature_x,temperature_y,temperature_z,tail_threshold,tail_represented_observations_m-2,tail_positive_x_fraction,tail_negative_x_fraction,tail_directional_population_imbalance,tail_mean_velocity_x,tail_longitudinal_energy_fraction\\n");',
        "velocity-space output header",
    )
    source = replace_once(
        source,
        """            fprintf(moments,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g\\n",
                phase, phase_fraction, region, phase_eedf_region_name[region],
                phase_eedf_region_min[region], phase_eedf_region_max[region],
                macro, represented, macro > 0
                    ? static_cast<double>(phase_eedf_overflow[phase][region]) /
                      static_cast<double>(macro) : 0.0,
                mean, sqrt(variance));
""",
        """            const double count = static_cast<double>(macro);
            const double ux = macro > 0 ? phase_eedf_vx_sum[phase][region] / count : 0.0;
            const double uy = macro > 0 ? phase_eedf_vy_sum[phase][region] / count : 0.0;
            const double uz = macro > 0 ? phase_eedf_vz_sum[phase][region] / count : 0.0;
            const double vx2 = macro > 0 ? phase_eedf_vx2_sum[phase][region] / count : 0.0;
            const double vy2 = macro > 0 ? phase_eedf_vy2_sum[phase][region] / count : 0.0;
            const double vz2 = macro > 0 ? phase_eedf_vz2_sum[phase][region] / count : 0.0;
            const double drift_energy = 0.5 * E_MASS *
                (ux*ux + uy*uy + uz*uz) / EV_TO_J;
            const Ullong tail = phase_eedf_tail_macro[phase][region];
            const double tail_count = static_cast<double>(tail);
            const double tail_v2 = phase_eedf_tail_vx2_sum[phase][region] +
                phase_eedf_tail_vt2_sum[phase][region];
            fprintf(moments,
                "%d,%.17g,%d,%s,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,15.8,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\\n",
                phase, phase_fraction, region, phase_eedf_region_name[region],
                phase_eedf_region_min[region], phase_eedf_region_max[region],
                macro, represented, macro > 0
                    ? static_cast<double>(phase_eedf_overflow[phase][region]) /
                      count : 0.0, mean, sqrt(variance), ux, uy, uz,
                2.0 / 3.0 * max(0.0, mean - drift_energy),
                E_MASS / EV_TO_J * max(0.0, vx2 - ux*ux),
                E_MASS / EV_TO_J * max(0.0, vy2 - uy*uy),
                E_MASS / EV_TO_J * max(0.0, vz2 - uz*uz),
                tail_count * represented_per_macro,
                tail > 0 ? static_cast<double>(phase_eedf_tail_positive_x[phase][region]) / tail_count : 0.0,
                tail > 0 ? static_cast<double>(phase_eedf_tail_negative_x[phase][region]) / tail_count : 0.0,
                tail > 0 ? (static_cast<double>(phase_eedf_tail_positive_x[phase][region]) -
                    static_cast<double>(phase_eedf_tail_negative_x[phase][region])) / tail_count : 0.0,
                tail > 0 ? phase_eedf_tail_vx_sum[phase][region] / tail_count : 0.0,
                tail_v2 > 0.0 ? phase_eedf_tail_vx2_sum[phase][region] / tail_v2 : 0.0);
""",
        "velocity-space output rows",
    )
    source = replace_once(
        source,
        "if ((t % 2) == 0) measure_phase_eedf(x_e[k], energy, t_index);",
        "if ((t % 2) == 0) measure_phase_eedf(x_e[k], mean_v, vy_e[k], vz_e[k], energy, t_index);",
        "synchronized velocity-space diagnostic call",
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
        "transform": "passive_edupic_c_regional_phase_anisotropy_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "phase_bins": 200,
        "regions": ["x000_010", "x010_020", "x020_040", "x040_060",
                    "x060_080", "x080_090", "x090_100"],
        "sample_every_timesteps": 2,
        "sampling_order": "pre_collision_synchronized_velocity_at_old_position",
        "tail_threshold_eV": TAIL_THRESHOLD_EV,
        "velocity_moments": ["vx", "vy", "vz", "vx2", "vy2", "vz2"],
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
