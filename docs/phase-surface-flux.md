# Phase-resolved internal-surface flux

AuroraPIC can measure directional particle and kinetic-energy transport across
fixed internal surfaces in a 1D simulation. The diagnostic is generic: it is
configured for one species and any strictly increasing set of internal
positions. It records left-to-right and right-to-left crossings separately by
RF phase and particle energy.

## Configuration

The diagnostic uses the RF window and phase bins of `spatial_average`:

```ini
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 1
spatial_average_end_step = 4000
spatial_average_rf_frequency = 13560000
spatial_average_rf_cycles = 1
spatial_average_phase_bins = 200

phase_surface_flux = true
phase_surface_flux_species = electrons
phase_surface_flux_positions = 0.005,0.015
phase_surface_flux_energy_bins = 320
phase_surface_flux_energy_max = 80.0
```

In SI configurations, positions are metres and energy bounds are eV. In
normalized configurations they use the configured normalized units. Positions
must satisfy `0 < x_0 < ... < x_n < length`. The diagnostic currently rejects
periodic boundaries because wrapped trajectories require an explicit unwrapped
path convention.

Crossings are detected during particle drift and are recorded on every species
timestep inside the configured spatial-average window. They are not
subsampled by `spatial_average_interval`. A particle crossing several configured
surfaces in one drift contributes once at every crossed surface. Crossing
energy uses the leapfrog drift velocity (`v_half`) and includes transverse
velocity components in 1D3V.

## Output

`phase_surface_flux.csv` contains the represented crossing histogram for every
phase, surface, direction, and energy bin. `probability_density` is normalized
by all represented crossings in that phase/surface/direction, including any
overflow population.

`phase_surface_flux_summary.csv` reports:

- macro and represented crossing counts;
- the fraction above `phase_surface_flux_energy_max`;
- represented particle flux; and
- represented kinetic-energy flux.

Fluxes assume the 1D planar unit-area convention. Direction rows contain
positive magnitudes. Signed net rightward particle or energy flux is obtained
as `left_to_right - right_to_left`.

## Restart contract

Checkpoint v17 stores the complete phase, surface, direction, and energy-bin
accumulators. A restart requires an identical diagnostic contract unless
`phase_surface_flux_reset_on_restart = true`, in which case stored accumulators
are consumed and discarded and a fresh window begins. Enabling the diagnostic
from an older checkpoint that already contains spatial-average samples likewise
requires the reset flag.

The diagnostic measures transport directly; it does not infer a flux from a
moment balance. It is currently 1D-only. Multidimensional internal-surface
crossings require geometry-aware intersection and area weighting and remain a
separate extension.
