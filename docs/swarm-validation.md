# Electron-swarm validation

`aurorapic_swarm` exercises the same three-velocity null-collision kernel used
by imported-geometry simulations in a homogeneous prescribed electric field.
It isolates gas transport from mesh, wall, and space-charge effects and scans
one or more reduced electric fields:

```bash
aurorapic_swarm examples/synthetic_swarm.swarm
```

The checked-in example and its cross sections are synthetic interface data,
not argon or another material. For physical work, set `gas_data_file` to a
locally imported version-2 SI gas package and use the collision-frequency
recommendation from its `audit.json`.

## Configuration

```ini
swarm_config_version = 1
gas_data_file = local-gases/argon/Ar.gas
neutral_density = 2.4e20
reduced_fields_td = 1, 2, 5, 10, 20, 50, 100
max_frequency = 2.0e8
timestep = 2.5e-10
steps = 20000
burn_in_steps = 10000
particles = 10000
uncertainty_blocks = 10
work_item_limit = 1500000000
initial_mean_energy_ev = 0.1
max_energy_ev = 100
seed = 271828
output_file = argon_swarm.csv
```

The runner requires `max_frequency * timestep <= 0.1`. `max_energy_ev` must
not exceed the upper energy of any collision table, and the run fails as soon
as an electron exceeds it. This turns energy coverage into an enforced
precondition rather than constant extrapolation beyond measured data. The
sampling-step count (`steps - burn_in_steps`) must divide evenly into
`uncertainty_blocks`. A scan is also bounded to 10 million allocated
particles and, by default, 100 million particle-step-field work items. Raising
`work_item_limit` is an explicit acknowledgement that a larger run belongs on
an appropriately monitored compute host.

The CSV embeds the dataset identity, version, citation, provenance, retrieval
date, license text, gas-manifest path, numerical controls, and per-field seed.
For each E/N it reports:

- signed mean electron velocity along the electric field and the conventional
  electron drift velocity opposite the field;
- reduced mobility `N * drift_velocity / electric_field`;
- mean electron energy;
- longitudinal and transverse endpoint diffusion estimates;
- maximum observed energy;
- null-collision candidate counts and per-channel rates;
- block standard errors for drift and mean energy, plus Poisson counting
  errors for channel rates.

The reduced-mobility unit is `1 / (V m s)`. One Townsend is
`1e-21 V m^2`.

## Model boundary

This first swarm implementation intentionally uses a fixed electron
population. Ionization divides the available excess energy using the engine's
current equal-sharing model and increments the ionization rate, but the
secondary electron is not added to the ensemble. The reported rate is
therefore not a Townsend avalanche coefficient. The model also assumes
stationary zero-temperature neutrals and isotropic scattering; the configured
neutral mass is active in elastic recoil.

Consequently, this benchmark can validate the current MCC implementation's
drift, mean energy, diffusion trend, and collision rates. It cannot yet claim
high-accuracy transport for datasets requiring differential angular
scattering, thermal neutral motion, non-equal ionization energy sharing, or
electron multiplication.

## Production study checklist

1. Use one complete and internally consistent collision set.
2. Compare against independent measured or evaluated swarm coefficients.
3. Repeat with smaller timesteps and higher particle counts.
4. Extend burn-in until block means no longer show a transient trend.
5. Confirm the maximum observed energy remains comfortably below
   `max_energy_ev`.
6. Record the imported package and `audit.json` alongside the CSV.
7. Do not tune cross sections against device results before the homogeneous
   swarm comparison is understood.
