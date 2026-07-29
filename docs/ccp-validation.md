# Turner helium CCP validation target

AuroraPIC targets the four one-dimensional helium capacitively coupled plasma
(CCP) cases published by Turner et al., *Physics of Plasmas* 20, 013507
(2013), DOI
[10.1063/1.4775084](https://doi.org/10.1063/1.4775084). These are
code-to-code physics benchmarks with quantified statistical agreement, not
experimental validation of a particular reactor.

## Published contract

The model is a 6.7 cm planar gap containing uniform 300 K helium. One
electrode is grounded and the other is driven by

```text
phi(t) = V*sin(2*pi*13.56e6*t),
```

so the applied voltage is zero at time zero. Electrons and singly charged
helium monomer ions are absorbed at both electrodes, with no secondary
emission. The initial electron temperature is 30,000 K and the initial ion
temperature is 300 K.

| Case | Neutral density (m^-3) | Voltage amplitude (V) | Initial plasma density (m^-3) | Cells | Particles/cell/species | Steps/RF period | RF periods |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `9.64e20` | `450` | `2.56e14` | 128 | 512 | 400 | 1280 |
| 2 | `3.21e21` | `200` | `5.12e14` | 256 | 256 | 800 | 5120 |
| 3 | `9.64e21` | `150` | `5.12e14` | 512 | 128 | 1600 | 5120 |
| 4 | `3.21e22` | `120` | `3.84e14` | 512 | 64 | 3200 | 15360 |

The table uses the paper's cell count. AuroraPIC's 1D Dirichlet `nx` is the
number of endpoint-inclusive nodes, so the corresponding configurations use
`nx = cells + 1` (for example, `nx = 129` for Case 1).

The final 32 RF periods are averaged. The primary comparison is the
time-averaged ion-density profile, evaluated using the published mean,
standard-deviation, and chi-squared data. Useful secondary quantities include
mid-plane ion density and electron temperature, electron/ion power per area,
and electrode ion current.

The prescribed collision set is Biagi 7.1 electron-He elastic momentum
transfer, two excitation channels, and ionization, all with isotropic
center-of-mass scattering. Ionization divides residual energy equally between
the primary and secondary electrons. He+-He scattering uses separate
isotropic and backward components. Intermediate cross sections are linearly
interpolated and values above the table range use the final tabulated value.

## AuroraPIC readiness

The first two bounded prerequisites are complete:

- 1D Dirichlet electrodes accept independent static offsets plus sinusoidal
  amplitude, frequency, and phase;
- the field solve applies the voltage at the new field time level;
- restart reconstructs the waveform phase from the stored simulation time;
- `scalars.csv` records the actual Dirichlet `phi_left` and `phi_right`;
- `examples/rf_electrode_1d.cfg` checks the zero, quarter-cycle, and
  half-cycle values with a bounded normalized 1D3V run;
- `velocity_dimensions = 3` retains transverse velocity through initialization,
  energy diagnostics, BGK, isotropic elastic/excitation MCC, and deterministic
  velocity-aware checkpoint/restart while preserving 1D1V as the default;
- `examples/mcc_relaxation.cfg` exercises the 1D3V MCC command-line path.

The configuration is:

```ini
boundary = dirichlet
phi_left = 0
phi_right = 0
phi_right_amplitude = 450
phi_right_frequency = 13.56e6
phi_right_phase = 0
```

In SI mode, voltage is in volts, frequency is in hertz, time is in seconds,
and phase is in radians. In normalized mode, frequency is cycles per
normalized time unit. A nonzero amplitude requires a positive frequency and
Dirichlet boundaries. Driven cases require `mode = transient`; the existing
instantaneous energy-window steady-state test is intentionally rejected until
cycle-averaged convergence is implemented.

## Blocking capabilities

AuroraPIC must not claim a Turner result until all of these are complete:

1. support simultaneous electron-neutral and ion-neutral MCC models in one
   discharge;
2. create equal-weight electron/He+ products during 1D ionization while
   preserving charge, energy accounting, bounded storage, and restart;
3. implement the benchmark's two-component He+-He scattering;
4. import and pin the exact benchmark collision tables with their permitted
   redistribution terms and provenance;
5. report species-resolved electrode current, deposited power, spatial
   density, ionization source, and phase/time averages;
6. implement statistically bounded campaign and chi-squared comparison
   tooling;
7. implement whole-RF-cycle convergence and phase/time averaging;
8. run the full case only through an explicit production profile. Case 1
   alone requires 512,000 steps at the published resolution, so it must never
   become an ordinary laptop/CI regression.

Open-access implementations such as the
[WarpX capacitive-discharge example](https://warpx.readthedocs.io/en/latest/usage/examples/capacitive_discharge/README.html)
are useful independent integration references. LXCat is an open-access
platform, but access does not erase dataset-specific attribution, version,
or redistribution conditions. AuroraPIC's existing local import workflow
therefore remains the route for user-supplied data until an exact
redistributable benchmark package is identified and reviewed.
