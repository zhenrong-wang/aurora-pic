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

For Case 1, the published chi-squared acceptance ranges are `55–303` at 95%
and `48–405` at 99%. Passing one stochastic realization is evidence of
statistical consistency, not a guarantee for arbitrary discharges; repeated
runs should leave those ranges at approximately their declared frequencies.

The prescribed collision set is Biagi 7.1 electron-He elastic momentum
transfer, two excitation channels, and ionization, all with isotropic
center-of-mass scattering. Ionization divides residual energy equally between
the primary and secondary electrons. He+-He scattering uses separate
isotropic and backward components. Intermediate cross sections are linearly
interpolated and values above the table range use the final tabulated value.

## AuroraPIC readiness

The first five bounded prerequisites are complete:

- 1D Dirichlet electrodes accept independent static offsets plus sinusoidal
  amplitude, frequency, and phase;
- the field solve applies the voltage at the new field time level;
- restart reconstructs the waveform phase from the stored simulation time;
- `scalars.csv` records the actual Dirichlet `phi_left` and `phi_right` plus
  live macro-particle counts for each species;
- `examples/rf_electrode_1d.cfg` checks the zero, quarter-cycle, and
  half-cycle values with a bounded normalized 1D3V run;
- `velocity_dimensions = 3` retains transverse velocity through initialization,
  energy diagnostics, BGK, isotropic elastic/excitation MCC, and deterministic
  velocity-aware checkpoint/restart while preserving 1D1V as the default;
- `examples/mcc_relaxation.cfg` exercises the 1D3V MCC command-line path.
- named 1D collision models can target electrons and ions simultaneously with
  qualified diagnostics and a combined restart fingerprint;
- 1D3V ionization creates equal-weight electron/ion pairs, removes the
  configured threshold energy, preflights bounded storage, and defers newborn
  collisions until the next timestep;
- `examples/mcc_ionization_1d.cfg` exercises simultaneous electron/ion MCC,
  reactive products, diagnostics, capacity, and checkpoint integration.
- elastic tables can explicitly select projectile or center-of-mass lookup
  energy, and 1D3V supports separate finite-mass isotropic and exactly backward
  channels. For equal-mass He+-He, the backward channel exchanges projectile
  and neutral velocities while conserving momentum and energy. Both choices
  are restart-fingerprinted.

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

1. retain and audit the locally normalized, checksum-pinned collision tables
   without changing their interpolation contract, while keeping the
   publisher-derived files local unless redistribution permission is
   established;
2. report species-resolved electrode current, deposited power, and ionization
   source; restart-safe spatial-density time averaging is complete;
3. connect the implemented chi-squared comparator to a statistically bounded,
   run-contract-validated campaign;
4. implement whole-RF-cycle convergence and phase/time averaging;
5. run the full case only through an explicit production profile. Case 1
   alone requires 512,000 steps at the published resolution, so it must never
   become an ordinary laptop/CI regression.

Open-access implementations such as the
[WarpX capacitive-discharge example](https://warpx.readthedocs.io/en/latest/usage/examples/capacitive_discharge/README.html)
are useful independent integration references. LXCat is an open-access
platform, but access does not erase dataset-specific attribution, version,
or redistribution conditions. AuroraPIC's existing local import workflow
therefore remains the route for user-supplied data until an exact
redistributable benchmark package is identified and reviewed.

## Publisher supplement lock

The exact publisher supplement was manually acquired on 30 July 2026. The
committed `examples/turner_ccp.sources` registry pins the 133,969-byte archive
as SHA-256
`a0a5fe93900d7d7b213157f1eab664e06aab6e718f2189910e65f23bd699d661`
and separately pins all four members. It contains the Biagi 7.1 electron
table, the 101-row centre-of-mass He+-He table, original results for all four
published grids, and refined results. Case 1 has the required 129 nodes over
the 6.7 cm gap.

The electron table carries an embedded `All rights reserved` notice.
Consequently the archive and extracted tables remain under ignored `tmp/`;
AuroraPIC commits only identity and semantic metadata. Verify a local copy
without modifying or redistributing it with:

```sh
python3 scripts/verify_turner_source.py \
  examples/turner_ccp.sources \
  tmp/013507_1_supplements.zip \
  --output tmp/turner-source-verification.json
```

Normalize the verified bytes into local AuroraPIC manifests, channel tables,
reference CSV files, and a transformation audit with:

```sh
python3 scripts/normalize_turner_source.py \
  examples/turner_ccp.sources \
  tmp/013507_1_supplements.zip \
  --output-dir tmp/turner-normalized-v3
```

Electron thresholds are converted from eV to joules using the 2006 CODATA
elementary charge required by the paper. Normalization audit version 2 records
that constant edition explicitly. Ion energy remains tabulated in centre-of-mass eV and ion
cross sections remain in their source `1e-20 m2` scale; the gas manifest
performs both SI conversions at load time. No table is resampled. Named 1D
collision models load either generated manifest through `gas_data_file`, with
only ionization product-species mappings supplied by the simulation deck.

## Exact Case 1 campaign preflight

[`examples/turner_helium_ccp_case1.case`](../examples/turner_helium_ccp_case1.case)
locks the published Case 1 physics and numerical contract independently of
the restricted local tables. After normalization, generate the production
deck and its machine-readable preflight report with:

```sh
python3 scripts/prepare_turner_case.py \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --output tmp/turner-case1-campaign/turner_case1.cfg \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN
```

The preparer verifies the normalization audit and every normalized file,
derives the exact 65,536 macro-particles per species and weight, locks the
512,000-step/1,280-cycle duration and final 12,800-sample window, and audits
thermal-neutral collision majorants through the declared 10 keV projectile
envelopes. The scan follows the production kernel's eight-standard-deviation
thermal bound, applies a safety factor, and records both sampled peaks and
configured majorants. During execution, the collision kernel independently
checks every encountered particle and fails rather than biasing the result if
a majorant is exceeded.

Generation never launches the calculation and the report explicitly carries
`physics_claim = none`. The exact campaign has a floor of 67,108,864,000
initial-particle updates before collision work and particle creation, so it
must first receive a bounded, single-core runtime qualification. Use
`aurorapic_cli --validate-only` to parse the generated deck without allocating
or advancing its particles.

Run that qualification explicitly with:

```sh
python3 scripts/qualify_turner_runtime.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --steps 4 \
  --max-initial-updates 1000000 \
  --timeout-seconds 60 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_PROBE \
  --report tmp/turner-case1-runtime-qualification.json
```

The qualifier regenerates and audits the exact production deck, changes only
duration, output, checkpoint, averaging, and serial-runtime controls for the
probe, then requests lower process priority and one CPU. It will not exceed
its built-in two-million-initial-update ceiling, refuses to overwrite a
report, and never forwards the CLI's large-run authorization. Its projection
is an initial-population-only planning estimate, not a promise: ionization,
electrode losses, checkpoints, storage, and long-run hardware behavior can
all change the production cost.

The next bounded rung executes exactly one RF cycle in two half-cycle stages:

```sh
python3 scripts/run_turner_startup.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --work-dir tmp/turner-case1-startup \
  --report tmp/turner-case1-startup/report.json \
  --max-initial-updates 60000000 \
  --timeout-seconds 120 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_ONE_CYCLE_TURNER_STARTUP
```

This stage advances 52,428,800 initial-particle updates, checkpoints at the
half-cycle, restarts to the full cycle, and retains both checkpoints and both
diagnostic directories under ignored `tmp/`. Its report requires exact
integer population and collision-counter continuity, floating diagnostic
continuity to roundoff, the prescribed electrode waveform, finite energy and
charge histories, charge-paired ionization/population balance, and a
129-node final field. It reports an early boundary-to-bulk field indicator
but explicitly does not interpret it as a stationary sheath.

Extend that retained one-cycle checkpoint through cycles two to four with:

```sh
python3 scripts/extend_turner_horizon.py \
  build/aurorapic_cli \
  examples/turner_helium_ccp_case1.case \
  tmp/turner-normalized-v3 \
  --prior-work-dir tmp/turner-case1-startup \
  --prior-report tmp/turner-case1-startup/report.json \
  --work-dir tmp/turner-case1-horizon \
  --report tmp/turner-case1-horizon/report.json \
  --additional-cycles 3 \
  --max-initial-updates 160000000 \
  --timeout-seconds 120 \
  --acknowledge-cost I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_HORIZON
```

The horizon runner verifies the prior report and checkpoint hashes, advances
one whole cycle per low-priority serial process, and writes a distinct
checkpoint and phase-matched diagnostic set after every cycle. A completed
horizon report/work directory can be supplied as the next `--prior-*` pair,
so bounded blocks chain without repeating earlier cycles. Per-cycle
reports contain collision deltas, charge-paired population balances, inferred
electrode losses, energy and charge ranges, waveform error, restart
continuity, and boundary-to-bulk field structure. Four added cycles are the
built-in maximum; this is a startup trend screen, not a stationarity test.

Screen a contiguous hash-chained report history before authorizing any
published-profile comparison:

```sh
python3 scripts/analyze_turner_stationarity.py \
  tmp/turner-case1-startup/report.json \
  tmp/turner-case1-horizon/report.json \
  tmp/turner-case1-horizon-next/report.json \
  --window-cycles 4 \
  --max-population-change 0.005 \
  --max-observable-span 0.05 \
  --output tmp/turner-case1-stationarity.json
```

The default engineering screen requires both species to change by no more
than 0.5% per cycle throughout the final four cycles and requires the
ionization count, phase-zero boundary field, and phase-zero total energy each
to span no more than 5% across that window. These are conservative AuroraPIC
startup gates, not thresholds published by Turner et al. Passing only permits
a longer confirmation window; failing explicitly keeps the published
chi-squared comparison inapplicable.

The case manifest separately pins the published Case 1 summary from Table III:
mid-plane ion density `1.40e14 m^-3`, electron temperature `9.36 eV`,
electron and ion power `34.3` and `90.6 W m^-2`, ion current
`0.219 A m^-2`, and approximately 31,900 total macro-particles. The
stationarity report includes the current-to-reported macro-particle ratio as
context, but does not turn that reported stochastic population into an
acceptance gate.

## Restart-safe density averaging

The primary Turner observable is now produced by a generic 1D post-step
spatial averager. It deposits each species' represented particle number with
the same linear particle-grid shape and node-volume convention used for
charge, then accumulates nodal density independently of the ordinary output
interval. For Case 1, the exact final-32-cycle contract is:

```ini
spatial_average = true
spatial_average_interval = 1
spatial_average_start_step = 499201
spatial_average_end_step = 512000
spatial_average_rf_frequency = 13.56e6
spatial_average_rf_cycles = 32
```

This selects 12,800 post-step samples. Configuration validation requires an
integer number of timesteps per RF cycle, a whole-cycle window ending at the
time-zero drive phase, an interval that divides the cycle, and agreement with
every active electrode-drive frequency. `spatial_average.csv` uses long-form
species/node rows; `spatial_average_metadata.json` records the window,
timestep, sample count, species, and a `complete` gate.

1D checkpoint v5 stores the averaging contract, sample count, and every nodal
sum. A changed averaging window is rejected on restart, and a legacy v1-v4
checkpoint can restart only when spatial averaging is disabled. A bounded
regression proves byte-identical continuous and checkpoint-split profiles and
represented-number conservation.

The baseline statistical comparison uses the original-grid reference, not the
numerically refined profile. For every mesh node it computes

```text
X^2 = sum((candidate_ion_mean - reference_ion_mean)^2
          / reference_ion_population_stddev^2).
```

The population standard deviation is the denominator; the standard deviation
of the reference mean is not. Once a candidate final-32-period profile exists,
run:

```sh
python3 scripts/compare_turner.py \
  --case 1 \
  --reference tmp/turner-normalized-v3/turner_case1_benchmark.csv \
  --candidate tmp/turner-case1-run/spatial_average.csv \
  --candidate-metadata \
    tmp/turner-case1-run/spatial_average_metadata.json \
  --species ions \
  --normalization-audit tmp/turner-normalized-v3/audit.json \
  --output tmp/turner-case1-comparison.json
```

The tool requires exact node coordinates and does not interpolate. It rejects
partial data unless metadata proves the exact case timestep and final 32 RF
cycles. A report still makes no overall physics claim because density
metadata alone cannot prove the prescribed particle count, collision model,
or prior stationary state.

## Bounded execution ladder

Case 1 remains the smallest whole-discharge target, but shortening it changes
the statistical benchmark. Work therefore proceeds without relabeling a
reduced run as a pass:

1. **C0, collision-law verification:** exact center-of-mass energy lookup,
   finite-mass isotropic recoil, backward velocity exchange, electron
   equal-sharing ionization, and restart fingerprints. This stage is complete.
2. **C1, source and diagnostic lock:** the electronic supplement is acquired
   outside the repository and its archive/member hashes, structure, and usage
   constraint are pinned. Exact tables are locally normalized without
   resampling and load through provenance-bearing gas manifests. Restart-safe
   whole-cycle spatial-density averaging and the statistical comparator are
   complete. Implement the remaining power/current/source observables to
   complete the broader C1 diagnostic set.
3. **C2, bounded startup screening:** run the exact Case 1 grid, timestep,
   particle weight, waveform, and collision model for a small declared number
   of RF periods. Check invariants, collision balance, sheath formation,
   restart, and resource behavior. This cannot use the steady-state
   chi-squared acceptance range. The exact-population, few-step runtime
   qualification and its hard resource gates are complete. A checkpoint-split
   whole-RF-period startup now proves restart continuity, species/collision
   balance, waveform timing, finite diagnostics, and early boundary-field
   formation. The bounded four-cycle horizon adds phase-matched population,
   collision, and field trends; longer blockwise stationarity remains.
4. **C3, checkpointed steady-state comparison:** continue only through
   explicit low-priority blocks, establish whole-cycle stationarity, average
   the final 32 periods, and apply the published chi-squared test.

This ladder gives useful bounded evidence early while preserving one honest
definition of a Turner benchmark pass.
