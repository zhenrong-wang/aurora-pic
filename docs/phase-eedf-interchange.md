# Phase-EEDF cross-code interchange

AuroraPIC uses a solver-neutral, versioned interchange for independent
phase-resolved regional electron-energy distribution function comparisons.
This avoids treating AuroraPIC's native CSV layout, particle weights, or
diagnostic implementation as part of the benchmark definition.

An interchange directory contains `manifest.json`, `distributions.csv`, and
`moments.csv`. Version 1 records probability mass per energy bin rather than a
probability density. Each distribution is identified by a periodic RF phase
bin and a named spatial region with SI bounds. Exact mean energy,
drift-separated temperature, observation count, and out-of-range probability
are carried separately. Energy is in eV and spatial bounds are in metres.

The comparison contract requires the same case identifier, species, phase-bin
centres, region names, and region bounds. Energy grids may differ: the
comparator treats each histogram as piecewise uniform and integrates total
variation over the union of both grids. It also reports relative differences
in exact mean energy and drift-separated temperature, overflow differences,
and optional energetic-tail fractions.

Export an AuroraPIC diagnostic window with:

```sh
python3 scripts/export_phase_eedf.py diagnostic-output cross-code/aurorapic \
  --code-version COMMIT --case-id turner-case1 --species electrons
```

Perform a descriptive comparison with:

```sh
python3 scripts/compare_phase_eedf.py cross-code/reference cross-code/aurorapic \
  --tail-eV 11.55 --tail-eV 19.82 --json comparison.json
```

Descriptive mode deliberately reports `acceptance.passes = null`. A formal
pass/fail claim is made only when all three acceptance limits were declared
before examining the candidate result:

```sh
python3 scripts/compare_phase_eedf.py cross-code/reference cross-code/aurorapic \
  --max-tv TV_LIMIT \
  --max-mean-energy-relative MEAN_LIMIT \
  --max-temperature-relative TEMPERATURE_LIMIT
```

The format does not by itself make two simulations comparable. A campaign
must additionally lock voltage phase/sign, gas state, collision cross sections
and scattering conventions, electrode and particle boundary models, initial
state, averaging horizon, region bounds, and phase sampling. Independent code
agreement is verification evidence; agreement with measurements is still
required for physical validation.

## Native velocity-space diagnostics

AuroraPIC's `phase_eedf_moments.csv` also exposes component temperatures and
energetic-tail directionality. Set `phase_eedf_tail_threshold` to an energy in
eV that is non-negative and below `phase_eedf_energy_max`. A value of zero
includes the complete sampled distribution in the tail columns.

The component temperatures remove the mean drift independently in each
velocity component. The tail columns report represented population, positive
and negative longitudinal fractions, their signed population imbalance, mean
longitudinal velocity, and the fraction of tail kinetic energy carried by the
longitudinal component. These observables distinguish an energy-distribution
difference caused by longitudinal sheath acceleration from one caused by
transverse scattering or energetic-particle residence. They are accumulated
with the EEDF and preserved by checkpoints; changing the threshold on restart
therefore requires resetting spatial-average diagnostics.

These extra native columns are diagnostic evidence, not part of interchange
version 1. A cross-code claim using them must lock the velocity sampling time,
velocity-coordinate convention, and energetic-tail threshold in its
prospective rule.

Set `phase_eedf_history = true` to add passive, per-particle history for the
selected EEDF species. History begins at `spatial_average_start_step` and is
updated every species timestep through `spatial_average_end_step`, independent
of the histogram sampling interval. Tail samples then report particle age,
total and consecutive steps above the configured tail threshold, energetic
duty fraction, threshold-entry count, collision exposure by process class,
and the fraction born during the measurement window. Particles already present
at the start are left-censored; the born-during-window fraction makes that
limitation explicit.

History counters are aligned with reusable particle slots, reset when a slot
is reused by an ionization product, and preserved by checkpoint version 20.
Restarting an older checkpoint with history enabled requires
`spatial_average_reset_on_restart = true`. The diagnostic changes neither
particle state nor the random stream, but it does add memory proportional to
the selected species' particle-storage size.

The same option writes `phase_eedf_threshold_crossings.csv`. Unlike the tail
history means, this ledger is unconditional: its denominator is every live
selected-species particle timestep in each phase/region. It records energetic
occupancy, below-to-above and above-to-below transitions between consecutive
synchronized pre-collision states, accepted collision transitions by process,
and energetic versus subthreshold births. Rates per million electron steps are
included for the interstep transitions.

An interstep transition is deliberately not labelled as field-only. It spans
the previous collision stage and the following push/drift/synchronization, so
it is the exact net change visible to the pre-collision sampler. The separate
collision before/after columns localize accepted collision crossings without
claiming that the two event sets form an additive decomposition. Threshold-
crossing accumulators share the phase-EEDF restart/reset contract.
