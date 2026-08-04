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
