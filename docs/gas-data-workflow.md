# Local real-gas data workflow

AuroraPIC does not download or redistribute measured cross sections. Users
obtain a complete, appropriately licensed LXCat/BOLSIG+ text export and run
the converter locally:

```bash
python3 scripts/import_lxcat.py ~/Downloads/argon.txt \
  --output-dir local-gases/argon \
  --gas Ar \
  --neutral-mass 6.6335209e-26 \
  --dataset-id provider.argon \
  --dataset-version 2026-01 \
  --provenance "provider, dataset, source URL" \
  --citation "citation requested by the contributor" \
  --retrieved 2026-01-15 \
  --license "applicable dataset terms" \
  --neutral-density 2.4e20
```

The input remains in the user's local environment. The output directory
contains:

- one strict two-column table per supported collision channel;
- a version-2 SI `.gas` manifest with eV-to-joule scaling;
- `audit.json`, containing the source SHA-256, channel coverage, thresholds,
  energy ranges, point counts, peak cross sections, elastic mass-ratio checks,
  and a collision-frequency recommendation over the audited energy range.

The importer supports BOLSIG+ `ELASTIC`, `EXCITATION`, and `IONIZATION`
blocks. It rejects malformed tables, nonzero inelastic data below threshold,
ambiguous units, mismatched elastic mass ratios, incomplete electron sets, and
unsupported processes by default. `--allow-partial`,
`--skip-unsupported`, and `--ignore-mass-ratio` are explicit acknowledgements,
not silent fallbacks. AuroraPIC does not yet convert `EFFECTIVE`, `MOMENTUM`,
`ATTACHMENT`, or `ROTATION` blocks.

LXCat/BOLSIG+ exports use eV for energy and square metres for cross sections.
The generated manifest therefore declares `units = si` and cannot be loaded
by a normalized simulation. Gas manifest version 2 makes this unit contract
mandatory. Version 1 remains readable as legacy normalized data.

The recommended `max_frequency` is

```text
safety_factor * neutral_density * max_E(sum_i(sigma_i(E)) * speed(E))
```

over the requested audit interval. The default upper energy is the largest
tabulated energy. Set `--max-energy-ev` to the full energy envelope expected
in the simulation. This calculation is a preflight bound for stationary
neutrals and the imported tables; it does not replace runtime enforcement or
an energy-coverage study.

Before publication or device simulation:

1. Review the original contributor's references, comments, and license.
2. Keep a complete internally consistent set from one source when possible.
3. Retain the source file and compare its SHA-256 with `audit.json`.
4. Run swarm validation over the intended reduced-field range.
5. Use the same validated gas package unchanged in the geometry simulation.

The checked-in LXCat fixture is synthetic parser-validation data only.
