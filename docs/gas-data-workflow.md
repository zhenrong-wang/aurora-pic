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

The importer supports BOLSIG+ `ELASTIC`, `EXCITATION`, `IONIZATION`, and
`ATTACHMENT` blocks. Attachment onset is retained in its tabulated
cross-section values because that BOLSIG+ block has no separate threshold
parameter. It rejects malformed tables, nonzero inelastic data below threshold,
ambiguous units, mismatched elastic mass ratios, incomplete electron sets, and
unsupported processes by default. `--allow-partial`,
`--skip-unsupported`, and `--ignore-mass-ratio` are explicit acknowledgements,
not silent fallbacks. AuroraPIC does not yet convert `EFFECTIVE`, `MOMENTUM`,
or `ROTATION` blocks.

Imported elastic channels remain explicitly isotropic. The converter does not
guess a differential phase function or reinterpret momentum-transfer data as
a total elastic cross section. Energy-dependent anisotropic scattering can be
added only when the source provides enough documented information to build
the separate mean-cosine table described in
[`collisions.md`](collisions.md).

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
neutrals and the imported tables. Positive-temperature SI runs additionally
enforce a conservative bound over the reachable bounded-Maxwellian relative
speeds, so the stationary recommendation may need to be raised. The importer
audit does not replace runtime enforcement or an energy-coverage study.

Before publication or device simulation:

1. Review the original contributor's references, comments, and license.
2. Keep a complete internally consistent set from one source when possible.
3. Retain the source file and compare its SHA-256 with `audit.json`.
4. Run the documented
   [electron-swarm validation](swarm-validation.md) over the intended
   reduced-field range.
5. Use the same validated gas package unchanged in the geometry simulation.

The checked-in LXCat fixture is synthetic parser-validation data only.

## Pinned eduPIC reference package

The independent eduPIC 1.0 argon target has a separate strict importer. The
upstream code contains a dormant `test_cross_sections()` routine that writes
its five analytic channels to a six-column, one-million-row table. Enable that
routine only in a local checkout of pinned commit
`32050728c961a317d6d6acd6bc86d026da403326`, build it locally, and then run:

```bash
python3 scripts/import_edupic_cross_sections.py \
  /local/edupic-run/cross_sections.dat \
  --output-dir /local/gases/edupic-argon \
  --source-sha256 <sha256-of-cross_sections.dat> \
  --retrieved 2026-08-10
```

The importer does not contain or reproduce the GPL cross-section formulas. It
requires exactly 1,000,000 rows on the `0.001 eV` source grid, validates all
six finite non-negative columns and both inelastic thresholds, and writes
separate electron and ion manifests. The electron package declares elastic,
excitation, and Opal-style ionization with a `10 eV` ejected-energy scale. The
two inelastic electron channels declare the finite-mass center-of-mass
transform used by the reference. The ion package declares center-of-mass
isotropic and backward elastic channels.
All five channels declare lower-bin cross-section evaluation, matching
eduPIC's integer selection on the `0.001 eV` source grid. Every generated
artifact is hashed in `audit.json`.

These generated tables remain local and untracked. The workflow records their
GPL-derived provenance but makes no redistribution or relicensing claim. Use
`scripts/prepare_edupic_argon_case.py` with
`examples/edupic_argon_ccp_reference.case` to validate the exact local hashes
and generate only a bounded, non-production AuroraPIC preflight deck.
