# Data layout

- `external/` — vendored external input required for reproducibility. The LepTraits snapshot is byte-identified and hash-guarded.
- `frozen/s1_species_manifest.json` — exact 339-species S1 manifest used by the resource reconstruction.
- `frozen/figure_sources/` — byte-exact inputs consumed by the audited five-figure renderer.

Larger HOSTS, WCVP, GBIF and CHELSA inputs are not vendored wholesale. Their pinned identities, reconstruction rules and provenance are documented in `docs/REPRODUCIBILITY.md`.
