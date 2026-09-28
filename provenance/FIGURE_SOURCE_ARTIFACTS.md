# Manuscript figure source artifacts

The five submitted manuscript figures were rendered from four hash-verified source files produced by frozen TTF workflow runs. The historical rendering workflow is preserved at `provenance/workflows/butterfly-specialization-manuscript-figures-v01.yml`.

The original GitHub Actions copies are temporary, but the exact source bytes are now vendored under `data/frozen/figure_sources/`. The vendored copies are the durable release inputs and are hash-checked by `tests/test_reproducibility_snapshot_integrity.py`.

| Input | TTF workflow run | Artifact ID | Original file inside artifact | Durable chocho path | SHA-256 | Actions expiry (UTC) |
|---|---:|---:|---|---|---|---|
| S1 resource descriptors | 36222306369 | 10898504566 | `butterfly-resource-envelope-pilot-v0.1/s1_resource_descriptors.csv` | `data/frozen/figure_sources/s1_resource_descriptors.csv` | `894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc` | 2026-10-26 06:04 |
| Anthropogenic species metrics | 36243549495 | 10906299578 | `butterfly-anthropogenic-resource-expansion-v0.1/species_metrics.csv` | `data/frozen/figure_sources/anthropogenic_species_metrics.csv` | `7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30` | 2026-10-26 13:05 |
| Host-contribution metrics | 36270981110 | 10915562949 | `species_mechanism_metrics.csv` | `data/frozen/figure_sources/host_contribution_metrics.csv` | `b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa` | 2026-10-26 21:07 |
| Independent climate result | 36270581743 | 10919985332 | `primary_result.json` | `data/frozen/figure_sources/independent_climate_primary_result.json` | `a73dca6e8b669f721d8f2745f27198d9b47a5cbd05dde6146cc8d4f1ddfbf79b` | 2026-10-27 00:42 |

The hashes above are the identities consumed by the audited figure run recorded in `manuscript/butterfly_specialization_figures_v0.1.json`. No analysis rerun or result change was used when the files were moved into the ecology repository.

The CSV sources retain their original CRLF byte representation; this matters because the release gate verifies SHA-256 on raw bytes, not merely parsed table equality.

The figure renderer itself is active source code in `scripts/paper/render_butterfly_specialization_manuscript_figures.py`. The historical workflow file is stored only as provenance because the ecology repository was intentionally separated from the broader TTF execution surface.
