# Reproducibility environment

This repository separates the lightweight offline paper test environment from the external-data reconstruction environment.

## Python

- Supported Python: 3.11+
- Reference CI: Python 3.12 on Ubuntu
- Core dependency: `numpy>=1.26`
- Paper test dependency: `shapely>=2,<3`
- Spatial rebuild dependency: `rasterio>=1.4`
- Figure dependency: `matplotlib>=3.9,<4`

A clean paper test environment is:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
pytest -q
```

The active paper CI validates the current v0.2 manuscript/submission package and current claim boundaries rather than the superseded v0.1 manuscript.

## R host-distribution sidecars

The pinned WCVP/HOSTS reconstructions use base-R sidecar scripts. Current post-hoc concentration/null reconstruction uses `scripts/build_wcvp_hosts_null_sidecars.R`; the historical native/contemporary builders are retained as provenance.

The pinned upstream identities are:

- rWCVPdata / WCVP v13 commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`
- HOSTS mirror commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`

The sidecar scripts use **base R only**. The historical CI used `r-lib/actions/setup-r@v2` with `r-version: "release"`. Because that label moves, the archival guarantee is the pinned source identity plus git-tracked/frozen outputs and receipts, not an assumption that future R releases are bitwise identical.

## Current v0.2 figure rebuild

Install the figure extra and run:

```bash
python -m pip install -e ".[figure]"
python scripts/render_butterfly_specialization_v02_figures.py \
  --anthropogenic-csv data/frozen/figure_sources/anthropogenic_species_metrics.csv \
  --matched-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_matched_null_v0.2.json \
  --hostbias-null-json provenance/reviewer_defenses/results/butterfly_resource_expansion_hostbias_null_v0.1.json \
  --plant-prominence-json provenance/reviewer_defenses/results/butterfly_host_plant_prominence_expansion_v0.1.json \
  --host-concentration-json provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json \
  --occurrence-csv data/frozen/figure_sources/occurrence_resource_validation_species.csv \
  --occurrence-null-json provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json \
  --occurrence-species-robustness-json provenance/reviewer_defenses/results/butterfly_occurrence_species_robustness_v0.1.json \
  --ceiling-json provenance/reviewer_defenses/results/butterfly_expansion_ceiling_sensitivity_v0.1.json \
  --regional-json provenance/reviewer_defenses/results/butterfly_regional_robustness_v0.1.json \
  --climate-csv data/frozen/figure_sources/climate_distance_sensitivity_species.csv \
  --climate-effect-json provenance/reviewer_defenses/results/butterfly_climate_effect_size_v0.1.json \
  --output-dir results/butterfly-specialization-v02-figures
```

The expected output is four PDFs and four PNGs: Figure 1, Figure 2, Figure 3 and Supplementary Figure S1. The current `.github/workflows/paper-ci.yml`, blinded-DOCX workflow and anonymous-bundle workflow all use this same renderer and source surface.

## Historical reconstruction source hashes

The earlier reconstruction/figure-source files remain byte-exact and hash-tested because they anchor the route to the current results. Their identities are documented in `provenance/FIGURE_SOURCE_ARTIFACTS.md` and checked by `tests/test_reproducibility_snapshot_integrity.py`.

The current v0.2 submission should be reproduced from the current git-tracked renderer, figure-source tables and result receipts listed above; the historical five-figure renderer is retained only for provenance.
