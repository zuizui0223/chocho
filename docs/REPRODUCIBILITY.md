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

The active GitHub Actions paper CI installs the test extra, so Shapely-dependent tests are executed rather than silently skipped.

## R host-distribution sidecars

Two reconstruction scripts use R:

- `scripts/resource/build_wcvp_hosts_sidecar.R`
- `scripts/resource/build_wcvp_hosts_contemporary_sidecar.R`

They use **base R only**. There are no CRAN/Bioconductor package dependencies and no `library()` or `require()` calls.

The historical frozen WCVP/HOSTS sidecar workflow used:

- Ubuntu GitHub Actions runner
- `r-lib/actions/setup-r@v2`
- `r-version: "release"`

The scripts read the pinned upstream files directly:

- rWCVPdata / WCVP v13 commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`
- HOSTS mirror commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`

Example rebuild:

```bash
Rscript scripts/resource/build_wcvp_hosts_sidecar.R \
  /path/to/rWCVPdata \
  /path/to/HOSTS \
  results/wcvp_native_sidecar

Rscript scripts/resource/build_wcvp_hosts_contemporary_sidecar.R \
  /path/to/rWCVPdata \
  /path/to/HOSTS \
  results/wcvp_contemporary_sidecar
```

Because the historical workflow selected the moving R label `release` rather than a semantic R version, the archival guarantee is the frozen input identity plus output SHA-256, not an assumption that a future R release is bitwise identical. Release provenance should therefore retain the exact frozen sidecar/source outputs and their hashes.

## Figure-source inputs

The exact four inputs consumed by the audited manuscript-figure workflow are vendored in `data/frozen/figure_sources/` and recorded in `provenance/FIGURE_SOURCE_ARTIFACTS.md`. Their raw-byte SHA-256 values are part of the release boundary.

A figure rebuild no longer requires historical Actions artifacts:

```bash
python -m pip install -e ".[figure]"
python scripts/paper/render_butterfly_specialization_manuscript_figures.py \
  --descriptors-csv data/frozen/figure_sources/s1_resource_descriptors.csv \
  --anthropogenic-csv data/frozen/figure_sources/anthropogenic_species_metrics.csv \
  --mechanism-csv data/frozen/figure_sources/host_contribution_metrics.csv \
  --climate-primary-json data/frozen/figure_sources/independent_climate_primary_result.json \
  --output-dir results/butterfly-specialization-manuscript-figures-v0.1
```

The active paper CI executes this rebuild and requires five PDF plus five PNG outputs.
