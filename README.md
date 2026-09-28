# chocho — Butterfly specialization ecology

Reproducible code, frozen analysis receipts, and manuscript materials for:

> **Anthropogenic host redistribution expands butterfly resource geography through contrasting host portfolios in specialists and generalists**

Target journal: *Global Ecology and Biogeography*.

## Scientific result

The paper treats butterfly specialization as a hierarchy of partially coupled, non-interchangeable ecological dimensions rather than a single specialist–generalist axis.

The frozen manuscript claim map supports six claims, presented here in manuscript order:

1. Introduced host ranges expand reconstructed resource opportunity for most butterflies, but proportional expansion is not concentrated in broad family-level generalists.
2. Similar aggregate expansion can be assembled through concentrated host contributions in specialists or distributed host portfolios in generalists.
3. Family-level host breadth only partly tracks geographic larval-resource breadth.
4. Species-level host richness reveals strong specialization structure even within a fixed host-family breadth category.
5. Climate commonly filters realized butterfly geography within reconstructed contemporary larval-resource opportunity.
6. The independent prediction that broader host-family diets weaken climate filtering after controlling contemporary resource breadth was **not supported** and is retained without retuning.

See `manuscript/butterfly_specialization_claim_map_v0.1.json` for the exact claim boundaries and source receipts.

## Repository boundary

This repository is self-contained for the butterfly ecology paper. Method-development history and unrelated response domains are not part of the scientific argument presented here.

Migration provenance is retained separately in `provenance/migration/MIGRATION_PROVENANCE.md` and `provenance/migration/SOURCE_SNAPSHOT.txt` so the origin of code and frozen artifacts remains auditable without entering the manuscript narrative.

The scientific question/hypothesis lineage is documented separately in `provenance/SCIENTIFIC_ORIGIN_AND_HYPOTHESIS_LINEAGE.md`. It distinguishes discovery chronology from manuscript presentation order and is not part of the manuscript argument.

## Main evidence path

```text
fixed external sources
  -> response-blind S1 species manifest
  -> host-resource reconstruction
  -> ecological analyses
  -> frozen result receipts
  -> claim map
  -> manuscript
```

Important fixed external identities include:

- HOSTS mirror commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`
- rWCVPdata / WCVP v13 snapshot `65bed76bae9d644ccb6ad200c05f9f5071d89e05`
- exact LepTraits snapshot SHA-256 `6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b`

The repository stores the frozen identities and reconstruction logic rather than treating changing upstream resources as interchangeable.

## Repository map

Each major directory has its own short navigation index.

- `manuscript/README.md` — submission-facing manuscript, blinded source, title page, cover letter, claim map and readiness metadata.
- `benchmarks/README.md` — frozen result receipts and which ones support the manuscript.
- `docs/README.md` — reproducibility/release instructions and frozen ecological protocols.
- `scripts/README.md` — normal paper entry points versus diagnostics/history.
- `src/butterfly_specialization_ecology/` — reusable analysis modules used by the paper scripts and tests.
- `tests/README.md` — scientific, manuscript, release and layout test map.
- `data/README.md` — vendored external input, frozen S1 manifest and byte-exact figure sources.
- `provenance/README.md` — migration, workflow, figure-source and hypothesis-lineage audit history.

## Reproducing the offline test suite

Python 3.11+ is required.

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
pytest -q
```

The GitHub Actions workflow `.github/workflows/paper-ci.yml` runs this suite without network-dependent tests. The ecology-only pruning state has passed this CI.

The `analysis` optional dependency group additionally provides Shapely and rasterio for spatial analyses; `figure` provides Matplotlib.

For full reconstruction requirements, including the two base-R WCVP/HOSTS sidecar builders and pinned upstream commits, see `docs/REPRODUCIBILITY.md`. The R scripts have no external R-package dependencies; the historical frozen workflow used `r-lib/actions/setup-r@v2` with `r-version: "release"`.

## Interpretation boundaries

The reconstructed host envelopes represent **potential regional larval-resource opportunity**, not confirmed butterfly occupancy.

Introduced host distributions are not interpreted as proof that host introduction caused butterfly range expansion. GBIF non-observation is not treated as true absence. Climate filtering is an association within the declared resource-opportunity design, not proof of physiological causation.

The independent host-breadth climate-release prediction was not supported (`partial Spearman rho = -0.166`, one-sided `p = 0.2237`, `n = 24`). No response-driven replacement predictor search is used to rescue that result.

## Submission state

The scientific story is frozen, and both review-package builders are complete: the privacy-scrubbed line-numbered DOCX and the de-identified reviewer code bundle are generated and CI-audited. Remaining release tasks are administrative/external:

- choose and add the software license;
- complete final author/affiliation/CRediT/funding/conflict metadata;
- create/tag the final `v1.0.0-butterfly` release and mint the permanent archive DOI;
- upload the anonymous reviewer bundle behind a non-identifying reviewer-access URL and insert that URL before the final DOCX rebuild.

Until the archival release is minted, cite the eventual versioned release rather than a moving branch.

The v1 release identifiers are fixed as package/CITATION version `1.0.0` and intended GitHub tag `v1.0.0-butterfly`. To show only the remaining administrative blockers, run:

```bash
python scripts/paper/release_preflight.py
```

This preflight deliberately fails until the software license, final author/CRediT metadata, public archive DOI, and anonymized reviewer-access URL are supplied. It does not request or authorize additional ecological analyses.

For the exact DOI-first release order, use `docs/RELEASE_PROCEDURE.md`. The v1 archive should use a manual Zenodo draft so the DOI can be reserved and written into the exact green release commit before `v1.0.0-butterfly` is tagged.
