# chocho — Butterfly specialization ecology

Reproducible code, frozen analysis receipts, and manuscript materials for:

> **Butterfly specialization is hierarchical: host portfolios structure resource opportunity while climate filters realized geography**

Target journal: *Global Ecology and Biogeography*.

## Scientific result

The paper treats butterfly specialization as a hierarchy of partly independent ecological dimensions rather than a single specialist–generalist axis.

The frozen manuscript claim map currently supports six claims:

1. Family-level host breadth only partly tracks geographic larval-resource breadth.
2. Introduced host ranges expand reconstructed resource opportunity for most butterflies, but proportional expansion is not concentrated in broad family-level generalists.
3. Similar aggregate expansion can be assembled through concentrated host contributions in specialists or distributed host portfolios in generalists.
4. Species-level host richness reveals strong specialization structure even within a fixed host-family breadth category.
5. Climate commonly filters realized butterfly geography within reconstructed contemporary larval-resource opportunity.
6. The independent prediction that broader host-family diets weaken climate filtering after controlling contemporary resource breadth was **not supported** and is retained without retuning.

See `manuscript/butterfly_specialization_claim_map_v0.1.json` for the exact claim boundaries and source receipts.

## Repository boundary

This is the ecology-paper repository. It is intentionally separate from `zuizui0223/TTF`, which contains the transferability/qualification methodology.

The split is pinned to TTF commit:

`1a112334cca2f2ef5e234c3ae1fc1a80b8266956`

Migration provenance is recorded in `MIGRATION_PROVENANCE.md` and `SOURCE_SNAPSHOT.txt`.

Genetic-response, phylogatR, and generic transferability-development code have been removed from this repository after the paper-specific test suite demonstrated that they are not required by the GEB analysis.

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

- `manuscript/` — main, blinded, title-page, claim-map, and submission files.
- `benchmarks/exploratory/` — frozen ecological result receipts supporting the manuscript.
- `docs/exploratory/` — frozen protocols and analysis rules.
- `scripts/` — acquisition, reconstruction, analysis, gate, test, and rendering entry points.
- `src/ttf/` — the small dependency surface inherited from the precursor repository and still required by this paper.
- `tests/` — offline paper-specific tests.
- `data/external/leptraits_consensus_v1.0.csv` — the exact LepTraits snapshot used by the reconstruction.
- `provenance/` — preserved workflow/provenance material needed to audit frozen executions.

## Reproducing the offline test suite

Python 3.11+ is required.

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
pytest -q
```

The GitHub Actions workflow `.github/workflows/paper-ci.yml` runs this suite without network-dependent tests. The ecology-only pruning state has passed this CI.

The `analysis` optional dependency group additionally provides SciPy, Shapely, and rasterio for spatial analyses; `figure` provides Matplotlib.

## Interpretation boundaries

The reconstructed host envelopes represent **potential regional larval-resource opportunity**, not confirmed butterfly occupancy.

Introduced host distributions are not interpreted as proof that host introduction caused butterfly range expansion. GBIF non-observation is not treated as true absence. Climate filtering is an association within the declared resource-opportunity design, not proof of physiological causation.

The independent host-breadth climate-release prediction was not supported (`partial Spearman rho = -0.166`, one-sided `p = 0.2237`, `n = 24`). No response-driven replacement predictor search is used to rescue that result.

## Submission state

The scientific story and double-anonymous manuscript package are complete. Remaining release tasks are administrative:

- choose and add the software license;
- freeze the archival release;
- mint the permanent archive DOI and insert it into citation/title-page metadata;
- complete final author/affiliation/CRediT/funding/conflict metadata.

Until the archival release is minted, cite the eventual versioned release rather than a moving branch.
