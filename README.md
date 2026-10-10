# chocho — Butterfly resource geography under plant globalization

Reproducible code, frozen analysis receipts, and manuscript materials for:

> **Plant globalization expands and homogenizes butterfly larval-resource geography**

Target journal: *Global Ecology and Biogeography*.

## Scientific result

The paper asks whether plant globalization changes only the **amount** of butterfly larval-resource opportunity or also its **community structure**. Known butterfly–host identities are held fixed while host geography changes from native-only to contemporary distributions.

The current v0.2 evidence supports four linked ecological results:

1. Introduced host ranges expand reconstructed resource opportunity for **206/239 butterflies**, increasing aggregate butterfly × WGSRPD3 coverage by **54.9%**.
2. Plant redistribution **homogenizes resource geography**: mean Jaccard similarity among regional butterfly resource assemblages rises from **0.277 to 0.462 (+66.8%)**, and mean overlap among butterfly resource envelopes rises from **0.222 to 0.328 (+47.4%)**. Both exceed conservative fixed-margin null expectations (**p = 0.002**).
3. Exact shared-host geography expands from **62,473 to 141,885 butterfly-pair × region units (+127.1%)**; **920/986** host-sharing butterfly pairs gain new shared-resource regions. This is potential resource co-use, not observed competition.
4. A secondary occurrence analysis recovers **66/115** outside-native butterfly × region observations, linking some reconstructed resource opportunity to contemporary presence.

Added opportunity remains concentrated among host plants (38/670 species account for half), yet **58.9%** of added butterfly × region units are supported by only one contributing introduced host. Diet breadth is retained as a secondary modifier: family-level host breadth has little relationship to proportional expansion (rho = 0.008; bootstrap 95% CI -0.111 to 0.128).

See `manuscript/butterfly_specialization_claim_map_v0.2.json` for exact claim boundaries, especially the distinction between resource sharing and realized competition.


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

The full-panel diet-breadth association is near zero but is not presented as proven equivalence: a post-hoc ±0.10 TOST is narrowly inconclusive (`p = 0.078`). The independent host-breadth climate-release prediction was also not supported (`partial Spearman rho = -0.166`, one-sided `p = 0.2237`, `n = 24`). No response-driven replacement predictor search is used to rescue either result.

## Submission state

The scientific v0.2 manuscript is frozen for initial submission. The privacy-scrubbed line-numbered DOCX and de-identified data-and-code review bundle are generated and CI-audited.

For **initial GEB peer review**, the anonymous review bundle can be uploaded directly as supplementary review material. An external anonymous reviewer URL, public archive DOI and public GitHub release are therefore not prerequisites for initial submission. The remaining initial-submission blockers are author-specific metadata and declarations:

- final author list, affiliations, emails and ORCIDs, with exactly one corresponding author;
- CRediT contributions;
- acknowledgements, funding and conflict-of-interest statements;
- corresponding-author signature in the cover letter;
- exact-commit green rebuild of the blinded DOCX and anonymous review-bundle ZIP.

Run the initial-submission gate with:

```bash
python scripts/submission_preflight.py
```

See `docs/GEB_INITIAL_SUBMISSION_PROCEDURE.md` and `manuscript/geb_initial_submission_manifest_v0.1.json` for the upload mapping.

Public archival tasks are tracked separately and may be completed before publication: choose/apply the archive licence, deposit a stable public data/code archive, record its DOI, and create/tag `v1.0.0-butterfly`. Those later tasks are checked by:

```bash
python scripts/release_preflight.py
```

Neither preflight requests nor authorizes additional ecological analyses.
