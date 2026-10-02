# Release manifest — v1.0.0-butterfly candidate

This manifest defines the archival boundary for the current butterfly resource-geography paper. It is a release checklist, not a new scientific analysis.

Release identifiers are fixed as package/CFF version `1.0.0` and intended GitHub tag `v1.0.0-butterfly`. The DOI and release date remain unset until an actual archive/release exists.

## Current scientific state

- Claim map: `manuscript/butterfly_specialization_claim_map_v0.2.json`
- Main manuscript: `manuscript/butterfly_specialization_ecology_v0.2.md`
- Double-anonymous source: `manuscript/generated/butterfly_specialization_ecology_blinded_v0.2.md`
- Supplement: `manuscript/butterfly_specialization_supplement_v0.2.md`
- Submission readiness: `manuscript/butterfly_specialization_submission_readiness_v0.2.json`
- GEB checklist: `manuscript/butterfly_specialization_geb_submission_checklist_v0.2.json`

The current v0.2 inference boundary is frozen for release: introduced host distributions add 14,553 butterfly × WGSRPD3 units (+54.9% aggregate coverage); added opportunity is strongly concentrated across host plants; host-family breadth has a near-zero point estimate but does not pass the strict post-hoc ±0.10 equivalence test; occurrence validation is secondary; network prominence remains exploratory/SI only.

## Frozen input identities

- LepTraits snapshot: `data/external/leptraits_consensus_v1.0.csv`
  - SHA-256: `6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b`
- HOSTS mirror commit: `808e0b869f9ec1adf8efff87cf6a395adda103e0`
- rWCVPdata / WCVP v13 commit: `65bed76bae9d644ccb6ad200c05f9f5071d89e05`

The historical reconstruction inputs remain hash-pinned under `data/frozen/figure_sources/` and documented in `provenance/FIGURE_SOURCE_ARTIFACTS.md`.

## Current v0.2 figure/review evidence surface

The current renderer is `scripts/render_butterfly_specialization_v02_figures.py` and produces three main figures plus Supplementary Figure S1. It consumes git-tracked frozen sources and result receipts, including:

- `data/frozen/figure_sources/anthropogenic_species_metrics.csv`
- `data/frozen/figure_sources/occurrence_resource_validation_species.csv`
- `data/frozen/figure_sources/climate_distance_sensitivity_species.csv`
- `provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_expansion_equivalence_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_occurrence_species_robustness_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_expansion_ceiling_sensitivity_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_regional_robustness_v0.1.json`
- `provenance/reviewer_defenses/results/butterfly_climate_effect_size_v0.1.json`

The exact release commit in Git is the archival identity for these tracked files; current paper CI regenerates all four PDF and PNG display files from them.

## Reproducibility environment

See `docs/REPRODUCIBILITY.md`.

- Python: 3.11+; reference CI uses 3.12.
- Clean paper tests install `.[test]`.
- Figure CI installs `.[figure]` and regenerates the current v0.2 figure set.
- WCVP/HOSTS sidecar reconstruction uses the pinned upstream commits and base-R scripts documented in the reproducibility guide.

## Review-package preparation completed

- `scripts/build_blinded_review_docx.py` and `.github/workflows/build-blinded-review-docx.yml` reproducibly build the editable double-anonymous review DOCX from the current blinded v0.2 Markdown plus the four regenerated display files.
- The DOCX builder checks continuous visible line numbering, embedded figures and scrubbed creator/lastModifiedBy metadata.
- `scripts/build_anonymous_review_bundle.py` and `.github/workflows/build-anonymous-review-bundle.yml` build the de-identified reviewer code bundle.
- Internal response/decision memos are excluded from the anonymous review bundle.
- The anonymous bundle runs the scientific tests and regenerates the current v0.2 figures.

## Required release actions still outside the scientific package

Follow `docs/RELEASE_PROCEDURE.md`.

1. Choose the software license and add the corresponding `LICENSE` file.
2. Confirm final authors, affiliations, corresponding author, ORCIDs, CRediT roles, acknowledgements, funding and conflicts.
3. Create a manual Zenodo draft and reserve its DOI before publication.
4. Add that reserved DOI to `CITATION.cff` and the identifying title-page/final-public metadata.
5. Keep the peer-review bundle separate from the public archive: for initial GEB review, upload the de-identified reviewer ZIP directly as supplementary review material.
6. Run the current v0.2 paper/review-package CI on the exact public release commit and require `scripts/release_preflight.py` to report READY.
7. Tag that exact green commit as `v1.0.0-butterfly`.
8. Upload the exact tag archive to the existing Zenodo draft and publish it, activating the reserved DOI.

## Release invariant

Release administration may change author metadata, licensing, archive identifiers and submission formatting. It must not change the frozen ecological statistics, species panels, thresholds, inference boundaries or result receipts without reopening the scientific audit.
