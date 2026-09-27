# Release manifest — v1.0.0-butterfly candidate

This manifest defines the archival boundary for the butterfly-specialization ecology paper. It is a release checklist, not a new scientific analysis.

## Scientific state

- Claim map: `manuscript/butterfly_specialization_claim_map_v0.1.json`
- Main manuscript: `manuscript/butterfly_specialization_ecology_v0.1.md`
- Double-anonymous source: `manuscript/generated/butterfly_specialization_ecology_blinded_v0.1.md`
- Submission readiness: `manuscript/butterfly_specialization_submission_readiness_v0.1.json`
- Scientific results are frozen; no response-driven analysis is authorized for the release gate.

## Frozen input identities

- LepTraits snapshot: `data/external/leptraits_consensus_v1.0.csv`
  - SHA-256: `6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b`
- HOSTS mirror commit: `808e0b869f9ec1adf8efff87cf6a395adda103e0`
- rWCVPdata / WCVP v13 commit: `65bed76bae9d644ccb6ad200c05f9f5071d89e05`

## Vendored manuscript-figure source inputs

All four audited figure inputs are stored under `data/frozen/figure_sources/` and are verified as raw bytes in CI.

| File | SHA-256 |
|---|---|
| `s1_resource_descriptors.csv` | `894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc` |
| `anthropogenic_species_metrics.csv` | `7b2a20387d656dbfcd7b3c38ec6a7fe2505474ad51f783dd1f201eb6b7eabd30` |
| `host_contribution_metrics.csv` | `b0f16c5fa9a5b4a0842d6d23f69de7a1f5e938a4a96fea426c97df2dd73e63aa` |
| `independent_climate_primary_result.json` | `a73dca6e8b669f721d8f2745f27198d9b47a5cbd05dde6146cc8d4f1ddfbf79b` |

Historical artifact IDs and expiry dates remain recorded in `provenance/FIGURE_SOURCE_ARTIFACTS.md`, but the release no longer depends on their continued availability.

## Reproducibility environment

See `docs/REPRODUCIBILITY.md`.

- Python: 3.11+; reference CI uses 3.12.
- Clean paper tests install `.[test]`, including `shapely>=2,<3`.
- Figure CI installs `.[figure]` and regenerates all five figures from the vendored source inputs.
- WCVP/HOSTS sidecar scripts use base R only. Historical CI used `r-lib/actions/setup-r@v2` with `r-version: "release"`.

## Required release actions still outside the frozen scientific package

1. Choose the software license and add the corresponding `LICENSE` file.
2. Confirm final authors, affiliations, corresponding author, ORCIDs, CRediT roles, acknowledgements, funding and conflicts.
3. Create the final release commit and tag `v1.0.0-butterfly`.
4. Archive that exact tag in Zenodo and mint the public DOI.
5. Add the DOI/version to `CITATION.cff` and the identifying title-page/final-public metadata.
6. Supply an anonymized reviewer-access link in the blinded review manuscript; do not expose the public author-identifying archive there during double-anonymous review.
7. Export the editable review document with continuous line numbering and embedded figures, and remove author-identifying document properties.
8. Run final CI on the exact release commit.

## Release invariant

The release process may change administrative metadata, licensing, archive identifiers and submission formatting. It must not change frozen ecological statistics, species panels, thresholds, hypotheses, claim boundaries or result receipts.
