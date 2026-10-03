# GEB initial-submission procedure

Checked against the current *Global Ecology and Biogeography* author guidelines on 2026-10-02.

## What is required for initial peer review

GEB uses double-anonymous review. The initial submission therefore separates the identifying title page from the blinded main manuscript.

The journal requires data and code supporting the paper to be accessible during peer review. A stable public repository is required for publication, but peer-review access may be provided through supplementary materials. This repository therefore uses the already-tested anonymous review bundle as a file uploaded directly with the submission, rather than requiring an external anonymous reviewer URL.

## GBIF occurrence DOI — complete

The 53,434 occurrence records used in the secondary occurrence validation were originally retrieved through the GBIF occurrence search API, so the historical API retrieval has no automatic DOI. The exact-record GBIF occurrence download has been created: **https://doi.org/10.15468/dl.pp5nc9** (download key `0008693-260928105237408`).

This DOI is recorded in both manuscript surfaces and cites the external occurrence data actually used; it is distinct from the public DOI for this study's code/reproducibility archive, which remains deferred until publication.

## Upload order

1. **Blinded main manuscript** — use the line-numbered DOCX produced by `.github/workflows/build-blinded-review-docx.yml`. Confirm that the document properties and visible text contain no author identity.
2. **Identifying title page** — fill `manuscript/butterfly_specialization_geb_title_page_template_v0.2.md` with the final author list, affiliations, emails, ORCIDs, one corresponding author, CRediT contributions, acknowledgements, funding and conflict-of-interest statement.
3. **Cover letter** — fill the corresponding-author signature in `manuscript/butterfly_specialization_geb_cover_letter_v0.2.md` and upload it separately as a PDF.
4. **Supporting Information** — upload `manuscript/butterfly_specialization_supplement_v0.2.md` in the portal's supporting-information slot after rendering to the desired submission format.
5. **Anonymous data-and-code review supplement** — upload the CI artifact `butterfly_specialization_anonymous_review_bundle.zip` produced by `.github/workflows/build-anonymous-review-bundle.yml`. Treat this as supplementary review material/data-code review archive.

Do **not** expose the public GitHub repository in the blinded manuscript or anonymous bundle.

## Not required before the initial submission

The following can be completed after peer review but before publication/public release:

- choosing and applying the public archive licence;
- depositing the stable public data/code archive;
- inserting the persistent public archive DOI;
- creating/tagging the public `v1.0.0-butterfly` release.

These remain requirements of the public/release workflow and are checked separately by `scripts/release_preflight.py`.

## Initial-submission gate

Run:

```bash
python scripts/submission_preflight.py
```

This gate intentionally checks the journal-facing initial submission rather than the later public archival release. With the current repository state, the expected remaining blockers are the exact-record GBIF occurrence-download DOI plus author-specific metadata and declarations.

After filling those items, rebuild the blinded DOCX and anonymous review bundle on the exact commit to be submitted and require all submission workflows to pass.

## Scientific freeze

No additional ecological analysis is required for submission. Do not reopen the unsupported climate-release hypothesis, the within-butterfly portfolio architecture, or response-driven predictor searches merely to clear submission metadata.
