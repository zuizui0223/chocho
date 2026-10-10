# GEB initial-submission procedure

Checked against the current [Global Ecology and Biogeography author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html) on **2026-10-10**, and aligned with `manuscript/geb_initial_submission_manifest_v0.1.json`.

## Journal-facing requirements

GEB offers free-format initial submission, but requires a **structured abstract**, **double-anonymous** review materials, a **separate identifying title page**, a **separately uploaded PDF cover letter**, and **supporting information in separate files**. The blinded main text must contain a Data and Code Availability Statement and continuous line numbering. Data and code must be available to reviewers; a stable public archive is required before publication. An anonymous source/data/code ZIP uploaded directly as supplementary review material is the chosen peer-review access route here.

Do not put author names, identifying document properties, public repository URLs or identifying GitHub Actions URLs in files provided to anonymous reviewers. The separate title page and cover letter are identifying files handled by the journal editorial process.

## GBIF occurrence DOI — complete

The 53,434 frozen occurrence records were retrieved through the GBIF search API. An exact-ID archival download was requested and assigned [DOI 10.15468/dl.pp5nc9](https://doi.org/10.15468/dl.pp5nc9), download key `0008693-260928105237408`. GBIF returned 53,144 records on archival download, **290 fewer than the frozen analysis set**. This is disclosed in the blinded and unblinded Data and Code Availability sections. The DOI is **not an outstanding initial-submission blocker**.

## Six files to upload at initial submission

The authoritative six-item file inventory is `manuscript/geb_initial_submission_manifest_v0.1.json`.

1. **Blinded main manuscript DOCX.** Use `butterfly_specialization_GEB_blinded_review.docx` from the latest successful `build-blinded-review-docx.yml` job for the final PR head. It has continuous line numbers and must be checked for authorship metadata.
2. **Identifying title page.** Complete `manuscript/butterfly_specialization_geb_title_page_template_v0.2.md` with approved authors, affiliations, emails, relevant ORCIDs, exactly one corresponding author, CRediT, acknowledgements, funding and conflict-of-interest disclosures. Convert into the submission system's accepted editable file format.
3. **Identifying cover letter PDF.** Complete `manuscript/butterfly_specialization_geb_cover_letter_v0.2.md`, especially the corresponding-author signature block; export separately to PDF. GEB asks for a journal-interest paragraph under 250 words.
4. **Editable Supporting Information DOCX.** Upload `butterfly_specialization_GEB_supporting_information.docx` from the **same exact-head blinded review job**. It contains Tables S1–S9, rendered as 12 native Word tables; do not upload Markdown source in its place. The build checks its anonymous properties and page-extracted table text.
5. **Primary anonymized reproducibility ZIP.** Upload `butterfly_specialization_anonymous_review_bundle.zip` from the successful exact-head `build-anonymous-review-bundle.yml` workflow as supplementary peer-review material.
6. **Separate Table S9 anonymous source-addendum ZIP.** Upload `bce_clarke_anonymous_s9_reproducibility_v01.zip` from that **same workflow artifact**. It contains the candidate BCE/Clarke source links, WCVP taxonomic crosswalk and Level1 null-model sensitivity, with checksums. Do not merge or silently omit it.

**Figures:** The manuscript review builder embeds the main figures; the figure workflow also outputs PDF/PNG versions of three main figures and one supplementary figure. Upload the independent figure files only as required by the live submission portal; this does not replace any of the six manifest items.

## Exact-head quality gate

1. Complete author-specific information and obtain all coauthors' approval before submission.
2. Run `python scripts/submission_preflight.py`. With the frozen scientific inputs, expected remaining blockers concern author list/affiliations/emails/ORCIDs, corresponding-author signature, CRediT roles, acknowledgements, funding and conflicts; the GBIF DOI and S9 supporting ZIP are already complete.
3. Verify the final PR **head SHA** and require success at **that SHA** for the paper tests, blinded DOCX build, manuscript figures and anonymous-review-bundle workflows. If any file on the branch changes, re-check new exact-head artifacts rather than relying on an older green run.
4. Download and **retain locally** the review DOCX, Supporting DOCX, figure files and **both** anonymous ZIPs with checksums. GitHub Actions artifacts are configured to expire after **30 days**. Upload files themselves rather than identifying links to GitHub.
5. Inspect visible anonymity, accessibility, figure numbering and rendered Word tables; confirm upload slots and the journal portal's current prompts manually. CI cannot approve authorship or submit to the journal.

## Scientific freeze and later publication

The paper addresses reconstructed *potential* butterfly larval-resource geography, not observed competition, demographic rescue or butterfly community homogenization. Do not convert post-hoc exploratory link-identity, temporal or larval-performance results into confirmatory main-paper claims just to strengthen the submission. No additional scientific test is required to clear the existing submission gate.

Public repository deposit, stable public archive DOI, licence and `v1.0.0-butterfly` release/tag can be completed before publication, as documented in `docs/RELEASE_PROCEDURE.md` and `scripts/release_preflight.py`. This later release gate is separate from initial peer-review submission.
