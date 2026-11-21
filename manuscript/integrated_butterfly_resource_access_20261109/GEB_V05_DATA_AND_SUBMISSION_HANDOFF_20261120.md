# GEB v0.5 — data/code submission handoff, 20 November 2026

**Status: STAGE A PREPARATION; not ready for coauthor review or journal submission.** Updated 21 November 2026. Applies only to *Ports, the sea and two layers of island floras reshape butterfly larval-host opportunities* (Draft PR #179, `paper/geb-single-island-port-endemics-20261118`), not the earlier *Butterfly specialization* v0.2 Research Article. Complete non-coauthor source, permission and anonymous reviewer-access preparation before requesting coauthor review.

## Confirmed / current outputs

- Single source current paper: `ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md`; Discussion 4.2 corrected at commit `a8a97e19`: EuDiS 22-trait-island morphology is not a native-flora marine-isolation model.
- New **v0.5** journal-facing manuscript and separate SI: `scripts/build_GEB_v05_blinded_main_and_SI.py` creates four-figure anonymized main DOCX and S1–S8 independent SI DOCX with two embedded figures, and an honest source-gate preflight receipt. Verify the actual latest build artifact before sending.
- **November 21 citation revision:** source-verified APA6 export has 4,987 main-body words, 292 abstract words, 24 cited main references and nine SI references. All 39 historical references and scientific source files remain unchanged. The curated export manifest is `GEB_REFERENCE_COPYEDIT_V01.json`; source/visual QA is recorded in `GEB_CITATION_COPYEDIT_QA_20261121.md`.
- The source-backed main 239-butterfly resource result, 126 botanical islands, 104-island undirected marine result and independent 155/224 butterfly records are distinct cohorts. Do **not** pretend one archived GBIF DOI covers all of them.
- Journal-scope editorial cover letter for coauthor feedback: `GEB_V05_COVER_LETTER_DRAFT_20261120.md`; its 191-word GEB relevance paragraph is drafted, author signature/assurances unverified.
- **Critical legacy trap:** `.github/workflows/build-anonymous-review-bundle.yml`, `scripts/build_anonymous_review_bundle.py`, `scripts/submission_preflight.py` and `docs/GEB_INITIAL_SUBMISSION_PROCEDURE.md` explicitly target **v0.2 butterfly_specialization**, NOT the current island-port-sea paper. Green status from these jobs does **not** certify a complete present-v0.5 reviewer data/code archive.

## Exact data-citation ledger

| Input / output family | Confirmed provenance | Remaining gate |
|---|---|---|
| Larval-host plant links | NHM HOSTS DOI `10.5519/havt50xw`; frozen mirror commit `808e0b8` | Respect source terms / exact network transform |
| Plant range & origin | WCVP v13, 21 May 2024 DOI `10.34885/0yex-xv26` (not misleading v10 DOI); pinned rWCVPdata commit `65bed76` | Include required attribution and licence notice |
| Original GBIF download for 155-species forecast and 224-species community analysis | `10.15468/dl.umqmsf`; GBIF.org-recommended citation dated 29 October 2026; original forecast README and pinned community annual-derivation source | Same 14,637,152-record download, distinct outcomes/cohorts; package permitted derived inputs without conflating their model replays or the separate historical/port acquisitions |
| Original 155-species forecast model workspace | 54 source-verified output assets and 11 model/panel scripts; shared GBIF download identified | Recover the 842,964-row `work/temporal/forecast/species_cell_period_panel.parquet` and 12 associated inputs before claiming an independent refit |
| Previous 31-butterfly / 53,434 occurrence validation | `10.15468/dl.pp5nc9`, records returned 53,144 (290 missing from the archived re-download) | **Out of scope** for the integrated 224 and 225/51 samples; never substitute this DOI |
| 225 candidate butterflies × 51 landmasses (first-index source) | Query/shard and checksum receipts exist; the source record contributing datasetKeys/DOIs are not complete | Reconstruct source-specific publisher attribution and register a distinct derived dataset/download; record drift if reconstruction differs |
| Historical butterfly–island–year source (83,573 rows) | Privately recovered, source-checksum verified and independently decoded in the PR #189 source audit; original bytes remain outside public Git | Reviewer redistribution and original GBIF contributor attribution remain uncleared; distinguish recovered private bytes from missing 155-species inputs |
| 1,596-port plant queries | Completed port-buffer acquisition receipts; 25,742 derived candidate plant×port pairs recorded | Recover original GBIF record-level or datasetKey attribution for **this** port cohort; do not label the 224-butterfly DOI a port DOI |
| WPI 2019 / GSHHG 2.3.7 | Original publication/version/source checksums identified | Reconcile exact WPI original→converted GeoJSON chain; preserve rights/licence notices |

The derived-dataset DOI may be created **only** from an actual source attribution table linked to a persistent deposit, using an authenticated GBIF.org account. Never fabricate or retroactively mislabel a DOI.

## Reviewer access and stable publication — a clean route

1. Build a **new current-v0.5** anonymous data+code package, with selected original/derived input tables (where redistribution rights allow), transformations, version-locked scripts, output receipts, environment instructions and source-to-figure instructions. At minimum it must support all main analyses and verify that each omitted original file has a specific legal or data-provenance explanation; old v0.2 ZIP fails this requirement.
2. Scan both source content **and file metadata** for names, repository accounts, URLs and identifying commit IDs. Verify on independent extraction that all current data source receipts, re-run commands and Figure 1–4 / SI output dependencies resolve. Do not claim full reproducibility if restricted raw data or scripts are absent.
3. For initial GEB review, upload the anonymous source package as separate supplementary files **or** use a properly anonymous Dryad private-for-peer-review deposit; use only the link/actual uploaded package that exists. GEB permits supplementary access at review, but requires a stable public repository for eventual publication. Native GitHub is **not** itself a stable publication archive.
4. Obtain exact separate GBIF citations for the historical butterfly and port plant sources; confirm third-party media/port-data rights. Make the final Data and Code Availability statement identify accurately which inputs are externally cited, what is shared for review, and what remains legitimately restricted.
5. Collect the complete identifying title page from the approved authors (author order, full affiliations, one corresponding email/ORCID, funding, CRediT, conflicts, acknowledgments), coauthor approval and cover-letter signature. The November 21 journal-export reference copyedit is complete; exact attribution for the separate historical acquisition families remains step 4. No placeholders may be uploaded as approved declarations.

## Source and policy cross-check

- GEB 2026 author rules: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html — double-anonymous manuscript + separate title page; structured abstract ≤300; approximately 5,000 main words; cover-letter PDF with <250-word relevance paragraph; review-accessible data and code; stable archival deposit for publication; free-format references permitted at initial submission if consistent.
- GBIF citation and derived datasets: https://www.gbif.org/citation-guidelines and https://www.gbif.org/derived-dataset/about — authenticated derived registration needs contributing dataset identities/counts and a permanent dataset URL.
- Current notebook/legacy manuscripts are NOT authoritative evidence of new reviewer-archive completeness. The new journal-facing preflight must refer to **v0.5**, not the old `submission_preflight.py`.
