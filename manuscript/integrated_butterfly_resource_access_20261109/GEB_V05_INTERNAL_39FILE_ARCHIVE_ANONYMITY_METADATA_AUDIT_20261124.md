# GEB v0.5 — independent limited internal ZIP metadata and direct-identifier check (24 November 2026)

**Evidence scope: internal archive only. This is NOT double-anonymous reviewer clearance, third-party redistribution permission, or a public repository release.**

## Byte-exact inspected source

- GitHub Actions run: https://github.com/zuizui0223/chocho/actions/runs/44970859551
- Actions wrapper artifact ID: `14733606979`
- Outer downloaded artifact: **563,382 bytes**, SHA-256 `3a4beb33c72093aaf146ecca8da28650e2264fa221f7f9f8ab12bf8c9abbd09d`
- Contained internal ZIP `GEB_V05_INTERNAL_COMPOSITE_REPRO_NOT_FOR_JOURNAL_REVIEW.zip`: **560,519 bytes**, SHA-256 `40369ba5820f0e40ed56758f64165f951c82c6d3d64e5613175533d3c2055c95`
- **39 files** with no ZIP archive comment, no individual file comment, no per-file extra metadata, and uniform fixed timestamp **2026-11-21 00:00:00**.

The ZIP and its 39 members were inspected as bytes with the Python standard library's `zipfile` read interface. **No bundled code was executed**, no archive contents were redistributed, and no raw GBIF or private model panels were introduced.

## Explicit pattern checks on every member

| Inspected feature | Matches |
|---|---:|
| Repository account token from the project owner | 0 |
| Known direct author-name forms | 0 |
| Email-address pattern | 0 |
| Local `/home/`, `/Users/`, `/workspace/`, `/mnt/data/` paths | 0 |
| `github.com/` URLs | 0 |
| Common git/author-signature field patterns | 0 |
| Embedded external HTTP(S) URLs | 2 |
| DOI-like references | 18 |

**Interpretation:** This independent, explicitly scoped check found no *listed direct-identifier patterns*. It is not proof that an adversarial reader cannot recover identity from the paper's cited original-download DOI, data provenance, numbers, authorship of specific source code, institutional metadata outside this ZIP, or external web searches. Two external URLs and 18 DOI references remain and must be reviewed for indirect disclosure before any double-anonymous sharing.

## Identified release blockers

1. The archive's manifest explicitly says `INTERNAL_ONLY_NOT_ANONYMITY_OR_LICENCE_APPROVED`, `indirect_author_anonymity_certified = false`, `reviewer_access_URL_exists = false`, `permanent_archive_DOI_exists = false`, and `source_licensing_and_redistribution_clearance = false`. These boundaries are retained.
2. The `source_qa_receipts/GBIF155_original_source_and_missing_RF_panel_gate.json` is a **historical pre-recovery receipt** stating the forecast panel was missing, while separate subsequent 22-Nov source recovery and model refits have superseded this old current-status claim. The internal README labels it historical, but a reviewer-facing package must either omit it from active QA or separate it clearly from the latest source-status receipt.
3. The 39-file package replays limited **derived inputs** and **saved scores**, *not* a complete source-exact original GBIF→HOSTS/WCVP→all figures/SI pipeline. Original WPI conversion, full SI coverage and source-specific reuse rights are still missing.
4. This completed check did **not** examine exported DOCX/PDF, any future Dryad deposit, platform-side access-control settings, or reviewer-facing links; their anonymity/data-rights verification remains separate.

**Practical next release gate:** Assemble a replacement manuscript-current, version-locked source package with current audited source-status receipts, the exact component-to-figure/SI coverage, licence-by-file decisions, the approved omitted-original-inputs explanation and zero direct-identifying metadata. Check its actual uploaded bytes and access URL, not the draft's source tree. Do not publish or supply any private data until source-specific authorization is resolved.

The scientific main paper, model evidence and the 39-file reference ZIP were not changed by this inspection. PR #179 and #189 remain Draft.
