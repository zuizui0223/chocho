# GEB citation copyedit implementation plan

> **For agentic workers:** Use superpowers:executing-plans. The user's action-first continuation authorizes reversible preparation and verification without another design-approval round.

**Goal:** Deliver an internally consistent, source-verified APA6 journal-facing main article and separate Supporting Information while retaining the complete historical bibliography and scientific source evidence.

**Baseline:** Draft PR #179 head `cdefceca781e22b6a8864ad06a9e31a751c81460`; all 14 Actions runs at that head succeeded, including the 104-island HOSTS-masked flora audit. Work proceeds in an isolated shallow clone on local branch `paper/geb-citation-copyedit-20261121`.

## Design and constraints

- Apply guarded citation and reference corrections only in the journal exporter. Leave the source article, archived v0.3 references/appendices, numeric receipts and scientific figures byte-identical.
- Record full author metadata, DOI/version evidence, reference placement and each exact prose substitution in a versioned manifest. Keep unused historical references in the source archive; do not invent scientific uses to retain them in the journal bibliography.
- Correct EuDiS version 5 to its issued year 2023, distinguish the dataset from its accompanying article, and attribute the actual GBIF and geographic data sources.
- Preserve the WPI converted-file edition uncertainty, distinct historical GBIF acquisition gaps, unavailable 155-species model inputs, signed adverse results and other Stage A limits.
- APA6: three-to-five-author first citations name all authors; up to seven reference authors are printed, eight or more use first six plus last. Journal and volume are italic; titles use sentence case; reference entries have hanging indents and no bullets.
- Main body stays at or below 5,000 words; abstract stays at 292 words. Four main figures and two SI figures remain embedded, with eight SI appendices and anonymous document properties.
- No model refit, contributor messaging, merge or journal submission. Initial GEB Free Format permits consistent references; this copyedit is preparation, not an additional submission authorization.

## 1. Source and use audit

- [x] Confirm all prior CI runs succeeded and recover a valid checkout without altering prior work files.
- [x] Check official GEB guidance and source author/DOI metadata; independently audit main/SI citation usage.
- [x] Resolve actual OSCAR product identity and the GBIF download's prescribed citation; record any remaining source ambiguity explicitly.

## 2. Journal export implementation

- [x] Add a curated reference/citation manifest, exact source guards and a citation-coverage receipt.
- [x] Apply export-only text corrections and independent main/SI reference lists; retain all historical source entries.
- [x] Connect source guards, hanging reference formatting and workflow triggers to the existing exporter.

## 3. Verification and delivery

- [x] Run the immutable scientific/editorial validator, current document builder, citation matching and scientific-number preservation checks.
- [x] Render and inspect all main/SI pages; resolve any material source, layout or independent-review findings.
- [x] Update the readiness audit, deliver verified Word/PDF revisions, and prepare the reviewed change for the existing Draft branch; publication uses a normal, expected-head-guarded update.

## Execution record

- The exporter verifies 24 main references and nine SI references against actual callouts. All 39 current historical source entries and the original archive remain untouched.
- Main body 4,979 words; abstract 292. All 381 non-citation numerical tokens remain in order. The immutable scientific/figure validator passes.
- Final local PDFs: 21 main pages, 11 SI pages. Every page was viewed; the Figure S5.1 lead-in pagination finding was fixed. Re-review covered the four changed SI pages; the remaining 28 page images were byte-identical.
- Independent source/code and visual reviews are closed. The four delivered PDF/DOCX files were updated successfully as version 1 of their existing file identities.
- Stage A remains active. Source-access and author gates remain explicit; no coauthor review, merge or journal submission is authorized by this copyedit.

Verification uses the existing source checks plus meaningful citation/data guards and rendered-page review. A presentation-only copyedit does not require implementation-mirroring unit tests.
