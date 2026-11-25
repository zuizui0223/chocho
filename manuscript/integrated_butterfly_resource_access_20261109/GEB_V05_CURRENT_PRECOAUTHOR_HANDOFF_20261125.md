# GEB v0.5 — current preparation-first handoff, 25 November 2026

**Stage A: technical preparation; NOT ready for coauthor review, external release or submission.** This dated record supersedes the older manuscript title, word counts, archive sizes and draft letter in the 20–24 November handoff files, but does not supersede their unresolved provenance/rights conditions. PR #179 remains unchanged; PR #192 is the current manuscript candidate, PR #193 its observation-source audit, and PR #189 the internal archive audit. No PR is marked Ready, merged or submitted.

## 1. Freeze the exact current candidate

| Item | Current verified identity |
| --- | --- |
| Paper candidate | PR #192, commit 7e182df03dc3f031688159e4ae0ac2914acd0735 |
| Full title | Including introduced plants weakens the island isolation gradient in mapped butterfly host resources |
| Export | [CI #45568373174, artifact #15000204579](https://github.com/zuizui0223/chocho/actions/runs/45568373174/artifacts/15000204579), ZIP SHA-256 508d6494d5155ab5f5e1645fd848228a902a3e4229759231d426a3af6eaba447 |
| Article/SI | 4,318 body words, 280 abstract words; 19/16 pages; 4/2 embedded figures; 27/9 cited references; 39 historical references preserved. Six exact-paper-head workflows passed. Overall contact-sheet scan, not a detailed manual page audit. |
| Source delivery | [CI #45572084747, artifact #15003375426](https://github.com/zuizui0223/chocho/actions/runs/45572084747/artifacts/15003375426), outer SHA-256 3a9bce3a344068da68e293d66676aa94f038c38c0efad0935c5bc4568368a92c; inner SHA-256 d9af86284a8cec0b82257c096094ad21c91724c8b2ae0abeeffef0d004a8a577 |
| Exact internal archive | 240 members, all previous 223 preserved byte-for-byte; new 17 members, 26 source families. All 240 per-file hashes/lengths and authorization flags verified; 0 authorized. |
| New statistical checks | 1,516 S3.6 and 1,284 S5.5 all224 extracted derived-panel checks passed. Raw-source rebuild and without-monarch replay are not repeated by this archive test; original CI separately retained. |

Use the current artifact, not the obsolete v0.2 review-bundle workflow or old 39/192/223-member package. The existing 240-member root wrapper/README is historical and does not execute added S3.5/S4.1/S3.6/S5.5 components automatically. Separate replay instructions and immutable receipts accompany those components. A complete reviewer-facing runner and current release documentation must be prepared only after the file-specific route is settled.

## 2. New scoped anonymity audit

Local downloaded exact archive and paper ZIP hashes match the identities above. Scanned all 240 inner members as text, including decompressed contents of four gzip files, for the known author repository account string, known author-name variants, absolute personal-home path patterns, GitHub Actions run URLs, ORCID and email-like strings. Zero files match the first four or ORCID patterns; six email-like matches in two Python files are matrix-multiplication expressions rather than contact addresses. Twelve files contain 14 GitHub URLs to upstream third-party owners; 43 files contain 70 DOI-pattern matches. ZIP archive and member comments/extra fields absent. Both exact DOCX core creator and lastModifiedBy fields empty. This is a limited pattern and metadata inspection, not a proof of double anonymity. Public manuscript history, source DOI chains, cited records, internal filenames and contextual clues may identify authors; outer audit/register material and links are not reviewer-ready. Do not alter source provenance, strip necessary attribution or change repository visibility automatically.

## 3. File-level decision queue

Use EXACT_240_FILE_SOURCE_DELIVERY_DECISIONS.csv and SOURCE_REGISTRY.json from artifact #15003375426, not a blanket release licence. All 240 members remain HOLD.

| Current register group | Files | Next action |
| --- | ---: | --- |
| Author approval/source summary review | 90 | Verify content and dependencies; confirm author authority to distribute. |
| Author code approval/provenance review | 58 | Confirm code licence and third-party notices; preserve exact original scripts separately from reviewer wrappers. |
| Unresolved source dependency | 38 | Resolve exact source-specific conditions or propose permitted external retrieval/access. |
| Mixed/restricted derivative review | 10 | Assess material and applicable terms; do not assume CC0 compatibility. |
| Historical attribution | 6 | Preserve original query scope and limitation; obtain advice on incomplete publisher attribution. |
| Added S3.5/S4.1/S3.6/S5.5 review | 38 | Inspect the 15 + 6 + 17 new members and dependency mapping; retain original evidence, verify reviewer anonymity and rights. |

Totals: 93 metadata/summary, 68 code, 60 tables, 19 validation; 240 files. Dependency mappings overlap: AUTHOR 240, TDWG 236, WCVP 220, HOSTS 190, WPI 103, LEP 100, GBIF 70, HIST 39, CIRAD 9. Counts describe conservative register dependencies, not legal findings that every derivative inherits upstream terms.

Carry forward the 24-Nov repository-route and source-delivery decisions: exact HOSTS/WCVP source identities; original shared 155/224 GBIF DOI 10.15468/dl.umqmsf with 1,716 contributors; incomplete original historical/1596-port provider attribution; exact WPI converted layer lineage; source-specific permissions (including mixed GBIF record licences, CIRAD CC BY-NC-SA and GloNAF/TDWG source issues); original155 independent geographic numerical reproduction and two sensitivity discrepancies. No invented DOI, no substitutions, no raw third-party release, no automatic relicensing. A stable anonymized reviewer-access route must be actually created and tested before Stage A is closed. The earlier draft editorial and GBIF clarification requests remain unsent and need author authorization.

## 4. Separate identifying title-page template (NOT part of the blinded files)

**Full title:** Including introduced plants weakens the island isolation gradient in mapped butterfly host resources

**Running title:** Island butterfly host resources

**Authors and order / affiliations / corresponding author and email / ORCID:** [Author-supplied and coauthor-approved]

**Funding / acknowledgments / CRediT contributions / conflicts / relevant ethics or permits:** [Actual author declarations required; do not infer]

**Data and code availability:** [Insert actual verified anonymized peer-review route and approved archival/rights statements when available; no invented URL or DOI]

## 5. Revised GEB cover letter — draft only, not sent

Dear Editors of *Global Ecology and Biogeography*,

Please consider our manuscript, “Including introduced plants weakens the island isolation gradient in mapped butterfly host resources,” for publication as a Research Article.

Island isolation shapes colonization, but consumer resources depend on the geography of interaction partners as well as distance. Holding documented larval-host identities fixed for 239 butterflies across 369 botanical regions, we show that native-host opportunities decline with island isolation more strongly than host-family presence does. Including introduced plants makes the contemporary resource gradient nearly flat, although the paired slope contrast remains uncertain under conservative geographic inference. More than half of mapped island host opportunities depend on introduced plants alone. Large proportional gains reflect a small native baseline, and just 38 of 670 contributing plants account for half of new global opportunities. In a shared 224-species panel, independent butterfly records include 193 species–island-region combinations with mapped introduced hosts but no mapped native host, of which 143 recur across years. Recorded introduced-only shares also increase with isolation, though stricter recording-effort subsets give uncertain estimates. Together with independently documented island-endemic larval feeding and the monarch's historical Pacific expansion, these results connect plant redistribution to consumer geography without equating mapped resources or occurrence with establishment. This interaction-based perspective complements classical island biogeography and identifies why the geography of familiar food plants matters for consumers in a human-modified world.

The main advance is a consumer-resource geographic perspective based on fixed documented host identities, rather than inference of realized butterfly range expansion from plant distributions alone.

Sincerely,

[Approved corresponding author, affiliation, email]

Internal checklist (remove from any eventual submitted letter): confirm originality, no simultaneous submission, all-author approval, complete declarations, data/code access, source permissions and any suggested-reviewer conflict statement. The significance paragraph is 194 whitespace-delimited words, below the previously checked 250-word limit. No signature or compliance assertion is made. Only export the separate letter PDF once approved.

## 6. Stage A exit gates, in order

1. Resolve source/derivative permissions, author code/data authority and incomplete original-provider attribution; retain unresolved cases as explicit editorial questions rather than claiming clearance.
2. Construct a file-specific permitted stable archive and current independent replay instructions; test access and indirect anonymity from a reviewer perspective.
3. Close or explicitly document remaining original155 numerical discrepancies and exact source limitations; no new biological hypothesis is necessary.
4. Obtain the author's identifying information/declarations and decision to proceed to coauthor review. Prepare the exact manuscript/title/cover letter package from the agreed candidate.
5. Only after actual coauthor approval and the author's submission decision may the manuscript be integrated for submission and transmitted.

Current automated export receipt still explicitly says DOCUMENT_STRUCTURE_VALID; JOURNAL_SUBMISSION_BLOCKED. This handoff does not change manuscript numerical claims or grant approval, access, merge or publication.
