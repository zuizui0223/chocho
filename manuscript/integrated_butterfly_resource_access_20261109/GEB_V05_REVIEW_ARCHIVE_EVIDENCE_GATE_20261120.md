# GEB v0.5 reviewer-data source gate — current paper only

**20 November 2026. Internal coauthor handoff; NOT a peer-review-ready anonymized deposit.**
The current scientific manuscript remains Draft PR #179 at snapshot
`99b876421af2b4f51e58569862994622a02717d7`.
This gate belongs to new Draft PR #189 and changes **no manuscript result, figure, or claim**.

## Executed, inspectable audit

- Independent review-pack audit code: `scripts/audit_GEB_v05_review_source_bundle_20261120.py`.
- **Final source-ZIP extraction/replay CI (SUCCESS):** https://github.com/zuizui0223/chocho/actions/runs/44678903555 ; [artifact #14603227495](https://github.com/zuizui0223/chocho/actions/runs/44678903555/artifacts/14603227495) (90-day retention, internal/non-anonymous).
- Initial inventory-only CI: https://github.com/zuizui0223/chocho/actions/runs/44678394840
- Previous zipped internal preview: https://github.com/zuizui0223/chocho/actions/runs/44678230957/artifacts/14603351378 (90-day GitHub retention, not archival DOI).
- Two Git checkouts are essential: **paper branch** and pinned original numerical source commit
  `3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250`. Selected original figure model receipts are **not** at their original paths in the paper checkout alone. A single-branch ZIP could conceal this gap.
- **35 input/code/numeric files** were present, hashed SHA256 and copied into an INTERNAL evidence preview; four current Figure SVGs and two actual SI SVGs parsed, and the eight figure numeric source paths were found at the frozen source commit. Fixed 239×126 and 123-active-island isolation panel values were cross-checked.
- The preview contains selected code, figure receipts, SVGs, and derived summary tables. It **does not** contain all third-party original records or claim end-to-end raw-data reproduction. ZIP hash and machine-readable per-file inventory are emitted by the Action.
- **Subanalysis actually replayed from the extracted ZIP:** the frozen 126-island fixed-239 source model and 2,000-km Moran/MEM10 analysis, after downloading and SHA256-verifying the pinned official TDWG Level3 boundary GeoJSON. Recomputed port/flora, native and current host-slope contrast, and introduced-only-resource share coefficients matched the article's frozen numeric receipt to 1e−4. This independently tests **one** model family from packaged source tables, not the historical plant/GBIF record extraction or the 155/224 butterfly analyses. CI label: `ZIP_EXTRACTED_FIXED239_AND_MEM10_SUBANALYSIS_REPLAY_SOURCE_EXACT`.
- Three explicitly referenced **original input paths were absent from both inspected Git trees** (paper snapshot `99b8764` and pinned source `3dbf038`): `work/barriers/ports_2019.geojson`, `work/barriers/global_original_species_island_year.parquet`, and `data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json`. They could exist in private archives or old Action artifacts, but are not reproducibly checked out from either Git commit.

## Explicit author-anonymity findings

The unmodified internal candidate contains the owner/repository account token in **five** source files:

1. `manuscript/integrated_butterfly_resource_access_20261109/ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md`
2. `manuscript/integrated_butterfly_resource_access_20261109/GEB_SUPPORTING_INFORMATION_V05.md`
3. `manuscript/integrated_butterfly_resource_access_20261109/FIGURE_CAPTIONS_V03.md`
4. `figures/integrated_island_port_butterfly_v03/FIGURE_V04_SOURCE_NUMBERS.json`
5. `scripts/build_GEB_v05_blinded_main_and_SI.py`

These are identification risks **inside the internal source preview**, not a claim that author identities leaked in the separately generated **blinded journal-facing DOCX/PDF**. That builder removes visible repository URLs and has its own independent metadata/page QA. For source-code sharing, preserve immutable provenance **internally** while producing a **separate verifiably anonymous release copy** with accurate citation substitutions, not fabricated anonymous Git history. Re-run content and ZIP/Word metadata checks on actual extracted files.

## Data family citation ledger

| Cohort / source | Already supported source | Still needed for independent review |
|---|---|---|
| 239-butterfly fixed known hosts | NHM HOSTS DOI `10.5519/havt50xw`; pinned HOSTS mirror | Source terms and exact original extracted host-edge transformation |
| Native/introduced vascular flora | Kew WCVP **v13, 21 May 2024** DOI `10.34885/0yex-xv26` (not the wrapper's older v10 DOI) | Dataset licence attribution and a permitted reproducible extract/retrieval route |
| Observed community, separate 224 butterflies | Original GBIF occurrence DOI `10.15468/dl.umqmsf`, pinned annual aggregation | Original download and exact processing scripts/review-accessible permitted derivatives |
| Original 155-species later-record prediction | Source-pinned output receipts, model code | Confirm which original occurrence downloads and source datasets support this separate cohort, then give exact citation and reproduction route |
| 225-candidate butterflies × 51 landmasses (historical/auxiliary) | Frozen aggregate query/shard and checksums in source notes | Record-level original GBIF publisher `datasetKey`/DOI attribution **not complete**; new snapshot must be labelled as new, not retroactive proof |
| 1,596 port botanical-buffer queries | 25,742 derived candidate plant×port combinations and query receipts noted | Exact GBIF contributing datasets / original occurrence provenance **not complete** |
| Historical butterfly–island–year panel | Original private **83,573-row** Parquet referenced | Licensed original or permissible derived source and explicit restricted-access explanation |
| WPI ports, GSHHG | WPI 2019 Publication 150 (27th ed.), GSHHG 2.3.7 DOI `10.5281/zenodo.7007502` | WPI original input → converted GeoJSON identity/rights chain remains open; GSHHG pin must be preserved |

**DOI scope rule:** `10.15468/dl.umqmsf` is exclusively the 224-butterfly GBIF acquisition; `10.15468/dl.pp5nc9` belongs to an older 31-butterfly pilot. Neither is a citation for the 155, 225×51, or 1,596-port datasets without an actual source match. No GBIF derived-dataset DOI has been registered in this audit.

## November 20 — complete recovered original year-table quality audit

The source ZIP's three previously missing files were recovered **privately**, but remain intentionally absent from this public Git tree. Original WPI coordinates independently reproduce all 1,275 archived distances for the 51 island groups. This establishes port geometry identity, not actual ship transport.

**Important update to the earlier footer-only inspection:** The original SHA256-pinned private Parquet was fully decoded: **83,573 distinct species×landmass×year rows, 224 species, 3,384 landmass IDs, 2001–2024, and 14,425,536 indexed record counts**, six required fields with no nulls or duplicate keys. `landmass_id=0` is an actual source ID for the largest-area landmass (50,654,050.6945 km²), not a missing value: it includes **2,636 rows and 4,386,686 record counts (30.41%)**. Rejecting 0 would corrupt the full source, but **the original explicitly selected 51-island roster excludes it**, so 51-island results are not affected by that invalid general filter. Existing source scripts already preserve 0 outside their explicit island selection. An offline private-source validator was added as `scripts/audit_chocho_historical_parquet_20261120.py`, with a GitHub Actions **self-test permitting 0 and rejecting negative IDs**.

**Independent historical biology-index replay:** The original 51-island year-table slice is **9,739 source annual rows / 129 observed species**. Rechecking the 14 port-representative-distance pairs under 200 km reproduces **41** source-defined shared first-index combinations in 2011–2022 and **21** earlier-island reindexing-before-later-initial-index patterns. All 14 pair-specific counts match the original pinned pair CSV. The five previously selected host-associated links reproduce **5/5** original year lists and time ordering, not proof of resource-dependent butterfly colonization or true dispersal. `precise1km_records` equals `records` for all 83,573 rows, an *explicitly expected condition in the original analysis scripts*, not a newly verified original occurrence-coordinate filter. The restricted source rows remain unshared.

**Source-supported public aggregate result:** [row-by-row private source audit receipt](../../../results/geb_v05_private_source_recovery_20261120/ORIGINAL_YEAR_ROWWISE_AUDIT.json) and [full evidence and remaining limits](GEB_ORIGINAL_SOURCE_ZIP_RECOVERY_AUDIT_20261120.md). The restored historical Parquet still has **no `gbifID` or `datasetKey` columns**. It is a different source selection from the published observed 224-species GBIF DOI and separate 225-candidate GBIF reindex. Publisher attribution for historical butterfly and 1,596-port original records cannot be recovered from this aggregate alone. Do **not** invent or transfer the other cohort's DOI. Anonymous, redistribution-licensed reviewer source data and coauthor authorization remain blocked.

## Exact release blockers after this audit

1. **Dataset provenance:** recover actual original GBIF contributing dataset identifiers/counts for each distinct published query cohort, with licensed original data or clearly documented restrictions. Do not mint or cite a nonexistent derived DOI.
2. **Port conversion:** identify the original WPI 2019 download, transformation tool and original-to-derived port-ID audit matching the frozen source file; establish redistribution rights for included derived layers.
3. **True reviewer archive:** include reproducible source and dependencies for each main result (including 155/224 consumer endpoints), original or permitted derived input data, and explicit legal omissions. Verify it after extracting its ZIP. A GitHub Actions artifact does not meet the journal's persistent archive rule, and the current ZIP is neither complete nor anonymous.
4. **Double anonymity:** separately scan document properties, source records, comments, ZIP paths, filenames, code snippets and active reviewer URLs. Do not upload the internal preview as a double-anonymous deposit.
5. **Human-required:** coauthor approval, approved one-corresponding-author title page and declarations, permissions, and chosen actual reviewer-access URL (e.g. a Dryad private-for-peer-review link if deposited).

## Decision and source-status firewall

**Current scope is SOURCE-LINEAGE AUDITED / INTERNALLY PREVIEWABLE / NOT JOURNAL DATA-RELEASE READY.**
The primary manuscript and separate S1–S8 are already coauthor-reviewable; the data/code archive and authorship/repository gates are still blocked. This report does not reinterpret the ecological science, erase the South American negative transfer, change the 9-Level1 cluster paired-slope p=0.068, or treat old v0.2 anonymous-review bundle checks as validation of v0.5.

Source and journal rules: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html ;
https://www.gbif.org/citation-guidelines ;
https://www.gbif.org/derived-dataset/about .
