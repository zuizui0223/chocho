# GEB one-paper v0.3 — dataset citation, provenance and journal preflight

**Audit date:** 2026-11-10. **Article:** [single integrated GEB manuscript v0.3](ISLAND_PORT_BUTTERFLY_RESEARCH_ARTICLE_V03.md). **Status:** revised coauthor-review draft, NOT formally submitted, NOT final release/rights approval.

## 1. Material verified from frozen original source plus publisher records

| Component | Source used by this study | Version / licence / citation status | Submission treatment |
|---|---|---|---|
| Butterfly–host links | HOSTS via the `globalbioticinteractions/HOSTS` mirror, commit `808e0b869f9ec1adf8efff87cf6a395adda103e0` (2023-08-29) | NHM archived HOSTS dataset DOI **10.5519/havt50xw**, **CC0**. Source-specific mirror `globi.json` also identifies NHM query DOI `10.5519/qd.bsucrxdz`; both are distinct. | Original dataset DOI in References, commit in reproducibility docs; source-specific query DOI may be listed in SI if its exact input is documented |
| **Native/introduced host-plant geography** | `matildabrown/rWCVPdata` commit **`65bed76bae9d644ccb6ad200c05f9f5071d89e05`**; `docs/REPRODUCIBILITY.md` and `docs/exploratory/butterfly_resource_envelope_protocol_v0.1.json` agree | **WCVP version 13, 21 May 2024; Kew DOI 10.34885/0yex-xv26; dataset CC BY 3.0.** The wrapper R package `rWCVPdata 0.5.0` has package license GPL ≥ 3 (distinct from Kew dataset license). | **CLOSED citation-version gate**: v13 dataset citation added to manuscript References and Methods; acknowledge Kew and preserve version |
| WCVP **citation error in the upstream package** | The exact pinned `rWCVPdata/inst/CITATION` claims v13 but embeds DOI `10.34885/nswv-8994` | Kew's permanent landing record for **10.34885/nswv-8994** says **WCVP v10 (2022)**, not v13. The correct Kew v13 landing record is **10.34885/0yex-xv26**. | **Do not copy the misleading package DOI** into the article; cite Kew's authoritative version-13 record |
| Islands / shorelines | Original GSHHG **2.3.7** and source SHA-256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41` | 1996 methods and 2.3.7 Zenodo dataset DOI `10.5281/zenodo.7007502` already cited; source terms in prior v0.2 audit | SOURCE-IDENTIFIED, keep exact version and limitations; no new route claim |
| Ports | NGA World Port Index **2019**, Publication 150, **27th edition** | NGA original 2019 publication states no US-government copyright claim under Title 17. Original port GeoJSON had SHA-256 `d580e8520dc5bd77c0ab464df3bdc9ea95653c47b00c42a510dc7c1519feaa78`; actual downloaded/conversion bytes of the 2019 source must be reconciled with the original paper's derived port file. | Publication edition confirmed, **source-file provenance/crosswalk OPEN** |
| Butterfly traits | LepTraits 1.0 `RiesLabGU/LepTraits` consensus commit identified in source protocol; 2022 Scientific Data article | Scientific Data article DOI `10.1038/s41597-022-01473-5`; original data-file distribution license requires separate source review | Source identity frozen, exact raw table rights **OPEN** |

### Primary official records / source identity

- Kew official **v13** archive: https://kew.iro.bl.uk/entities/product/5d36f109-7e49-461e-9541-6a85aac339fa , DOI https://doi.org/10.34885/0yex-xv26 , CC BY 3.0
- Kew official **v10** archive (the incorrect v13 package DOI): https://kew.iro.bl.uk/entities/product/a867652b-666f-4587-8ce2-7fb2c233708a , DOI https://doi.org/10.34885/nswv-8994
- Frozen R source: https://github.com/matildabrown/rWCVPdata/blob/65bed76bae9d644ccb6ad200c05f9f5071d89e05/inst/CITATION
- Kew terms: https://www.kew.org/science/collections-and-resources/data-and-digital/terms-of-use
- NHM source: https://data.nhm.ac.uk/dataset/hosts
- US NGA 2019 WPI edition: https://msi.nga.mil/Publications/WPI

## 2. GBIF citation gate remains genuinely OPEN: no simulated DOI

The new **225 butterflies × 51 landmasses = 11,475 butterfly×landmass index cells** source is an aggregate of GBIF species-key period-index queries. Its six 2026-11-09 shards and source-hash/GBIF-query-date receipts are in `data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json`. **A butterfly×landmass count is not a GBIF occurrence ID or a `datasetKey` source list.** The original 2026 snapshot cannot be assumed to be exactly reproducible from later live searches; even a re-query needs a new snapshot and drift ledger.

The global port source had **1,596 completed queries** out of 3,630 WPI candidates; **25,742 host-candidate plant×port combinations** are derived presence inventories, not original GBIF record IDs or source-publisher attribution. The additional 17-island, 24-port/nonport pilot is also a distinct GBIF acquisition.

The existing **10.15468/dl.pp5nc9** applies **only** to the earlier 53,434-occurrence/31-butterfly validation (313 contributing GBIF datasets) documented in `docs/GBIF_OCCURRENCE_DOI.md`. **It cannot be passed off as the DOI for the 225×51 index, the global port plant inventory or the later port/nonport matched queries.**

GBIF's official instructions for search-API outputs require attribution to underlying dataset publishers and recommend registering a **derived dataset** using the contributing `datasetKey`/DOI list (ideally source record counts), an accessible permanent URL and a clear description of how the aggregates were produced. See https://www.gbif.org/citation-guidelines and https://www.gbif.org/derived-dataset/about .

**Executable route, retaining valid scientific claims:** (1) inventory the exact original queries and available original record IDs/`datasetKey` from all three acquisition families without assuming aggregates contain occurrence-level provenance; (2) if unavailable, reconstruct equivalent occurrence-level provider-attributed records into a *new* documented snapshot and quantify drift, without replacing the frozen analyses or retroactively declaring unchanged inputs; (3) archive a persistent derived dataset with publisher counts and query rules; (4) use the actual authenticated GBIF derived-dataset DOI and GBIF-provided citation in manuscript/SI. **Account authentication and DOI creation were NOT performed.** Do not invent any DOI.

For the restricted older **83,573-row private historical Parquet**, and any unreleased original host-edge/point subsets, preserve the restricted-file statement rather than treating public index aggregates as complete replay of those analyses.

## 3. Official GEB research-article formatting audit

GEB official author guidelines: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html (checked 2026-11-10).

| Requirement | Current state | Action |
|---|---|---|
| Global patterns/relevance | Island resource amplification, port flora and island-community record geography are framed as global macro-biogeography | Retain explicit generality, do not claim shipborne colonization |
| Research Article **structured abstract ≤300 words** | Canonical v0.3 abstract 291 words by source-audit tokenizer at the 2026-11-09 freeze; verify after any edit | PASS, preserve seven required sections |
| Main body **roughly 5,000 words** typical | Existing source ~5.7k words depending on tokenization; five figures, below journal's typical six–eight display pieces, may offset some length | **EDITORIAL COMPRESSION STILL ADVISED**; important limitation is not optional |
| Short running title **<40 characters** | Added `Island butterfly host resources` | PASS |
| **6–10 alphabetized keywords** | Eight keywords added to canonical source | PASS |
| Continuous **line numbering** | Not yet certified in the clean coauthor Word/PDF | Open pre-submission formatting |
| Double-anonymous main text and **separate title page** | Earlier anonymous export was for an older resource-only manuscript; the local coauthor version is not an accepted anonymous GEB upload | Submission-grade new anonymous clean manuscript and separate real author/affiliation/funding/COI page still required |
| Main text file, separate Supporting Information, figure captions/figures in review document | Source-verified five figures and a clean coauthor Word/PDF exist, but earlier review version is not a journal submission package | Separate SI provenance and raw source-access limits; check panel labels (a–d), caption/page joins and file metadata |
| Mandatory **Data and Code Availability** | Explicit draft statement exists; no new GBIF source DOI or immutable archive DOI yet | **OPEN — blocks claims of a final submission** |
| Cover letter | GEB requires a separately uploaded PDF and a <250-word significance paragraph | Prepare a truthful editor-focused draft, but do not send without the corresponding author |

**Journal submission is not authorized or implied by this audit.** Do not anonymize a document merely by deleting names from visible text; inspect Word properties, author metadata and hidden review information. The actual title-page authorship, affiliations, corresponding author and funding must come from the collaborators themselves.

## 4. Updated submission gate decision

- **Completed now:** Correct WCVP v13 source and citation found and independently verified, 2024 pinned source identified, base dataset rights distinguished from R package software license; main manuscript updated, eight GEB keywords and a running title added; official GBIF, GEB, Kew, NHM source requirements reconfirmed.
- **Remain BLOCKING for final submission:** New GBIF occurrence provenance and DOI(s); raw original/source citation and distribution rights, especially restricted historical records and exact WPI source provenance; distinct anonymous manuscript/title page and fully checked main/supporting files; corresponding author and coauthor approval.
- **Optional editorial improvement:** Reduce approximately 0.5–0.8k main-text words without removing the historical null, first-index limits, and AIS/wind falsification controls. Main figures remain source-verified, **not** publisher-approved artwork.

**Next focus: GBIF provenance, not a new ecological hypothesis.**
