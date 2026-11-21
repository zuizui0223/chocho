# GEB v0.5 Fig. 2D: limited anonymous-identifier-screened replay ZIP

**2026-11-21 | source-only, coauthor-review QA, not a full anonymous GEB review-data deposit.**

## Actual completed verification

- Immutable scientific draft: PR #179 at `99b876421af2b4f51e58569862994622a02717d7`. This PR **was not modified or merged**.
- The independent source audit PR #189 supplies `scripts/build_GEB_v05_blind_fixed239_subanalysis_20261121.py` and `.github/workflows/verify-geb-v05-limited-blind-fixed239-review-source.yml`.
- [Passed source package build/extract/replay CI #44769070926](https://github.com/zuizui0223/chocho/actions/runs/44769070926). The reviewer-screened **limited** source artifact is [Actions artifact #14645249184](https://github.com/zuizui0223/chocho/actions/runs/44769070926/artifacts/14645249184) (90-day temporary Actions retention, not a stable DOI).
- Exactly **8 files** in `GEB_V05_SCREENED_LIMITED_FIXED239_NOT_FULL_REVIEW_ARCHIVE.zip`; exact SHA256 **`c4a9f05f593f59924bdd188c88b7eb036b5d5a89804d665369ab0f4621ba9783`**:
  - `README.md`: correctly scoped source citations (HOSTS, WCVP v13, WPI2019, TDWG) and explicit limitations/rights gate
  - `REPRODUCE.sh`: installs dependencies, downloads source-hashed TDWG geography and calls both published scripts
  - `VERIFY_REPLAY.py`: compares eight source-pinned original and MEM10 coefficients at tolerance **1e−4**
  - `EXPECTED_FIXED239_SUBMODEL.json`: numerical expectations extracted from the pinned original Fig. 2D receipt **without an author-identifying source URL**
  - `data/region_metrics.csv`: 369 regional fixed-239 HOSTS/WCVP candidate host-coverage entries, source Git blob `4aba89565354edc6c44d566523864cec66a30cbe`
  - `data/whole_flora_wpi_126_islands.csv`: 126 botanical islands with port counts and whole native/introduced floral richness, source Git blob `53b618fa51a1a418ccc6e4363f4c866c389c90cc`
  - `scripts/audit_island_port_isolation_239_20261120.py` / `scripts/audit_island_port_spatial_239_20261120.py`: the two source-locked island distance and Moran/MEM10 Python programs.
- Build-time checks: both code blobs and both source table blobs match immutable originals, all eight ZIP file names/contents and fixed ZIP timestamp/mode are directly scanned for identifying account/name/email/user paths, no record-level `gbifID` or `datasetKey` source columns, no raw GBIF sighting files or privately recovered WPI/annual butterfly Parquet.
- **After extracting the ZIP**, the GitHub runner installed dependencies, fetched the fixed TDWG Level3 GeoJSON and confirmed SHA256 `c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840`, recomputed the original 126-island models, then checked Fig. 2D coefficients and top ten Moran-eigenvector residual sensitivity with absolute error below 1e−4.
- Console check: `GEB_LIMITED_FIXED239_INDEPENDENT_REVIEW_SOURCE_REPLAY_PASS` and `ZIP_EXTRACTED_ANON_IDENTIFIER_SCREENED_FIXED239_MEM10_SOURCE_REPLAY_PASS`. Existing full manuscript, blinded DOCX, internal evidence preview and legacy anonymous v0.2 CI checks also passed for the audit source branch after a scanner-source false positive was fixed by splitting a regex literal without weakening the scanner.

## Exact limited replay numerical scope

The internal source-only expected JSON identifies: **126** insular botanical regions, **239** fixed candidate butterflies, **123** currently host-covered islands, native known-host distance **β−0.16756537**, introduced-only resource-share distance **β+0.42509951**, paired current-minus-native distance **β+0.18675537**, whole introduced flora port **β+0.35045012**; Moran/MEM10 native distance **β−0.13302902**, introduced-only distance **β+0.39489401**, paired distance **β+0.20390851** and introduced-whole-flora port **β+0.38770982**. All are source-indexed *potential larval-host geography*, not actual host consumption or observed butterfly arrival. Clustered uncertainty in the *paired slope change* remained open in the paper (nine-Level1 p=0.068).

## Remaining release and interpretation blockers

**This package does not complete a reviewer-deposit gate.** It reproduces only one source-processed Figure 2D subanalysis, **not** full WCVP/HOSTS original records, the whole manuscript's four figures and additional SI, historical 83,573-row butterfly history and full 225×51 index, 1,596-port historical original GBIF provider identities, two independent 155/224 recorded butterfly outcomes, individual larval images, or historical cargo/plant introductions.

Filename/content scans reduce direct identity leakage but do not prove perfect indirect anonymity; source code and data may already have been present in publicly searchable Git history. An actual editor-approved review route needs a separate anonymous account/workspace/deposit and metadata checking **plus author rights approval** of aggregate table redistribution, coauthor-approved title page/disclosures and a durable archive plan. This temporary GitHub Actions artifact is not a stable citation. No GBIF derived DOI or journal upload was performed.

**Practical next release decision:** retain this ZIP for confidential coauthor QA, then package a complete journal-specific independent source and provenance inventory with source-specific access (original public DOI, aggregate derived table or openly labelled restricted input) and a real reviewer-access URL; do not claim raw data are available when they are not.
