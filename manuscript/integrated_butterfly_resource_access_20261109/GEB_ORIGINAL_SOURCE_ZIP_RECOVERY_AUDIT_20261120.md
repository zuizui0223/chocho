# Original chocho analysis ZIP: recovered-source audit (internal upload, no raw-data publication)

**2026-11-20 | verified from user-supplied ZIP bytes in the analysis environment.**
Applies to [Draft PR #189](https://github.com/zuizui0223/chocho/pull/189) only. The manuscript [Draft PR #179](https://github.com/zuizui0223/chocho/pull/179) is unchanged. The original ZIP is **NOT** included in this Git repository or public archive: it contains source-machine metadata and potentially redistribution-restricted WPI/GBIF-derived datasets.

## Recovery decision

The three paths formerly absent from both inspected Git commits **are now present in a separately provided private recovery ZIP**. They remain *absent from the published Git trees*. Recovering file bytes does not establish that all original occurrence-level inputs, dataset publisher attribution or rights to redistribute those bytes are available.

| Recovered path | Actual bytes | Exact SHA256 from inspected bytes | Verified scope |
|---|---:|---|---|
| `work/barriers/ports_2019.geojson` | 859,406 | `d580e8520dc5bd77c0ab464df3bdc9ea95653c47b00c42a510dc7c1519feaa78` | 3,630 unique WPI-source-style port OBJECTID point features; **matches previously source-pinned original converted GeoJSON SHA256** |
| `work/barriers/global_original_species_island_year.parquet` | 411,552 | `8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f` | Valid `PAR1` head/tail, Thrift footer specifies **83,573 rows**, one row group, six named columns; same SHA exists in prior source-audit register |
| `data/exploratory/public51_butterfly_index_20261109/full51_source_validation.json` | 6,202 | `58b697b86925fb6cf753e2f8e68fe646076abd28038283958accdfb63a926741` | Source receipt for a separate **225 candidate butterflies × 51 landmasses = 11,475 cells**, six November 9, 2026 asynchronous GBIF query shards; **not the original historical Parquet** |

ZIP CRC validation and all seven data-file byte sizes and SHA256 values matched its supplied manifest. The ZIP additionally contains the five source tables/receipts needed for original port-to-island identity and geodesic reconstruction, the manifest, and `validation/port_coordinates_receipt.json` (nine archive entries total).

## Original 51-island geometry independently re-executed

Using the **actual recovered WPI port coordinates**, the source `SINGLE_LANDMASS` 5-km port-to-island mapping, and the `query_complete` status, the archived spherical port-coordinate averaging method was independently implemented and tested against every original island-pair distance:

- **3,630** source port OBJECTIDs, all unique.
- **1,596** botanically queried/complete ports across the wider inventory; **194** eligible complete original port coordinates within **51** eligible landmasses.
- **1,275 = 51 × 50 / 2** unique island pairs were evaluated.
- Largest absolute difference from the original spherical pair distance field: **6.61 × 10⁻¹¹ km**, vastly below the source reconstruction tolerance of 0.01 km. Archived source script/receipt reports **2.996 × 10⁻¹¹ km** using its own floating-point evaluation; differences are roundoff only.
- Corrected ecological scope: a representative **spherical mean of selected mapped WPI2019 ports**, not an island polygon centroid, historical ship trajectory or direct plant arrival.

This establishes source-geometry **identity and reproducibility**, not actual transport.

## Historical Parquet: important new evidence and unresolved defect

The recovered file is a real Parquet payload with six footer-defined columns: `speciesKey`, `landmass_id`, `year`, `records`, `precise1km_records`, `landmass_area_km2`. The footer independently confirms **83,573 rows**, not raw occurrence events. Neither `gbifID` nor `datasetKey` is present. Consequently **this file cannot alone reconstruct the contributing GBIF publisher/occurrence DOI citation**, even though the aggregate counts are now physically recovered.

The supplied source manifest additionally reports **224 species**, **3,384 landmass groups**, **2001–2024**, **14,425,536 indexed source records in aggregate**, and **2,636 rows with `landmass_id=0`**. Those five distributional/count summaries were *not independently recomputed from column values* in the present environment; row count and column names were independently read from footer metadata. The manifest records that a legacy importer **rejects landmass_id=0**. Treat that as a confirmed **replay-interface blocker to investigate**, not license to silently drop 2,636 rows or rewrite the original dataset.

A distinct 224-butterfly observed-community GBIF download DOI `10.15468/dl.umqmsf` and 14,637,152 original occurrence rows are described elsewhere in the paper. Do not equate that original cohort to this separately recovered Parquet's purported 14,425,536 records, despite the same broad butterfly count.

## Public51 JSON: restored but not provenance-complete

The recovered JSON source receipt records 6 shard hashes and timestamps from **2026-11-09**, plus a comparison to the older October 31 sample showing query-index drift. Its explicitly declared units are **GBIF index events**, not colonizations or proven newly established butterflies. The JSON contains no original `gbifID`, `datasetKey`, or complete contributing-dataset records. Thus the original scientific caveats still apply. It does not support minting a retrospective GBIF derived DOI without recovered provider attribution or newly attributed, clearly versioned record-level retrieval.

## Remaining access/rights gates

1. The converted WPI2019 `ports_2019.geojson` **is now checksum-identified**, and the 51-island port linkage independently repeats. However, the original NGA input file and conversion recipe were not included, and permission to redistribute the derived port geometry is not established.
2. The historical Parquet is privately recovered, but an importer/landmass-ID handling issue, full original GBIF provider citations, the 155-/224-species independent models' complete source datasets, and legal sharing terms remain open.
3. The standalone public51 merged cell-level CSV and GBIF six-shard raw queries are **not** in this nine-entry ZIP; only their hash/summary receipt is present.
4. A separate **reviewer-anonymous, source-licensed stable archive** and coauthor-approved title-page/approval inputs remain necessary. The full original uploaded ZIP, source-machine absolute paths, and identifying internal receipt must **not** be deposited into a public/non-anonymous PR or GEB reviewer upload without rights and identity review.

**Status update:** previously **FILE MISSING** three-file input gate → **RECOVERED PRIVATE BYTES / CHECKSUM VERIFIED**. WPI-derived port geometry for the original 51-island test → **independently reproduced**. Original GBIF publisher attribution and complete reusable anonymous review dataset → **not yet cleared**. No change to the paper's conclusions, its nine-cluster p=0.068 caution, or the original 239-butterfly fixed known-host population.
