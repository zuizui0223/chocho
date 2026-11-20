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

## Original historical Parquet — full row-level audit and zero-ID repair

**2026-11-20: all six SNAPPY/RLE_DICTIONARY Parquet columns, all five data pages per column, and all 83,573 original rows were independently decoded from the private bytes.** This supersedes the earlier footer-only/manifest-only diagnosis. The original Parquet SHA256 matches the pinned 2026-10-30 source receipt.

| Source-original measurement | Independently computed value |
|---|---:|
| Unique species × landmass × year rows | 83,573 / 83,573, no duplicates |
| Butterfly species / landmass IDs | 224 / 3,384 |
| Year range | 2001–2024 |
| Sum of indexed source-record counts | **14,425,536** |
| Required-field missing cells | **0** in all six fields |
| Rows with `landmass_id="0"` | **2,636** |
| Indexed records associated with ID 0 | **4,386,686 (30.41% of source count)** |
| Distinct butterfly species for ID 0 | 126 |
| Constant landmass area for ID 0 | **50,654,050.6945 km²**, the largest area in this source |
| Records identical to `precise1km_records` | **83,573 / 83,573** |

**Critical finding:** `landmass_id=0` is a **populated, largest-area landmass identifier**, not a null/sentinel for this source. A hypothetical full-panel importer that rejects all zero IDs or silently filters with `landmass_id > 0` would lose 30.4% of source-indexed record counts. Preserve the original zero ID when validating the global source. The source-selected **51 island groups explicitly do not contain ID 0**; their original Parquet subset is **9,739 annual rows / 129 butterflies with records**. The new [private full-row audit CLI](../../../scripts/audit_chocho_historical_parquet_20261120.py) and self-test distinguish legitimate key zero from truly invalid negative or missing IDs. The preexisting source `recompute_near_island_relay_holdout_20261031.py` already casts landmass IDs to integers without a blanket `>0` filter; do **not** rewrite working analyses or reclassify the island cohort.

**Precision interpretation:** `precise1km_records` equals `records` on every original row, and source-pinned 2026-10-30 reconstruction scripts **explicitly require this equality**. It is not an independently validated precision-filter mechanism and cannot be used as a separate covariate describing georeferencing uncertainty. The original record-level geographic-precision filters would require original GBIF occurrence-level data to audit.

**Temporal replay from actual private Parquet:** Using 51 source-selected island IDs and the archive's original 1,275 island pairs, all **14 pairs within 200 km** produced exactly the earlier **41** butterfly×island-pair combinations where both islands were first indexed during 2011–2022. Of these, **21** had a re-indexing on the earlier island within 5 years and before the other island's first indexing. All fourteen pair-specific source counts matched the archived `shared_new_species` field. The five *a priori* focal species/island links also reproduced their complete year lists and original **5/5** temporal ordering. These source-selected, dependent sequences remain **observational recording patterns, not proof of relay, colonisation or shipping**.

**Published inspectable aggregate-only receipt:** `results/geb_v05_private_source_recovery_20261120/ORIGINAL_YEAR_ROWWISE_AUDIT.json` (contains counts/limited previously public five-link years, **not** the private 83,573-row Parquet). The original bytes remain private and were not pushed to GitHub. Local reproduction (requires `pandas` and `pyarrow`):
```bash
python scripts/audit_chocho_historical_parquet_20261120.py \
  --private-zip /path/to/chocho_original_analysis_requested_files_20261120.zip \
  --out /path/to/private-aggregate-audit.json
```
A public no-input `--self-test` confirms that ID 0 is allowed but negative IDs are rejected.

**GBIF citation remains unresolved:** This original 224-species 14,425,536-record **aggregate** has six fields (`speciesKey`, `landmass_id`, `year`, `records`, `precise1km_records`, `landmass_area_km2`), without `gbifID` or `datasetKey`. It cannot independently recover contributing dataset publishers or mint a source-accurate derived GBIF DOI. The **separate** original 224-butterfly observed-community study has a published 14,637,152-occurrence-source GBIF DOI and must not be substituted for this aggregate. The more recent 225×51 public reindex is a third, distinct source and snapshot.

## Public51 JSON: restored but not provenance-complete

The recovered JSON source receipt records 6 shard hashes and timestamps from **2026-11-09**, plus a comparison to the older October 31 sample showing query-index drift. Its explicitly declared units are **GBIF index events**, not colonizations or proven newly established butterflies. The JSON contains no original `gbifID`, `datasetKey`, or complete contributing-dataset records. Thus the original scientific caveats still apply. It does not support minting a retrospective GBIF derived DOI without recovered provider attribution or newly attributed, clearly versioned record-level retrieval.

## Remaining access/rights gates

1. The converted WPI2019 `ports_2019.geojson` **is now checksum-identified**, and the 51-island port linkage independently repeats. However, the original NGA input file and conversion recipe were not included, and permission to redistribute the derived port geometry is not established.
2. The historical Parquet is privately recovered and ID 0 is validated, but full original GBIF provider citations, the 155-/224-species independent models' complete source datasets, and legal sharing terms remain open.
3. The standalone public51 merged cell-level CSV and GBIF six-shard raw queries are **not** in this nine-entry ZIP; only their hash/summary receipt is present.
4. A separate **reviewer-anonymous, source-licensed stable archive** and coauthor-approved title-page/approval inputs remain necessary. The full original uploaded ZIP, source-machine absolute paths, and identifying internal receipt must **not** be deposited into a public/non-anonymous PR or GEB reviewer upload without rights and identity review.

**Status update:** previously **FILE MISSING** three-file input gate → **RECOVERED PRIVATE BYTES / CHECKSUM VERIFIED**. WPI-derived port geometry for the original 51-island test → **independently reproduced**. Original year-table 83,573 rows and island-pair chronology (14/41/21; five focal links) → **full-row independently replayed**. Original GBIF publisher attribution and complete reusable anonymous review dataset → **not yet cleared**. No change to the paper's conclusions, its nine-cluster p=0.068 caution, or the original 239-butterfly fixed known-host population.
