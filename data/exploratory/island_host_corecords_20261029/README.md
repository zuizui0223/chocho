# Original26 source records joined to exact physical islands

## Scope and result

This is a **fresh source acquisition**, authorized and frozen before querying in `provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md`. It preserves the original26 butterfly–country cases,44 plant–country query definitions,46 candidate-host edges,16 butterfly taxa and24 plant taxa. The observation interval is2010–2025. Original numeric taxon keys are explicitly queried against GBIF Backbone checklist `d7dddbf4-2cf0-4f39-9b2a-bb099caae36c`; they are not replaced by the current COL XR island-panel keys.

**Primary frozen96 panel:** the acquired original26 butterfly–country queries contribute **0 butterfly records** inside these96 polygons. There are **13 plant records**, all *Robinia pseudoacacia* on the Isle of Wight (`gshhg237_3c58026612a9987c`):2 in2010–2017 and11 in2018–2025. Twelve records pass conservative full reported-uncertainty containment, also with reported uncertainty≤1km. Consequently the original-source co-recording count is0. This is a source-cohort geographic non-overlap result, **not absence of these16 butterfly taxa on the96 islands under other country queries**, and not absence of host plants. The96 roster is unchanged.

**Separate original8260-island inventory:**341 butterfly and3,324 unique plant–country record indices lie inside source islands≥5km². Same-island candidate-host co-recording within±3calendar years occurs for287 butterfly indices, in12 original butterfly–country cases comprising6 butterfly taxa, on3 physical islands. These are source-index counts, not independent populations or verified biological interactions.

| Physical source island | Exact source key | Butterfly indices with±3year host co-recording | Both full reported uncertainty envelopes contained | Both also≤1km |
| --- | --- | ---: | ---: | ---: |
| Cuba | `gshhg237_14f1fba44d402d46` |71|24|20|
| Hispaniola | `gshhg237_b8e2b5c8b125a3e0` |215|160|144|
| Great Britain | `gshhg237_0c2e36f13d43c4d7` |1|0|0|
| Total | |287|184|164|

Human-readable island names were checked against the exact source polygons using rounded public cartographic interior anchors at Santa Clara, Santiago de los Caballeros, Oxford and Newport, respectively. Biological observation coordinates were not used for naming. The hashed source keys remain the authoritative identities.

The12 cases are *Aguna asander*–CU/DO/HT; *Anteos clorinde*–CU/DO; *Epargyreus clarus*–GB; *Historis acheronta*–CU; *Papilio demoleus*–CU/DO/HT; and *Phoebis agarithe*–DO/HT. *P. demoleus* contributes219/287 co-recorded butterfly indices (76.3%). This strong concentration and the original positive-country ascertainment preclude treating this inventory as a representative global host-effect estimate.

## Time and uncertainty

Butterfly year assigns the output window. A candidate plant can be within±3years across the2017/2018 boundary, so cross-window-only matches are explicitly separate.

| Butterfly window | Butterfly indices on any source island | Plant indices on any source island |±3year co-recorded butterfly indices|Both uncertainty envelopes contained|Both also≤1km|Cross-window-only matches|
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
|2010–2017|61|1,264|25|12|9|12|
|2018–2025|280|2,060|262|172|155|0|
|Total|341|3,324|287|184|164|12|

There are275 same-window co-recorded indices and12 cross-window-only indices;160 have a same-year source-host record. Of341 island-assigned butterfly indices,252 pass the full reported-uncertainty containment screen. Unknown uncertainty stays unknown. Across all countries,94/467 butterfly and342/3,597 plant records lack usable uncertainty.

Any-point membership uses the exact unsimplified GSHHG2.3.7 high-resolution L1 polygon. Strict containment requires a conservative rectangle enclosing the reported WGS84 uncertainty disk to lie strictly inside that polygon; boundary points, missing uncertainty, and uncertified dateline/polar envelopes do not pass. This assumes reported uncertainty covers actual error and does not validate GSHHG shoreline accuracy, shoreline age or omitted lakes. Same-island co-recording is not proximity within5km, feeding, wild plant availability, oviposition, establishment, arrival or causality.

## Acquisition and provenance

First source call: **2026-10-29T13:48:34.628675Z**. Last source call: **2026-10-29T13:51:56.987118Z**. All70 source queries completed:467 butterfly plus3,597 plant rows,132 HTTP attempts,0 retries; minimum observed call-start spacing1.500092seconds. Every acquired row passed the source coordinate/year/taxon gate. The per-query5,000, total12,000 and200-HTTP-attempt caps were respected.

Two count changes distinguish this new retrieval from the older facet snapshot: *Robinia pseudoacacia*–GB3,203→3,204 and *Asclepias syriaca*–DE86→87. All other source counts were unchanged. Missing/capped/failed sources would be incomplete, not zero. A genuine complete zero means no indexed records under the fixed query, not ecological absence.

Raw exact biological coordinates, occurrence IDs and response bodies existed only in transient process memory. No such data were saved to repository, cache, logs or artifacts. Response hashes, request identities, timing, dataset×year counts and derived island aggregates are retained. This new local pipeline leaves the prior Actions-only pipeline unchanged.

## Files and verification

| File | Role |
| --- | --- |
| `source_island_query_lock.json` | Immutable pre-query plan, implementation/test, original-source, query-contract and geography identities. |
| `frozen_original_island_universe.csv` | The independent8260-island universe: derived keys, area, geometry hashes and frozen96 membership. |
| `source_query_receipts.json` |70 exhaustive-query receipts, source-body hashes, count changes and dataset×year diagnostics. |
| `source_island_year_counts.csv` |174 sparse source×island×year aggregate rows; no individual record identities. |
| `frozen96_case_island_window_corecords.csv` |All4,992 primary26×96×2 cells, including complete zeros. |
| `original8260_positive_source_case_inventory.csv` |52 case×island×window rows on islands with a queried butterfly or candidate-plant record. |
| `original46_edge_island_window_inventory.csv` |94 source-edge×island×window rows retaining exact candidate-host identity. |
| `host_corecord_manifest.json` |Acquisition/scientific summaries and hashes of acquisition outputs. |
| `verification_receipt.json` |Subsequent independent weighted-year recomputation of every co-recording cell, hash checks and field-level privacy checks. |

Sixteen focused new/original point tests passed. The independent weighted-year check reproduced all4,992 primary,52 diagnostic and94 edge cells, including timing and uncertainty sensitivities. This README and the subsequent verification receipt were added after acquisition and are not retroactively inserted into the immutable acquisition manifest.

## Replay

The script refuses to overwrite an existing live-snapshot manifest. A new run is another current-index retrieval, not reconstruction of unsaved historical coordinates. Use a fresh output directory and keep the fixed source cache/plan:

```bash
/tmp/chocho-geo-venv/bin/python scripts/acquire_island_host_corecords_20261029.py \
  --source-cache /tmp/chocho-island-sources \
  --plan provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md \
  --output /tmp/chocho-host-corecords-new-snapshot
```

Append `--prepare-only` to verify/freeze source and geometry identities without any GBIF source queries. Coordinate clients must run sequentially with other GBIF clients to preserve the declared source rate limit.
