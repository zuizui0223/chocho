# Dataset × observation-year audit of seven butterfly country detections

**20 October 2026; posthoc cohort, within-cohort query contract fixed before
the dataset–year cross-tabs were inspected.** This is a source-provenance
diagnostic, **not** a prospective ecological hypothesis or proof of colonization.

## Frozen source and exact scope

Use the already completed seven non-singleton species × country cells selected
with **zero country records in 2010–2017**, **>=2** in 2018–2025, and >=3
*other* focal butterfly taxa recorded in each period. Read these seven
unchanged identities from
`gbif_country_provenance_facet_qc_20261020.json`, which includes complete
datasetKey/year/basisOfRecord facet totals matching the original country data.
Use the six GBIF dataset UUIDs from that file. No live taxonomy or new cohort
selection. In total there must be **11 positive dataset × species × country
combinations**. The pre-registered original 229-taxon country model is
unmodified.

## Query and validation

For each of the eleven combinations, request a **compact indexed-count facet**:

- GET `/v1/occurrence/search` with the frozen taxonKey, country ISO2, GBIF
  datasetKey, `year=2018,2025`, `hasCoordinate=true`,
  `hasGeospatialIssue=false`, `occurrenceStatus=PRESENT`, the
  **same five** accepted basisOfRecord values as the source, `limit=0`,
  `facet=year`, `facetLimit=300`.
- Preserve exact source query URL, response SHA-256, retrieval UTC, HTTP status,
  and attempts, but do **not** acquire or publish individual occurrences,
  record IDs, coordinates or observer identifiers.
- Require count to equal the existing datasetKey marginal total and yearly
  facet values to be complete, unique, nonnegative integers in 2018–2025.
  Sum the dataset-conditioned year counts back to the frozen **country-level
  year marginals** and country total. On mismatch, record INDEX_DRIFT or
  FACET_INCOMPLETE, never substitute live counts into source-frozen effects.
- Preserve the original dataset-registry titles, GBIF dataset keys, DOIs, and
  publication dates, but distinguish a dataset **created in 2025** from
  observations **made in 2025**.

## Before looking at the conditional year table, fixed readouts

For `Papilio demoleus` in Israel (86 records, 43 iNaturalist, 42 INPA, 1
eButterfly; 61 total records dated 2025), report each publisher's fraction of
the 2025 total, annual counts, and 2024→2025 change. Recalculate the 2025
share excluding INPA **only as a provenance sensitivity**, not as a
population trend. For `Aguna asander` in the Dominican Republic, report
the 9-record iNaturalist series and compare its 6 separate years with the
Israel concentration. Report the other five cells, without cherry-picking
favorable trends.

Do not claim that source change explains year variation until the year ×
dataset marginals pass all checks. Even a large source-related jump does not
distinguish observer effort from biological abundance without standardized
survey design. Do not equate database independence with independent wild
butterflies or local larval recruitment.

## Output and interpretation

Save a compact 7-cell × dataset × year matrix, verification hashes, complete
and incomplete receipts, a human-readable discussion of source concentration,
standalone no-network tests, and an Actions run. All historical analysis
results and frozen GEB manuscript remain unchanged. Put these files outside
the anonymous frozen GEB v0.2 snapshot.