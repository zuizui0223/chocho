# Source-provenance audit for 2018–2025 GBIF late-only country detections

**Protocol frozen before inspecting the new live GBIF dataset/year facets.**

## Focal units and rationale

Use ONLY the seven already identified country × butterfly cells that have
**zero indexed records in 2010–2017** and at least **two indexed records in
2018–2025**, after the previously fixed both-window other-taxon effort gate
(>=3 recorded other butterfly taxa in each window). Do not select a new species
on a favorable effect sign. The seven are derived mechanically from the archived
2026-10-20 [record-turnover summary](../results/gbif_country_20261020/gbif_country_record_turnover_summary_20261020.json),
not selected from newly queried GBIF results. Keep all 7, including the
four non-added-resource cells. The two >=5 cells are *Aguna asander* in the
Dominican Republic and *Papilio demoleus* in Israel.

The fixed baseline is **23,433 country facet rows**, the original **239-name**
receipt (231 resolved / 229 non-colliding), and full-index country occurrence
filters used in the 924-query national archive. Names and accepted taxon keys
must come from that exact saved receipt. No live taxonomy resolution.

## New data acquisition and QC before any interpretation

For each of the **seven** pairs, issue one compact GBIF occurrence search:
taxonKey=<frozen accepted key>, country=<frozen ISO code>, year=2018,2025,
hasCoordinate=true, hasGeospatialIssue=false, occurrenceStatus=PRESENT,
the five original accepted basisOfRecord values, limit=0, facetLimit=300.
Request three independent aggregate facets: datasetKey, year, basisOfRecord.
Preserve URL, UTC retrieval timestamp, raw response SHA-256, HTTP status and
facet category count. Do **not** download or publish individual occurrence
records/observer identifiers; preserve the complete compact aggregate JSON
response only when validated.

Require all three facets to exist, report exactly one facet per requested field,
contain unique facet values with nonnegative integer counts and no category cap.
Facet totals must equal response count; years must remain in 2018–2025. A
provider/schema error is a DATA_UNAVAILABLE outcome, not a zero. The current
GBIF index may have drifted from the Oct 20 frozen snapshot: record both new
count and frozen count, but **never silently substitute** the new count into
the October 20 country-effect calculation.

Record for each pair the number of years represented, distinct publishing
dataset IDs, largest dataset fraction, largest year fraction and basis split.
These provide *indexing provenance heterogeneity*, not verified independent
observations or multiple wild generations. Multiple years and datasets can
still contain correlated copies of the same source observation.

## Decision rule and ecological limits

This audit may classify records as single-dataset-dominated, temporally
concentrated, multi-year/multi-dataset, or unverified due to provider failure.
It cannot classify colonization, larval host use, establishment, survival,
population growth or temporal plant arrival. In particular, years in an
occurrence index are **not** introduction dates for host plants. If no cell
meets a stronger independence criterion, the original country model stays
unresolved; do not upgrade a source-audit category into a biological claim.

Report each of the seven pairs, all partials/errors, and the exact GitHub CI
run/retrieval timing. Do not change the frozen GEB manuscript or its anonymous
supplemental archive. This is an exploratory diagnostic, not a confirmatory
study or positive-effect discovery exercise.
