# Subcountry resolution and juvenile-evidence audit — locked 48 positive species × country cells

**20 October 2026; exploratory diagnostic protocol fixed before retrieving
GADM level-1 or life-stage facets.** This is not an independent confirmatory
hypothesis, colonization test, population count, local host-use study,
or a replacement for the source-frozen 229-taxon country comparison.

## Motivation and source cohort

The primary GBIF country host-resource comparison retains exactly **48**
positive species × country cells, from **20 species / 19 genera**, after all
source, country, climate, native-host-guard and other-butterfly reporting
filters. They comprise **23 countries**. The conditional expected
added-resource overlap is 24.121933 among 26 observed, but only 21 informative
species × region strata remain vs a prospective minimum of 50.

Country-level co-coverage cannot determine whether the focal butterfly
and introduced host plant ever occurred in the same province, much less
the same larval habitat. The current test measures how coarsely 48 observed
butterfly presences are geographically described in the *GBIF index*.

Lock the exact 48 positive cells and species scientific names from
`provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json`;
derive taxonKey only from
`data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json`.
No new focal taxa, case-specific exception or live taxonomy re-resolution.

## Query contract

Pilot first on the two previously audited non-singleton candidates:
`Aguna asander` (DO, 9 indexed records in 2018–2025) and
`Papilio demoleus` (IL, 86 indexed records in 2018–2025).
Use fixed 2018–2025 years and compare GBIF reported counts to the source
receipt. If pilot GADM facet is unsupported or incomplete, stop; do not
claim successful localization.

Then acquire **all 48** positive species × country cells with the original
2010–2025 *pooled* window (both original 8-year windows combined).
For each cell use one compact /occurrence/search GET, exact locked taxonKey,
country, year=2010,2025, hasCoordinate=true,
hasGeospatialIssue=false, occurrenceStatus=PRESENT, and the five original
accepted basisOfRecord categories. Use limit=0 and facetLimit=300.

Request three independent facets: `gadmLevel1Gid`, `lifeStage`, and
`establishmentMeans`. The first is the compulsory geospatial subcountry
partition; other fields are evidence-quality diagnostics whose absence
is an explicit metadata coverage outcome. No individual occurrence JSON,
record keys, coordinates, observer IDs or point geometries are acquired.

For GADM, require unique valid GADM level-1 names, nonnegative integer
counts, fewer than 300 facets and sum ≤ reported total. If GBIF returns no
GADM facet or an unsupported field, mark NOT_EVALUABLE. Summation below
the returned count denotes **unmapped/unassigned** geodata, not absence.
Never re-normalize unknown GADM areas away. Administrative-unit counts
must not be interpreted as independent butterfly populations.

For lifeStage and establishmentMeans, preserve all raw category labels,
exact counts and unmapped remainder. **Do not** infer larvae, eggs, wild
reproduction or naturalization from a non-standardized category string
without inspecting record-level evidence. Report whether any records
are stage-annotated; report stage classification only as unverified
provider metadata. Invalid/truncated responses fail source completeness.

## Fixed descriptive outcomes

Across 48 positive species × country cells, report:
- Source-complete GADM facet fraction and fraction of *indexed records*
  geocoded into level-1 units; do not equate with spatial accuracy.
- Number of distinct mapped GADM1 units per positive cell, concentration
  of records in its largest unit (over **all** records, not just mapped),
  and how often a cell occupies exactly one unit.
- Cross-tab the fixed **added-host-country flag** vs spatial metadata
  concentration without implying local host overlap. These groups are
  source-conditioned and may differ in species, climate and reporting.
- Life-stage and establishmentMeans non-missing coverage; preserve exactly
  how many of the positive cells and records have interpretable labels.
- Report both pilot and full-source source totals; any API index drift from
  October 20 must remain explicit, never overwrite or recalibrate the
  original country-effect estimate.

This may help determine if finer-scale observational evidence exists.
It cannot add absent GADM1 candidate control cells, generate a fine-grain
host-resource overlap estimand, verify establishment, or repair the
21-informative-strata global country test.

## Validation, output and release boundary

Save safe aggregate query URLs, UTC timestamps, raw response SHA256, all
partial/failure receipts, exact source cohort SHA, per-country summary and
a reproducible Python script with unit tests. No author-only native-host
geometry or raw HOSTS/WCVP/CHELSA/TDWG inputs are redistributed.

All new content stays on the separate draft PR #51 analysis branch and
outside the frozen GEB v0.2 anonymous review snapshot and manuscript.
