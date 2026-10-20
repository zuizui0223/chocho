# Matched-dataset recording-effort denominators for the seven butterfly country cases

**2026-10-20; posthoc test.** The seven positive species × country cells
and 11 dataset × case groups are already selected. Freeze these *new*
denominator comparisons before obtaining the denominator data. This does
not revise the original 229-taxon host-plant macroecology protocol or manuscript.

## Scientific purpose and source limitations

The 2025 Israeli `Papilio demoleus` records increase across iNaturalist
and INPA. The within-dataset year counts are known, but counts cannot be
interpreted without each platform's overall reporting. We ask only:
**did focal-species indexed records rise relative to the platform's total
recording volume?** This is not an abundance or causal study. All counts
are GBIF indexed records, not independent insects or survey units.

Use the same 11 dataset × country cells from
`gbif_country_dataset_year_matrix_20261020.json`, all their frozen
country/species identifiers, 2018–2025 year window, all five accepted basis
classes, geospatial coordinate quality and PRESENT status. Retain 2010–2017
as fixed initial-zero reference only; do not re-query earlier occurrences.

## Two non-interchangeable reporting denominators

For each dataset × country group, query `limit=0&facet=year&facetLimit=300`,
`year=2018,2025`, `hasCoordinate=true`, `hasGeospatialIssue=false`,
`occurrenceStatus=PRESENT`, all original `basisOfRecord` values, and
the **same datasetKey/country**.

1. **Broad recording activity:** no taxonomic filter, all taxa within that
   dataset/country. This includes the focal species and is deliberately a
   general source-activity denominator.
2. **Lepidoptera recording activity:** identical filters plus
   `orderKey=797`, the GBIF 2025 legacy backbone order key. Because GBIF's
   taxonomy has been migrated in 2026, report a **compatibility check**:
   every result must be >= the focal-species count in the same dataset/country/
   year to be interpreted as a valid nested denominator. Otherwise mark
   TAXONOMY_KEY_NOT_COMPARABLE and do not estimate Lepidoptera-normalized rates.
   Never silently swap keys or combine taxonomic checklists. The all-taxa
   denominator is the fallback, still subject to its own exact nesting check.

Do not compare unfiltered all-GBIF taxa with a frozen butterfly species if
even the source year/country/dataset scope differs. Save query URLs, UTC,
raw response hashes and completeness statuses (no individual occurrences).
Any zero or missing denominator with a positive species numerator is
INVALID_SOURCE_NOT_ZERO. On provider drift or taxonomy mismatch, do not
replace the original source metrics.

## Descriptive readouts fixed up front

For Israel `Papilio demoleus`, compare **2024→2025**, separately for
iNaturalist and INPA, after dividing the source-locked 2024 and 2025 focal
counts by the exact same source's 2024 and 2025 country recording denominator.
Report each numerator, denominator, rate per 10,000 GBIF records, and
ratio-of-ratios from 2024 to 2025 for all-taxa and validated Lepidoptera
denominators. Do not interpret such rates as an abundance index unless
observation design is standardized, which has **not** been established.

For Dominican `Aguna asander`, report the six sampled 2019–2025 years
and the annual fraction against iNaturalist's country-recording denominator,
not a test of repeated wild generations. List all seven cell × provider
combinations and their year consistency, keeping cases with smaller counts.

No significance-test optimization, no cherry-picked denominator, and no
claims about larval success, naturalization or recent plant arrival.
Preserve the frozen GEB v0.2 anonymous review snapshot boundary.
