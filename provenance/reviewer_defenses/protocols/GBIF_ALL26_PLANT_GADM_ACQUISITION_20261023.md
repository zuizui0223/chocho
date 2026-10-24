# Frozen observational GBIF GADM1 query plan for all 26 source-added butterfly country cells

**23 October 2026**. This query protocol is fixed before retrieving
the **24 WCVP-accepted introduced candidate plant species** in
**44 unique country×plant pairs**, reconstructed in the original
source-only outcome-blind inventory
`provenance/reviewer_defenses/results/gbif_country_20261020/gbif_all26_source_host_country_inventory_20261023.json`. The 44 are grouped botanical observations, not 44
independent ecological experiments.

## Inputs and acquisition

- Use the exact source `unique_plant_country_query_basis` set
  containing 44 distinct accepted-ID×ISO2-country pairs; **no removal
  based on whether GBIF records exist**, butterfly juvenile outcome
  or taxonomic convenience.
- Look up each of 24 original single accepted WCVP source binomials
  against GBIF `/v1/species/match` using strict `Plantae`,
  `rank=SPECIES`, `strict=true`; only exact accepted canonical
  species may be assigned a GBIF taxonKey. No fuzzy/synonym matching
  or invented accepted taxon IDs. Save failed resolutions explicitly.
- For every taxonomically resolved pair, query exactly one compact
  country GADM1 facet, matching the original 2010–2025 year window,
  `hasCoordinate=true`, `hasGeospatialIssue=false`,
  `occurrenceStatus=PRESENT`, and the same five
  `basisOfRecord` categories. Require `limit=0`,
  `facet=gadmLevel1Gid`, `facetLimit=300`; do not retrieve
  individual occurrence records or coordinates. Preserve GBIF URLs,
  response SHA-256, retrieval timestamps and missing GADM count.
- Source-frozen **five-case 20 botanical country pairs** are
  included in 44; their previous 20 GBIF query totals, taxonomy
  IDs and GADM bins are a *historical drift diagnostic only*.
  Do not silently replace an older observation in the earlier
  five-case reference.
- The remaining 22 non-added butterfly positive country cells
  must remain unaltered and are NOT known host-absent controls.
- Rate limit at <=1.0 GBIF request/s, bounded retries for
  HTTP429/transient errors, and keep all provider-failure statuses
  separate from zero biological/GBIF count.

## Pre-fixed descriptive outcomes

1. Source-recognized botanical taxon identities resolved exactly
   out of 24; completed GADM1 facets out of 44 and positive GBIF
   plant-country source pairs out of 44.
2. Number and proportion of the **26** added-resource butterfly
   positive cells containing >=1 indexed candidate plant country
   record; and their candidate-host GADM1 coverage (any positive
   registered host inside butterfly-positive GADM units), without
   juvenile selection.
3. For all 26, compare observed butterfly GADM1 count and the
   count of those units where a listed host has GBIF records.
   Report the numbers of butterfly GADM1 records in observed
   host-positive administrative divisions and whether the complete
   five-case old results drifted.
4. If complete, calculate descriptive host links per country,
   candidate source richness and observed-taxonomic coverage.
   Report positive-query fractions by source candidate richness
   (1 vs >1), but DO NOT treat extra-host richness as an
   experimentally estimated biological benefit.
5. There is **no** valid 26-cell larval recruitment comparison
   simply because adult butterfly presences were initially
   required and there is no comparable unrecorded control effort
   background. No p-values or causal claims from this positive-only
   source sampling frame.

If GBIF indexing changes, report new snapshot separately. All
derivatives stay outside the frozen GEB manuscript and anonymous
review bundle.