# Fixed within-country localization of GBIF juvenile life-stage metadata — 21 October 2026

This exploratory follow-up separates **where a GBIF record is indexed** from
**where an identified butterfly host is observed**. It is a conditional
metadata audit of the original 48 source-positive species×country cells
only, not a colonization, breeding, or global host effect test.

## Source-lock before outcomes

Read the existing source-exact 48-country-cell file at
`provenance/reviewer_defenses/results/gbif_country_20261020/gbif_positive48_gadm_lifestage_facets_20261020.json`, with 756 original GBIF indexed occurrence records,
736 GADM1 mapped records, 500 life-stage labels, and 39 putatively
immature-source labels within 11 cells and three butterfly species.
The only included raw lifeStage labels are exactly
`Egg`, `Larva`, `Caterpillar`, `Pupa`, as present in the
**previously archived full 48-cell source**. All other cells remain
zero **labeled juvenile records** rather than biological non-occupancy.
No additional outcome-based selection. Do not change host exposure flags,
original plant-host sources, 2010–2025 time window, taxa, countries or
record basis filters. No 2026 partial year records.

## Indexed facet requests and exact quality checks

For every nonzero (species, country, exact raw lifeStage label), send
one GBIF `/v1/occurrence/search` count-only request:
`taxonKey` from the pinned source, `country` original ISO2,
`year=2010,2025`, `lifeStage` original exact label,
`hasCoordinate=true`, `hasGeospatialIssue=false`,
`occurrenceStatus=PRESENT`, five original accepted
`basisOfRecord` values, `limit=0`, `facet=gadmLevel1Gid`,
and `facetLimit=300`. Save query URL, UTC timestamp, HTTP status,
complete facet category counts and SHA-256 of raw compact JSON.

GBIF documents LIFE_STAGE and GADM_LEVEL_1_GID as occurrence search
parameters, but their indexed normalization can change. Count-only
stage filtering is **accepted only when GBIF query total exactly equals
the original raw stage-label marginal count**. A missing or unsupported
response, a changed total, or an unrecognized GADM facet is NOT_EVALUABLE,
not a biological zero. A stage-conditioned GADM1 count may be less than
the stage total due to missing locality; retain the unassigned part.
For every cell the number of juvenile records in each GADM unit may not
exceed the original all-stage count in that same GADM unit. An original
GADM1 unit absent from a stage-conditioned result does not establish
biological absence.

No individual occurrence JSON bodies, observer identities or point
coordinates are required. Do not infer larvae feeding on introduced
hosts from country/resource flags alone.

## Fixed ecological reporting boundary

Report (1) numbers and fractions of the 39 juvenile-*labeled* records
whose GADM1 subdivision can be resolved; (2) distinct GADM1 divisions
across 11 source cells; (3) unconditioned butterfly GADM1 units versus
stage-conditioned units; (4) stage-label concentration by the frozen
added-host-country flag and **within each focal butterfly species**;
(5) mismatch/partial failures, including cross-border GADM prefixes.

Even if a larval label occurs in multiple provinces, the source does
not include a verified plant host, dated plant arrival or survival to
adulthood. Do not use individual occurrence counts as independent
population or larval fitness replicates. This diagnostic cannot raise
the primary 21-informative-strata global ecological evidence threshold.

Keep this exploration outside the original double-anonymous GEB v0.2
snapshot and do not modify its manuscript. Store safe aggregate outputs,
offline tests, a GitHub Actions completion receipt and limitations.
