# Independent country-grain GBIF macro analysis — 20 October 2026

## Why a separate analysis is necessary

The original 135-species, WGSRPD3 ordinal-page result remains unchanged.
Recovery transport repeatedly truncates ordinary GBIF record bodies, including
75-, 15- and 10-record requests. A compact `limit=0&facet=country` pilot succeeds
and can describe the full live matching index without the 100,000-result paging
ceiling. This is a **new country-grain source**, not successful recovery of the
96 incomplete WGSRPD3 samples and not the same estimand at a higher sample size.

This protocol was written before inspecting country occurrence/host associations.
It follows the unsuccessful transport investigation and a source-feasibility
pilot, so it is an exploratory, prospectively specified continuation, not an
independent preregistered confirmatory study.

## Fixed taxonomic and occurrence source

Use the original 239-name source receipt, SHA256
`18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2`.
Keep its accepted GBIF keys; do not re-resolve live names to improve coverage.
Eight unresolved taxa stay unresolved. Oeneis bore and Oeneis polixenes both
map to accepted key 9804158; exclude **both** from ecological analyses before
looking at outcomes. Maximum non-colliding target/background cohort: **229**.
Retain all 239 source receipts and report totals by unique accepted key as well
as by source name; the latter otherwise double-counts a shared GBIF index.

Acquire country facets for 2010–2017 and 2018–2025, with coordinates,
`hasGeospatialIssue=false`, `occurrenceStatus=PRESENT`, and the following scopes:

- Primary accepted bases: HUMAN_OBSERVATION, PRESERVED_SPECIMEN, OBSERVATION,
  MACHINE_OBSERVATION, OCCURRENCE.
- Core sensitivity: HUMAN_OBSERVATION and PRESERVED_SPECIMEN.

Both windows must complete and validate within a basis scope for a species to
enter that scope. Check unique country codes, nonnegative integer counts,
facet totals no larger than the reported matching count, and no facet-limit
truncation. Preserve query URLs, timestamps, raw-response hashes, source keys,
missing-country counts and failed queries. A zero count is no indexed record
under these filters, not a biological absence. Index totals are occurrence
records, not deduplicated organisms, populations or independent observations.

## Geographic crosswalk fixed independently of occurrence outcomes

Intersect source-frozen WGSRPD3 polygons with Natural Earth 1:50m admin-0 **map
units**, pinned to `nvkelso/natural-earth-vector` commit
`117488dc884bad03366ff727eca013e434615127`. Use audited ISO_A2_EH country codes
and dissolve units sharing a code. Use equal-area EPSG:6933 areas. Preserve
source hashes, all nonzero intersections, and country-code/boundary anomalies.
Partition overlapping source polygons into exact disjoint membership areas;
count each selected geographic portion once when evaluating a resource or
climate subset. Unmapped land stays in the full country-area denominator.

A material link has intersection area / min(country area, WGSRPD3 area) >=0.05.
A wider guard link uses >=0.005. Primary countries require >=95% geometric
coverage by material links, one WGSRPD Level-1 continent across material links,
all material and guard units available in the 365-unit climatic source, and no
unresolved ISO or documented map-unit/country-definition ambiguity. Do not
silently assign disputed/unmapped map units to neighboring ISO countries.

For each butterfly-country pair:

1. Exclude the country from the outside-native-resource comparison if **any
   guard-linked** WGSRPD3 unit has a documented native host. This protects
   against declaring a country native-unsupported after dropping a small part.
2. Mark an added-host country when at least one **material-linked** unit has
   introduced-added host opportunity, provided step 1 does not exclude it.
3. Require at least 90% of its **total country land area** to fall in
   the butterfly's native-host BIO1/BIO12 5th–95th percentile analogue.
   This inherits the native-host climate proxy; it is not a butterfly niche
   estimated from occurrences. Keep the unmapped remainder in the denominator;
   do not renormalize it away. The independent methods review made this
   denominator clarification before any outcome association was inspected.

These are country-scale potential-resource labels. A butterfly record and a
host may occur in different parts of a country. Neither co-occurrence, local
host use nor breeding is established.

## Primary comparison and fixed sensitivity checks

Primary: all accepted bases; pooled windows; >=1 focal record in a country;
other-focal-butterfly richness >=3; q=0.05 climate bounds; country climate
coverage >=0.90; species × Level-1 strata. Targets and the effort background
are all source-complete, non-colliding taxa. Other-species richness excludes
the focal taxon. The background is an internal reporting proxy, not measured
survey effort.

Within each stratum let N be eligible countries, A added-host countries,
K positive countries, and O positive added-host countries. Retain resource
contrasts with 0<A<N. The conditional uniform-country expectation is AK/N.
The primary effect is 100 × sum(O−AK/N) / sum(K), in percentage points.
Use 9,999 conditional hypergeometric draws and 2,999 genus-cluster bootstrap
draws for the primary, deterministic seed 20261020. Report counts even if the
information gate fails; require >=50 strata with 0<A<N and 0<K<N, >=25 mixed-
resource target species and >=25 positive country-species cells for the
sample-size label. All results remain exploratory if this gate is passed.

Prespecified sensitivities:

- Both basis scopes × pooled/early/late windows × other-species richness
  gates 1/3/5 × q=0.05/0.10 × focal positive threshold 1/5 records (72 cells).
- Primary settings with all material units in the climate analogue;
  with a single material-linked WGSRPD3 unit; and with introduced-added
  coverage >=90% of total country area (exclude intermediate exposures
  rather than reclassifying them as no added resource).
- Primary with country area classes (<10,000, 10,000–100,000, >=100,000 km²)
  or other-species record intensity classes (<100, 100–9,999, >=10,000 records)
  added to species × Level-1 strata. These are bias diagnostics, not recovery
  of true effort or causal exchangeability.
- A continuous area-score diagnostic at primary settings: within each retained
  species × Level-1 stratum compare the mean fraction of total country area
  that is both introduced-added and climate-analogue in positive countries
  with the K-conditioned uniform-country expected mean. Preserve K, shuffle
  country labels without replacement, and use the same genus bootstrap. This
  avoids treating a small supported portion as equivalent to an entirely
  supported country. Report native, added, analogue and added-within-analogue
  fractions for audit. This remains a coarse opportunity score, not local use.
- Original-complete source taxa (after the same key-collision exclusion) as
  both target and effort cohort, and all available targets with that original
  effort background. These separate target expansion from effort expansion
  **within the new country data**. They cannot reproduce a WGSRPD3 result.

Use 999 null draws and 499 genus-bootstrap draws for non-primary settings.
Publish the full grid, including failures of the information gate and changes
of sign. Do not promote the largest association or smallest null-tail value.
Do not interpret temporal contrasts as colonization or population change.

## Interpretation and deliverables

The ecological question is whether observed country-level butterfly coverage
outside documented native-host resource geography is enriched in countries
with introduced-added resources, conditional on the stated geographic and
reporting filters. Shared plant ranges, country area, uneven GBIF contributions,
catalogue omissions, adult movement and coarse resource/climate labels can
produce or mask associations. The conditional null is a geographical benchmark,
not a biological or causal p-value; genus intervals do not capture all shared
geographic/host dependencies.

Save reproducible code, tests, public GBIF country aggregates and acquisition
receipts, derived crosswalk audits, and analysis summaries. Do not redistribute
the author-only HOSTS/WCVP/TDWG raw archive. Preserve the original 135-species
replay and report this new result beside it, with its changed spatial grain.
