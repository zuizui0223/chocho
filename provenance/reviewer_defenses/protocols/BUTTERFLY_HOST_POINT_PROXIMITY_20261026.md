# Pre-outcome protocol: exact butterfly and candidate-host GBIF point proximity

Date: 2026-10-26. STATUS: new exploratory observational **feasibility and falsification** branch stacked above chocho PR #72. Frozen *before* acquiring any point-level GBIF occurrence bodies or measuring point distances. This study does **not** test a realized climatic niche shift, butterfly population expansion, host-driven recruitment, or ecological competition.

## Fixed original observation and resource cohort

Inputs are already source-verified in the parent PR #51:
- \`gbif_all26_source_host_country_inventory_20261023.json\`: **26** existing positive butterfly×country cases with introduced WCVP/HOSTS/BCE accepted larval-host resource candidates, 16 butterflies / 11 countries, 24 accepted plant species, **46** butterfly×host×country links and **44** unique accepted plant×country queries. The original frozen 48 positive butterfly×country comparison also contains 22 non-added cases, which are **NOT true nonhost or absence controls**.
- \`gbif_positive48_gadm_lifestage_facets_20261020.json\`: exact accepted butterfly species keys, original 2010–2025 source facet filters and indexed counts; all 48 originally complete.
- \`gbif_all26_host_gadm_source_20261023.json\`: exact accepted plant species GBIF keys, original botanical source queries, **44/44** complete country facets and 3,595 indexed botanical records in the earlier snapshot.

Freeze these names/keys/query identities and the source cross-SHA validation **before** point records are accessed. Do not expand taxa based on successful host co-detection. All source names are candidate recorded larval hosts, not authenticated local oviposition.

## New live GBIF source snapshot

Request the standard documented \`GET /v1/occurrence/search\` with exact original filters: source accepted \`taxonKey\`, country ISO2, \`year=2010,2025\`, \`hasCoordinate=true\`, \`hasGeospatialIssue=false\`, \`occurrenceStatus=PRESENT\`, and precisely the source's five basisOfRecord values (HUMAN_OBSERVATION, PRESERVED_SPECIMEN, OBSERVATION, MACHINE_OBSERVATION, OCCURRENCE). **No stage or host-text filtering** at acquisition. Convert original \`limit=0&facet=...\` query URLs to paginated \`limit<=300&offset=0,300,...\`. The GBIF search API has a documented 100,000-record hard offset+limit maximum; this pilot instead caps each taxon×country query at 5,000 records and the complete newly fetched panel at 12,000, with nonzero over-limit queries explicitly marked unscorable (not truncated and not zero). Set >=1.05 seconds between queries and bounded retry for HTTP 429/5xx. Fail closed on transport and count/page inconsistency.

Keep per-query request URL template, current count, original source facet count, source response SHA hashes, retrieval timestamp, failed query reasons, taxon cross-checks and semantic exclusion totals. No raw occurrence coordinates, collection localities, observer identities, or occurrence-body JSON will be pushed to the repo or uploaded as an artifact. Both 2010–2025 windows are pooled at retrieval; source index can drift, so never silently claim this is the identical historical index.

## Primary spatial and temporal measurements

For every original source-positive butterfly×country case:
1. Validate original GBIF point coordinates and year 2010–2025; reject nonfinite/out-of-range lat/lon, (0,0), invalid years, and records whose speciesKey/acceptedTaxonKey resolves to a different accepted source species (source mismatch stays separately recorded).
2. Use only the exact source accepted host species in **that butterfly×country**. Match each butterfly record to the geographically nearest indexed candidate-host record from the SAME country whose index year is within ±3 calendar years.
3. Primary outcome: fraction of all coordinate-valid butterfly GBIF indices whose nearest year-matched source-listed plant is within **5 km great-circle distance**. No botanical indexed record means **no GBIF co-detection**, not proven plant absence. Also report denominator-positive eligibility: how many butterfly indices have ANY plant candidate recorded in the fixed window.
4. Predeclared spatial sensitivity at **1 km and 25 km**; index-year sensitivity at exactly same year (0) and all years (no year constraint). No search over new thresholds after results.
5. A strict **known coordinate uncertainty <=1 km on both butterfly and plant observations** subset must be reported separately and cannot be replaced with unknown-uncertainty data. The all-valid primary is an optimistic index-level upper-opportunity metric, not confirmed true local interaction.
6. Preserve all records including missing life stage, and separately report labels for Egg/Larva/Caterpillar/Pupa. For independently verifiable local host-use annotations require a juvenile label, exact accepted candidate-host binomial in original text, and an unambiguous usage/action phrase; otherwise report "host use unverified". Do NOT convert administrative/local proximity to larval feeding or establishment.
7. Report a species×country case table, uncertainty/missingness and dominated-case exclusions (e.g. leave out the top contributing butterfly country) for descriptive source concentration. No country/record-independent p-values; all 26 cases have butterfly presences selected by design.
8. A stricter proximity measure counts butterfly indices and separately distinct coarse coordinate bins (0.01° lon/lat bins), to show duplication/reporting concentration; 0.01° cells are **not** identical-area units.

## Information and scientific decision gates

Source-complete only when all original **26 butterfly queries and 44 unique plant queries** finish (including genuine zero-result queries) and all declared source taxon identities are maintained; if not, emit **INCOMPLETE_SOURCE_SNAPSHOT**, explicitly listing incomplete queries and never turn them into ecological zeroes.

A point-level spatial overlap number is descriptively reportable if at least 5 of the 26 butterfly×country cases have >=5 coordinate-valid butterflies and at least one source candidate plant record; otherwise label an underpowered feasibility pilot. **Even when complete, proximity in space/time is an observational co-indexing statistic, not host use or a causal realized-niche result.**

The next separate macro hypothesis remains the 73 climate-analogue versus climate-novel positive species×country records in PR #72, sampled without cherry-picking; do not extrapolate this preselected 26-case analogue-biased pilot to that panel.

Original GEB draft PR #38, main, and previously completed PR #72 must remain untouched.