# GloBI larval feeding bridge: outcome-blind public-source feasibility gate (2026-10-24)

**Scope:** investigate whether `chocho` can replace potential host-map overlap with geographically and life-stage-resolved *observed* butterfly–plant feeding, without reanalysing published single-species experiments. Exploratory acquisition, **no ecological effect estimated or inferred**.

## Current epistemic boundary
The frozen 239-butterfly HOSTS/WCVP reconstruction measures potential resources, not recorded caterpillar feeding or emergence. The 229-butterfly GBIF country index lacks usable explicit exact host-feeding tags in the audited 39 juvenile records. This is an October 2026 source-feasibility pilot, not a new observation season. Previously reanalysed monarch and Kyoto experiments are prior art, not independent global recruitment evidence.

GloBI aggregates recorded interaction links with provenance, possible specimen life stage, event date, and coordinates. Its live API is **not version-stable** and may include the original HOSTS literature or duplicated sources; a GloBI hit is not automatically independent of HOSTS. Its versioned Zenodo release (doi:10.5281/zenodo.3950589) is required for later publishable full-population analysis. GloBI references: https://www.globalbioticinteractions.org/data and https://www.globalbioticinteractions.org/how-to .

## Fixed sampling and acquisition (before results)
Use the 239 exact species names in frozen `data/frozen/figure_sources/anthropogenic_species_metrics.csv`; verify its expected Git blob SHA `9f83d2d9ac8e46f1cfd4ded572c7a023273503dd`. Select **20 focal butterflies** with the smallest unsigned SHA-256 hash of `20261024|<exact scientific name>`, independent of ecology. Add two predeclared query-system positive controls `Danaus plexippus` and `Pieris brassicae` when not already selected. Controls are **not** included in the random-denominator readiness gate. Save selected names and complete per-taxon queries.

One `https://api.globalbioticinteractions.org/interaction.csv` request per species with `sourceTaxon=<exact name>`, `includeObservations=true`, `type=csv`, `limit=1024` and specific fields:
`source_taxon_name,interaction_type,target_taxon_name,target_taxon_path,source_specimen_life_stage,event_date,latitude,longitude,study_source_id,study_citation,study_source_citation`.

Timeout 25 s per request; one attempt, no result-based retries. Record response status, byte length, sha256, row count, source IDs; failures/empty/non-exact replies remain **missing, not negative ecological evidence**. A source reaching exactly 1,024 rows is `PAGINATION_REQUIRED` and cannot be declared complete.

## Predeclared descriptive QC
- Exact source taxon name equality (binomial), evidence of plant target path (Plantae or Viridiplantae); unknown target kingdom is classified unknown rather than guessed from plant-like name.
- Interaction types categorize `eats`, `feedsOn`, `hasHost` or related explicit trophic actions vs pollination/visitation, and are reported separately.
- Strict immature-source stage matches `larva|larval|caterpillar|juvenile|instar` case insensitive; missing stage is not imputed.
- Strict source independence: keep all distinct `study_source_id` values; mark any containing `nhm`, `hosts`, `naturalhistorymuseum`, or `globalbioticinteractions/hosts` as **potential HOSTS-derived** (conservative, may overexclude). Shared citations cannot be assumed independent unless provenance is inspected manually.
- Georeference: both valid numerical latitude in [-90,90] and longitude [-180,180], excluding artificial missing 0/0. Calendar date: parse ISO-prefix `YYYY-MM-DD` or `YYYY` with year 1800-2026, but a year without day/month is only `year-known`. Keep counts of both.
- Candidate realized **feeding evidence** requires an exact focal butterfly, plant target, explicit trophic interaction, explicit immature stage, coordinates, and non-HOSTS-labeled source. This is **use**, not recruitment or local native/introduced host origin.
- Zero candidate rows with any missing source/transport/truncation is not evidence for global absence.

## Prospectively defined decision
If the full 20 response-blind sample is retrieved with no missing/truncated requests, at least **10/20** taxa have `>=1` high-specificity life-stage/georeferenced non-HOSTS-marked trophic record from at least **5 distinguishable source IDs total**, then proceed to a pinned stable-version GloBI full-239 acquisition and a continent/host/effort matched hypothesis. Otherwise this public source does **not yet** bridge the potential-to-realized gap at macro scale. No posthoc relaxation of strict stage/geography or source independence to claim a positive effect.

The gate requires a **source-verified** externally generated receipt and code. Preserve read-only GEB `main`, primary disturbance branch and original GBIF cohort; do not mix exploratory source-level evidence with adult recruitment.


## Source-transport amendment v0.2 — after v0.1 source gate (explicitly post-first-acquisition)

The first acquisition returned 20/20 fixed random taxa but **five** reached the API default/capped `limit=1024` and are missing source-complete records; both positive controls can also be capped. This violates the prospective completeness prerequisite. The following is a **transport-only correction** before the independent full replay, not a change in focal species, ecological labels, response or 10-taxon/5-source gate:

1. Reissue the **same 20 original fixed taxon queries** at `limit=1024`, `skip=0,1024,2048...` until a page is shorter than 1024. A fixed maximum of **12 pages per focal species** applies; cap overflow is missing, never a biological negative. Positive controls use only one page and remain excluded from the gate.
2. Archive each page's immutable read-time response SHA256, row count, offset and HTTP status; reject exact duplicate page-body SHA (indicating skipped pagination) and unknown headers. Any endpoint error makes a taxon incomplete; no result-based retry.
3. Do not silently deduplicate identical records or treat multiple GloBI source rows as independent organisms. Stage/plant/trophic/latitude/source definitions remain byte-for-byte consistent with v0.1. Field-level missing values remain missing.
4. Keep the **v0.1 receipt unchanged**. Store v0.2 page-corrected results separately and label live-source race/drift. Source provenance IDs or GloBI linkage still cannot verify independence from HOSTS without deeper citations/underlying data checks.
5. Maintain the original decision: full 20 focal taxa complete and >=10 with strict evidence across >=5 distinguishable source IDs. If not met, stop the geographic real-use macro claim and return to experimental or original source data rather than lowering the evidence threshold.
