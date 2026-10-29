# Expanded island recording panel — implementation and analysis plan

> **For agentic workers:** use the executing-plans workflow for the owned tasks below; root integrates and independently verifies the outputs.

**Goal:** extend the 24-island pilot to a fixed maximum of 96 exact physical islands, test whether recorded butterfly richness supports a natural-isolation association, and create a conservative taxonomy bridge to the original host panel.

**Architecture:** preserve the source-locked GSHHG/WPI measurements and all original pilot selections and receipts. Add geography-only maximin selections, acquire complete GBIF aggregate facets, then calculate prespecified record-count-standardized richness and design-support gates. Taxonomy matching is a separate source-backed join; no country host status is copied onto physical islands.

**Tech stack:** existing Python analysis environment; stdlib acquisition, numpy/scipy statistical calculations, source geopandas/shapely geometry; pytest for consequential numerical/source-identity logic.

**Spec and authorization:** the next biological-test scope in `ISLAND_PORT_BARRIER_RESULTS_20261029.md`, followed by the user's instruction to continue. This plan starts from repository commit `c3fb7bc3e437d1c5aca7fd95d9facaf3b37f216e`. It is an exploratory extension, not a registered confirmatory study. The original 24 outcomes were already known; the new selection and primary analysis decisions below are fixed before inspecting the additional 72 islands' outcomes.

## Global constraints

- Keep original GEB/main and original reviewer deliverables frozen; work on draft PR #102.
- Reuse exact GSHHG 2.3.7 high-resolution dissolved L1 polygons and original WPI snapshot pins. No rectangle or simplified-geometry substitute.
- Retain empty and failed selected cases, distinguishing complete zero records from missing queries. Do not replace islands based on outcomes.
- All results concern retrospective indexed records. They do not establish absence, abundance, arrival, establishment, causal port effects, historical exposure, or functional source populations.
- COL XR string keys and original numeric GBIF-backbone keys remain separate until the explicit audited bridge.

## 1. Fixed geographic expansion

The same 12 strata are used: area [20,100), [100,1000), [1000,20000) km²; continental coast gap certainly <100 or >=500 km; WPI proximity certainly <=50 or >100 km. Representative latitude remains within [-55,55], with Antarctic-ambiguous comparisons excluded.

Retain the original two selected islands in each stratum. Repeatedly select the remaining candidate whose representative point maximizes its minimum spherical distance to already selected islands in that stratum; source island key breaks ties. Stop at eight or candidate exhaustion. This is purposive geographic coverage, not probability sampling or proof of independent ecological replication. Representative coordinates are selection/clustering covariates, never substitutes for coastline distances.

Freeze the complete roster and source hashes before any new occurrence query. Reused pilot queries retain their original receipts and are identified as previously observed. Additional 2010–2017 and 2018–2025 exact-polygon queries use the same COL XR Papilionoidea concept and complete species, basis, dataset and year facets. Check margins, duplicate bins, pagination and an independently fetched first page. Retain errors and resume original successful bytes; no failed query becomes a zero.

## 2. Primary estimand and gates fixed before new outcomes

The primary response is **record-count-standardized recorded richness**, the exact expected number of COL species bins among 20 species-assigned records sampled without replacement from an island's 2018–2025 indexed record pool:

`E[S20] = sum_s (1 - choose(N - n_s, 20) / choose(N, 20))`.

This is a mathematical subsample of records, not individuals, visits, true community richness or sampling completeness. No random rarefaction draw is used. An island is primary-eligible only if its query is complete, at least 20 records have species bins, at least 90% of records have species bins, and at least three calendar years have positive records. All excluded islands remain in coverage tables with reasons.

Primary exploratory OLS model:

`log E[S20] ~ log10(area km²) + absolute latitude / 10 + log10(1 + continental coast gap km) + longitude-quadrant fixed effects`.

Longitude quadrants are fixed half-open bins [-180,-90), [-90,0), [0,90), [90,180], labeled only as longitude blocks; they are not biological realms. Use the sampled coast-gap upper estimate; its <=1 km sampling interval is retained. The isolation coefficient is reported as a recorded-richness ratio for a tenfold increase in (1 + gap), conditional on the included covariates. Broad block/latitude adjustment does not remove unmeasured habitat, climate, regional species-pool or recording differences; even a successful fit is an exploratory recorded-assemblage association.

Do not fit the primary model unless all conditions hold: >=30 eligible islands; >=10 in each continental-gap class; >=2 longitude quadrants each containing both gap classes; >=12 geographic clusters and >=5 clusters represented by each gap class; full-rank design with condition number <=100 after centering/scaling nonconstant columns. Geographic clusters are connected components of the selected roster with representative-point distances <=500 km, frozen without outcomes; they are dependence approximations, not certified independent archipelagos. Covariance uses cluster scores, finite-sample correction and a t reference with clusters minus one degrees of freedom.

Report gate failures rather than relaxing thresholds. The secondary port term adds the fixed within50 versus beyond100 WPI indicator only if all primary gates hold, all four gap × port cells contain >=3 eligible islands, each port class spans >=4 clusters, and the augmented design passes rank/conditioning. No port × isolation × host interaction is fitted. Port coefficient concerns present registry proximity and cannot be interpreted as a transport/arrival effect.

## 3. Prespecified sensitivities and geographic diagnostics

Use the same gates, without tuning to results, for m=10 and m=50; the earlier 2010–2017 window; >=5 positive years; >=90% human-observation records; and dominant dataset share <=90%. These last two are whole-island filters using marginal facets, not within-dataset or within-basis species rarefaction. Repeat the primary covariance with 1000 km roster clusters. Report leave-one-longitude-block and leave-one-cluster coefficient signs only where the reduced design is identifiable; these are influence diagnostics, not replacements for the gated primary model.

Measure nearest other >=5 km² and nearest other >=1000 km² landmass coast gaps for the fixed roster using the already tested spherical sampling/cap pruning. These are secondary geographic access diagnostics. Neither candidate universe defines occupied source populations. No additional biological model is selected after inspecting the resulting gaps.

## 4. Original-host taxonomy bridge

Acquire strict verbose matches for the 239 original names in the same current COL XR checklist. Preserve exact/synonym status, accepted keys, classification, alternatives, source bytes/hashes and retrieval times. Flag fuzzy, non-species, outside-Papilionoidea, unresolved or competing exact accepted concepts; quarantine old-backbone and current-COL many-to-one collisions. Inspect the overlap of conservatively eligible concepts with island species facets. Do not infer local host presence from national plant ranges or synonymy alone.

## Implementation tasks and checks

- [ ] Expansion owner: add `scripts/acquire_expanded_island_panel_20261029.py`, selection/source tests, and `data/exploratory/island_expansion_20261029/`; verify preservation, maximin ties, source identity, complete/zero/missing distinctions and actual acquisition.
- [ ] Root: add `scripts/analyze_expanded_island_recordings_20261029.py` and meaningful tests for analytic rarefaction against enumeration, eligibility, clustering and gated fitting; produce auditable metrics/support/model outputs.
- [ ] Root: add `scripts/measure_expanded_island_access_20261029.py` and tests for focal exclusion/candidate area thresholds; run measured coast diagnostics for the frozen roster.
- [ ] Taxonomy owner: add `scripts/bridge_island_host_taxonomy_20261029.py`, ambiguity/collision tests and `data/exploratory/island_taxonomy_20261029/`; acquire and audit source matches.
- [ ] Root/reviewer: independently reconcile source receipts and numerical results, review interpretation, exclude new exploration from the original reviewer bundle, run the required suite and CI, and save the continuation to PR #102.

## Review focus

Selection drift and previously known pilot outcomes; taxonomic namespaces and many-to-one joins; failed requests versus zeros; species-assignment residuals and unequal recording; geographic support, rank and cluster dependence. In particular, a model cannot create missing biological contrasts, and rarefaction cannot remove selective reporting.

## Method references

- [vegan authors' rarefaction documentation](https://vegandevs.github.io/vegan/reference/rarefy.html): the Hurlbert expectation for without-replacement subsamples. Here the finite pool is explicitly records, not inferred individuals.
- [GBIF taxonomy interpretation](https://techdocs.gbif.org/en/data-processing/taxonomy-interpretation) and [occurrence API](https://techdocs.gbif.org/en/openapi/v1/occurrence): namespace-specific concepts and retrieval.
