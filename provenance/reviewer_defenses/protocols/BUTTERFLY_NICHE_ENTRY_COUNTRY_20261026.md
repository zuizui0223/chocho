# Butterfly potential climate-niche entry: fixed source-index comparison

**Frozen before the first new niche-class analysis, 2026-10-26.** Status: new exploratory, posthoc to the source creation, **NOT a confirmatory ecological or butterfly colonisation test**. The fixed GEB manuscript PR #38, main, and prior #46/#51 results must remain unchanged.

## New question

Among countries outside every documented WCVP native-host resource polygon of a focal butterfly, are **indexed adult-or-unspecified butterfly presences** more frequently recorded in countries where the added potential host area is predominantly *outside* (new climate resource opportunity) versus *inside* (climatically analogous host-resource release) the focal butterfly's **native HOST-PLANT** climatic envelope?

This is **NOT** measured butterfly climate tolerance, niche evolution, range expansion in time, establishment, reproduction, larval host use, or ecological competition.

## Frozen source, cohort, spatial estimand

- Baseline: the original 239 butterfly HOSTS/WCVP graph with 1,027 evidence-filtered BCE added accepted links to 81 European-skewed taxa, as already represented by the same fixed #51 source-derived `country_resource_area_audit.csv.gz` (5th–95th native-host BIO1 mean temperature + BIO12 annual precipitation rectangle).
- Country geography: existing 140 source-qualified Natural Earth × WGSRPD3 country units with exact disjoint area atoms. An excluded or unmapped country cannot be recoded as butterfly/host absence. Exclude a species×country whenever `native_guard_present` is true, using the original WCVP guard buffer. All fractions use entire mapped country area as denominator, not equal WGSRPD region counts.
- Source butterfly observations: **unchanged** full-index GBIF country facets (2010–2017 and 2018–2025); select the 229 unambiguous accepted keys with both fixed-window queries complete. Never add overlapping `core` and `all_accepted_bases` basis groups. Records are indexed occurrences; zero records are **not** verified species absences.
- Do not retrieve new GBIF records, author-only WCVP/HOSTS/CHELSA/TDWG data, or arbitrary species from their absence/presence outcomes. The country area audit was constructed once from the existing author-only source; reproduce SHA, row identity, and numeric constraints.

## Pre-outcome main classification, identical for all taxa

For each source-qualified species×country cell, let `a` be the fraction of the **entire country** covered by added host resources, `b` the fraction covered by the intersection of those added resources with the native-host climate rectangle. Use `novel=a-b`, `analogue=b`.

Only outside-native-host countries with `a >= 0.05` enter the primary pairwise opportunity comparison. Define:
- **analogue-dominant (potential constraint release)**: `b / a >= 0.80`.
- **climate-novel-dominant (potential climate niche entry opportunity)**: `(a-b) / a >= 0.80`.
- All other positive-`a` countries = mixed and excluded **from the paired contrast**, but retained in source/coverage audits. `a<0.05`, exactly zero, native-host countries, unknown climate countries, and unqualified butterfly taxa receive explicit disposition counts.
- An indexed butterfly record of any stage cannot be called colonisation; missing host-use labels remain missing. The climate envelope describes *recorded native hosts*, not butterfly realized or fundamental climate niche.

## Main detection contrast and information gate

For all 229 resolved taxa and 140 eligible countries, use pooled 2010–2025 **all_accepted_bases** counts; focal taxon is excluded from its own effort indicator. Country must have >=3 *other* complete focal butterfly taxa with >=1 indexed GBIF record.

Within each **butterfly species × receiving WGSRPD Level-1**, compare only the two dominant types when both are present. Let N = eligible analogue+novel countries, A = novel countries, K = focal-positive countries (>=1 indexed record), O = focal-positive novel countries. Under a fixed-margin artificial geography null, O has hypergeometric expectation AK/N and distribution Hypergeom(N, A, K). Sum independent *mathematical draws* across strata to obtain an exact convolution **conditional reference**, not a biological randomization test. Estimate `100*(sum(O)−sum(AK/N))/sum(K)` percentage points; direction >0 means greater recording in climate-novel than analogue opportunities **conditional on source index, preselection and matching**.

Minimum information gate: >=50 strata with `0<K<N` and both opportunity classes, >=25 contributing butterfly species, >=25 positive species×country cells. Report an unqualified/no-conclusion status when any criterion fails. Emit counts and conditional null anyway as **nonconfirmatory diagnostics**, not proof of a zero effect. Species sharing the same introduced host, geography and GBIF contributors violate independent ecological replication.

## Frozen sensitivity panel (NOT outcome-based thresholds)

- Opportunity purity `0.80` (primary) and `0.90`; country added-host area minimum `0.05` (primary) and `0.20`.
- Other-observed-butterflies effort gate `3` (primary), `1` and `5`.
- GBIF basis `all_accepted_bases` (primary) and `core` (sensitivity, overlapping).
- Record threshold >=1 (primary), >=5 (indexing-depth diagnostic).
- Separate 2010–2017 and 2018–2025 windows (neither proves host arrival time). Identical source rows and classification across these windows.

Explicitly distinguish *negative detection contrast* from *lack of statistical identification*. Do not re-label failed threshold configurations as biological evidence; do not cherry-pick positive configurations.

## Source and inference red lines

1. Georeferenced coarse country presence is not an occurrence in the added-plant subarea; even a country with 90% novel host-resource area does not prove the observed butterfly occupied that resource.
2. No host-arrival dates, host availability at larval stages, land-use or standardized zero-survey sampling. Therefore this cannot test realized climatic niche shifts, dated range expansion, actual demographic benefit, or causal release from host limitation.
3. Species-level novel climate is the novel **source plant climate proxy**, not known physiological tolerance. Multidimensional climatic novelty, warm-vs-cold direction, trait adaptations, host quality, and butterfly interactions remain untested.
4. Primary source age/time and native HOSTS lists are fixed; no recoding after viewing observed butterflies. Supplemental BCE is partially European.
5. Do not change GEB submission bundle or upgrade a nominal mathematical upper/lower tail to an ecological p value. Source-only aggregate outputs may be archived; copyrighted original botanical distribution data may not.

## Deliverables

- Reproducible offline script and synthetic tests on source-independent fixtures.
- Stable summary JSON: exact source checks, cohort/area audits, primary fixed classification and conditioning, null and uncertainty, sensitivity panel, explicit inferential gate.
- Optional CI artifact and concise ecological verdict. **No scientific discovery status until the frozen full source run has succeeded and passed integrity checks.**
