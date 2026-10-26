# A coarse larval-host opportunity does not automatically produce nearby co-recorded host plants

**2026-10-26 — complete, prefiltered small-volume public GBIF source audit.** This document is an *early observational feasibility diagnostic*, **not** the all-26 cohort primary result and not butterfly climatic niche expansion. The fixed diagnostic selection was committed **before accessing point-level outcomes** in [protocol](../exploratory/BUTTERFLY_POINT_LOWVOLUME_DIAGNOSTIC_PROTOCOL_20261026.md).

## Source and realized measurement contrast

Original source: 26 butterfly-positive introduced-host opportunities, 44 exact WCVP accepted plant×country source pairs and verified GBIF species keys; choose *every* source butterfly×country case with prior butterfly source count <=20 and total candidate botanical source count <=20. **12/26 cases** qualify, without conditioning on whether the plant has an indexed record, or on butterfly/host distance. Exact original accepted GBIF taxonKey, ISO2 country, 2010–2025 year window, valid coordinates/geospatial flags, PRESENT and the same five original basisOfRecord categories were retained. All **30 unique butterfly/plant query identities** produced complete responses: **zero missing sources and zero indexed-count drifts** versus the earlier frozen facet archive.

- 46 complete butterfly indexed point records from 12 cases.
- 40 source-candidate botanical indexed records summed over the 12 butterfly cases (not 40 independent botanical populations or necessarily 40 unique source individuals, because host source may be reused).
- 8 cases have at least one candidate plant GBIF point record and 4 have no indexed candidate botanical point record; the latter is **not** verified plant absence.
- **Within 1 km, ±3 years:** 0/46 butterfly index records.
- **Within 5 km, ±3 years (prespecified main diagnostic): 0/46**.
- **Within 25 km, ±3 years:** 10/46 (21.7%). Three contributing cases are *Anteos clorinde* / Dominican Republic 3, *Colias philodice* / Honduras 1, *Papilio demoleus* / Haiti 5, plus *Phoebis agarithe* / Haiti 1.
- Only 26/46 butterfly point records document coordinate uncertainty ≤1km; zero of those have a known-uncertainty-≤1km botanical point within 5km and three years. Unknown source uncertainty never counts as high precision.

**This is a measurable cross-scale discrepancy, not an inferred ecological limitation.** The original country-level analyses counted butterfly and candidate-host plant resource availability within the same country or coarse botanical regions. Point-level GBIF source coordinates reveal that in this **small and inherently sparse** prefiltered sample, no observation has a source-documented candidate plant within the narrow 5km/3-year neighborhood. The most honest interpretation is the **failure of the existing observational coverage to verify local host access**, not a proven lack of host plants, feeding, colonization, fitness, or habitat quality.

### Individual case evidence, no raw locations

| Original butterfly source / country | Butterfly indexed points | Candidate plant indexed points | 5 km/3 years | 25 km/3 years |
|---|---:|---:|---:|---:|
| Aguna asander / CU | 2 | 2 | 0 | 0 |
| Aguna asander / DO | 9 | 4 | 0 | 0 |
| Aguna asander / HT | 3 | 1 | 0 | 0 |
| Anteos clorinde / DO | 15 | 2 | 0 | 3 |
| Anteos clorinde / HT | 1 | 0 | 0 | 0 |
| Battus philenor / GT | 1 | 6 | 0 | 0 |
| Colias philodice / HN | 1 | 5 | 0 | 1 |
| Kricogonia lyside / HN | 3 | 0 | 0 | 0 |
| Papilio demoleus / HT | 8 | 9 | 0 | 5 |
| Papilio machaon / NA | 1 | 0 | 0 | 0 |
| Phoebis agarithe / HT | 1 | 11 | 0 | 1 |
| Vanessa atalanta / GA | 1 | 0 | 0 | 0 |

## Integrity, caution, next test

This independent live public snapshot had **30/30 query responses complete**, no transport errors, 46/46 matched source butterfly counts, and current botanical indices consistent with all 12 original source case sums. Only aggregate case counts were saved to [machine-readable result](../../provenance/reviewer_defenses/results/butterfly_host_point_lowvolume_20261026.json); **no raw occurrence keys, coordinates, observer identifiers or exact locality are included**. Do not interpret the 46 indices as 46 independent wild butterflies, or compute binomial confidence intervals on them.

The all-26 original study stays in the source-frozen [PR #73](https://github.com/zuizui0223/chocho/pull/73), with its separate GitHub Actions [workflow](https://github.com/zuizui0223/chocho/actions/workflows/butterfly-host-point-proximity-20261026.yml). The source-only full-26 workflow has not yet produced an accepted complete-cohort result. **Do not pool this 12-case convenience subset with that main cohort.** 

For a true niche-entry conclusion, additional target evidence must link **the PR #72 climate-novel versus climate-analogued 73 positive source observations**, butterfly point climates and host plant range/local occurrence, and direct larval host-use/repeated wild breeding with effort controls. Neither the current 12-case co-indexing nor a future 26-case descriptive pilot can establish this on its own.
