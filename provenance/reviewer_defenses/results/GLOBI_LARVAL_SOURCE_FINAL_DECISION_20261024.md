# Final decision: independent realized larval host-use macro bridge remains unsupported

**24 October 2026. Source-verified exploratory acquisition.** Separate from the frozen 239-butterfly GEB manuscript and earlier disturbance trait exploratory branch.

## What was executed
- Pre-registered-as-workflow (not a study preregistration) response-blind species-name hash chose **20/239 butterfly taxa**, plus query controls. The source set was fixed in the [protocol](../protocols/GLOBI_LARVAL_HOST_SOURCE_FEASIBILITY_20261024.md) before retrieving outcomes.
- Original [source receipt v0.1](globi_larval_host_feasibility_20261024.json): **15/20 complete**, 5 reached row limit 1,024; 14 provisional strict metadata matches in 5 sampled species, drawn from two GloBI dataset source IDs. This was incomplete, not a species absence or a confirmed field feeding survey.
- The post-first-pilot transport-only amendment used GloBI `skip` pagination, same original 20 taxa, unchanged GloBI source-stage-plant-location definitions and unchanged 10/20 + five-source criterion. The [source receipt v0.2](globi_larval_host_feasibility_paginated_20261024.json) documents page offsets/SHAs and a [successful GitHub Actions job](https://github.com/zuizui0223/chocho/actions/runs/40078420682).

## Exact v0.2 outcome
| Source availability/statistic | Count |
|---|---:|
| Fixed response-blind sample | **20 butterfly species** |
| Original selected species with complete paginated retrieval | **19/20** |
| Transport error (not biological absence) | **1/20** — `Elymnias hypermnestra` |
| Number of species with one or more provisional strict metadata record in retrieved v0.2 pages | **5/19** |
| Provisional record rows passing mapped plant, named consumer, `eats`, larval-stage label, numeric coordinates, non-HOSTS-labelled data source | **26 source rows** |
| Distinct contributing dataset source IDs | **2**: GloBI iNaturalist and GloBI Menares 2025 herbivory |
| Maximum qualifying species if the one missing species is positive | **6/20** |
| Fixed prospectively defined minimum for scaling to a real-use macro study | **10/20 species and ≥5 qualifying data-source IDs** |

The prospective source-feasibility gate **cannot be passed even under best-case imputation for the one incomplete butterfly**. No further API retry or selection of another convenient butterfly cohort is necessary to reach that decision. A source-missing butterfly is not coded as zero in any ecological model.

## Important biological-status downgrade after inspecting original dataset mappings

The [source-provenance audit](GLOBI_LARVAL_SOURCE_PROVENANCE_VERDICT_20261024.md) confirms Menares's GloBI configuration assigns **region-wide centroid coordinates** and an historical time-range to literature-derived butterfly-larval host associations. Those source records are **not individually georeferenced field feeding events** even when the API reports numeric coordinates, a date-like string and the life stage `larva`. GloBI iNaturalist uses contributed interaction observation fields mapped to `eats`; source imagery, explicit plant identification and original observation IDs have not been independently audited. Therefore **26 source rows do not equal 26 confirmed caterpillar feeding encounters, let alone 26 independent breeding successes**. Multiple rows within one dataset are not biological replicates.

It would be invalid to interpret the 5/20 metadata coverage as the real-world prevalence of larval host use. The sample was selected independently of outcomes, but data availability is highly heterogeneous, the data source is opportunistic and adult colonization was never measured.

## Decision and conservation boundary

**No new supported general claim** that globally redistributed host plants actually increase butterfly recruitment, local coexistence or native butterfly displacement. Do not insert this exploratory dataset into PR #38, treat missing source links as negatives, collapse 19 individually surveyed taxa into independent countries, or publish the source-availability ratio as an ecological effect.

**What survives:** the 239-species geographic *potential* resource expansion, contribution concentration, resource homogenization, climate-conditioned small residuals and known-host support redundancy remain separately source-supported with their stated limits. The prior strong-disturbance trait association remains provisional because the adequate-host catalogue subset weakens its magnitude and CI includes zero.

**Next genuine advance is conditional on new independent evidence**, not further reshuffling HOSTS or GloBI filters: versioned source-verified individual local larval plant use with actual dates and locations, source image/observation IDs, duplicate detection, and survival-to-adult or independently verified establishment in a sufficiently large cross-species panel. Until that exists, prioritize the already source-supported GEB resource-geography manuscript rather than claiming an unverified demographic mechanism.

All GloBI pilot scripts, source receipts, protocol and provenance audits remain isolated in `analysis/globi-larval-host-source-pilot-20261024`.
