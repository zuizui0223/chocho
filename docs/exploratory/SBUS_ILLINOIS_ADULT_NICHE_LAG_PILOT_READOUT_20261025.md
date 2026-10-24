# US butterfly adult counts × documented larval host niches: Illinois preliminary readout
**25 October 2026 | exploratory original-data pilot | no evidence of causal competition or anthropogenic-host-driven population increase is established.**

## Source provenance and identifiability

Henry et al. (2025), *Twenty years (2000–2020) of butterfly monitoring data across the contiguous United States*, *Scientific Data*, doi:10.1038/s41597-025-05513-8. Archived Figshare dataset doi:10.6084/m9.figshare.27934602.v1 (CC BY 4.0); 55,066,096-byte ZIP, MD5 `3161b1c8802a0d67771b2044e18da928`, uncompressed `4_res/complete-share.csv` contains exactly **1,321,882 actual species-survey count records** and 19 separately labeled monitoring programmes, not GBIF observation-frequency proxies.

Only exact names among original 239 butterflies were retained, linked to source-verified accepted WCVP host IDs from pinned HOSTS `808e0b869f9ec1adf8efff87cf6a395adda103e0` and WCVP `65bed76bae9d644ccb6ad200c05f9f5071d89e05`. Documented host sharing is a trait-proxy and does not verify that a sampled caterpillar at a specific U.S. monitoring transect ate that plant. The source adult count survey does not include historical plant introduction dates, so **no inference of butterflies becoming more abundant because nonnative plants arrived** is licensed.

## Original pilot eligibility

Focal **Illinois Butterfly Monitoring Network**: 121 qualified site identities, 851 site-years across 2000–2020, 47 exactly matched butterflies, 91 recorded-host-sharing species pairs and 990 recorded-nonsharing pairs. Site-year screening required ≥5 unique standard transect events in ≥3 31-day seasonal windows; sites required ≥3 consecutive eligible years. This yields 9,178 species×site year-to-next-year transitions with both sharing and nonsharing peers, from 1,718 species×site groups. Survey count is divided by source duration; it is a **survey effort-indexed adult encounter rate**, not a census of all adults or larval recruitment.

**Critical ascertainment:** 7,501 site-year×butterfly panel cells were assigned zero from absence of a reported count for a butterfly previously observed at that same monitored site, whereas 8,953 had an original source record. These conditional zeros form ~45.6% of annual species×site panel cells and are **not independently validated absences**. They require explicit sensitivity analysis, not silent certainty.

Michigan (43 qualifying monitored sites) and Colorado (34) did not meet the **≥50 independent site** source gate, despite 1,610/1,619 candidate species×site transitions. Their models were **NOT estimated** under the original locked rule, so this dataset does not yet provide genuine independent programme replication.

## Original adjusted adult time-lag association (NOT a competition effect)

Observational response: `log(1 + focal adult count/recorded survey hour at t+1) - log(1 + count/hour at t)`.

Predictors: lagged focal own encounter rate, mean lagged encounter rate of other **host-sharing** butterflies at the same site, mean of **non-host-sharing** butterflies, change in total survey duration, aggregate same-site year mean adult abundance proxy; butterfly species×site and calendar-year fixed effects. Regressions use 121-site cluster sandwich intervals.

| Illinois original pilot | coefficient β | site-cluster normal 95% interval (when available) |
|---|---:|---|
| Focal own previous count | **−1.0048** | −1.072 to −0.937 (mathematically coupled to change-score response; regression-to-mean caution) |
| Other butterflies sharing exact known larval host | **+0.02194** (SE 0.01717) | −0.0117 to +0.0556 |
| Other butterflies not sharing recorded host | **+0.30545** (SE 0.12678) | +0.0570 to +0.5539 |
| **shared minus nonsharing association** | **−0.28351** (SE 0.11719) | **−0.5132 to −0.0538** |

The shared coefficient is **not negative**. This is a stronger covariate-adjusted *positive association for the nonsharing comparator*, not a measured negative density response to competitors. The disparity could be caused by species traits, shared responses to abiotic conditions, compositional adjustment, different adult flight phenologies or source-zero imputation.

**Biological decision:** neither adult population *increase attributable to introduced host plants* nor larval *interspecific competition* is supported by this analysis. The possible contrast is a new original-data **hypothesis generator**, not established ecological mechanism.

## Required frozen follow-ups
- Refit without the site-year community mean covariate (mechanically related to both peer abundances).
- Restrict focal response to butterfly species with positively observed counts in **both** successive source years; test whether zero-completion dominates, and label this survivor-only sensitivity.
- Use 199 source-locked butterfly/plant-family **bipartite degree-preserving** host-link null graphs; current observational comparison cannot be considered more than a source-conditioned network relation until its original β is evaluated against a graph-constrained reference.
- Independent site-specific habitat, weather, plant biomass and caterpillar observations remain absent; future manipulative study would require larval density×host leaf biomass×light microhabitat factorial treatments to identify biological competition.

See [source-verified aggregate JSON](../../provenance/reviewer_defenses/results/sbus_illinois_host_niche_abundance_pilot_20261025.json) and [fixed ecological model specification](../../provenance/reviewer_defenses/protocols/SBUS_ILLINOIS_HOST_SHARING_DYNAMIC_PILOT_20261025.md). Main GEB `main`, PR #38 and the submission draft PR #70 remain unchanged; this is an independent exploratory source branch.
