# How much climate alignment survives finer plant-climate and geographic matching?

**Analysis date: 2026-10-19. Status: original WCVP/HOSTS/BCE/CHELSA source verified; 3×3 and four-BIO-variable nulls QUALIFIED; local 1500km full matched-panel inference requires 399-draw extension.** This source study concerns **potential larval host-resource geography**, not measured adult butterflies, recruitment or local ecology.

## Scope

Starting with 239 butterflies, the original HOSTS accepted-butterfly→host plant graph is extended by **1,027 additional evidence-filtered accepted-host links in 81 butterfly species** from a European-biased partial BCE/Clarke list (679 added plant IDs, families resolved using 323 original HOSTS accepted IDs + 304 original single-family genera + 52 pinned GBIF botanical families). This is not an independent worldwide second sample. Plant WCVP range/family, native host exclusions, CHELSA climate and TDWG Level3 geometry come from an exact-SHA author-only source archive, and every source/MCMC result is regenerated from the source in GitHub Actions.

Observed augmented resource climate alignment: **9,949/14,961 = 66.4996%** of added butterfly×WGSRPD3 resource cells lie within the native **host-plant** BIO1/BIO12 5–95% climate rectangle; **7,370/14,961 = 49.2614%** lie in the 10–90% rectangle. Exactly 235 butterflies qualify. The same-original-eligibility 234 butterfly sensitivity yields **9,903/14,903 = 66.4497%** and **7,331/14,903 = 49.1914%**, excluding only newly climate-eligible *Satyrium pruni*. Two extra geographic cells (`ANT`, `HMD`) for *Coenonympha oedippus* lack regional climate and are excluded from climate denominators only.

All nulls preserve **fixed butterfly–host identities**, exact native plant ranges, plant family, each host's introduced number of regions in its climate/continent stratum, and each region's incoming introduced plant-family count. Plant redistribution changes the union count of added butterfly resource cells; therefore all nulls compare **aligned added / all newly added resource cells in each individual simulated graph**, not absolute aligned-cell counts.

## Independent qualified MCMC results (two seeds ×199 draws per mode)

| Source-conditioned null / complexity | 235 eligible: residual climate 5–95% | 235 eligible: residual climate 10–90% | Same original-eligibility 234: 5–95% | Same 234: 10–90% | Strict per-response quality gate |
| --- | ---: | ---: | ---: | ---: | --- |
| Previous plant family×continent×BIO1/BIO12 **2×2** | +0.972–0.977 pp | +1.571–1.599 pp | +0.982–0.993 pp | +1.566–1.592 pp | QUALIFIED ([39221327128](https://github.com/zuizui0223/chocho/actions/runs/39221327128)) |
| **Temperature/precipitation 3×3** | **+0.387–0.415 pp** | **+0.808–0.819 pp** | **+0.389–0.419 pp** | **+0.809–0.816 pp** | QUALIFIED, original swap viability 88.1% |
| **BIO1/BIO7/BIO12/BIO15 2×2×2×2** | **+0.560–0.572 pp** | **+1.029–1.038 pp** | **+0.576–0.586 pp** | **+1.028–1.043 pp** | QUALIFIED, original swap viability 81.9% |
| 1,500km geographic original-footprint tether + prior BIO1/BIO12 2×2 | +0.340–0.349 pp | +0.583–0.584 pp | **DO NOT PROMOTE**: ESS29.988 in one of two chains (numerical median +0.348–0.356 pp) | +0.580–0.588 pp | FULL 4-response gate **NOT EVALUABLE** at 199 draws; longer unchanged 399-draw check is separately scheduled |

These are **source- and nuisance-model-conditional exploratory magnitudes**, not a causal decomposition or tests of historical butterfly range expansions. The two climate models are **not a nested sequence**: BIO1/BIO12 3×3 partitions differ from the four-axis 2×2×2×2 strata, so the difference is a comparison of alternative controls, not an additive temperature/precipitation/seasonality variance allocation.

- All four **3×3** responses pass ESS≥30 in both chains, median concordance, ≥10% burn movement, ≥100 accepted swaps per sampling interval and viable-edge fraction≥75%; broad-climate ESS approximately 199 / 169, narrow ESS 123 / 199.
- All four **four-bioclimate** responses pass the same gates, with broad ESS approximately 194 / 164, narrow ESS 199 / 127.
- For the **spatially tethered model**, the original geometry is SHA-verified (`c172bcf6...`), 365 focal region polygon representative points are mapped and 1,500km neighbour distances control symmetric swaps. Full augmented **235** responses both pass ESS≥30, and matched 234 narrow climate also passes; the matched 234 broad climate ESS=29.9885423 is narrowly below the preset ≥30, so the *complete* two-panel spatial model is held as **NOT EVALUABLE** at 199 draws, not rounded up. See a transparent after-diagnostic extension to 399 draws in `docs/exploratory/BCE_FINER_CLIMATE_SPATIAL_FILTER_PROTOCOL_20261019.md`.

Source-verified run for all three [GitHub Actions 39227164974](https://github.com/zuizui0223/chocho/actions/runs/39227164974); independent derived JSON artifacts reside under that run. Numeric MCMC fractions/ESS can be reproduced without redistributing copyrighted WCVP/HOSTS/TDWG sources. Serial MCMC add-one right-tail fractions (where qualified, mostly the 0.005 minimum) **are not independently randomized-permutation p-values**.

## Ecological interpretation and stopping decision

1. A **large original climate bias** (+8.28 pp against a simple within-continent geographic expectation on the original HOSTS source) contracts to approximately **+1 pp** after preserving host-plant network and coarse plant climatic distributions, and to **~+0.4 pp** under a finer temperature/precipitation 3×3 climate-stratified plant network. This is strong evidence that much of the headline potential resource climate affinity is **inherited from plant-side geography**.
2. Positive conditional residual persists at 5–95% and 10–90% under two alternative stricter climatic nuisance models and with exactly matched original 234 butterfly eligibility. Its *magnitude*, however, is **sensitive to null-model climate resolution**; one should not elevate a sub-percentage-point residual to a generalized butterfly-specific evolutionary host preference.
3. Spatially local null medians suggest an even smaller residual, but the all-four-endpoint 199-draw qualification failed on one **near-threshold ESS**. No complete spatial generality claim is allowed before the independently tracked 399-draw unchanged-null run qualifies.
4. Host source remains Europe-biased under augmentation; recorded host-plant-native climate may not correspond to the insect's physiological climate niche, and large WGSRPD3 regions are not area- or abundance-standardised. There is no observed butterfly occurrence/establishment response in these outputs.

**Conclusion:** The more defensible ecological headline is the **dominant role of the redistributed plant climatic geography** in shaping spatial opportunity for specialist/generalist butterflies, with a small network-conditional residual whose size is highly dependent on what plant spatial/climatic constraints are conserved. This is NOT evidence that humans made butterflies reproduce more, shifted their ranges in historical time, or benefited endangered insects. Keep the frozen GEB PR #38 unchanged; these macro ecological experiments stay in exploratory [PR #46](https://github.com/zuizui0223/chocho/pull/46).
