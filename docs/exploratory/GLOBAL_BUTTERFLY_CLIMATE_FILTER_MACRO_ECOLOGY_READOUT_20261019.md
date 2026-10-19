# Host redistribution follows native-host climate analogues, but resource opportunity is uneven among butterflies

**Date:** 2026-10-19. **Status:** SOURCE-VERIFIED EXPLORATORY GLOBAL MACROECOLOGY; no measured butterfly recruitment or causal effect of plant introductions. **Do not merge into the frozen GEB v0.2 review manuscript PR #38 or main.** This is a forward macroecological analysis of the original HOSTS/WCVP/CHELSA and LepTraits-derived species panels, *not* a new extraction of individual published insect experiment results.

## Main nontrivial ecological question

Once already-recorded larval host plants are mapped beyond their native geographic ranges, do the new potential host-resource opportunities follow the **climate of a butterfly's established native host pool** more than the area/region composition of introductions would suggest? Does the share of climate-analog opportunity filled differ according to a butterfly's documented disturbance affiliation?

Simply adding introduced host range cells to native ranges monotonically increases the **set union** of potential resources. Neither **the climate of the additional cells** nor **the trait dependence of capacity filling** follows mechanically from that union: these are the empirically evaluated, potentially falsifiable ecological patterns.

## A. Independent source-verified climate macro result

Original source archive: the **author-only** original HOSTS/WCVP/CHELSA source reproducibility artefact, [workflow 39006554609](https://github.com/zuizui0223/chocho/actions/runs/39006554609), whose inner ZIP SHA256 and constituent source checksums were checked by the new analysis. Original 239 butterfly species, 365 WGSRPD3 units, **26,530 native butterfly×region cells**, **41,083 contemporary**, **14,553 introduced-added**. For native-host climate-envelope estimation, 234 species had at least 10 native regions with climate; the missing five were **excluded, not recoded as climate failures**, leaving 14,481 introduced-added cells (204 species have nonzero additions).

Within the 5th–95th central intervals of **CHELSA BIO1 mean temperature AND BIO12 annual precipitation across each butterfly's NATIVE-HOST regions** (a host-region climate proxy, NOT demonstrated insect climate tolerance), 9,451 of 14,481 added cells (**65.26%**) were climate analogues. If each butterfly's actual added cell count *in each WGSRPD Level1 area* were allocated uniformly across that area's currently native-unsupported WGSRPD3 units, the analytical expectation is 8,251.9 (**56.98%**), an observed excess of **+1,199.1 cells / +8.280 pp**.

Prespecified 10th–90th interval sensitivity: **6,948 of 14,481 (47.98%)** versus Level1 composition expectation **5,834.7 (40.29%)**, **+7.688 pp**. Genus bootstrap for added-cell-weighted excess, 5–95 interval: **+7.278 to +9.309 pp**; this captures genus-level resampling uncertainty, but **cannot account for resource cells shared through the same introduced plants or spatial autocorrelation**.

The sign is broad within the *studied resource panel*: **188 of 204** butterflies with added cells show positive 5–95 Level1-matched excess (15 negative, one exactly zero). Among species with at least 10 added regions, 175 of 187 are positive and added-cell-weighted excess is **+8.249 pp**. Removing the ten butterflies with the largest added counts leaves **+8.906 pp**; leave-one-butterfly-genus-out weighted excess ranges **+8.149 to +8.499 pp**. Species are **not** independent replications because they reuse host plants and geographic cells.

[Primary and climatic-threshold successful run 39209404496](https://github.com/zuizui0223/chocho/actions/runs/39209404496); [extended source-pinned same-data run 39211273830](https://github.com/zuizui0223/chocho/actions/runs/39211273830).

## B. Stronger, margin-preserving structural null — qualified, but still not host-network matched

A **posthoc** doubly constrained null now preserves **(1) the number of introduced-added regions of each butterfly in each WGSRPD Level1 area, (2) the total number of newly supported butterflies in each WGSRPD3 region, (3) every native-host-supported butterfly×region exclusion**, while performing continent-stratified double-edge swaps on the added butterfly×region incidence matrix. This does NOT preserve underlying individual introduced plant IDs, plant niche conservatism, shared host links or actual butterfly occupancy. It is a *more demanding geometric comparator*, not a plant-network-controlled causal experiment.

The initial 2×99-draw diagnostic was rejected as **NOT EVALUABLE** because effective sample size was below 30. Without changing either ecological response, observational data or model margins, the chain was extended to **two separate 199-draw chains with larger inter-draw spacing**. The new [source-exact CI 39211273830](https://github.com/zuizui0223/chocho/actions/runs/39211273830) passed its prespecified minimum burn movement (25.8–26.0% of originally added cells moved), spacing acceptances (>3,080 accepted), chain-median agreement and **ESS gates**:

| Native-host analogue | Observed matching added cells | Fixed-margin null medians (two independent chains) | Observed-minus-null median | ESS of two chains |
| --- | ---: | ---: | ---: | ---: |
| BIO1+BIO12 5–95% | **9,451** | 8,910 / 8,908 | **+541/+543 cells = +3.736/+3.750 pp** | 58.4 / 95.7 |
| BIO1+BIO12 10–90% | **6,948** | 6,407 / 6,408 | **+541/+540 cells = +3.736/+3.729 pp** | 119.8 / 51.0 |

The observed count exceeded all sampled 199 null draws in each chain (add-one, *serial* MCMC right-tail **0.005 per chain** for 199 draws). This is **NOT an exact, independent-permutation p value**, especially since the null was added after the weaker geographic reference's positive result was known. Do not call it a pre-registered discovery or causal host effect.

The stricter null accounts for much of the original +8.28 pp excess: compared with region-uniform expectations, the fixed region-degree/consumer-by-continent margins explain a further ~4.5 pp, **leaving ~3.74 pp source-conditional positive residual**. This distinction is essential for the paper's ecological credibility.

## C. Cross-biogeographic structure (recipient Level1, not butterfly populations)

The continent-uniform expectation from section A gives the following *observational source-panel* contrasts:

| WGSRPD Level1 receiving region | Added potential resource units | Analogue observed % | Level1 expectation % | Excess (percentage points) |
| --- | ---: | ---: | ---: | ---: |
| Europe | 1,044 | 65.42 | 57.90 | **+7.52** |
| Africa | 2,540 | 66.46 | 54.30 | **+12.16** |
| Temperate Asia | 1,258 | 38.31 | 29.94 | **+8.38** |
| Tropical Asia | 1,118 | 36.31 | 35.44 | **+0.88** |
| Australasia | 1,167 | 86.12 | 81.37 | **+4.75** |
| Pacific | 1,596 | 66.54 | 59.70 | **+6.84** |
| Northern America | 3,144 | 77.29 | 67.16 | **+10.13** |
| Southern America | 2,179 | 68.24 | 58.31 | **+9.93** |
| Antarctic/outer subantarctic region as coded in WGSRPD | 435 | 47.82 | 48.53 | **−0.71** |

These are pooled butterfly×region *potential plant-resource* counts, not area-weighted climate niches, regional abundance, butterflies surviving in Antarctica, or 9 independent test replicates. The dramatic recipient-area heterogeneity is an exploratory descriptive result, not a fitted causal continental climate effect. **Tropical Asia's 0.88 pp and the Antarctic −0.71 pp specifically argue against interpreting the geographic signal as universal.**

## D. Disturbance affinity: climate-analog capacity filling versus climate alignment surplus

From the original 339-species LepTraits-derived descriptors, the frozen HOSTS/WCVP 239-species panel retains **150/193 strong disturbance-associated** taxa but only **33/68 strong disturbance-avoidant** taxa. This is serious resource-data eligibility bias and over-representation of habitat generalists. **No balanced worldwide inference is justified.**

With native host-resource breadth (nonlinear), HOSTS species/family richness, wing size, voltinism, butterfly-family fixed effects, HOSTS recording adequacy, and **native-resource dominant WGSRPD Level1 area, native-host regional mean BIO1 and BIO12 climate** as covariates:

- Outcome = species-wise continuity-corrected logit of the fraction of *climatically analogous native-unsupported resource regions* newly filled by introduced host distributions, `log[(O+0.5)/(K−O+0.5)]`, K≥5. **179 butterflies** (146 strong associated, 33 strong avoidant), genus-cluster coefficient **+1.379**, approximate 95% **[+0.630,+2.129]**. This is an observational trait association, NOT host-induced consumer adaptation or population expansion.
- HOSTS-inventory adequacy restricted to 160 species: coefficient **+1.080**, 95% **[+0.423,+1.738]**. Adding number of analogue candidate regions K as an additional covariate: **+1.471**, 95% **[+0.740,+2.202]**.
- Tighter 10–90 climate criterion: geography-adjusted coefficient **+1.285**, 95% **[+0.479,+2.091]**. Held-out genus (fivefold repeated ten times) mean squared-error improvement from including the habitat label beyond covariates: **6.20%** at the primary climatic interval (5.09% tighter).
- **Negative control on a different estimand**: limiting analysis to 157 butterflies that obtained added cells, the difference between disturbance groups in species-wise `(O−E)/A` — *surplus climatic sorting among the added resources* — is **−0.0168**, 95% **[−0.0914,+0.0578]**, i.e. **no robust positive differentiation**. More climate-analog area filled does **not** imply more special climatic targeting **conditional on gaining a resource region**.

Positivity violation: among butterflies with strong habitat classes and climate support, native-resource-dominant Africa (27 species), Southern America (26), Tropical Asia (7) contain **zero strong disturbance-avoidant species**; comparative data with both groups exist primarily in Europe, Northern America and Temperate Asia. Including continent fixed effects does not manufacture missing treatment overlap. This is a central generalizability limit, not a minor caveat.

## E. What the global ecological story can and cannot be

**Defensible ecological claim, sample-conditional:** Redistribution of fixed known larval host plants is associated with an **ecologically structured supply of climate-analog potential resource regions**. Geographic/region-degree constraints explain some but not all of the climate alignment. Disturbance-affiliated butterflies, in this biased focal sample, potentially fill a larger fraction of their remaining host-resource space, yet do not show demonstrably greater *climatic targeting excess* once opportunities are conditioned on.

**Not established:** Actual butterflies colonised those places or bred on those plants; introduced hosts improved survival; local trophic interactions became monotonous; butterflies survived in high-latitude/Antarctic areas; species truly differ in climate adaptation; historical plant introductions led to butterfly population growth. Host catalogues (HOSTS and partially European BCE), climate envelopes, spatial resolutions and biological viability all limit ecological generality. The stricter butterfly×region null alone **does not preserve plant-host identities**, so we additionally completed two plant-family×continent and plant-family×continent×BIO1/BIO12 network-preserving nulls, both under the original and the partial BCE-augmented accepted-host inventory. The remaining climate-sorting effect is **small (~+1 pp with 5–95% and ~+1.4–1.6 pp with 10–90% native-host analogue widths), but remains positive and numerically MCMC-qualified** across both source catalogues. It is not evidence for butterfly demographic gains or fine-scale ecological interaction preference.

**Paper positioning:** Keep GEB PR #38 frozen as an independent resource-network paper. This is a promising **macroecological follow-up candidate**, but current ecological discovery is *potential host-resource niche geography*, not a realized consumer demographic law. Major next gates: (1) reclassify added-region identity under the alternate partial host catalogue (BCE) and re-test climate residuals; (2) source-preserving introduced-plant-identity null; (3) local occupancy validation using balanced global butterfly occurrences and climate-matched effort, if available. Do not make conservation claims about benefits of invasive host plants.

## Code and original evidence

- [Source-exact climate macro and climatology reconstruction](https://github.com/zuizui0223/chocho/blob/analysis/global-habitat-filter-20261019/scripts/analyze_butterfly_climate_analog_macro_20261019.py).
- [Dual-chain fixed-margin structural null](https://github.com/zuizui0223/chocho/blob/analysis/global-habitat-filter-20261019/scripts/analyze_climate_analog_fixedmargin_null_20261019.py).
- [Climatic biogeography, bias controls and holdout predictions](https://github.com/zuizui0223/chocho/blob/analysis/global-habitat-filter-20261019/scripts/analyze_butterfly_climate_analog_geography_sensitivity_20261019.py).
- [All three result JSON files, author-only derived (no redistribution of source plants), workflow run 39211273830](https://github.com/zuizui0223/chocho/actions/runs/39211273830).

**Result classification: Exploratory macroecological positive resource-structure pattern supported by multiple source-verified numerical contrasts and qualified posthoc geographic nulls, with habitat-based observational generalization explicitly restricted.** No novel field experiments or independent new butterfly observations were created.
