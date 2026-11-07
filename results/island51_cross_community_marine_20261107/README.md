# Historical port plants and butterfly communities: is maritime predictive information shared across taxa?

**Completed exploratory macroanalysis, 2026-11-07.** Based exclusively on the public original 194-WPI-port, 2001–2010 general *Plantae* speciesKey-by-port source snapshot archived in [PR138](https://github.com/zuizui0223/chocho/pull/138). **This analysis does not reconstruct the original restricted butterfly larval-host interaction data, actual ship itineraries, or true arrival/colonization dates.**

## Research question and main result

Can the same modern AIS-derived and mean-wind Circuitscape maritime features predict **historical (2001–2010) indexed general plant-community similarity, earlier (2001–2010) butterfly-community similarity, and later (2011–2022) butterfly-community similarity** when they are evaluated on the *same original island pairs*, with every pair's two island endpoints excluded from its model training folds?

**Result: a small AIS prediction increment is observed for later butterfly communities, but is not a generic property of the earlier indexed general flora or earlier butterfly communities.** This is observational and temporally asymmetric predictive association, **NOT proof of an historical change in marine traffic or causal dispersal**: the AIS activity source represents 2015–2021, which postdates the early plant/butterfly records.

## Exact identical botanical-indexed source cohort

- Original study: **51 island/landmass groups, 194 botanical WPI source ports, 182 marine-qualified ports**, 1,275 unique unordered original butterfly island pairs and the complete initial candidate-butterfly panel retained unchanged.
- Independently indexed Plantae in 2001–2010: **4,876 unique GBIF speciesKeys across 27,598 port×plant speciesKey rows**, recorded at least once in **45/51 islands**. The six original islands with zero indexed early GBIF plant species are **NOT ecologically plant-free** and were not excluded from the original PR138 analysis.
- This new direct historical-general-flora **response** is well-defined without zero-union imputation on the fixed **45 choose 2 = 990 source island pairs** where both islands have at least one indexed early plant speciesKey. Their source population was fixed by indexed botanical coverage, without looking at either butterfly outcome.
- **741 of the 990** have ≥10 originally indexed early general Plantae speciesKeys at **each** island and form a **source-only quality sensitivity** (no new model fitting on this subgroup). The 990-pair source population contains **188 within-longitude-block and 802 between-longitude-block pairs**, compared with 235/1,040 in the old unfiltered 1,275-pair cohort. The source-only selection changes population size and therefore **old PR138 MSE percentages are not directly interchangeable**.
- Indexed plant similarity is sparse: **440/990 pairs** share at least one recorded historical plant speciesKey, while **550/990 have Jaccard zero** in the retrieved GBIF record (not ecological absence). Median original indexed general plant speciesKey Jaccard among 990 eligible pairs is 0.

## Main original 990-pair double-unseen MSE result

| Community response/model panel | AIS added to the panel baseline (MSE reduction) | Wind only (MSE reduction) | AIS and wind added jointly (MSE reduction) | AIS added after wind (MSE reduction) |
|---|---:|---:|---:|---:|
| **2001–2010 indexed general port-plant community Jaccard** | **−1.2943%** | −0.9237% | −1.7152% | −0.7842% |
| **2001–2010 indexed butterfly community Jaccard** | **−8.5406%** | −2.8357% | −3.8859% | −1.0213% |
| **2011–2022 first-indexed butterfly community Jaccard** | **+0.8308%** | −4.9564% | −3.9745% | +0.9355% |
| **2011–2022 butterflies, also adjusting retrospective pooled 2001–2024 candidate-HOSTS Jaccard and host richness** | **+0.6776%** | −5.4258% | −4.6986% | +0.6898% |

Positive MSE reduction means better held-out pairwise prediction. Each row's percentage has a separate reference-model denominator; **do not add individual increments** or treat differences among panels as causal treatment effects. The retrospective host panel is labeled separately because the host Jaccard includes records across 2001–2024 and can *postdate* the later butterflies.

The AIS-only addition also improved later butterfly **MAE** by **+0.3109%** without retrospective host controls and **+0.1760%** with them. For indexed earlier general flora and earlier butterflies, AIS-only MAE changes were **−0.5448% and −3.5405%** respectively.

### Sensitivity: source geography and plant-indexing quality

**Between original longitude blocks (802 identical recorded-plant-eligible pairs):** adding AIS alone changes MSE **−2.6989%** for earlier general flora, **−11.6206%** for earlier butterflies, **+1.1942%** for later butterflies and **+1.1453%** for the later retrospective-HOSTS panel. These 802 paired observations are neither independent across island endpoints nor a new set of external regions.

**Both endpoints with at least 10 indexed general Plantae speciesKeys (741 pairs):** adding AIS alone changes MSE **−1.4991%** for early general flora, **−9.1118%** for early butterflies, **+0.5642%** for later butterflies and **+0.1837%** for the later retrospective-HOSTS panel. This additional source-only observation-quality restriction retains the same sign but attenuates the later butterfly increment. Do not claim statistical equivalence or high confidence from just these point estimates.

**Absolute accuracy matters:** on the 990-pair later butterfly source panel, the geography/recording-effort-only baseline MSE is **0.006644710365** and AIS-augmented MSE **0.006589506384**. On the recorded early general-plant Jaccard response, baseline MSE is **0.001713550320**, versus AIS-augmented **0.001735728843**; this ratio is negative. These outcomes have *different scales, recording opportunity and marginal distributions*, so predicting butterfly community similarity is not directly comparable to predicting historical general-plant similarity as a biological effect.

## Exact predictive design

The baseline for every panel includes exactly the original logged port great-circle distance and square, sum/absolute difference of log1p original source WPI port counts, fixed uniform ocean resistance, four previously frozen representative-latitude/hemisphere pair controls, two archived all-general-plant 2001–2010 WPI 5-km buffer **record-count** effort terms, two original **live re-queried** 2001–2010 GBIF 5-km buffer record-count terms, and two original historical general-plant indexed **speciesKey richness** pair terms. These are intended to reduce sampling and marginal opportunity confounding, **not prove detection-adjusted abundance or actual ecological absence**.

Earlier indexed butterfly Jaccard enters the **later butterfly** baseline, but never enters the earlier butterfly or earlier plant outcome baseline. The additional retrospective-HOSTS panel uses original pooled 2001–2024 candidate-host Jaccard and original marginal host richness, only for the later butterfly target; its historical invalidity is explicit. Historical plant Jaccard is the **response**, never a predictor in its own target model. The four models are BASE (without AIS or wind), AIS, WIND and BOTH (AIS + wind), with identical train-only standardized OLS, intercept and [0,1]-clipped predictions.

All models use the exact same **15 longitude-block-pair both-endpoints-unseen folds** from the earlier PR132–138 source: every test-pair source island is excluded from that fold's training, and training/testing restricted to the same 990 source-only indexed-plant-eligible pairs. Five prespecified incremental contrasts per panel: BASE→AIS, BASE→WIND, AIS→BOTH, WIND→BOTH and BASE→BOTH. Four source-only reporting scopes: all 990, within longitude, between longitude, and both endpoint islands with ≥10 indexed speciesKeys. No p-values or naïve confidence intervals are formed from the dependent island pairs/folds.

## Reproducibility, original-source checks and retained negative findings

The [frozen protocol](PROTOCOL.md) was committed BEFORE this cross-community fit, but the prior AIS and historical general-flora outcomes from PR138 were already known, so this is exploratory robustness and **not an independent preregistered study**.

- [Complete 990-pair and original 1,275-pair community source ledger](source_pair_community_panel.csv), retains every original source-pair index, two butterfly community responses, old plant speciesKey Jaccard, original 990 eligibility and 741 high-coverage flags, distances and longitude class.
- [All 80 model × community × geographic/quality comparisons](comparison_scope_summary.csv), full MSE, MAE and pairwise gain/worse counts.
- [All 15,840 individually held-out prediction records](heldout_predictions.csv.gz), [240 model fold errors](model_fold_losses.csv) and [300 comparison-fold records](comparison_fold_losses.csv).
- [Original 15 source-only fixed train/test holdout populations](source_eligible_geo_folds.csv), [execution source receipt](execution_receipt.json) and [standalone independent auditor](verify_cross_community.py).

A separately written Python-standard-library auditor (without importing fitting code) independently reconstructed the original **194 WPI ports, 27,598 port×plant speciesKeys, original 1,275 butterfly source pairs, 45 indexed-plant islands, exactly 990 matched source pairs, 741 quality sensitivity pairs, all 15 geographic double-unseen holdouts, 15,840 individual model predictions, 240 model-fold errors, 300 fold contrasts and 80 reported comparisons**. All original input tables, historical GBIF source snapshot and original PR138 plant/butterfly records remain unchanged. The source runner and entire artifact are reproducible without an internet or GBIF request.

## What this does and does not mean

**Supported:** on these 990 historically plant-indexed original island pairs, AIS-derived marine resistance contributes small positive prediction information to **later indexed butterfly communities** but not to **earlier general port plant community or earlier butterfly Jaccards** under the fixed baseline. Wind-speed resistance alone does not improve all-pair MSE in any panel. These sign differences survive the source-only ≥10 indexed plant speciesKey quality sensitivity. This is useful to frame an ecological *timing/community-specificity* hypothesis, not proof of transport mechanisms.

**Not supported:** AIS directly carried butterfly eggs/larvae, plant seeds or vegetation; AIS represented historical 2001–2010 shipping; the ship network preceded the 2001–2010 plant observations; plant Jaccard captures larval host plants, historical alien introductions or resident whole-island vegetation; butterfly first GBIF indexing equals colonization, feeding, survival or reproduction; maritime connectivity caused butterfly communities to expand. AIS represents modern activity-density-derived **undirected** circuit resistance, not actual port-to-port vessel itineraries. Both taxa have unequal GBIF index quality. The source cohorts and randomization tests are observational and geographically dependent.

The previously established GEB global potential larval-host resource opportunity increase of +54.9%, original 51-island later-butterfly AIS effects, and PR137 geography-constrained placebo (actual AIS is not distinguishable from nearby source-node swaps) remain valid with their **separate original comparison populations and inferential limits**. This result cannot be used to reverse the negative PR137 source-specificity diagnosis.
