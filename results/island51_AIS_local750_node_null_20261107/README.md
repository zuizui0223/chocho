# Locally constrained AIS circuit-vertex placebo on the original 51 islands

Exploratory analysis completed 2026-11-07, stacked on PR #136. It retains all original butterfly/host/circuit source data and applies a **strict geography-only 750-km maximum source-node relocation** to whole-graph AIS surrogate assignments.

**Main result:** AIS retains a modest positive **absolute held-pair prediction increment** for later butterflies (+1.0020% MSE), but **the original placement of AIS vertices is not unusually good relative to 199 whole-graph swaps only among original island pairs at most 750 km apart**. Actual AIS alone is beaten/equaled by 100 of 199 geographically local surrogates (descriptive upper-tail fraction 101/200 = **0.505**). This distinguishes predictive information from uniquely verified bilateral shipping routes.

## Later butterfly first-indexed community similarity, 2011–2022

| AIS model contrast | Evaluation | Original AIS MSE improvement | Median of 199 local swaps | Local swap 2.5–97.5% gain range | Local swaps at least as good as actual | Descriptive tail rank |
|---|---|---:|---:|---:|---:|---:|
| AIS alone All 1,275 +1.0020% +1.0600% -0.9340% to +3.9410% 100/199 0.505 |
| AIS alone Within 235 +0.3812% +1.2268% -0.9276% to +7.6281% 122/199 0.615 |
| AIS alone Between 1,040 +1.3347% +0.8056% -3.2769% to +4.3833% 75/199 0.380 |
| AIS after wind All 1,275 +0.5446% +0.9826% -1.2169% to +3.9605% 128/199 0.645 |
| AIS after wind Within 235 +0.6460% +1.1309% -1.0523% to +7.9872% 107/199 0.540 |
| AIS after wind Between 1,040 +0.4936% +0.5745% -2.9324% to +4.4764% 106/199 0.535 |

The later AIS-alone between-longitude-block result remains a positive +1.3347% MSE improvement on the fixed 1,040 source pairs, but **75/199 local swaps** do at least as well, compared with only **4/199** in the earlier broad latitude-only swaps. AIS after original wind is positive +0.5446% on all pairs, but **128/199** local swaps equal or exceed it. The observed AIS graph is therefore not distinctively superior under the tighter local source-geometry control.

**The changed placebo changes the interpretation, not the original positive score.** In PR136, the source AIS observed gain +1.0020% was contrasted against latitude-only placebo median −0.5042%, with 5/199 placebos at least as good. Here the same source AIS improves +1.0020% against local 750-km placebo median +1.0600% and 100/199 placebos at least as good. The smaller, more geographically plausible perturbation changes only AIS labels of 18 distance-matched pairs (36 of 51 islands), leaving 15 node labels fixed. Those unchanged labels and the limited 750-km exchange can leave much original geography aligned by construction, so this result does **not** prove the AIS information is entirely reducible to distance. It does show that the evidence for exact island-level AIS source specificity is not robust to this null specification.

The MAE result retains a weaker positive all-pair actual AIS-alone gain +0.8485%, compared with local surrogate median +0.5417%; 63/199 local surrogates have equal/better MAE gain (rank 0.320). For within-longitude-block pairs the observed AIS-alone MAE is adverse (−0.2060%) and 0.830 of the null-tail rank; between-block observed AIS-alone MAE is +1.2258%, local median +0.5963%, descriptive rank 0.225. Scores are descriptive, not independent-island p-values.

## Earlier butterfly first-indexed community similarity (2001–2010)

| AIS model contrast | Evaluation | Original AIS MSE improvement | Median of 199 local swaps | Local swap 2.5–97.5% gain range | Local swaps at least as good as actual | Descriptive tail rank |
|---|---:|---:|---:|---:|---:|---:|
| AIS alone All 1,275 -8.1571% -1.4864% -14.9214% to +0.5775% 176/199 0.885 |
| AIS alone Within 235 -1.1376% -1.0266% -4.3399% to +0.3240% 110/199 0.555 |
| AIS alone Between 1,040 -11.3743% -1.7266% -20.7103% to +1.1681% 176/199 0.885 |
| AIS after wind All 1,275 -5.9276% -1.7296% -11.6780% to +0.7165% 173/199 0.870 |
| AIS after wind Within 235 -1.3188% -1.1393% -4.2083% to +0.4826% 110/199 0.555 |
| AIS after wind Between 1,040 -8.0676% -2.1623% -16.2170% to +1.6332% 172/199 0.865 |

The earlier real-AIS source remains harmful: AIS alone changes all-pair heldout MSE **-8.1571%**, compared with local null median -1.4864%, with **176 of 199** geographically local surrogates at least as good. This should be reported alongside later positive scores. The source AIS activity layer is modern and cannot identify 2001–2010 vessels or colonization events.

## Fixed design and source-only geography

The original 51 source island groups, 194 original WPI botanical ports, 182 qualified marine terminals, all 1,275 unordered island pairs and 225-butterfly candidate list, pooled 2001–2024 candidate-host Jaccard, marginal host-richness pair features, original pre-2011 all-botanical observation-effort pair terms, geographical distance, port-count/latitude controls, uniform-ocean resistance and original mean-wind resistance are retained unchanged. The 15 original longitude-block double-endpoint-unseen source splits are unchanged, with 235 within- and 1,040 between-block island pairs tested exactly once per model/map. The later target also contains earlier butterfly similarity as a predictor, and the earlier target does not use its own outcome as a predictor.

The [source-only local pairing table](local_swap_pairs.csv) was constructed with a 6371.0088-km Earth-radius haversine using the same source original WPI representative centroids, restricted to original island pairs in the same signed hemisphere with distance at most 750 km. Eligible pairs are sorted by increasing original great-circle distance and numerical IDs, then disjointly matched by deterministic greedy selection. Exactly **18 nonoverlapping source pairs (36 original islands)** result, leaving 15 original island vertices unpaired and immutable. Maximum chosen pair separation is about 730.5 km. None of the matching uses the butterfly response.

Using NumPy PCG64 seed 20261107, each null globally swaps the two AIS node labels in a randomly selected subset of the original fixed 18 pairs, rejecting only repeated or all-unchanged configurations. There are exactly **199 unique nonidentity** maps and a saved real original identity map 0, reused unchanged for both butterfly targets and both fixed comparisons. Each map permutes the complete **symmetric AIS circuit graph**, not pairwise butterfly similarities, sampled pairs or other source covariates. The map's geographical swap cap does NOT mean original AIS edge routes were constrained to 750 km; actual whole-graph source connections can extend far beyond that.

The two fixed models assess AIS alone beyond original candidate-host similarity, marginal richness, older botanical observation effort, geographical distance, original WPI port counts, four source latitude/hemisphere terms and uniform ocean resistance, and AIS conditional on additionally keeping the original mean-wind-speed circuit. Train-only standardized OLS with intercept and clipped [0,1] predictions was refitted under each fixed original both-islands-unseen split. The observed source graph must reproduce all 12 earlier PR135 full-training MSE+MAE source comparisons before accepting the null.

## Independent audit and complete outputs

The [frozen before-fit protocol](PROTOCOL.md), [source-locked evaluation runner](../../scripts/run_island51_AIS_local750_node_null.py), [separate independent Python-standard-library auditor](verify_AIS_local_null.py), [published comparison table](observed_vs_latitude_node_null.csv), [source pair distances](local_swap_pairs.csv) and [GitHub Actions complete output](https://github.com/zuizui0223/chocho/actions/runs/42338910247) provide reproducibility. Scientific CI artifacts have 30-day retention and can be regenerated by rerunning the source-pinned workflow.

The independent auditor reconstructed **18 geographical source pairings, 15 immutable island vertices, 10,200 source-to-donor labels, 1,020,000 augmented AIS model predictions and 5,100 baseline forecasts**, all 2,400 comparison/scope MSE+MAE scores, and all 12 empirical rank statistics. Original source butterfly target values, 15 geographic heldout training populations and no test-island endpoint leakage were separately checked. The real graph reproduces previous PR135 scores to maximum absolute difference 4.79×10^-13. Published-table consistency is checked by dedicated CI.

## Ecological limits and positioning

The correct supported claim is: original AIS-derived maritime connectivity contains modest *additional held-out prediction information*, but its very precise island-to-circuit mapping **cannot be distinguished from 199 local swaps under this fixed source-distance-capped placebo**. A rank of 101/200 cannot be interpreted as an equivalence test or proof of no particular AIS relationship, since 15 vertices are fixed, 36 move only locally and the null may have low power to generate meaningfully different predictors.

Neither this placebo nor the earlier latitude-only one is a formal causal or exchangeable permutation null for real shipping. AIS is activity-density-derived resistance, not actual bilateral voyages or cargo routes. The original source host Jaccard still pools plants indexed 2001–2024 and may postdate 2011–2022 butterfly observations; older all-plant counts are not historical host species identity. Butterfly first indexing is not verified arrival, reproduction or colonization. The island observations are correlated and reused across folds; do not assign independent-sample p-values. Preserve the prior global GEB +54.9% potential larval-host portfolio expansion, older island/port correlations, and earlier positive and negative observations under their original estimands.
