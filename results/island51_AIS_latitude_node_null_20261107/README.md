# Original AIS island-network alignment exceeds most latitude-band relabeling placebos

**Exploratory source-locked result, 2026-11-07.** The original AIS maritime-circuit assignment to island pairs contains extra predictive information for the later (2011–2022) butterfly community similarities relative to 199 same-hemisphere, similar-latitude graph-node relabelings. This is a quantitative predictive specificity check, not proof of actual shipping or ecological causality.

## Main results: later first-indexed butterfly communities, 2011–2022

| AIS predictor comparison | Evaluation | True AIS MSE reduction | Median of 199 vertex-reassigned AIS networks | Approximate 2.5–97.5% relabeling range | Relabelings at least as good as true AIS |
|---|---|---:|---:|---:|---:|
| AIS alone beyond host composition, candidate richness, original geography, latitude and older botanical effort | All 1,275 pairs | **+1.0020%** | **−0.5042%** | −2.8221% to +0.9756% | **5/199** |
| Same AIS-alone contrast | Between longitude blocks, 1,040 pairs | **+1.3347%** | −0.7078% | −3.7986% to +1.0011% | **4/199** |
| Same AIS-alone contrast | Within longitude blocks, 235 pairs | +0.3812% | −0.2667% | −1.9818% to +1.8640% | 50/199 |
| AIS added beyond the original wind circuit | All 1,275 pairs | +0.5446% | −0.8815% | −3.1254% to +0.9717% | 14/199 |
| Same wind-conditioned AIS contrast | Between longitude blocks, 1,040 pairs | +0.4936% | −1.2734% | −4.3616% to +1.2338% | 19/199 |
| Same wind-conditioned AIS contrast | Within longitude blocks, 235 pairs | +0.6460% | −0.2658% | −2.0310% to +1.8494% | 33/199 |

Positive MSE reduction means lower held-out prediction error. The observed AIS-alone full-sample improvement exceeds 194 of 199 finite latitude-band-reassigned circuits. Its descriptive upper-tail rank fraction is **6/200 = 0.030**. Between longitude blocks, the descriptive upper-tail rank fraction is **5/200 = 0.025**. This is NOT an independent-source or formally calibrated causal p-value: the null deliberately breaks the AIS graph's island labels and does not preserve the precise geographical, historical and ecological exchangeability conditions required for such an interpretation. The conditional AIS-after-wind result is positive but has a weaker rank (15/200 = 0.075).

The separate MAE check supports the later AIS-alone all-pair result directionally: true AIS improves MAE by **+0.8485%**, compared with a latitude-reassigned median of **−0.2562%**, with **13/199** shuffled source networks at least as good. In between-block pairs, true AIS improves MAE +1.2258% against null median −0.3350% (rank 0.035). In within-longitude pairs, AIS-alone MAE worsens −0.2060%, and the null rank is 0.545. All adverse metrics are retained in the complete numerical source table.

## Earlier butterfly communities, 2001–2010: essential unfavorable result

For the original earlier indexed butterfly communities, the actual AIS assignments are **worse** than the surrogate networks on the all-pair MSE criterion. AIS added without wind gives **−8.1571%** held-out MSE change compared to a median **−0.4382%** among the 199 relabelings; **199/199** surrogates have MSE gains at least as large as actual AIS. With original wind already present, real AIS addition yields **−5.9276%**, versus the null median **−0.8501%**; again **199/199** null networks are at least as good for that all-pair criterion. No negative evidence was removed or replaced.

These differences are NOT a formal test of historical maritime traffic, causal butterfly range expansion, or changing transport mechanisms. The source AIS layer is modern (2015–2021), whereas the earlier butterfly community target is 2001–2010. Butterfly first indexing is a GBIF observation event, not a verified colonization date.

## Exactly what was tested

This new branch stacks on the fully audited [PR135](https://github.com/zuizui0223/chocho/pull/135). The candidate-host Jaccard pooled over 2001–2024, original 51 landmass groups, 225 butterfly candidate species, 194 original botanical ports, 182 marine-qualified terminal ports, original distance/port-count/latitude/hemisphere controls, 2001–2010 **all-plant** recording-effort proxies, uniform resistance, native AIS and wind Circuitscape resistance matrices, and both butterfly response tables are unchanged. All **1,275 original pairs** are predicted once for each model/null, with **both island endpoints absent from training** in each of the same 15 fixed longitude-block-pair folds. Of the original pairs, 235 are within and 1,040 between longitude blocks.

The original latitude data were sorted separately in the north and south hemispheres, then placed in exactly **10 source-only latitude bands**, with 5+5 southern and 6+5+5+5+5+5+5+5 northern islands. The [source-island band table](latitude_bands.csv) is fixed independent of the butterfly outcomes. We generated exactly **199 unique nonidentity node relabelings** using NumPy PCG64 seed 20261107 and tested the original AIS identity assignment as replicate 0. Each relabeling exchanges the *vertex labels of the entire symmetric AIS resistance graph* within these narrow-ish fixed latitude strata; this preserves the graph up to isomorphism and the complete edge-value distribution but destroys which particular original island receives which AIS connectivity values. This is NOT independent shuffle of individual island pairs, not a valid directional shipping-route model, and may still alter distance or historical covariation.

Two comparisons were prespecified before the new null scores:
1. AIS alone beyond HRE_LAT, which already contains host candidate composition, richness, historic all-botanical effort and full existing geographic latitude controls;
2. AIS added beyond original HRE_WIND_LAT, controlling the fixed wind resistance field as well.

The baseline, target values, original covariates, train-only standardization, OLS with intercept, prediction clip to [0,1], source folds, and 199 maps were identical for both target years and fixed model contrasts. **Only the AIS feature in the augmented predictor was substituted with its globally vertex-reassigned copy.** We report both MSE and MAE for all/within/between-longitude pair scopes without selecting the favorable split.

## Reproduction, source verification and archived tables

The original real AIS forecast models reproduce all **12 relevant latitude-adjusted PR135 FULL comparison MSE and MAE values**, with maximum absolute difference approximately **4.8×10⁻¹³**. An independently written Python standard-library auditor reconstructed:
- all 51 source island identities, 10 signed-latitude bands, the original identity map and all 199 unique same-band nonidentity node mappings (10,200 mapping rows);
- all 1,275 original unordered butterfly pairs, both original response values and 15 longitude-block holdout populations, with no test pair endpoint allowed in training;
- **5,100 fixed reference predictions + 1,020,000 AIS-augmented predictions**, every corresponding squared/absolute error;
- **2,400 original/permuted source/model/scope losses and all 12 observed-to-null rank comparisons**.

The numerical outcome is permanently retained as the [12-row observed-versus-placebo table](observed_vs_latitude_node_null.csv), [51-island original latitude-band table](latitude_bands.csv), the [frozen before-fit protocol](PROTOCOL.md), the [full model runner](../../scripts/run_island51_AIS_node_null.py), the [independent raw-prediction and rank verifier](verify_AIS_node_null.py), and the [published-table consistency checker](verify_published_tables.py).

The [GitHub Actions scientific reproduction and artifact](https://github.com/zuizui0223/chocho/actions/runs/42336810716) includes full-precision 2,400 score rows, all 199 null maps, 5,100 baseline forecasts, 1,020,000 gzip-compressed augmented forecasts, numeric/source checksums and environment. The workflow artifact has a retention limit of 30 days; the branch retains the fixed analysis code and source-locked protocol so that the evidence can be regenerated after expiry.

## Inference boundary and paper positioning

This is **additional observational prediction information arising from the observed AIS maritime network arrangement versus this ONE precisely specified graph-label placebo**. The corresponding finite 199-shuffle rank fractions are descriptive; latitude bands do not preserve every geographical or biological property of the island vertices and are not an exact null of vessel transport. AIS itself represents indexed vessel activity intensity and source marine resistance, not bilateral ship trips, cargo carriage, butterfly dispersal, resource-mediated survival or local population establishment. Source historical plant-species identity is missing: recorded candidate hosts pool 2001–2024 and can **postdate** the later 2011–2022 butterfly observations. The early all-plant observation-effort controls do not make the candidate-host predictor historically valid.

The positive source-specificity check should be paired with the adverse earlier-butterfly findings and all previous sensitivity limits. The global GEB larval host portfolio reconstruction (+54.9% potential resource opportunity), original port/island correlations, AIS/wind source comparisons and their earlier negative scores remain distinct and unchanged.
