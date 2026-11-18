# Fixed post-hoc source-factor decomposition on longitude block holdouts

## Context
Previously inspected results: PR132 preserved the positive later-butterfly combined native AIS+wind increment +4.4744% held-pair MSE under a strict both-endpoints-unseen, geographically grouped evaluation. Earlier-butterfly combined circuit increment was −10.2742%. These results were known before this source-factor split was specified. This is EXPLORATORY robustness on the SAME 51 islands, not independent causal or time-valid evidence.

## Identical geography, population, and method
The existing target-free five longitude blocks (11,10,10,10,10 islands), all 15 block-pair folds, all 1,275 source island pairs, all 194 original botanical ports, the previous 182 marine-qualified ports, original candidate-host Jaccard (pooled 2001–2024), earlier butterfly similarity as later-target baseline, original log-distance and port counts, candidate-host marginal richness and old 2001–2010 all-botanical record-count pair controls are unchanged. The audited uniform ocean resistance remains in every model. Source matrices are the original direct-native AIS and original mean wind-speed undirected Circuitscape graphs, reused without recalibration. Their physical correspondence with vessel movement or directional dispersal is unverified.

## Prespecified 2-by-2 marine source comparison
All models retain common HRE covariates including candidate-host Jaccard and all source-count controls:
- HRE: no additional native-AIS or wind resistance;
- HRE_AIS: HRE + direct-native AIS resistance only;
- HRE_WIND: HRE + original wind resistance only;
- HCRE: HRE + AIS + wind (exactly reproducing PR132 HCRE).

Fixed contrasts, both targets and all pair classes:
1. HRE→HRE_AIS, AIS main incremental information without wind;
2. HRE→HRE_WIND, wind incremental information without AIS;
3. HRE_AIS→HCRE, **wind** incremental information conditional on AIS;
4. HRE_WIND→HCRE, **AIS** incremental information conditional on wind;
5. HRE→HCRE, original combined marine increment (must reproduce PR132 both target records and all three scopes to 1e-9).

No outcome-selected source scaling, group selection, interaction modeling, model pruning, or parameter tuning is allowed. Predictive differences are conditional, not unique causal contributions, and different comparison denominators must not be added as if an exact partition. Retain all adverse comparison results.

## Double-unseen evaluation and checks
For every unordered block pair, train on pairs involving neither held-out block and predict pairs with both endpoints from the selected held-out block(s). The original all-pair sample is partitioned exactly once: 235 within, 1,040 between. Train-only z-scoring, OLS intercept, [0,1] clipped predictions. Test both later (2011–2022) and earlier (2001–2010) butterfly Jaccard separately. Export all 10,200 predictions, 120 model/fold losses, 150 contrast/fold results and 30 target/contrast/geography summaries. Independently reconstruct original response identity, block and train/test support, prediction errors, all fold and summary results. No pair- or fold-independent p-values/intervals.

## Caveats
Pooled candidate-host list can include plant records later than butterflies; pre-2011 controls measure general botanical recording effort, not historic host-species identities. AIS density is not observed actual routes and the wind field is undirected mean speed. "Marine circuit improves prediction" cannot be translated into "plants/butterflies were transported by ships" without independent route, arrival, feeding and colonization data. Longitude blocks are target-free sampled strata, not independent regions. Older target is not a negative control for modern environment because its observation effort, composition and climate differ.
