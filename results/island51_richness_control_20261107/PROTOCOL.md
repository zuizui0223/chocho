# Fixed exploratory sensitivity: host identity versus host-list richness

## Goal and timing
This model specification is frozen in the new commit before the sensitivity fit runs. The parent host/circuit results and atlas signs were already seen. This is **post-hoc exploratory robustness**, not preregistration before any relevant results, not a new population, and not independent evidence.

## Unchanged evidence and support
Use the original 51 island landmass IDs and 1,275 unordered pairs from the 225-butterfly panel. Use all 194 original biological ports for baseline port-count controls and the already frozen native-AIS, wind, and uniform-ocean resistances from PR128. The later butterfly target is first indexed 2011–2022; the earlier target, first indexed 2001–2010, is secondary context. No island selection, response replacement, source fetching, plant species addition, or circuit recalculation is allowed.

## Prespecified controls, models and comparisons
For island i, define r_i = log(1 + alien_host_species_i), using the original eligible-island roster. Append two **response-free** pair predictors: r_i+r_j and abs(r_i-r_j). These are marginal candidate-host richness controls. The existing models remain:
- C: original geographic reference with native AIS and wind circuits;
- H: original geographic reference with candidate-host Jaccard;
- HC: geographic reference with both host Jaccard and AIS/wind.
CR, HR, HCR append exactly the same two marginal controls to C, H, HC respectively. No other features or interactions are tested.

Primary comparison: CR -> HCR (host Jaccard beyond marginal richness and circuits).
Secondary comparison: HR -> HCR (circuits beyond host Jaccard and marginal richness).
Context: original C -> HC, H -> HC, C -> CR, H -> HR, HC -> HCR. Neither the ordering, richness transform, nor choice of metric is to change after seeing the scores.

## Estimation and evaluation
For each target and 51 folds, train on 1,225 pairs among 50 islands and predict the 50 pairs incident to the excluded island. Use only training rows to standardize, ordinary least squares with intercept and prediction clipping to [0,1]. Preserve both endpoint-held-out predictions for each unique pair. Compute individual squared and absolute errors before averaging by held-out island; average the 51 folds with equal weight. For all comparisons report full MSE, MAE, relative reductions, improved/worsened/tied island counts without p-values or treating folds as independent replicates.

Reproduce and validate all 306 original C/H/HC target/fold MSE and MAE values from the archived parent before accepting results. Verify the parent input lock and the atlas-pinned hash of the parent fold table. Preserve every poor result rather than selecting favorable islands. Source and parent artifacts must not be rewritten.

## Interpretation limits
A candidate-host richness adjustment does not remove plant recording effort, host quality, suitability, shared ancestry, climate, or colonization timing. Host union data pool 2001–2024, potentially later than the target. AIS density and undirected wind resistance do not give actual vessel routes or windborne butterfly tracks. First-indexed butterfly composition is neither verified arrival nor establishment. Positive predictive increments are conditional observational associations, not transport or resource-mediated causal effects.
