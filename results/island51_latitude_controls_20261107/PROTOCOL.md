# Locked exploratory latitude robustness test of the original 51-island AIS/wind prediction

## Known before specification
PR133 had positive later-butterfly double-unseen predictive gains: AIS beyond wind +3.1581%, wind beyond AIS +2.3030%, combined AIS+wind +4.4744% mean squared error (MSE) reduction over 1,275 island pairs. For earlier 2001–2010 butterflies, combined circuits worsened MSE by 10.2742%. These were known before this new latitude sensitivity was specified. This is post-hoc exploratory robustness in the SAME island panel, not independent preregistration or a new validation sample.

## Source and geography
Use exactly the same 51 source island IDs, 1,275 unique pairs, original 194 biological ports, 182 qualified marine ports, original pooled 2001–2024 candidate-host Jaccard, fixed candidate-host richness controls, pre-2011 all-botanical recording-effort controls, earlier butterfly baseline (only for the later response), great-circle distance and original port-count controls, uniform ocean circuit and existing native AIS and wind resistance matrices.

Both island-pair endpoints must remain absent from each of the same 15 held-out longitude-block training sets (five longitude blocks of 11,10,10,10,10). Evaluate every source pair once, with 235 within-block and 1,040 between-block test pairs. No longitude block definition or model scaling changes.

## Exactly four new response-free controls
Use the original 51 port-group representative latitudes from the audited atlas file results/island51_host_circuit_atlas_20261106/island_geography.csv, verified by Git blob f0c5b1139e54b4872d28682f0991f89f2368ab43. For every unordered island pair (i,j) append:
1. lat_i + lat_j (signed degrees);
2. abs(lat_i) + abs(lat_j);
3. abs(lat_i - lat_j);
4. an indicator that lat_i and lat_j have strictly opposite signs (zero latitude is not opposite).

No learned climate proxies, nonlinear selection, outcome-driven cutoff, or additional control substitutions. These terms represent coarse latitude and hemispheric geometry, NOT complete climate or geographic barrier control.

## Eight models and twelve comparisons (all predefined)
Retain exactly PR133's HRE, HRE_AIS, HRE_WIND and HCRE, and create latitude-augmented versions HRE_LAT, HRE_AIS_LAT, HRE_WIND_LAT and HCRE_LAT. No other model changes.

For each target window (later 2011–2022 primary, earlier 2001–2010 context) and each scope (all pairs, within longitude block, between longitude blocks), report:
- all five original PR133 AIS-only, wind-only, wind conditional on AIS, AIS conditional on wind, and combined marine circuit contrasts;
- those same five contrasts with latitude controls present in both compared models (PRIMARY);
- HRE to HRE_LAT and HCRE to HCRE_LAT (additional latitude information in the baseline and joint models).

All 30 original source-factor comparison rows must reproduce with source pair identity, both MSE and MAE reference/augmented values (tolerance 1e-9) before accepting new results. Train-only standardization, OLS intercept, [0,1] clipping and pair-weighted MSE and MAE unchanged. Preserve every adverse comparison and all 72 target/comparison/scope summary rows. Target 20,400 source-pair predictions, 240 model-fold rows, 360 comparison-fold rows.

Implement an independent standard-library verifier checking the original 1,275 pair identities, two target values each, all 15 folds, that no test endpoint occurs in training, every prediction and squared/absolute error, all 240 model folds, 360 contrast folds and 72 summaries. No naive p-values or independent-fold confidence intervals.

## Inference limitations
AIS reflects vessel activity index not direct vessel routes; wind is UNDIRECTED mean wind-speed resistance not a flight vector or drift path. Both may share geography, climate and observation-effort covariance beyond latitude. Candidate-host Jaccard still pools 2001–2024 and can postdate the butterfly response; early botanical effort is not dated host taxon identity. First-indexed butterflies are not verified arrivals, larval feeding, reproduction, colonization or abundance. Longitude-stratified double-unseen CV is a stricter prediction task in the same already selected set, not an independent geographic source or a causal transport test.
