# Original 51-island double-unseen port-roster leverage sensitivity (exploratory)

## Prior information and purpose
PR134's longitude-block both-endpoints-unseen results were already inspected: among later-butterfly records, latitude-adjusted AIS-only improved MSE +1.0020%, AIS conditional on wind +0.5446%, latitude-adjusted wind alone worsened MSE -4.6618%, and AIS+wind combined worsened -4.0918%. No new source population is claimed. This new sensitivity is exploratory and frozen before the following re-fits, not independent preregistration.

## Source-only leverage selection
Retain the exact original 51 source landmass islands, all 1,275 pairs, original 225 candidate butterflies, fixed 194 source botanical ports and 182 qualified marine circuits, AIS/wind/uniform original matrices, pooled 2001-2024 candidate-host Jaccard, marginal host richness controls, 2001-2010 all-plant port buffer record-effort controls, original distance, port-count and latitude/hemisphere controls, earlier butterfly similarity (only for later target), and exactly the same 15 longitude block-pair heldout splits from PR134.
Determine port-heavy original landmasses from the source-only original port-count roster, sorting by descending count then ascending numeric ID. The expected first two are landmass **15 (39 original ports)** and **9 (21 original ports)**. Fail source validation if these counts, identities or original total 194 differ. This selection is NOT based on any butterfly score. These are landmass IDs, not claims of unique named islands, historical port operation, or complete country-level coverage.

## Fixed four training populations (identical test pairs)
1. FULL: original PR134 train pairs.
2. DROP_15: exclude any training pair incident to landmass ID 15.
3. DROP_9: exclude any training pair incident to landmass ID 9.
4. DROP_15_9: exclude any training pair incident to either ID 15 or ID 9.

For each scenario and each of the **15** fixed unordered longitude-block-pair folds, maintain exactly the original heldout test pair set (one of 45, 55, 100 or 110 original pairs). Training contains no endpoints belonging to the two held-out longitude blocks and none of the additional scenario-excluded IDs. Some scenario exclusions are redundant in folds where the focal island was already held out: DO NOT force a change to training size in those cases and do not drop that fold or its test pairs. Each of all **1,275 original pairs** is tested ONCE per scenario and model; within-block 235, between-block 1,040. Overlapping island test folds are not statistical replicates.

## Exactly four latitude-adjusted models, five source contrasts
Use the four **unchanged PR134 latitude-controlled models** HRE_LAT (host candidate Jaccard + marginal host richness, early botanical effort, geography including 4 latitude terms and uniform-ocean circuit), HRE_AIS_LAT (+ AIS), HRE_WIND_LAT (+ wind), and HCRE_LAT (+ both). Fit train-only-standardized OLS with intercept and predictions clipped to [0,1].

Fixed five contrasts: HRE→AIS only, HRE→wind only, AIS→AIS+wind (wind conditional on AIS), wind→wind+AIS (AIS conditional on wind), and HRE→AIS+wind combined, each using LAT models. Both response windows (later 2011–2022 primary, earlier 2001–2010 context), three pair classes (all 1,275, within-longitude 235, between-longitude 1,040), and all four scenarios must be output. Both MSE (primary) and MAE (secondary) with raw scored folds, improvements AND deteriorations. No model/threshold selection after seeing scores, and no independent p-values.

Expected records: 4 training scenarios × 2 targets × 4 models × 1,275 predictions = **40,800 individual predictions**; **480** model/fold records; 4 × 2 × 5 × 15 = **600 contrast/fold records**; 4 × 2 × 5 × 3 = **120 scenario/target/contrast/scope summaries**; all 15 original source-only geographic folds and all source exclusions must be preserved. Parent FULL scenario's original 30 PR134 latitude-specific contrasts must independently reproduce the committed original 72-row target/contrast/scope audit for each MSE and MAE reference/augmented value to absolute 1e-9. Other scenarios change training only; all saved test observations, pairs, statuses and source predictor values are unchanged.

An independent standard-library audit must re-derive the source block folds from the original pair table and original representative longitude sources, independently check source-only ranking of the original ports, every scenario training-pair exclusion and every test pair assignment, source pair observation identities, all individual prediction losses, fold and scope summaries. Publish full exact result CSVs, selected port table, matching source/scenario verification receipts, and a standalone summary that does not cherry-pick outcomes.

## Interpretation boundary
This directly tests whether the previously observed latitude-adjusted AIS/wind gains are sensitive to TRAINING data associated with port-heavy island groups. It does NOT estimate the biological effect of the ports themselves, establish shipping routes, disentangle wind direction or latitude/climate, or reconstruct early host species identity. Pooled candidate-host Jaccard (2001-2024) can postdate the later butterfly outcome. Port counts and historical all-botanical record counts represent heterogeneous potential sampling exposure, not total transport volume or direct plant absence. The same 51 islands and model scope are used, so results remain exploratory and observational, not independent new-region replication.
