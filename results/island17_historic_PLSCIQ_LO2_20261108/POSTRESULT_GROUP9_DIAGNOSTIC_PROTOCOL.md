# Post-result diagnostic: distinguish group9 TEST population from group9 TRAIN inclusion

Exploratory post-hoc diagnostic conducted AFTER viewing the already frozen 17-island source model:
- Historical 2006–2010 median island PLSCIQ SUM: later 2011–2022 butterfly test 136 pairs, double-endpoint-unseen MSE gain **−2.923251534344%**.
- The precommitted group9-excluded 16-island re-fit predicts 120 original source pairs with MSE gain **+8.146126807524%**.

These original contrasts and their cohorts remain fixed; a source reproducibility auditor has passed on the original 2048 model forecasts.

Original island group 9 (Borneo landmass) bundles multiple early-indexed WPI ports from at least two countries (ID, MY); it may have unusual source geography and observation opportunity, but the above direction flip does NOT establish geographic aggregation as a causal reason.

## Exact diagnostic from SAVED predictions only — NO refitting/tuning
Read source-fixed results/island17_historic_PLSCIQ_LO2_20261108/LO2_heldout_predictions.csv.gz, old 136-pair manifest, and original 20 contrast summary. Preserve every original 2048 record unchanged.

For BOTH butterfly targets and all 4 source-fixed model families, define 4 exact score panels:
A. ALL17_ALL136: original 136 held-out pairs from the full17 model, frozen source result.
B. FULL17_TRAINED_NO9_TEST120: **the same original full17-trained model forecasts**, restricted to the 120 original test pairs that do NOT contain original group9. In the full17 source LO2 folds, group9 source records can still enter training; neither heldout pair endpoint ever enters fitting.
C. REFIT_NO9_TEST120: same 120 test pairs, but original PRECOMMITTED full no-group9 refit forecasts, so group9 was absent from both training and testing.
D. FULL17_TRAINED_GROUP9_TEST16: the remaining 16 original group9-containing held-out pairs, using original full17 model forecasts. Small 16 pairs is a descriptive decomposition, NOT independent inferential replicate.

Exactly 2 original target periods × 5 original model contrasts × 4 panels = **40 source prediction-loss summary comparisons**, with MSE, MAE, signed relative improvement and pairwise better/worse counts. **No optimization, no model recalibration, no new validation folds, no invented confidence intervals or p-values.**

Independent standard-library verifier recomputes all 40 rows from original 2048 archived individual per-pair squared/absolute errors, verifies exact 136/120/16 coverage, parent 20 contrast score preservation, and that full17 and refit-no9 row labels / original source pair indexes are equal on the 120 common pairs.

Interpretation must distinguish:
- test-composition effect: A vs B (removing the 16 group9 test pairs WITHOUT changing fitted models),
- training-composition effect: B vs C (same 120 test pairs, refitting without original group9 training pairs),
- nonindependence and different source train populations (105 vs 91 source pairs/fold).
This is not causal proof that the Borneo multi-country source grouping explains all differences. Do not drop other islands based on scores, do not combine source AIS modelled edges as historical vessels, and do not change or overwrite earlier GEB / later AIS positive and adverse results. No journal submission.
