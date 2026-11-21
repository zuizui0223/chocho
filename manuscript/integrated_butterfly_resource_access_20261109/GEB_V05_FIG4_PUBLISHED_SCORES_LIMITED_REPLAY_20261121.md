# GEB v0.5 — independent GBIF Fig.4 **evaluation-score** source replay, NOT raw model refit

**2026-11-21. Coauthor analytical QA on Draft PR #189; manuscript PR #179 unchanged.**

## Scope of executed analysis

- The 155-species retrospective butterfly-record forecast has **eight original saved evaluation modes × seven model scores** in the original source-locked `forecast_validation_all_metrics.csv` (Git blob `3e2dcf109b48d3c3766e34ff7603d5061cb850b2`). Source `M3_both` is background climate + land cover, while `M5_host_composition` adds the published known-host resource/composition information. The score-only replay recalculates the species-equal MSE percentage reduction (M3−M5)/M3 × 100 on every original cohort.
- The 224-butterfly recorded community has **169 heldout botanical regions × four models = 676 original region-output score rows** in `gbif_global224_observed_beta_region_holdout_v1.csv` (Git blob `be2018b452e6c44eb3c60282e9d4af6414d9a70c`). The replay independently averages 169 regional MSE values for each model; it compares `M2_other_aliens` to `M3_both` to obtain the incremental known-host H after other-alien O improvement.
- Both original model outputs use a **shared initial GBIF original download `10.15468/dl.umqmsf`** (14,637,152 downloaded occurrence records; publisher list of 1,716 underlying datasetKeys independently verified in prior data-provenance audit). They are *different model response/eligibility samples using the same source download*, not independent GBIF acquisitions.
- The source score tables do not contain individual original `gbifID` record points or the complete GBIF/WCVP/HOSTS model feature and train/test splits. **This validation reruns the arithmetic on original model outputs, not full source-to-prediction models.**

## Source-exact values (recomputed, not inferred)

| Evaluation | Source population | Recomputed species/region-equal MSE improvement |
|---|---|---:|
| Retrospective 155-species temporal evaluation | 155 species, 208,482 butterfly×cell tests | **+4.1201756%** |
| 155-species spatial block evaluation | Same 155 species, 208,482 source test combinations | **+2.9952060%** |
| 155 species, 100-km spatial buffer scope | 155 species, 102,728 source test combinations | **+2.7356579%** |
| 155 species excluded from fitted training | 155 source-evaluated butterflies (same row count), distinct training geometry/model constraint | **−0.2085202%** |
| 155-species effort≥20 temporal sensitivity | 159,084 test combinations | **+2.9757493%** |
| 155-species 250-km spatial limit | 180,654 test combinations | **+3.7626497%** |
| 155 species, largest GBIF source provider removed | 191,004 test combinations | **+5.5455643%** |
| 155 species, original coordinate precision≥100m subset | 135,532 test combinations | **+5.3449401%** |
| 224-species recorded-community region holdouts | 169 regions, 14,196 dependent region pairs; 169 region tests per model | **+0.09129078%** incremental H after O (original, not effort-thinned) |

The **negative unseen-species result is retained**; it is not an outlier to delete. The modest 224-species unthinned held-region increment is not mischaracterized as a large independent effect. All temporal/spatial subsets are not new biological populations; all are model evaluation units using original source records.

## Reproducibility implementation

- [Pure-stdlib recomputation script](../../../scripts/replay_GEB_v05_fig4_stored_model_scores_20261121.py)
- [Successful CI #44838019372](https://github.com/zuizui0223/chocho/actions/runs/44838019372), [5-file ZIP in Actions artifact #14679506760](https://github.com/zuizui0223/chocho/actions/runs/44838019372/artifacts/14679506760)
- ZIP `GEB_V05_FIG4_155_224_ORIGINAL_EVALUATION_SCORES_ONLY_NOT_MODEL_REPLAY.zip`: SHA256 **`aa87dcf6b58f68d2533557e2940a1b294b372ddd3c096257dc6e0de26bf74d25`**.
- The CI pins each original derived score table's Git blob, verifies results against original forecast temporal/spatial/buffer JSON and original 224-species model result JSON, scans five ZIP file names/contents and metadata for direct author identity/paths, extracts and reruns the ZIP in a separate directory. The independent extracted ZIP also ran successfully in the working environment. Execution labels: `GEB_FIG4_155_8_FORECAST_MODES_SCORE_SOURCE_RECALCULATED_PASS`, `GEB_FIG4_224_169_HELD_REGION_H_AFTER_O_SCORE_SOURCE_RECALCULATED_PASS`, `GEB_FIG4_RAW_GBIF_MODEL_RETRAINING_AND_THINNING_NOT_CLAIMED`.

## Not closed by this package

Full GBIF original occurrence-to-forecast-panel construction, 155-species random-forest model training, 224 species community Jaccard response modelling, its observed provider/day effort data, the 29-draw thinning simulations, entire world holdout re-fits and South American dyad calibration, full source WCVP/HOSTS host-edge transforms, and third-party data reuse/anonymous reviewer repository approvals remain outside this five-file bundle.

The published 224-species raw observation score does **not** substitute for the 29-draw effort-standardized source claim; those results have separate source receipts. No full journal data archive was deposited, no new GBIF derived DOI was minted, and no manuscript scientific numerical results were changed.
