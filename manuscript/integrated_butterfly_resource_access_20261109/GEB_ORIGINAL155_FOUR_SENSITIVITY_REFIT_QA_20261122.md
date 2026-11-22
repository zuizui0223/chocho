# Original155 four-sensitivity source refit — evidence audit (22 November 2026)

**Stage A, Draft source audit only.** All four original sensitivity models were freshly retrained from checksum-verified original derived panels (radius250 required one targeted retry). **Two numerically reproduced** all original predictions and metrics under frozen tolerances; two retained positive host-predictor MSE reductions but failed exact numerical reproduction. This is not raw-observation-to-panel rebuilding. No claim is made that the original observations-to-panel chain was rerun.

| Sensitivity | Original M5 vs M3 MSE reduction | Fresh MSE reduction | Numerical result |
|---|---:|---:|---|
| Observation-effort threshold 20 provider days | +2.975749% | +2.975749% | **PASS 193/193 groups**, zero out-of-tolerance model predictions |
| Radius 250 km | +3.762650% | +3.762650% | **PASS 193/193 groups** on targeted rerun, zero model-prediction errors beyond tolerance |
| Exclude largest provider | +5.545564% | +5.538677% | **NUMERIC MISMATCH** 78/193 groups, despite completed fitting |
| Position precision ≤100 m | +5.344940% | +5.354669% | **NUMERIC MISMATCH** 78/193 groups, despite completed fitting |

Method: the original derived panel and archived predictions for each condition were extracted from the immutable private draft archive with manifest and per-member SHA-256 checks. Models were freshly fitted with the archived SHA-verified source script, the original RF settings and the same Python/replay environment as the already verified original temporal fit. Each comparison used keyed row predictions, all seven score groups, and per-species errors. RTOL = 1e−6, ATOL = 1e−9. No tolerance or original manuscript estimate was changed.

**Failure scope:** For the largest-provider-exclusion condition, out-of-tolerance model-prediction comparisons by M0–M5 were 180,170, 175,309, 160,023, 19,879, 3,286 and 41,299; for precision100 they were 115,507, 115,714, 128,963, 115,220, 77,778 and 3,540. These are overlapping per-model comparisons, not unique butterfly occurrences or a direct estimate of field occupancy. Both newly fitted M5-vs-M3 aggregate results remained positive. Neither is a numeric reproduction PASS.

**Radius250:** The first matrix job stopped during private source recovery, *before* model fitting; its underlying failure cause remains unproven. Source-manifest checks established the exact original panel/prediction paths. A [targeted radius-only rerun #44985167678](https://github.com/zuizui0223/chocho/actions/runs/44985167678) then **succeeded**, verifying the original 727,932-row restricted panel and reproducing all 193 group, prediction and species comparisons using 167,192 training and 180,654 test transitions. [Successful aggregate artifact #14739916146](https://github.com/zuizui0223/chocho/actions/runs/44985167678/artifacts/14739916146).

**Boundaries:** The original full archive, original panels, individual predictions, individual source paths and third-party datasets are not uploaded in these aggregate receipts. A source manifest lists an original fitting-script candidate whose hash matches the saved public script, but it does not prove the historical executable/runtime actually used for each sensitivity or explain the held-species fold-3 discrepancy. Thus no cause of the two numeric mismatches is asserted. Source provenance, the separate 155-species unseen-species discrepancy, 224-species cohort, and reviewer-access/right-to-redistribute gates remain distinct.

[Four-condition refit workflow #44984767232](https://github.com/zuizui0223/chocho/actions/runs/44984767232) · [Source-manifest checks #44985042656](https://github.com/zuizui0223/chocho/actions/runs/44985042656) · [Species-fold-3 M0 probe #44974105799](https://github.com/zuizui0223/chocho/actions/runs/44974105799)

No manuscript PR #179 scientific results are altered, no Draft PR is merged, no coauthor review is initiated, and no journal submission or premature reviewer data release has occurred.
