# Host information improves known-species temporal and spatial forecasts

This is a verified research milestone, not completion of the entire analysis.

Read [日本語の統合結果](forecast_validation_results_ja.md) and [comparison figure](forecast_validation_transfer.png). Host information reduces MSE by approximately 4.1% in temporal validation, 3.0% in held spatial blocks, and 2.7% at test sites at least 100 km from every training cell. The >=100 km result is a post-hoc geographic subset of the same spatial predictions, not a separate refit or prospectively fixed buffered validation. The held-species evaluation does not show this gain. Report both MSE and MAE, all model comparisons, and uncertainty; do not treat these predictive gains as causal butterfly range expansion.

The response is a rarefied recording probability, not confirmed occupancy, colonisation or breeding. Static host resources do not measure host arrival dates. Distance to historical records is a baseline predictor; geographic connectivity and sea barriers have not been separately tested. Distance from training sites is a validation device, not a dispersal estimate. The interaction analysis is exploratory and did not support the nine tested interactions after multiplicity adjustment.

The manifest identifies the completed validation modes in this snapshot; the Japanese report explicitly lists pending sensitivity analyses. Later verified increments will update this PR. The evaluation protocol and model lock retain their archived pre-fit execution-status statements; consult the integrated report and manifest for current status. No raw GBIF records, observer information, restricted WCVP/HOSTS source tables, credentials or fitted models are included.

## Verification

Run `python verify_snapshot.py` from this directory (standard library only). It reconciles published model scores, paired contrasts, census coverage, and snapshot checksums. This is a snapshot consistency check, not independent re-execution of the original scientific model. The source workspace separately recomputed scores from saved row-level predictions and checked source/geometry/rarefaction inputs.

`analysis_source/` contains archived analysis scripts for inspection. They require the original local analysis layout, derived panels and access-controlled host source archive; they are not turnkey from this public snapshot. Their original workspace-relative paths are retained. The live queued model process used frozen settings; no holdout-based tuning was performed.

Sources: [GBIF occurrence DOI](https://doi.org/10.15468/dl.umqmsf), [annual ESA LCCS aggregation](https://doi.org/10.5281/zenodo.11407287), NASA POWER monthly climate. The original host catalogue combines accepted HOSTS/WCVP identities with the pinned BCE supplement. Source and scope receipts are included without redistributing source tables.
