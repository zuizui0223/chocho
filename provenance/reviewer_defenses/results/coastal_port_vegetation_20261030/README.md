# Coastal vegetation comparison around ports

This retrospective 2021 WorldCover comparison uses coastal samples at 1 or 5 km from ports and remote shoreline controls on the same island. It does not establish host-plant identity, introduced status, temporal precedence or human-mediated butterfly arrival.

All 8,342 validated controls were resolved: 8,338 raster samples and four known absent-source locations. Requiring at least three valid points per side yields 208 islands / 625 ports at 1 km and 56 islands / 106 ports at 5 km. Each island receives equal weight after averaging its port comparisons. Shore proximity is approximate and differently measured between groups. Islands unable to provide remote same-island controls are excluded from this comparison, not from the butterfly candidate scope.

The Japanese report and figure show point-frequency differences, not area cover or plant abundance. Bootstrap intervals are descriptive; terrain, disturbance and other covariates are not controlled. Prior 15 km-ring estimates use different support and should not be interpreted as a before/after causal adjustment.

Run `python verify_snapshot.py` to verify hashes and independently recompute every island-weighted mean and bootstrap interval from paired aggregates. Point-level source checks are recorded in the validation receipt; raster pixels were not independently resampled. Included analysis scripts document construction but require the original local workspace inputs.

Sources: ESA WorldCover 2021 v200 (CC BY 4.0; https://worldcover2021.esa.int/), GSHHG 2.3.7 high-resolution L1 (LGPL), NOAA-distributed 2019 World Port Index catalogue. No raw GBIF records or restricted host-network edges are included.
