# Clean 51-island ocean-host model source replay (2026-11-18)

Source snapshots: [original ocean circuit PR #127](https://github.com/zuizui0223/chocho/pull/127) and [original joint-host PR #128](https://github.com/zuizui0223/chocho/pull/128), PR #128 head `17a4768cbb6cddf249449bd67e3ea1124e0db2f5`. This **clean branch starts from main** (`f9ae63509867c362d6e67fcb7142d366a0ae4672`) rather than merging the original stacked PR ancestry.

This branch restores **26** exact Git blobs missing from main, including 16 original small frozen source files with source-attribution/license notices, the original bundle manifest, original four tests, derived resistance matrices and the original indispensable joint-model source helper. No existing main file is overwritten. Frozen `input_lock.json` validates twelve scientific inputs, and the original source bundle manifest validates sixteen geographic/environmental inputs by recorded bytes/SHA256.

## Source rights and provenance
- World Bank / IMF Shipping Traffic Density transformed overview: **CC BY 4.0**, original attribution and transformation documented in `data/exploratory/island51_observed_ocean_20261106/README.md`.
- GSHHG 2.3.7 shorelines and derived network: **GNU LGPL v3+**, original three source notices restored with the data.
- NOAA wind/WPI: public source and exact source receipt documented in the original source bundle README. These files were already published at immutable original Git blob hashes. No hidden raw GBIF positions, unlicensed HOSTS original interaction table, or external unverified datasets are being added.

## Actual test (run on this branch)
`python -m pip install -e ".[test,analysis]"`, original four pytest modules, then `python scripts/run_island51_host_circuit_joint.py --out-dir "$RUNNER_TEMP/joint"` and `python results/island51_host_circuit_joint_20261106/verify_saved_predictions.py --results-dir "$RUNNER_TEMP/joint"`. The output must contain 25,500 endpoint predictions, 510 model-fold losses and original prior G/C Circuitscape numeric reproduction.

This is a **retrospective prediction of GBIF index-based community similarity**, not measured adult movement, actual plant seed drift, historically directed vessel voyages or newly established butterflies. Conditional host×marine interaction was not independently robust. Keep source rights/notice and side-by-side results with caveats if merged to main.
