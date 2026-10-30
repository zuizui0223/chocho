# Dated plant evidence for the fixed island host comparison

The target scope is the existing 28 islands and the early 2010–2017 window. The host source is the recovered original catalogue for the 93 butterflies recorded early in the fixed 96-island panel: 1,440 links to 1,062 source plant identities. Current taxonomic resolution and island plant records are distinct evidentiary steps.

## Which files to use

**Use `early_plants_resolved_manifest.json` and `early_plants_resolved_*_facets.csv` for the completed 28-island plant input.** The original `early_plants_manifest.json` intentionally retains its first-run 25 complete / 3 missing state. The separately frozen pagination repair completes those three cases without replacing any of the original 25 complete results. All 28 resolved queries contain 728,040 indexed plant records in total. No failure is interpreted as absence.

- `early_plants_run_lock.json`, the original manifest and raw archive retain first acquisition evidence.
- `early_plants_repair_run_lock.json`, repair manifest and repair raw archive retain the bounded repair, terminal probes and rechecks.
- `plant_taxonomy_bridge.csv` retains all 1,062 source IDs, query status, current concept evidence and explicit conservative eligibility. `SUCCESS` alone does not make a concept eligible.
- `plant_taxonomy_manifest.json`, collisions and raw responses preserve taxonomy provenance. The original acquisition script is retained separately from the corrected conservative processing code.
- `dated_host_primary_risk_cells.csv` retains all original 2,385 candidate cells and 66 later-positive cells, with recorded-host, unresolved-catalogue and missing-query states.
- `dated_host_butterfly_catalogue_coverage.csv`, `dated_host_target_support.csv` and `dated_host_recording_summary.json` retain denominators and descriptive summaries.

For a fully resolved catalogue, no matched plant record is still not ecological absence. A matched current plant concept does not establish larval use, resource abundance, plant origin, wild rather than cultivated presence, or a causal host effect. Species and year facets are separate margins: the time rule is within 2010–2017, not event-level simultaneity or a ±3-year point comparison.

## Reproduction

After verified taxonomy and plant evidence are available, run:

```bash
python scripts/analyze_island_dated_hosts_20261030.py --plants-manifest data/exploratory/island_dated_host_20261030/early_plants_resolved_manifest.json
```

The acquisition and repair scripts preserve original failures and raw request/response/receipt evidence. Read their fixed plans before any new acquisition. An API retrieval today may differ from the retained historical-index snapshot; do not silently replace the saved evidence.

See the [integrated result](../../../provenance/archive/host_hindcast/ISLAND_DATED_HOST_LATITUDE_RESULTS_20261030.md) and [continuation plan](../../../provenance/archive/host_hindcast/ISLAND_DATED_HOST_LATITUDE_PLAN_20261030.md).
