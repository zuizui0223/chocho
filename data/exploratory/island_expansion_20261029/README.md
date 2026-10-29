# Fixed 96-island butterfly recording expansion

The authoritative interpretation is [the results report](../../../provenance/archive/host_hindcast/ISLAND_EXPANSION_RESULTS_20261029.md); the [analysis plan](../../../provenance/archive/host_hindcast/ISLAND_EXPANSION_PLAN_20261029.md) was frozen before additional outcomes.

This directory contains the geography-only fixed roster, exact-polygon acquisition receipts/facets, finite-record richness metrics, support gates, and nearest-other-landmass measurements. All 192 island-period queries completed. The 48 original pilot queries retain original source receipts and are explicitly distinguished from 144 new queries.

The primary late-period recording screen retains 23/96 islands, and every prespecified model fails its support gate. No port, isolation, colonization or host effect estimate is supplied. Complete zero records are retained and are not treated as butterfly absence or zero true richness.

Key files:

- `expanded_island_panel_selection.json` / `.csv`: locked identities, geometry digests, selection order, geography and representative coordinates (selection/clustering only).
- `expanded_island_panel_manifest.json`: complete source query results, facets, timestamps, request/response hashes and reuse labels.
- `expanded_island_panel_counts.csv` and four `*_facets.csv`: inspectable recording pools; species/basis unassigned residuals remain explicit.
- `expanded_island_recording_metrics.csv`: all selected islands in both periods, finite-record expected species bins and exclusion reasons.
- `expanded_island_analysis_support.csv` / `expanded_island_recording_analysis.json`: primary/sensitivity support and no-fit outcomes; whole-roster geographic components.
- `expanded_island_access.csv` / `.json`: nearest other >=5 and >=1000 km² coast intervals, candidate identities/ambiguity and pruning checks. Neighbors are not verified biological sources.
- `expanded_island_recording_support.png` / `.svg`: descriptive support figure, visually inspected.
- `expanded_island_verification.json`: independent receipt, consistency, test and output-identity checks.

Source geometries and original occurrence API bytes remain in the bounded source caches; published facets and receipts define this snapshot. A future uncached retrieval can differ. The taxonomy bridge and original official taxonomy response archive are in `../island_taxonomy_20261029/`.
