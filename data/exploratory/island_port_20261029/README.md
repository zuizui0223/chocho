# Island, port and barrier measurements

The [result and interpretation](../../../provenance/archive/host_hindcast/ISLAND_PORT_BARRIER_RESULTS_20261029.md) is the entry point. The [measurement plan](../../../provenance/archive/host_hindcast/ISLAND_PORT_MEASUREMENT_PLAN_20261029.md) records the pre-calculation choices and geography-only pilot selection.

- `island_geography_exposures.csv`: all 8,260 ≥5 km² source islands, area, continental coastline-gap intervals, WPI proximity/count bounds and Antarctic ambiguity flags. One WPI (0,0) placeholder is excluded; 3,688 registered entries contribute.
- `island_geography_manifest.json` and `island_geography_support.csv`: source hashes, independent prior WPI digest match, named seeds and geographic comparison support.
- `named_island_neighbor_gaps.csv` / `.json`: eight named islands' nearest-other-landmass versus direct-continental gaps, including sampling ambiguity and candidate-pruning evidence.
- `island_host_support_summary.json`, `island_country_support.csv`, `island_stratum_support.csv`: preserved country-panel contrasts and zero positive support in the intended island host comparisons.
- `island_butterfly_pilot_selection.json`: 24 physical islands fixed by geography before response queries, plus source geometry hashes.
- `island_butterfly_pilot_manifest.json`, `island_butterfly_pilot_counts.csv` and the four facet CSVs: all 48 exact-polygon queries, response receipts and complete retained aggregate observations.
- `island_butterfly_pilot_design_support.csv` / `.json`: coverage after acquisition, without dropping the original selected cases or fitting biological effects.
- `island_port_geographic_support.png` / `.svg`: geographic comparison figure; no ecological-effect estimate.

Distances are spherical source-coast measurements with discretisation bounds, not butterfly dispersal cost. A nearby WPI entry may be on another landmass and does not establish historical shipping. No indexed butterfly records means no qualifying records, not absence. COL species keys in the new panel must not be silently merged with the older GBIF numeric-backbone keys.
