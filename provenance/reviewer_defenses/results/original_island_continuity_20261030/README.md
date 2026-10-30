# Original-coordinate island record continuity

This snapshot replaces interpretation of cell-centroid geography as island chronology. It retains the 225-species candidate scope; the 2011–2022 first-island-record event table has 205 species, 2,656 islands, and 12,028 events. The published tables contain aggregates, not raw occurrence coordinates or restricted host edges.

The Japanese report explains weighting and observation-coverage sensitivity. Distances link records strictly across years; they are not observed dispersal routes or evidence of island stepping stones. Host/port mechanisms remain untested by this figure.

Run `python verify_snapshot.py` in this directory. It verifies file hashes and independently recomputes distance-group summaries from species counts. The source validation receipt documents the earlier endpoint/year/distance checks but is not a substitute for independently rerunning original geography. `plot_original_island_continuity.py` documents figure construction against the local workspace source table; it requires that non-bundled table and its original workspace layout.

Source butterfly download: https://doi.org/10.15468/dl.umqmsf ; coastline: GSHHG 2.3.7 high-resolution L1 (LGPL). The source occurrence archive is not redistributed here.
