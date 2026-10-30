# Exploratory island / port associations

This snapshot follows the user's request to examine correlations before predictive performance. It does not replace the original 225-candidate / 155-eligible cohort or establish a causal role of islands, ports, or host introductions.

The descriptive island subset contains 137 species, 40 landmasses and 993 cells. Within-species / within-landmass rank correlations with >=10 cells cover at most 5 landmasses. The adjusted fixed-effect comparison requires >=2 cells per species-landmass and covers 93 species, 22 landmasses and 975 cells. These are different supports, not interchangeable estimates.

The adjusted response is 2019–2022 recording probability conditional on 2015–2018 baseline and recorded effort/environment controls. The exposure is distance to a 2019 port layer; this is a contemporaneous retrospective association, not an historically issued forecast. Host introduced proportion is static species-region data, not plant abundance near a port.

The port-distance by introduced-host interaction has intervals crossing zero in all three exploratory comparisons. This is neither evidence of absence nor support for transport or mediation. Species × landmass fixed effects absorb the invariant host main effect. Two-way species/spatial-block covariance is approximate, with 18, 15 and 9 spatial blocks; the covariance matrices have small negative eigenvalues, so no joint Wald tests are reported. The additional t-interval sensitivity preserves original normal-approximation intervals.

## Reproduction and limits

`verify_snapshot.py` runs with Python's standard library and checks hashes, arithmetic and stated support. It does not rerun models. The archived analysis scripts require the local derived parquet inputs under `work/temporal/forecast` and `work/barriers`, and resolve paths relative to their original workspace location; copying this folder alone is not a turnkey source-data pipeline. No raw occurrence records or restricted plant-host edge sources are included.

Source and upstream model context: preceding broad_landmass_validation_20261030 snapshot and PR #101. Port source: NOAA ArcGIS PIRO Digital Atlas MapServer layer 30, World Port Index 2019. Landmass source: GSHHG 2.3.7 high-resolution L1; source licence LGPL. Interpretation limits and local audits remain part of the Japanese results.

`stable_provider` filters absolute change in provider dominance <=0.1, not stable provider identity or composition. Equally dominant providers can turn over.
