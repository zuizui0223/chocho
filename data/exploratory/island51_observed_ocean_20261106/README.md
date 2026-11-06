# Frozen sources for the island51 observed marine circuit analysis

This bundle contains the external geography and environmental inputs for the fixed **51-island, 1,275-pair** analysis. Its manifest records the original source receipts, byte sizes, SHA256 hashes, and the protocol hash frozen before the new biological outcome models were fitted. The original biological pair table remains in `provenance/reviewer_defenses/results/island_host_network_20261030/`.

The analysis uses **182 marine-access ports from the original 194-port ledger**, retaining all 51 islands. Every original port remains in the eligibility audit. A qualifying port must lie within the fixed 25 km attachment limit of the largest coastline-defined marine component. The 12 exclusions include inland or river ports and other locations without qualifying access at the declared resolution; their original harbor types, coordinates and failure reasons are retained individually. Successful attachment distances are recorded; excluded rows have no accepted attachment distance. These exclusions do not remove biological island pairs.

## Fixed geography and observation scales

| Input | Preserved source and analysis representation |
|---|---|
| Coastline | GSHHG 2.3.7 intermediate level-1 polygons; **0.0625°** land/sea center mask with four-neighbor connections and longitude wrap. Separate connected water patches within each **0.5°** block remain distinct graph nodes. |
| Commercial AIS | World Bank/IMF January 2015–February 2021 position-report density. The checksum-verified **0.01° external overview** is averaged directly onto the exact **0.5°** grid; this is not interpolation of a previously derived 1° grid. |
| Wind | NOAA NCEP/NCAR Reanalysis 1, 74 monthly 10 m U/V fields for the same period. The bundle contains the derived **1°** mean monthly-vector speed; the runner repeats each value over its 0.5° children without adding observational resolution. |
| Ports | The NOAA-hosted NGA WPI2019 source GeoJSON and the original 194-port assignments. The runner reconstructs the original port-group centroids and verifies all 1,275 archived geographic distances. |

The 1° centroid-only topology failed a source-geography check. The first finer 0.125° mask still assigned two islands identical terminal distributions. The final 0.0625° refinement and explicit marine-port filter were fixed using coastline geometry, source port locations and harbor metadata before evaluating new biological predictions. Earlier failed geometry audits are retained under `geometry_audits/`. The final graph preserves the connectivity of the fine center mask; aggregation is not an exact electrical reduction of a fine-grid circuit.

The complete specification is [frozen protocol](../../../results/island51_observed_ocean_circuit_20261106/PROTOCOL.md). The frozen input inventory is [source_bundle_manifest.json](source_bundle_manifest.json).

## Reproduce from the repository root

The frozen bundle supports the circuit and held-out-island run without downloading the original large raster or NetCDF archives. Install the declared analysis environment:

```bash
python -m pip install -e ".[test,analysis]"
```

Run the source, geometry and numerical checks:

```bash
python -m pytest -q \
  tests/test_island51_ocean_circuit.py \
  tests/test_run_island51_observed_ocean_circuits.py \
  tests/test_build_island51_coastal_ocean_graph.py \
  tests/test_build_island51_1degree_commercial_AIS_source.py \
  tests/test_preflight_island51_real_GSHHG_ocean_circuit.py \
  tests/test_exploratory_island51_circuit_network_holdout.py
```

Run the frozen analysis:

```bash
python -u scripts/run_island51_observed_ocean_circuits.py \
  --input-root data/exploratory/island51_observed_ocean_20261106 \
  --coast-dir data/exploratory/island51_observed_ocean_20261106/coast_0p0625 \
  --marine-ports-only \
  --out-dir _ci_ocean_results \
  --network-dir /tmp/chocho-island51-observed-ocean-networks
```

The runner checks the frozen source hashes and protocol, preserves the complete port eligibility ledger, constructs uniform/AIS/wind circuit matrices on the same water graph, and evaluates the fixed island-held-out comparisons. Small receipts, pair tables, 51-fold losses and resistance matrices are written to the output directory. Large exported network edges stay in the separate temporary network directory. The scoped [GitHub Actions workflow](../../../.github/workflows/chocho-island51-observed-ocean-circuit.yml) runs these same checks and analysis, recording installed package versions with the small outputs.

## Sources, attribution and reuse

- **GSHHG 2.3.7:** [Original shapefile archive](https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip), developed by Paul Wessel and Walter H. F. Smith. The source identifies its license as **GNU LGPL version 3 or later**. Original [license notice](source_notices/GSHHG_LICENSE.TXT), [LGPL text](source_notices/GSHHG_COPYING.LESSERv3), and [source README](source_notices/GSHHG_README.TXT) accompany the derived graph. The graph receipt identifies the precise source member and archive checksum.
- **World Bank/IMF Global Shipping Traffic Density:** [Dataset catalogue](https://datacatalog.worldbank.org/search/dataset/0037580/global-shipping-traffic-density), **CC BY 4.0**. The original commercial-density ZIP URL and SHA256, overview CRC32, georeference and averaging method are recorded in `ais/commercial_AIS_mean_0p5deg_2015_2021_receipt.json`. This repository distributes the derived 0.5° grid, with changes from the original source explicitly documented; it does not imply source-provider endorsement.
- **NOAA NCEP/NCAR Reanalysis 1:** Data provided by NOAA PSL, Boulder, Colorado, USA; [source description](https://www.psl.noaa.gov/data/gridded/data.ncep.reanalysis.html), [monthly 10 m U](https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/uwnd.10m.mon.mean.nc) and [monthly 10 m V](https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/vwnd.10m.mon.mean.nc). The wind receipt preserves the source checksums, period, native Gaussian-grid resolution and transformation.
- **NGA World Port Index 2019, hosted by NOAA:** [PIRO DigitalAtlas Ports layer](https://gis.ngdc.noaa.gov/arcgis/rest/services/nccos/PIRO_DigitalAtlas/MapServer/30). The WPI receipt records the recovered 3,630 source points, paginated query receipts and GeoJSON hash; the coordinate-recovery receipt identifies the 194 original focal-port assignments.

## Interpretation of the reproduced outputs

The electrical measure uses unit current distributed equally over each island's qualifying ports; ports on the same island are not short-circuited. It measures circuit energy under this fixed allocation, rather than observed vessel movement between ports. AIS is a position-report density, includes stationary vessels, and its coarse source average includes valid land-zero pixels. The native-to-overview resampling method is unspecified. Wind is an undirected monthly 10 m speed proxy, not directional transport or butterfly flight-height wind.

The held-out comparisons concern retrospective community geography. The modern environmental fields do not identify historical plant introduction routes or prospective butterfly establishment. Pairwise observations and island folds share endpoints; their losses must not be treated as independent replicates. Biological prediction scores are produced by the runner and recorded with its outputs, separately from this frozen source inventory.

## Native-pixel source-quality sensitivity

The frozen inventory and primary runner above use the original declared external-overview average. A subsequent source-only audit found nearest-neighbor behavior in all159,744 sampled overview pixels and quantified discrepancies from native-pixel averages. The separate `ais_native_quality_check/` derivative is reserved for the fixed [post hoc native-source sensitivity plan](../../../results/island51_observed_ocean_circuit_20261106/NATIVE_SOURCE_SENSITIVITY_PLAN.md). It is not substituted into the original manifest or primary analysis. Its own receipt and the sensitivity outputs identify the new source-grid checksum and exact scope.
