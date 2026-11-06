# Direct-native AIS source-quality sensitivity input

This is the separate **post hoc source-quality sensitivity** specified in
[`NATIVE_SOURCE_SENSITIVITY_PLAN.md`](../../../../results/island51_observed_ocean_circuit_20261106/NATIVE_SOURCE_SENSITIVITY_PLAN.md).
The primary analysis had already been seen when a source audit found that the
published first overview sampled native pixels. The sensitivity plan was fixed
before deriving this grid or fitting its biological outcomes. The primary
overview-based grid, frozen manifest, protocol and results remain unchanged.

## Contents and method

- `commercial_AIS_native_mean_0p5deg_2015_2021.npz`: the direct native mean, using the same
  `AIS_mean`, `source_valid_mask`, `raster_transform`, `lat_centers` and `lon_centers`
  keys as the primary input. Shape is 340 × 720, with exact affine
  `[0.5, 0, -180, 0, -0.5, 85]` in EPSG:4326.
- `commercial_AIS_native_mean_0p5deg_2015_2021_receipt.json`: archive and native TIFF
  identity, exact source/target grids, source-only summaries, output hashes, and
  hashes of the plan lock and independent audit.
- `native_overlap_validation_39_cells.json`: every independently computed native
  overlap mean and the matching derived value, with row/column and errors. The
  full independent audit remains preserved in the results directory and is
  pinned by its SHA256.

The reproducer verifies the official archive SHA256 and the extracted native
TIFF's byte count and ZIP CRC32, then records its SHA256. It rejects changed
dimensions, transform, NoData, scale/offset, or any detected overview. External
overview discovery is disabled; the base-band `WarpedVRT` uses geographic
area-overlap averaging directly from the original 0.005° pixels. Source NoData
is excluded, valid zeros are included, and scale 1 / offset 0 are retained. No
unverified unit conversion is applied. Values remain the **published AIS
density index**, whose physical calibration has not been independently
established.

The source-only summaries compare all common valid cells, including land zeros.
Their unweighted grid mean is not a global vessel count or marine-only mean.
This derivation does not read biological outcomes or fit models.

## Reproduction

From the repository root, install `pip install -e '.[test,analysis]'`. Download the
official [World Bank commercial AIS archive](https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0045405/shipdensity_commercial_.zip)
and extract its native TIFF with a CRC-validating ZIP reader. The extracted TIFF
requires approximately 9.82 GB of temporary disk space; the averaging uses a
128 MB GDAL cache and a 128 MB warp-memory limit.

```bash
mkdir -p /tmp/chocho-native-ais/source
curl -fL \
  'https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0045405/shipdensity_commercial_.zip' \
  -o /tmp/chocho-native-ais/shipdensity_commercial_.zip
unzip -q /tmp/chocho-native-ais/shipdensity_commercial_.zip \
  ShipDensity_Commercial1.tif -d /tmp/chocho-native-ais/source
python scripts/verification/derive_island51_native_AIS_quality_grid.py \
  --source-zip /tmp/chocho-native-ais/shipdensity_commercial_.zip \
  --native-tiff /tmp/chocho-native-ais/source/ShipDensity_Commercial1.tif \
  --out-dir /tmp/chocho-native-ais/reproduced
```

The default primary grid, independent audit and locked plan are resolved from
the repository. Reproduction writes to a separate temporary directory in this
example. Large original TIFFs are not committed. The original source is licensed
under CC BY 4.0; attribution and additional source context are recorded in the
[bundle README](../README.md).
