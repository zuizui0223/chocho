# Global H3-r4 × WGSRPD3 source-verified geometry bridge — 2026-10-29

**Source-geography implementation COMPLETED.** This is not a historical host-effect model, confirmed local host occurrence, or butterfly climate-tracking inference.

The two original data Actions both succeeded: [first geographic source run #40815269681](https://github.com/zuizui0223/chocho/actions/runs/40815269681) and [overlap-corrected final source run #40815812202](https://github.com/zuizui0223/chocho/actions/runs/40815812202). All original source SHA256s, 369 original WGSRPD3 codes and full H3 index count verified.

## Exact measured results
- All **288,122** global H3 resolution-4 indexes visited; **89,161 cells** intersect nonzero botanical land; all original **369/369** WGSRPD3 polygons included. Source mapped land is **146,472,793.87 km²** in the original project equal-area CRS EPSG:6933.
- **98,529** H3×WGSRPD3 raw polygon intersection rows, **8,937** cells intersecting ≥2 botanical regions. **7,922** H3 cells have centroid outside source-mapped botanical land yet positive land overlap (these would be missing in a centre-only join). **80** cell centroids match multiple original polygons.
- Source WGSRPD polygons **genuinely overlap** in **3,864** grid cells. The largest source-region area overlap excess is **25.198% of a whole H3 cell**. Do not add raw region fractions as if all boundary polygons were disjoint.
- A new **102,660-row disjoint membership atom table**, with **4,282 multi-code membership atoms** in **3,990** cells, partitions source geometry into additive, nonoverlapping combinations (`AAA`, `BBB`, `AAA|BBB`). Each atom's fraction belongs to exactly one membership set.
- A separate independent streaming audit compared the sum of atomic fractions against mapped-land polygon union for **every one of the 89,161 cells**: **0** mismatch above 1e−6, largest absolute discrepancy **5.153×10⁻¹⁰** of cell area.
- H3 average hexagon area **1770.347654 km²**, average edge **26.071760 km**, original H3 Python/C version **4.5.0**.

## Full reusable output and SHA checksums
The completed GitHub Actions [derived artifact #12872072613](https://github.com/zuizui0223/chocho/actions/runs/40815812202/artifacts/12872072613) includes **only** the following original-free outputs, no redistributed raw botanical regions:

| Derived filename | Meaning | SHA256 |
|---|---|---|
| `global_h3r4_wgsrpd3_cells_20261029.csv.gz` | One H3 per row, centre and dominant region, coastal and boundary flags | `9f6299358d836e0f09b7e47766f7b698b33f9d3303d95e1764b3d170f1bdf98d` |
| `global_h3r4_wgsrpd3_overlap_20261029.csv.gz` | One cell × individual region, raw area overlap (not safely additive when source overlaps) | `edfe8906ed8c7543b6a0b42f12cba5dac5830b19a6a487d7b834a4ff033b4fd8` |
| `global_h3r4_wgsrpd3_atoms_20261029.csv.gz` | One cell × exclusive region membership set, **safely additive** | `783b3b3a52f2fb27cf61d2b1be1b56bbe7e806ccdb6289e3f5090b0090d770ba` |
| `global_h3r4_wgsrpd3_crosswalk_manifest_20261029.json` | Source metadata, versions, global QC and all output digests | Included |

Source exact archived WGSRPD3 `inputs/level3.geojson` SHA256 `c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840`, identical to earlier independently verified original country-crosswalk. Original 369 geographic regions differ from the older **365** climate-covered CHELSA regions: missing climate coverage must be treated as **missing**, not fabricated as a botanical absence.

## Next join contract
- **`h3_r4`** is the fixed spatial key, **`year`** the time key; species×H3×year panels must be built from independently audited butterfly monitoring/reporting opportunities.
- WCVP×HOSTS species×WGSRPD3 native/introduced statuses map onto the **exclusive membership atoms**, summing each disjoint land polygon once. These are static regional *potential resources*, NOT proof the host lives in that hexagon, nor of introduction/planting date or human maintenance.
- NDVI/land-cover from MODIS are available from approximately 2000/2001, and must not be fabricated as 1976–1995 historical predictors. NASA POWER source weather may be coarser than these H3 cells.
- Use historic, independently dated larval host records and an observational-effort-matched nonhost/plant background before inferring cultivated-host availability. GBIF zero records ≠ absence; actual original UKBMS survey visits supply a separate audited, qualified nondetection basis.

**Scientific status:** source-verified geographic crosswalk complete. Original host-mediated climate-tracking hindcast and CMIP6 forecast have not been inferred from these outputs. The frozen original GEB manuscript and `main` were not changed.
