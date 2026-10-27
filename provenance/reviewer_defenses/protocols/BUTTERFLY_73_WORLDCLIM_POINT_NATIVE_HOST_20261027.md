# Fixed source-consistent butterfly point-climate audit: WorldClim 2.1

**2026-10-27. Frozen before calculating observed point BIO1/BIO12 overlap.** Independent sensitivity study stacked above the successful exact-WCVP/CHELSA region-map work in PR #74. Existing source regional outcome is unchanged.

## Why this is needed

The verified 6,076 points in 73 butterfly species×country cases were previously categorized using the CHELSA climate of broad WGSRPD3 regions. That is **not climate at each actual butterfly point**. Furthermore, using native-host CHELSA percentiles directly against independently sourced WorldClim values would mix climate products/time windows and create artificial niche novelty. This test **rebuilds the native-host baseline and extracts each observation point in the SAME WorldClim product**.

## Versioned climate authority

- WorldClim version 2.1 public historical climatology **1970–2000**, official UC Davis archive `https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip`, global 10-arcminute bioclimatic layers. Read the 19 global BIO rasters from the official bundle; retain only `wc2.1_10m_bio_1.tif` (annual mean temperature, °C) and `wc2.1_10m_bio_12.tif` (annual precipitation, mm). Verify the TIFF grid, CRS, dimensions, missingness, source ZIP SHA and per-raster SHA. Do not put full licensed rasters into public Git or CI artifacts. Downloaded archive must match the official data schema and contain no path traversal. Source provenance: https://www.worldclim.org/data/worldclim21.html
- Reuse exact source-validated author-only original WCVP/HOSTS/BCE polygons and 365 accepted WGSRPD3 Level3 source regions, original archive SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`. Retain unchanged original 239 butterflies, 235 eligible source-native-host climate species, 1,027 supplement edges.
- **Both the baseline and focal butterfly points use WorldClim v2.1 10-arcminute pixels**. Never compare WorldClim butterfly points to CHELSA source native-host percentile cutoffs.

## Source-constrained spatial baseline

Rasterize the exact original 365 accepted WCVP/TDWG Level3 polygons onto the global WorldClim grid at cell *centres* (all_touched=False), marking regions with ambiguous coverage as missing where detectable. For each focal butterfly, form its native-HOST-PLANT supported region union according to the original frozen source graph. Use **all eligible WorldClim land-climate pixels inside those documented native host resource regions** to construct a species-specific BIO1 and BIO12 5–95% rectangular baseline (main), secondary 10–90%. Primary percentiles are weighted by cosine of pixel latitude (approximate surface-area weighting on the latitude/longitude climate grid). Equal-pixel 5–95% is fixed sensitivity; minimum valid native baseline pixel count >=50 and >=10 documented native host regions.

This baseline is the *climatic opportunity of the original native HOST PLANT resource geography*, **not the native butterfly's physiological or realized climatic niche**. Broad botanical regions include habitats never used by this insect and are not actual local plant-presence masks.

## Point mapping and outcomes

- Take **every** accepted GBIF source observation point in the 73 source-positive species×country cases, preserving each butterfly/country original accepted key, 2010–2025 and five original record-basis filters, all prior source completeness/5k-per-query and 40k-total gates. Coordinates and exact GBIF occurrence IDs are temporary, never uploaded.
- Sample WorldClim raster BIO1 and BIO12 at the pixel containing each point. Mark out-of-grid/nodata/boundary/unknown coordinate explicitly, not 'climate-novel'. Coordinate uncertainty <=1 km separately; 10-arcminute raster pixels limit site precision even when record coordinate is more accurate.
- Define the independent `outside_native_host_WorldClim_rect` indicator = BIO1 or BIO12 point value outside the same-product native-host 5–95% rectangle; retain 10–90% and equal-pixel baseline sensitivity. Do not infer directional climatic adaptation (warm vs cold, wet vs dry) without separate preregistered analyses.
- Report per original species×country case counts, 5–95 point novelty frequency, baseline source sample size, year and QA tallies, **equal-case vs record-weighted** results, and original novel-host-country vs analogue-host-country contrast. Two dominant Taiwan cases are independently preidentified from the previous indexed source roster; show fixed exclusion separately.
- A favorable geometric result is insufficient to claim host use or realized butterfly niche colonization. Full proof needs independent native-butterfly historical occupancy, dated plant arrival, effort-controlled detections, juvenile host use and repeated reproduction. If this source product is inaccessible/incoherent, do not substitute NASA/Open-Meteo values or invent missing climate numbers.

## Information gate and conservation

Claim only **source-conditioned observed butterfly point climate relative to baseline potential native-host resource climate** if all 73 accepted point source queries complete, all accepted 20 butterfly native-host baselines evaluable, and all 73 species×country cases present (unmatched climate remains separately missing). No biological preference p-value, physiological niche change, population expansion, competition, or invasive plant management prescription.

Maintain original GEB manuscript and #72/#73 sources without merge. Save only derived counts and exact source hashes, no raw point coordinates, full WorldClim grids, author-only polygons or identifying observation details.
