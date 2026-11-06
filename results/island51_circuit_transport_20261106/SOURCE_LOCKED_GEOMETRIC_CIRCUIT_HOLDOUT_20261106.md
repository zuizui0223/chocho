# Island51 circuit-theory steppingstone benchmark: distance baseline versus electrical resistance

**Source-locked 6 November 2026.** Research status: geographic network baseline completed; actual shipping/wind resistance raster NOT fitted.

## Decision

**A geometric steppingstone graph constructed solely from great-circle distances does NOT improve independent-island-out prediction over a nonlinear great-circle distance baseline.** This conclusion is about a **distance-only undirected graph**, not about AIS shipping, wind, marine currents or historical organism movement.

[Successful 51-fold benchmark source and logs](https://github.com/zuizui0223/chocho/actions/runs/42169528412).

## Source 51 landmasses / 1,275 dependent pairs

Frozen source under provenance/reviewer_defenses/results/island_host_network_20261030:
- island_host_network_islands.csv: 51 eligible original landmasses selected by existing candidate alien host records and earlier indexed butterfly observations.
- island_host_network_pairs.csv: all 1,275 source island pair records. Columns are island IDs, representative-port great-circle distance, candidate introduced-host Jaccard, butterfly Jaccard in 2001–2010 (earlier) and 2011–2022 newly indexed records (later).
- Previous source observations: host-versus-later butterfly Spearman +0.299806; host-versus-earlier butterfly +0.417826, and descriptive partial rank later controlling distance/earlier butterfly +0.215354.

**Data location gap:** The paired source table contains distances but no **absolute longitude/latitude** of the representative port positions, and landmass IDs are numerical source keys, not a georeferenced island catalog. Original raw port-to-landmass mapping and 2019 source port GeoJSON are named with checksums in island_host_network_audit.json but absent from this Git branch. Until exactly restored and verified, AIS or wind rasters CANNOT be sampled at the 51 true positions. Do not assign arbitrary island names or locations.

## Circuit-theoretic graph

A discrete 51-node graph implements the same effective-resistance Laplacian equations used in undirected **Circuitscape network mode**; the actual Julia package and any raster-based circuits were NOT run.

- Original source representative-port distance matrix, 51×51 from 1,275 pairs.
- Predeclared primary: undirected union of each node's **3 geographically closest islands**, plus minimum-spanning-tree bridges wherever disconnected; edge resistance proportional to kilometers.
- Primary topology: 95 k-nearest undirected connections plus three MST bridges = 98 edges. Maximum hypothetical edge **6,540 km**, hence NOT a shipping or feasible organism dispersal link.
- Laplacian eigen/pseudoinverse yields effective resistances. The complete test suite reproduces known series (two 1-ohm resistors in series -> 2-ohm) and parallel (opposite points of four-unit-resistor ring -> 1-ohm) circuit results.
- Secondary sensitivity parameters *fixed before observing response*: 2, 5, 8 neighbors with linear distance resistance; 3 neighbors with squared distance resistance.
- Predictor baseline: log(1+distance) and its square. For the **later butterfly** outcome, earlier-butterfly Jaccard also enters BOTH baseline and augmented model. Augmented model adds log(1+effective resistance).
- Strict held-out endpoint: each of 51 original islands left out with ALL its 50 pairs; training uses the other 1,225 pairs only. This gives 51 folds ×50 predictions, each original pair predicted twice withholding either endpoint. Fold losses are descriptive and *not* independent p-value replications.

## Source-verified primary results (51 island-omission folds)

| Target | Distance baseline MSE | + geometric circuit resistance MSE | MSE reduction (positive is better) | Island folds improved |
|---|---:|---:|---:|---:|
| Alien host Jaccard | 0.00371248 | 0.00371334 | **−0.023%** | 25/51 |
| Later butterfly Jaccard | 0.00671956 | 0.00673963 | **−0.299%** | 24/51 |
| Earlier butterfly Jaccard (negative comparator) | 0.01090385 | 0.01095206 | **−0.442%** | 23/51 |

Geometric resistance and the original great-circle distance themselves correlate (source Spearman **+0.769**), yet the extra predictor brings no incremental information to the primary fold tests.

## Five predetermined topology / resistance sensitivities (MSE reduction %)

| Graph | Host Jaccard | Later butterfly Jaccard | Earlier butterfly Jaccard |
|---|---:|---:|---:|
| k2 + spanning tree | −0.307 | −0.126 | −0.356 |
| **k3 + spanning tree (primary)** | **−0.023** | **−0.299** | **−0.442** |
| k5 + spanning tree | +0.035 | −0.489 | +0.018 |
| k8 + spanning tree | +0.287 | +0.504 | +0.157 |
| k3, squared edge distance | −0.204 | −0.479 | −0.418 |

The largest observed nominal improvement (0.50% for butterfly; 0.29% for hosts) depends on graph configuration and was not the primary. Do not tune k on this response then present it as a confirmed maritime or wind mechanism.

## Actual independent data required for ship-versus-wind model

1. **Ship traffic:** [World Bank/IMF Global Shipping Traffic Density](https://datacatalog.worldbank.org/search/dataset/0037580/global-shipping-traffic-density) (public CC BY 4.0) has a **commercial-ship AIS 0.005° raster ~458 MB ZIP** for Jan 2015–Feb 2021. These are ship AIS *observation density including stationary broadcasts*, not known plant carriage events, historic trade routes, or point-to-point flows. A coarse water-only raster needs tested aggregation and missing coverage gates; never use density as known bilateral ship freight.
2. **Atmospheric wind:** [NASA MERRA-2](https://registry.opendata.aws/nasa-merra-2/) requires Earthdata login for the AWS access path; [Copernicus ERA5 U/V climatology](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means) needs CDS registration/token plus licence acceptance for programmatic retrieval. Wind SPEED can enter an undirected circuit only as a directional information-losing proxy. Wind VECTOR transport must use a directed time-resolved advection model instead.
3. **Location source:** The published 51-node matrix has island distance but not georeferenced port positions. Recover the original indexed WPI GeoJSON and port→landmass crosswalk; validate every one of the 1,275 original pair distances BEFORE pixel sampling, because a distance matrix alone determines locations only up to global rigid rotation.
4. **Outcome validation:** The 51 islands were chosen using earlier butterfly records and candidate host data; later first GBIF indexing is not true colonization. Pooled 2001–2024 plant data can include plant observations recorded later than the butterfly records. Distinguish historical current-route proxy from a genuinely prospective cause.

## Manuscript decision

Preserve island +118.9% / mainland +40.9% potential host-resource gains, the port-number versus whole-introduced-flora and introduced host-taxon positive findings, and the +0.300 host–later butterfly assemblage comparison. **Do not replace JBI Fig. 4** with a supposed ship/wind mechanism on these geographic-only circuits. Actual separate AIS and wind resistance surfaces are required before any stronger transport claim.

Run source: scripts/exploratory_island51_circuit_network_holdout.py and tests/test_exploratory_island51_circuit_network_holdout.py. Archive artifact contains original SHA256 source hashes, all five resistance graph specifications, source graph diagnostics and complete per-island held-out prediction losses.
