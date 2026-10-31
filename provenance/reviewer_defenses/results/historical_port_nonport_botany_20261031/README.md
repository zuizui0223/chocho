# Historical island ports versus non-port coasts — 2001–2010 botanical control (2026-10-31)

## Biological purpose and inference boundary

The butterfly-first hypothesis asks whether **maritime introduction of compatible larval-host plants at island ports** can facilitate future butterfly records after ocean crossings. A necessary but insufficient intermediate prediction is that *independently earlier indexed* suitable plants are concentrated around ports beyond comparable non-port coasts.

This analysis records botanical indices only. It **does not demonstrate that the plants were introduced to the target island, that ships transported them, that they were used by local larvae, or that butterflies established**. It uses ten globally important potential larval-host resource plants from the frozen GEB host-resource reconstruction, **not** a complete locally verified alien-host inventory.

## Source-verified completed analysis

- GSHHG 2.3.7 original high-resolution level-one island membership and 2019 World Port Index source IDs were used to fix the geometry before examining botanical results.
- 17 islands and all **38 historical port candidates**, with **2 independently selected off-port coastline points per port**, yielded **114 sites** (38 ports + 76 controls). WPI sites within 5 km of controls were excluded. Both site types were queried using equal **1 km** radii.
- **All 114 GBIF plant queries for 2001–2010 completed**, 70 sites had at least one resolved speciesKey. The off-port controls in this pipeline are newly sampled coastpoints, not a claimed recreation of the older WorldCover 2021 landscape controls.
- An initial 17-port pilot mistakenly interpreted the canonical GBIF facet field `SPECIES_KEY` as camel-case `speciesKey`, incorrectly producing zero species even when occurrence records existed. The code was repaired to normalize canonical uppercase response fields and reject populated responses lacking a species facet. The **corrected full38 source** replaces that false zero diagnostic.
- Original 38-port publication push encountered a fast-forward refusal after concurrent PR commits; a non-destructive git-fetch/rebase retry workflow successfully republished the complete source. No force push or source overwrite.

## Results: top-ten candidate host index fraction, matched port minus non-port

Each site's fraction is **number of top-ten resource plants GBIF-indexed / total plant speciesKeys GBIF-indexed**. For one port pair to be eligible, the harbor and **both** controls require at least the specified number of resolved plant species. Eligible ports are averaged within island before averaging islands (independent island weighting).

| Minimum plant speciesKeys at every matched site | Eligible ports | Eligible islands | Island-equal port minus off-port fraction | Island bootstrap 95% interval |
|---|---:|---:|---:|---:|
| 1 | 14 | 10 | +1.02 percentage points | −0.47 to +3.05 pp |
| 5 | 9 | 5 | +2.05 pp | −0.93 to +5.09 pp |
| 10 | 8 | 4 | +2.56 pp | −1.13 to +6.24 pp |
| 20 | 5 | 2 | −1.13 pp | Not estimable with fewer than 3 islands |
| 50 | 3 | 1 | +0.65 pp | Not estimable |

At the 1-species cutoff only **3 of 10 island means are positive, one is negative, and six have exact zero index difference**. At the 10-species cutoff 3/4 island means are positive but the uncertainty is broad and still crosses zero. Restricting the weakest botanical-index opportunities changes both support and direction; there is **no geographically robust, independently attributable harbor enrichment for these ten candidate plants**.

The interval is a deterministic **5,000-draw island-level bootstrap of fixed source-derived island means**. It does not include GBIF reindexing, botanical detection/process error, shared region/phylogenetic structure or retrospective site-selection uncertainty. It is not a causal confidence interval. Using the low bits of an LCG would cause artificial degenerate resamples at 4 islands; the source audit instead uses high-bit scaled selection and explicitly tests against that defect.

## Reproduction and tests

- Matched input and geometry: `fixed_geometric_sites.csv`, `geometry_admissions.csv`
- All public aggregated botanical source facets: `historical_plant_site_facets.csv` (no original GBIF occurrence IDs, no private butterfly–host edges)
- Original contrasts and counts: `port_offport_historical_pairs.csv`, `historical_port_offport_results.json`
- Raw acquisition/reconstruction: `acquire_matched_botany.py` with SHA-pinned GSHHG 2.3.7
- Fixed post-processing and 5,000 island bootstrap: `audit_island_host_index_uncertainty.py`, `historical_38port_island_uncertainty.json`
- Offline regression tests: `tests/test_historical_port_nonport_botany_20261031.py`, `tests/test_historical_38port_island_uncertainty_20261031.py`

Main next test remains a **butterfly species × island × year at-risk panel** joined to actually prior, regionally alien-classified host observations, plus a non-port urban control and independent mobility data. A species' presence in contemporary GEB worldwide host-resource expansion cannot replace that historical field.