# Source-uniform historical GBIF index panel: 225 butterflies × 51 island landmasses

**Execution date: 14 November 2026. Complete, verified source; not an ecological colonization dataset.**

## Result

[GitHub Actions run 43586278558](https://github.com/zuizui0223/chocho/actions/runs/43586278558) successfully ran **six original-geometry shards (92 GSHHG polygon-fragment GBIF searches, zero failures)**, merged the 225 butterflies × 51 selected GSHHG 2.3.7 landmasses and checked the original 11,475 unique species×island identities. The same historical **1700–2000 indexed-GBIF record query** was applied to **all** 11,475 risk cells, rather than only subsequent positive first-index cases. [Source-locked full 11,475-row audit and summary artifact](https://github.com/zuizui0223/chocho/actions/runs/43586278558/artifacts/14110149052).

| Source index period | Original 2001–2010 reference | After source-uniform 1700–2010 GBIF reference |
|---|---:|---:|
| 2011–2017 apparent first index | **221** | **99** |
| 2018–2022 apparent first index | **96** | **73** |
| Eligible non-index species×island cells for 2011–2017 | 10,739 | 10,437 |
| Eligible non-index species×island cells for 2018–2022 | 10,518 | 10,338 |

Across all 11,475 combinations, **894 had at least one 1700–2000 GBIF index**. Of the original 317 later first-index events, **145 (45.7%)** were already indexed in this older period (122/221 early, 23/96 late). The earlier independent positive-outcome-only 317-case audit gave identical tallies (no observed same-day source-index drift). All geography and GBIF source response SHA-256 receipts remain in the source batch artifacts; the complete merged source table and model are archived separately. Original GSHHG ZIP SHA-256 was exactly **8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41**.

**Inference:** A 2011–2022 first index relative to an arbitrary 2001–2010 reference is demonstrably **not first known indexed presence** for many butterfly×island cells. The longer historical GBIF reference is temporally more consistent, but a non-index (even before 2011) is **not a verified biological absence or new colonization**.

## Reproducibility

- Source-acquisition code: `scripts/screen_pre2001_all_11475_risk_cells_20261114.py`. Fixed original six butterfly-family taxon keys plus complete `speciesKey` facet pagination were necessary. A direct 225-species single-URL query returned HTTP 400 on the Santa Catalina smoke test; this failure was isolated and corrected, then a species-key identity source test passed. The final full-source run used the same successful family-based query as the original 225×51 GBIF pipeline.
- Merge validation: `scripts/merge_pre2001_all_11475_risk_cells_20261114.py`.
- Source-locked workflow: `.github/workflows/chocho-pre2001-full11475-eligibility.yml`.
- The six source shards, original 11,475 table, source hashes, full period-specific risk flags and source limitations are retained under [the finished workflow](https://github.com/zuizui0223/chocho/actions/runs/43586278558).
- Independent model re-analysis: [the source-consistent port-host × sea-access audit](../port_host_ocean_long_history_20261114/README.md), [run 43586514578](https://github.com/zuizui0223/chocho/actions/runs/43586514578).

## Boundaries for the integrated GEB paper

The positive **globally reconstructed host-resource opportunity** results (239 butterflies, +54.9%; island relative +118.9%, mainland +40.9%), the **whole non-native flora at ports**, and the **retrospective 51-island Circuitscape group prediction** (+5.082% relative held-island MSE improvement) concern *different biological outcomes* and are **not invalidated** by this independent historical observation-index test.

However, **the term first GBIF index must not be read as new first island arrival, successful breeding, or dispersal relay**. The original older-window 2001–2010 reference selected which 2011–2022 case records appear novel; in the new same-source risk table all eligible non-index observations are source-defined and remain presence-only. Neither the 99 nor 73 revised source-index outcomes are established colonization events.

**Scientific gate:** Genuine introduction or establishment requires independently dated island occurrence history and reproductive evidence. Port flora may describe an opportunity, but historical alien status and actually utilized larval host identity are still unverified in most cases.
