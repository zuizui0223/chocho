# Historical GBIF verification of all 317 butterfly first-index cells (14 November 2026)

## Executed decision

**SUCCESS:** [Original GSHHG 2.3.7 and GBIF 1700–2000 historical audit, six independent shards plus verified merge](https://github.com/zuizui0223/chocho/actions/runs/43585614802) (Actions run **43585614802**). This is an **index-history diagnosis, not a survey of actual butterfly colonization**. None of the old GEB resource or Circuitscape community-prediction results is overwritten.

The frozen original population is **225 butterflies × 51 landmass groups (11,475 cells)**. Under the prior narrow 2001–2010 GBIF reference, there were **221** apparent 2011–2017 first-index species–island cells and **96** apparent 2018–2022 first-index cells.

| Original subsequent-index interval | Cases indexed after the 2001–2010 negative window | Historical same-species, same-island GBIF indices in 1700–2000 | Pre-2001 GBIF nonindex (genuine historical presence unknown) |
|---|---:|---:|---:|
| 2011–2017 | 221 | **122 (55.2%)** | 99 |
| 2018–2022 | 96 | **23 (24.0%)** | 73 |
| **Total** | **317** | **145 (45.7%)** | **172** |

Across the separate *port-introduced-host-indexed* population, 2011–2017 had **12/15** pre-2001 records, while 2018–2022 had **4/9**. These 15/9 subsets should not be treated as independent validation: they reuse the selected original ports and host catalogue.

**Why this matters:** A contemporary GBIF record from 2011–2022 is not necessarily a first arrival. A gap of records in 2001–2010 does not erase historical records from before that decade. Positive pre-2001 same-species/same-landmass indices preclude treating a later first index as the first documented source occurrence under the expanded search. Pre-2001 GBIF *nonindex* still **does not prove true prior absence** or novel colonization.

## Original-source integrity

- Each of the **six** original butterfly source shards used source-locked GSHHG 2.3.7 polygons verified against SHA-256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`.
- The source query is GBIF occurrence search restricted to indexed records 1700–2000, PRESENT occurrence status, existing georeferenced records without reported geographic issues, original speciesKey and original GSHHG landmass polygons.
- All **91 island-polygon-fragment searches** succeeded across the 6 shards (0 errors). Per-island responses, geometry fragment hashes, original-case membership and query response digests are archived.
- The six merged source CSVs preserve **317 unique species × original-landmass × post-index-period identities**; the original 15 introduced-host cases reproduce exactly **12 pre-2001 positives**.
- The sample is selected through port/host/butterfly indexing; percentages are rates of earlier source indexes **within this selected retrospective source**, not global butterfly-colonization error estimates.

**Reproduce:** `scripts/screen_pre2001_all_island_first_indices_20261114.py` (six original geographic batches), then `scripts/merge_pre2001_all_island_first_indices_20261114.py` using the [frozen GitHub Actions workflow](../../.github/workflows/chocho-pre2001-all-317-island-first-index.yml). The source complete 317-row CSV and computed summary are in the [source-verified artifact](https://github.com/zuizui0223/chocho/actions/runs/43585614802/artifacts/14110156536).

## Inferential boundary

- This historical check evaluated **positive subsequent index cells only**. It is **not** yet a fair recalculation of a new establishment model, because candidate unindexed species–island cells were not all screened with the same older source window.
- Even after full historical GBIF screening, the source remains presence-only and effort-heterogeneous. A scientific claim of colonization requires independent natural history, breeding, or other source-supported arrival chronology.
- 1700–2000 specimens can themselves have taxonomic/date/georeferencing errors; original records should be validated if they drive a biological claim.
- The central GEB findings on **resource opportunity** (global +54.9%, insular +118.9%, mainland +40.9%), **port-introduced flora**, and **marine community similarity** (uniform sea field held-island MSE improvement +5.082%) concern different outcomes and remain intact.
- Historical source outcome validity is a **specific new methodological finding**, not evidence the host resource opportunity is absent.

**Next identification gate:** screen the full frozen **11,475** species–island risk cells under the same historical year window before attempting a uniform longitudinal prediction; separate *re-indexing probability* from actual establishment throughout.
