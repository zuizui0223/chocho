# Refit port-indexed introduced larval hosts × marine access after uniform 1700–2010 GBIF history

**14 November 2026 | Completed source-backed exploratory model | not an island-colonization or establishment analysis.**

## Interpretation first

The original 225-butterfly×51-landmass first-index cohort was source-conditioned on **no GBIF record in 2001–2010**. A full independent historical source check of **every 11,475 species×island cell** found records before 2001 for **122 of 221** supposed 2011–2017 first indexes and **23 of 96** supposed 2018–2022 first indexes. The revised index-only eligibility excludes those pre-2001-positive cells, and **also upgrades historical source availability when the same species was indexed on another sampled island**. This produces a consistent **1700–2010 GBIF-index risk set**. Even this longer reference is not proof of actual butterfly absence, first arrival, breeding or establishment.

**Result:** The proposed positive `historical indexed introduced host × modeled uniform-sea access` interaction was **not confirmed** in either revised period. The source-verified prospective *index order* is an improvement over the former arbitrarily narrow historical reference, but the number of new butterfly indexes in the actual introduced-host-recorded group is just **3 or 5**.

## Direct model data

The 2001–2010 original-port host map is source-locked, with known larval-host relationships and *contemporary regional WCVP* introduced/native status. The marine predictor is the residual of `log(1+uniform Circuitscape resistance)` against quadratic `log(1+distance)` across the 1,275 original island pairs; smaller-than-expected modeled resistance means greater access. The nearest known prior conspecific GBIF-index island is selected *among the same 51 sampled landmasses*, now considering 1700–2000 historical indices, 2001–2010 records and (for the later response) 2011–2017 records.

The fixed-effect linear probability model includes butterfly-species and destination-island effects, known earlier native and introduced host records, distance, sea-access residual, and the introduced-host × access term. Approximate 95% intervals are two-way clustered by butterfly species and island, are exploratory, and are fragile for sparse exposed events.

| Interval | Destination island botanical source | Indexed first-event risk cells | Positive new GBIF indices | Indexes with introduced candidate host | Host × sea-access difference per 1 SD (percentage points; approximate 95%) |
|---|---|---:|---:|---:|---|
| 2011–2017 | All historical port plant records present | 4,730 | 76 | 3 | −1.335 [−2.526, −0.144] |
| **2011–2017** | **≥10 pre-2011 port plant records** | **4,125** | **65** | **3** | **−1.249 [−2.530, +0.032]** |
| 2011–2017 | Same-island and same-species exposed/unindexed support | 2,835 | 42 | 3 | −1.130 [−2.419, +0.160] |
| 2018–2022 | All historical port plant records present | 4,854 | 61 | 5 | −0.481 [−1.842, +0.880] |
| **2018–2022** | **≥10 pre-2011 port plant records** | **4,235** | **49** | **5** | **−0.311 [−1.685, +1.063]** |
| 2018–2022 | Same-island and same-species exposed/unindexed support | 2,879 | 26 | 5 | +0.190 [−1.099, +1.478] |

These **six analyses are sensitivity populations within one underlying source, not six independently replicated tests**. The main ≥10-record panel is deliberately conservative about botanical survey effort but not a random matched effort design. The all-positive-host event count of 3 and 5 makes any interaction-based mechanism inference weak; negative point estimates cannot establish ecological harm of alien plants. Coefficients differ from previous models **partly because both historical source distances and risk-set membership changed**.

## Full original source and reproduction

- Original historical source index: **225 butterflies × 51 landmass groups = 11,475**, screened with verified GSHHG 2.3.7 and indexed 1700–2000 GBIF records, 92/92 queries successful in six original geographic shards. Source runs: [full risk reconstruction](https://github.com/zuizui0223/chocho/actions/runs/43586278558), [11,475-cell source table](https://github.com/zuizui0223/chocho/actions/runs/43586278558/artifacts/14110149052).
- Conditional refit: [Actions 43586514578 **SUCCESS**](https://github.com/zuizui0223/chocho/actions/runs/43586514578); [source-hash-locked output with full six model details](https://github.com/zuizui0223/chocho/actions/runs/43586514578/artifacts/14110592733).
- Code: `scripts/analyze_port_host_ocean_long_history_20261114.py` and workflow `.github/workflows/chocho-long-history-port-host-ocean.yml`.
- Compared with original 2001–2010 reference, the early/late source-eligible GBIF index events became **99/73** among all 11,475 cells. In the more selective model source there were 76/61 indexed events, and after ≥10 historical plant records **65/49**. The other cases were excluded by exposure support and/or no sampled earlier island with an indexed conspecific; do not conflate denominators.

## Retained positive findings and future identification boundary

**Keep these results separate and intact:** Global potential larval-host resource geography across 239 butterflies (+54.9%); larger relative gains on islands (+118.9%) than mainland (+40.9%), because islands had fewer native host-covered combinations at baseline; positive port number–introduced whole-flora and known-host candidate correlations (the host-specific residual after whole flora is not independently positive); and 51-island retrospective Circuitscape community similarity prediction improvement (+5.082% uniform sea geography alone, +6.488% AIS/wind beyond hosts under represented held-island folds, with adverse strict geographic transfer).

**Not shown:** ships carrying host plant seeds or butterfly eggs; adult wind-assisted transoceanic journeys; the source of a new island butterfly population; locally proven larval feeding on introduced plants; port-enabled establishment; onward island stepping-stone dispersal.

**GEB main paper remains unchanged.** This source-corrected test is a transparent **exploratory limit of the proposed mechanism**, preserving the genuine resource and island-geographic findings. A field or independent dated-history study of host utilization and verified multi-year recruitment is needed for a true establishment mechanism.
