# LGM: exact climate suitability versus ice-free physical-domain spatial accounting

**2026-10-20; source-verified exploratory ecological calculation.** Main question: How much of the modelled 0 BP → 21 ka change in potential native larval-host suitable *area* falls on **physically unchanged ice-free land**, and how much comes from modelled suitable area lost or gained when PMIP4 land/ice status itself changes?

This is an **exact spatial bookkeeping identity**, not a factorial causal identification of physical and climate main effects, not an occupied 21 ka host landscape, not a dispersal calculation and not a new genotype test. The earlier TTF alternative-plant negative specificity result (#149) is unchanged.

## Frozen inputs

- **416** original native-WCVP modern-trained host plant 0.5° CHELSA/PMIP4 potential habitat binary masks in exactly identical accepted taxon order at 0 BP and 21 ka: main source SHA256 `169b9da084d0043abba5c89f008137d6867b554dfb7c7dbe283b750142d220af`.
- Separate **416-source-model-verified** modern-nonanalog climate-extrapolation sensitivity masks (same original main habitat as bitwise subset): SHA256 `67c0d99852b8b8896559dbc0da07504bcb07a6ad7c017afa0eb86e31d40acc2a`.
- Source PMIP4 ICE-6G_C 0ka and21ka source SHA-bound 0.5° ice-free masks; source ZIP and coarsened NPZ as documented in the parent PRs; NPZ SHA256 `d8283a8ef177b56fd4554dc5bcc115e62adc5556e445bf5b7492c2c636a2a9db`.
- Exact TTF original 641 insect species, 1037 accepted plant host links, SHA256 `a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044`.
- **513/641** insects have complete models for all 1–3 *recorded* hosts, corresponding to **349 distinct known-host combinations**. The remaining 128 species are explicitly unscored, not backfilled.
- No refitting, selective geographic subset, retuning, novel climate values, occupancy or genetic response. Same 360×720 global 0.5° grid and latitude-weighted cell areas.

## Exact modelled-area decomposition

Let `P0,P21` be the sets of PMIP4 ice-free land pixels and `H0,H21` be original or extrapolated host climate-suitable binary *and* physically passable pixels. Let common passable domain `S=P0 ∩ P21`, `L=P0 \ P21` (formerly available) and `G=P21 \ P0` (newly available). With source spherical grid-cell areas `area(...)`,

**area(H21) − area(H0) = [area(H21∩S) − area(H0∩S)] + area(H21∩G) − area(H0∩L).**

This partitions the actual **source model output** into climate suitability turnover on land present in BOTH eras plus site gains/losses from a **changing eligible physical domain**.

The land itself: **48,830** grid cells remain ice-free at both ages, **9,935** are passable only at0ka, **8,356** passable only at21ka. These are spatial cells at 0.5° majority >=5/9 original PMIP4 10' land/ice resolution. Newly passable includes newly emerged paleo-shelf and ice-loss; no ancient plant occupation is inferred.

Every one of **416 original plants**, **513 eligible insect host unions**, **349 unique host combinations** and their separate nonanalog-extrapolated masks satisfied the exact per-row identity. **Maximum numerical remainder 2.05×10⁻⁸ km²**.

## Main native-climate-analog-censored response

| Separate source unit | Original 21ka area declines | Climate-suitable area declines on *unchanged* passable land | Net terrain-domain replacement negative | Both terms negative |
|---|---:|---:|---:|---:|
| 416 plant models | 289 | **305** | 220 | 188 |
| 513 insect host unions | **399** | **404** | **359** | **314** |
| 349 distinct host identity combinations | 258 | **263** | **223** | **189** |

For 513 insects, **species-level** medians:
- Original full suitable area change: **−4,921,313 km²**.
- Difference in suitable area restricted to *unchanged physically passable land*: **−2,607,874 km²**.
- Net new-passable suitable area minus lost-passable suitable area: **−1,595,139 km²**.
- Median 21ka suitable area on newly passable land **+2,018,649 km²** and median present-day suitable area on soon-to-be blocked land **3,732,071 km²**.

**Do not add these independently calculated medians together.** The exact identity holds separately for each species or distinct source host combination; median does not distribute over sums.

Species-level source pattern: 314/513 have negative common-land climate term AND negative terrain replacement term. For **90/513**, a negative common-land climate term is partially counteracted by newly available terrain (positive physical-domain term); for **45/513**, a positive common-land climate term is opposed by net terrain loss. **304/513** have greater absolute common-land climate term than the signed net physical-domain replacement term, but this is not causal factor importance.

## Previously prespecified finite-nonanalog climate-extrapolation sensitivity

| Unit | Full area declines | Common-land climate suitability declines | Net terrain replacement negative | Both negative |
|---|---:|---:|---:|---:|
| 416 plants | 243 | 255 | 190 | 145 |
| 513 insects | **322** | **318** | **291** | **208** |
| 349 unique host sets | 201 | 214 | 174 | 125 |

Including finite climates outside modelled contemporary native-plant training support changes the magnitude and prevalence of estimated potential suitable-area decline, which is precisely why this sensitivity is retained. It **does not validate** those extrapolated environmental conditions as genuinely suitable or occupied by plants in the LGM.

## Source and inference firewall

- New standalone source-locked script `scripts/land_vs_climate_common_domain_accounting_v07.py`; focused synthetic tests and source ZIP/file SHA and full 416/641/349 identity checks PASS.
- [GitHub Actions #39407170764](https://github.com/zuizui0223/chocho/actions/runs/39407170764) **SUCCESS**, including reproducible compact CSV, source receipt, all original host IDs retained and explicit 128 incomplete insect statuses.
- Results in `results/lgm_land_climate_common_domain_accounting_v07/`: `summary.json`, `all_416_plants_land_climate_partition_v07.csv`, `all_641_insects_land_climate_partition_v07.csv` and `349_unique_hostsets_land_climate_partition_v07.csv`.

**Scientific conclusion:** The original modelled 21ka suitable-area contraction is not accounted for solely by disappearing ice-free land. Climate model suitability changes on geography common to both eras are also modelled to be negative for many accepted host sets. However, the decomposition cannot identify a unique causal attribution: the modelled suitability that would occur at 21ka climate on land which was physically blocked at21ka (and vice versa) is not observed in these source binary masks. A true `2×2 climate epoch × terrain epoch` experiment requires retrieving/validating the **ungated source plant climate scores on those counterfactual blocked cells**, and cannot be claimed from this accounting alone.

Neither these potential climatic resource landscapes nor matched-modern-geometry patch splits demonstrate occupied LGM host-plant refugia or historical genetic differentiation. Do not treat the negative TTF v0.5 known-host versus modern-climate-matched alternative-plant result as positive.
