# 2026-10-28 — Is the apparent introduced-host climate signal replicated within butterflies?

**Verdict: source-verified descriptive climate-recording contrast; NOT realized butterfly niche expansion, host-use benefit, or causation.** This is a **post-outcome exploratory diagnostic** of the original PR #74 same-WorldClim result, not a retrospectively preregistered test. Keep original GEB manuscript and `main` unchanged.

## Input / provenance

- Reused **all 73 original positive butterfly species × country cases**, **20 species**, exactly **6,076 accepted GBIF point-index records**, with original source-country climate classes: **32 climate-novel-dominant** and **41 climate-analogue-dominant**. No cases or points were selected on this new outcome.
- Source [successful WorldClim 2.1 1970–2000 analysis run #40435313884](https://github.com/zuizui0223/chocho/actions/runs/40435313884) and aggregate-only artifact **12703686785** (artifact ZIP SHA256 `a91cfd6ae4a1ef4c17c1a6e35b0d40eda078d312a69aa67cc6c1a17d164f9fa0`; internal JSON SHA256 `5d3ca4137821e60aaca7e1bf965c509b1cddee76c643fb00ee15df44ae775937`). The host baseline and butterfly points use **the same WorldClim product**, BIO1 and BIO12, with a geographical native-larval-host region reference. The source country labels were based on a separate CHELSA host-climate reference.
- Retained a minimal **73-row aggregate-only CSV**, SHA256 `e5d229bea32f93055d486d59a174bcca73ac8c7d6638e253d068238e5af764dc`. No exact GBIF points, observer identity, raw botanical source polygon, or climate-pixel data enter this branch.
- All 73 case sizes and climate class outcomes are exactly replayed from the original successful source run; 18 positive cases have **one** butterfly index record. The two largest cases (*Pseudozizeeria maha* and *Catopsilia pomona*, both Taiwan, climate-analogue class) contribute **4,418/6,076 = 72.7%** of point indices. They are GBIF records, not independent individual insects.

## Estimand and result

For each species represented in **both** climate-novel and climate-analogue source-country classes, first calculate the fraction of point records outside that species' **native-host-climate 5–95% rectangle** for each species × country. Average country fractions **equally within each species × class**, then take *novel − analogue*. Finally give each paired species equal weight. This removes the most obvious pooled record-count domination, but **does not match geography or observer effort** and does not solve conditioning on countries with observed butterflies.

| Sensitivity | Retained species × country cases | Paired species | Novel − analogue, equal-species average | Species direction +/−/= |
|---|---:|---:|---:|---:|
| Primary posthoc, area-weighted 5–95% source reference, >=1 record/case | 73 | **10** | **+42.42 percentage points** | 7 / 1 / 2 |
| At least 3 butterfly records per case | 47 | **2** | +33.05 pp | 2 / 0 / 0 |
| At least 5 butterfly records per case | 39 | **2** | +33.81 pp | 2 / 0 / 0 |
| At least 10 butterfly records per case | 30 | **2** | +39.24 pp | 2 / 0 / 0 |
| Narrower 10–90% climate reference, >=1 | 73 | 10 | +63.45 pp | 9 / 1 / 0 |
| Equally weighted rather than area-weighted 5–95% climate pixels, >=1 | 73 | 10 | +46.39 pp | 7 / 1 / 2 |
| Documented coordinate uncertainty <=1 km, >=1 record/case | 64 | **9** | +36.57 pp | 5 / 1 / 3 |
| Documented coordinate uncertainty <=1 km, >=5 records/case | 30 | **2** | +47.28 pp | 2 / 0 / 0 |

A full descriptive species-resampling interval for the first line is **+15.41 to +69.26 pp**, but it is **not** a valid causal or population confidence interval; species/climates are not randomly sampled, their country-class assignment is constructed from plant climate geography, and multiple cases share plant sources and geographical structure. **Do not report an exchangeable-label permutation p-value.**

Of ten nominally paired species, only **three** have at least **two different countries in each** class before applying record-number restrictions; the +100 pp species contrasts partly depend on **single-record novel-country cases**. At >=5 recorded points per species×country, only *Mechanitis lysimnia* and *Zizula hylax* retain both class comparisons. Two species cannot support a general comparative ecological effect, even if both show the same direction.

## Fundamental inference limits (not solved by statistical significance)

1. **Reference != butterfly physiology:** exceeding the mapped **native host plant** geographical climate support is NOT exceeding the butterfly's own native breeding-climate niche. A butterfly may use a local native host omitted by WCVP, be an immigrant adult, or occur in a cultivated garden outside the mapped wild plants.
2. **Geography is partially built into the classification:** the host-country source class is defined using climate availability, so a difference in *point climate* is compatible with trivial geographical sorting. A same-species, locally effort-matched null is required before treating the contrast as interaction-mediated niche release.
3. **Observation is not establishment:** the 73 cases were drawn from **GBIF-positive** species × country cells; they supply neither credible absences nor population denominators. Original source metadata found **0 juvenile/immature-index labels in climate-novel cases** versus 87 in analogue cases, which is missing evidence, not proof of no recruitment.
4. **Historical mechanism untested:** the point dates (2010–2025) are not introduction dates, and annual BIO1/BIO12 do not measure the feeding season or plant leaf persistence. The analysis does not infer additional generations, life-cycle completion, increased fitness or competition.

## Scientific gate and next independent test

**H1, host resource opportunity**, is established by the separate frozen GEB reconstruction. **H2, butterflies *occupy climate conditions unavailable on native hosts***, remains only exploratory positive-point co-occurrence. **H3, an introduced host supports breeding and new seasonal/climate niches**, is **NOT IDENTIFIED** by this dataset.

The next study should predefine candidate butterflies **before** selecting locations on positive records, build geographically and temporally comparable wild-butterfly effort and **non-host botanical** controls, map native vs introduced hosts at fine locality scales, and validate adult, egg/larva, multi-year recruitment separately. The most discriminating contrast is: **does an introduced-host locality exhibit reproducible larval recruitment outside the butterfly's *own native breeding-climate niche* beyond an effort-matched, other-plant geographic null?** A secondary phenology endpoint should compare host leaf availability and larval generation timing; monthly adult records alone cannot establish multivoltinism.

Recompute the aggregate-only diagnostic (offline, no GBIF/network):

```bash
python scripts/analyze_butterfly_73_paired_point_climate_20261028.py \
  --output /tmp/butterfly_73_paired_point_climate_20261028.json
pytest -q tests/test_butterfly_73_paired_point_climate_20261028.py
```

Machine-readable output: `provenance/reviewer_defenses/results/butterfly_73_paired_point_climate_20261028.json`. The analysis and claim firewall are separate from GEB's frozen submission results.
