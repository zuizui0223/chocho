# Native-screened host plants, 0 BP to 21 ka: first ecological climate-suitability evidence

**Date:** 2026-10-20. **Scientific status:** completed source-pinned exploratory climate-suitability hindcast; NOT a paleoplant occurrence map, time-specific dispersal network, larval host use probability or genetic structure test.

## Main progress and provenance

The biology code and results reside in **chocho**; TTF retains previously frozen inference/methods. This new analysis did **not** reopen the old TTF genetic response or revise the negative modern-climate matched alternative-plant control (TTF draft PR #149).

- Original WCVP native L3 geography (2024-05-21 version, rWCVPdata pinned commit 65bed76), 518 accepted host plant names/IDs and the original 641 insect-to-host associations.
- 296,697 source-bound modern plant observations checked against original native WCVP L3 regions; **220,137** had spatial positions inside documented native regions.
- Original pre-LGM source selection: 10 km presence thinning, cap200, ≥30 accepted native presences and ≥30 regionally native climate background sample sites; **416/518 taxa qualified**, with **73,386** selected plant observations.
- Frozen **4 CHELSA-TraCE21K variables**: bio01, bio07, bio12, bio15; 0 BP and 21 ka climate from the exact eight original SHA-bound raster files. 6,343 contemporary WCVP representative sites climate-finite and **6,337 common finite 0BP+21ka sites**.
- Pre-outcome model: quadratic degree2 logistic classifier, L2 C=0.3, liblinear, 1:1 class-balanced weights, native regional pseudo-background (NOT actual absence), one geography-blocked internal holdout per taxon, 10th percentile of fitted modern native presence scores to define *relative suitability*. 5th and 15th percentile threshold sensitivity saved.

**GitHub Actions model run:** [#39354040917](https://github.com/zuizui0223/chocho/actions/runs/39354040917) **SUCCESS**, 11 focused synthetic/source tests passing. Exact results artifact ID **12234343435**, full artifact digest SHA256 `8ac33d6f1c98b09fc057c73eaa3fb60278af89780b870864b22b44a9b494ecad`; `native_host_climate_surface_scores.npz` SHA256 `8559ec905cddd1a0b9579530551c9980ab3519daa92902df150a04687bf81159`.

## Results: first plant-level climate potential patterns

- **416/416 source-qualified plants fitted and projected**; the other source statuses remain 98 `INSUFFICIENT_NATIVE_OBSERVATIONS`, three `INSUFFICIENT_NATIVE_BACKGROUND` and one `NO_DOCUMENTED_NATIVE_WCVP`. They were not replaced.
- **366/416** modeled plants have a smaller fraction of the fixed reference sites passing the 0BP-fit suitability threshold under 21ka climate, **50/416** larger.
- **Median within-plant 21ka minus 0BP suitable-site fraction = -0.0743254** (−7.43 percentage points). The median of the *separate present site fraction* = **0.2081**, LGM = **0.1360**, which must not be subtracted as a species-paired estimator.
- Direction of decline remains under the frozen 5th or 15th percentile training-presence thresholds: **379/416** and **359/416** declined respectively. Under the fixed, per-plant 21ka-analogue subset, **365/416** declined, median paired change **−0.08659**. That support set differs among plants but is the same for each plant's paired time comparison.
- **Median share of 21ka sites outside the modern training climate predictor envelope: 0.4214** (42.1%). **148/416** exceed 50% nonanalog. Hence any global 21ka fraction is substantially extrapolative.
- Geography-blocked *presence versus opportunistic native-area background* AUC was computable for 406 plant models, median **0.794**. This tests sampling discrimination, not verified plant occupancy forecast, and cannot eliminate spatial observation bias.

The apparent majority contraction is **not automatically independent evidence of historical range contraction**: the classifier was trained on modern plants against a regional pseudo-background and then transferred to another climate period. The 6,337 support sites are not equal-area representation of the globe.

## Results: insect-specific aggregate potential host climate

The 641 frozen insect species each link to one to three **recorded** accepted plant host identities. Score actual modelled hosts with a binary OR over plant predicted 10th-percentile climate suitable points; **do not score any insect with an unmodelled accepted host**.

- **513/641** insects have models for *all recorded hosts*, 83 are partial (not scored), and 45 have no modelable known host. This is not independent evidence that any host catalogue is complete.
- Among 513 fully modelled, **478/513** have fewer potential climate-resource suitable sites at 21ka than 0 BP, **35/513** have more. **Median within-insect LGM−modern suitable-site fraction −0.0918416** (−9.18 percentage points), on the identical 6,337 fixed reference sites.
- Median separate suitable-site fractions: modern **0.24444**, 21ka **0.15039**. These separate medians differ from the paired median delta.
- The binary suitable-site sets rearrange: across insects, median **Jaccard overlap = 0.3411** between the fixed modern and 21ka sites. Every one of the 513 insects has *both* some sites gained and sites lost, even when the net site fraction declines. This is an **environmental suitability turnover**, not occupied species range relocation.
- Recorded host taxa per insect among 513: 323 with one, 125 with two, 65 with three. Median site fraction changes: −0.0841, −0.1013, −0.1215 respectively. A species with two/three *recorded* hosts is not necessarily a true generalist; the OR rule inflates the present climatic opportunity mechanically, so these are descriptive only.

**Source-verified linked analysis:** [#39355100124](https://github.com/zuizui0223/chocho/actions/runs/39355100124) **SUCCESS**; exact original WCVP plant host set, first climate model NPZ and climate input SHA hashes verified. Durable outputs:
`results/lgm_native_plant_hindcast_v01/summary.json`,
`results/lgm_native_plant_hindcast_v01/all_518_native_host_hindcast.csv`,
`results/lgm_insect_native_host_climate_opportunity_v01/summary.json`,
`results/lgm_insect_native_host_climate_opportunity_v01/all_641_insect_climatic_resource_opportunity.csv`,
plus source receipts. The full numeric 416×6,337 0BP and 21ka score matrix remains in the named original GitHub Actions artifact and is supplied in a separate download; it is not silently represented as a commit.

## The inferential boundary and genuinely next task

**The main question remains whether historical physical host-resource connectivity predicts independently measured insect genetic breaks across species.** No such genetic effect has been estimated in these data. The source-verified modern native WCVP polygon only establishes an area where a plant is currently listed as native; it does not establish the plant occupied that area at 21ka.

We can now construct *potential* 21ka suitability maps, but **must first assemble a valid global 21ka shoreline and ice-sheet/land barrier mask**, resolve sample-point connectivity versus 2D raster passability and model plausible taxon-specific dispersal before any claimed host resource route. The existing TTF ETOPO2022 Costa Rica native-tile pilot is **not** an authoritative global paleoland reconstruction. Explicit coastline/ice-sheet and historical occupancy uncertainty must enter path sensitivities.

Model qualification/fit and modern-to-21ka comparative suitability are complete; **historical plant occupation, plant-mediated gene flow and observed genetic barriers are not**. The prior TTF v0.5 modern-climate matched alternative-plant negative result stays negative and frozen.
