# PR #179: positive island–port–butterfly evidence on main (2026-11-19)

**Scope:** Source-only archive from [Draft PR #179](https://github.com/zuizui0223/chocho/pull/179). The unreviewed GEB article, main figures, image permissions and coauthor approvals are **not merged**. This source archive is not manuscript submission or causal transport identification.

## Source-qualified positives

| Dataset / estimator | Positive result | Limitation |
|---|---|---|
| Frozen 239 known-host butterfly species × 369 WCVP botanical regions | +14,553 potential consumer-resource cells (+54.9% globally); island +118.9%, mainland +40.9%. Native-only host coverage island 16.0% vs mainland 37.5%; contrast 21.4→17.7 percentage points | HOSTS×WCVP potential access; not butterfly abundance or new colonization |
| 126 matched botanical islands, contemporary whole-flora association with WPI ports | Introduced richness β +0.2854 (HC3 +0.0622,+0.5087) controlling botanical island area and TDWG Level1. Native richness β +0.0954 (−0.0774,+0.2683) | Direct introduced−native β +0.1900 (−0.0187,+0.3987), **uncertain**; port presence is not historical shipment |
| Port geographic robustness | Nine overlapping leave-one-Level1-out introduced β +0.2016 to +0.3797; 91 port-positive islands β +0.3110 (+0.0416,+0.5804) | Overlapping re-fits, not independent worldwide replication |
| 104 botanical island pairs: archived uniform-ocean circuit | Introduced plant-composition held-island MSE gain +4.239% vs native −0.110% | Model-defined undirected water resistance, not demonstrated currents/ships |
| 104 island plant-flora Jaccard → 239-species potential butterfly-resource Jaccard, 1,088 doubly-unseen dyads | E model: +26.363% MSE gain; retrospective S (target current-resource richness included): +33.768% | Shared WCVP/HOSTS input; **not** independent butterfly community geography |
| Constrained plant/resource gain-margin null (225 source-conditional states) | Excess above null median +7.389 pp (E) and +9.261 pp (S) | Exploratory source-conditional structure, not new biological experiments |
| Independent 155 known-host butterfly species later-index record prediction | Host geography improves GBIF-indexed MSE +4.12%; spatially blocked +3.00% | Indexed occurrences, not establishment |
| Independently reported Canary-endemic larvae | *Pieris cheiranthi* larvae used introduced *Tropaeolum majus* in Tenerife wild and La Palma garden | One direct use endpoint; no measured population benefit. Original 239 butterfly panel includes 0 of 18 endemic species on independent 41-species Macaronesian list |

## Adverse checks that must accompany the positives

- Strict 104-botanical-island two-endpoints-unseen **circuit**: alien +2.735% (descriptive fivefold −1.649,+5.508), native −0.044% (−0.984,+0.431). Alien/native contrast interval crosses 0. Different shortest-water-path model improves **both** native (+2.234%) and introduced (+1.775%) flora.
- Zero/nonzero Jaccard and overlap-index sensitivity: native zero Jaccard 16.51% vs introduced 9.79%; Sørensen and Simpson can reverse origin differences. Ocean-circuit positive is not universal across metrics.
- Potential butterfly-resource composition: island-island ΔJaccard +0.2367 vs mainland-mainland +0.1504; Level1×islandness resource-margin null leaves only +0.00119 residual, limiting independent island-specific ecological sorting.
- Independent 224-species GBIF community-composition increment from host geography is only +0.157%. Selected 51-island network retrospective +5.082% circuit and +6.488% AIS+wind are separate from stricter 990-pair double-unseen AIS +0.6776% and AIS+wind −4.6986%; do not subtract across non-identical baselines.
- Tonga 2011 native broad water dispersal 137/364 vs post-European alien 137/309 (native−alien −6.70 percentage points). It is **not** equivalent to EuDiS seawater thalassochory and does not establish a universal trait contrast.
- Neither alien-flora/port association, potential HOSTS resource availability, nor electrical sea resistance identifies historical port deliveries, seed establishment, butterfly colonization, abundance increases or competition.

## Provenance and next test

The linked `results/` directories preserve source-pinned protocols, JSON results and adverse diagnostics; matching `scripts/` and workflow specs are transferred without original third-party species-level raw rows. The original draft article and figures remain solely in PR #179 for review.

**Next identification gate (completed on PR #179; see verified result below):** reconstruct original **all 5,355** WCVP flora pairs exactly from accepted plant names, then evaluate **5,353** resource-defined dyads with **1,088** strict double-unseen test pairs. Remove plants in the frozen 239-butterfly HOSTS catalogue from introduced WCVP flora, forming `H` registered-host flora and `O` non-listed introduced flora. Compare whole `W`, `H`, `O` and their fixed-ridge increments, with empty-union indicators and both E/S baselines. **O is not proven non-host flora.** Even a positive H-only increment remains structurally source-linked.

**Unresolved publication gates:** independent GBIF provenance/DOI, tropical 2020 raw species crosswalk, media rights, independent endemic coverage and coauthor approval. Main source archive ≠ manuscript sign-off.


## 19 Nov 2026: HOSTS-registered vs non-registered introduced plant negative control (completed)

Original **5,355 matched WCVP coastal flora pairs all source-reproduced exactly**; 5,353 resource-defined pairs, five frozen spatial folds, **1,088 doubly-unseen** test dyads. Of 7,887 distinct alien accepted WCVP plant IDs on 104 islands, **668** are known hosts of the original fixed 239 HOSTS butterflies and **7,219** are not listed. Full accepted source masks contain 7,726 host vs 29,439 other alien species×region units; no undefined H/O union in model dyads.

| Conditional introduced-flora Jaccard predictor | E MSE improvement vs common base | S retrospective MSE improvement |
|---|---:|---:|
| W: whole introduced flora | +26.363% | +33.768% |
| **H: registered known-host introduced flora** | **+33.821%** | **+41.804%** |
| O: non-listed introduced flora | +19.539% | +26.112% |
| H+O | +33.707% | +41.584% |

After O, adding H improves prediction **+17.609%** [fivefold resampled descriptive 95% +8.539,+28.224] under E, **+20.939%** [+11.121,+31.503] under S. After H, adding O instead gives **−0.173%** [−0.481,+0.058] under E and **−0.379%** [−0.862,+0.133] under S. W→W+H is +9.740% E / +11.672% S, but W→W+O also improves (+7.583% E / +6.926% S). H alone improves over common base in all 5 E and all 5 S held-out folds; O helps in 4/5 E and 5/5 S folds.

**Interpretation:** host-name masking carries more information for a *resource response constructed using those very host names*. It discriminates potential-resource arithmetic from a completely nonspecific whole-flora correlation, but **does not identify independent host specificity, use, survival or colonization**. Shared WCVP/HOSTS source records and richness/plant-family structure remain possible explanations; the S baseline includes current target-derived richness. The five source-fold resampling intervals are descriptive, not iid p-values. A genuine consumer validation still needs independent occurrence/larval performance data. Existing independent global 224-butterfly assemblage host-geography addition is small (+0.157%).

Source receipt: [successful host-mask CI #44445904218](https://github.com/zuizui0223/chocho/actions/runs/44445904218), [full source result artifact #14493153837](https://github.com/zuizui0223/chocho/actions/runs/44445904218/artifacts/14493153837). Exact artifact JSON SHA256: `3d9531806d25c0d08022ec63b99b7c7224d89303383139abe964b31aa85ceb62`. [Condensed result summary](../island104_host_masked_alien_flora_20261119/VERIFIED_HOST_MASKED_SUMMARY_20261119.json).

## 19 Nov 2026: family-occupancy-matched pseudo-HOSTS null (199 randomisations)

Follow-up to the verified H vs O source prediction. **Original 5,355 WCVP flora dyads exactly matched**, 5,353 potential-resource dyads, same 1,088 double-unseen test pairs. Matched 668 known-host alien plants against WCVP accepted **family** and **104-island occupancy bin**, preserving the host count in each bin and conditioning on each island's H/O richness. Of 668 host taxa, 614 can exchange labels and 54 are fixed in host-only strata.

**E (primary):** observed addition of H after O lowers relative held-fold MSE **+16.365%**, versus matched-label null median **+7.699%** (central 95% range +1.708 to +17.597%). **9/199** null scores are at least as large, descriptive add-one upper rank **0.050**, and observed value is **inside** the null central range. H-alone improvement versus the matched null has rank 0.145. **Do not claim robust identified host-specific ecology.**

**S (retrospective only):** H given O +21.647% versus null median +10.380% (central range +3.366 to +19.256%), **1/199** null scores at least as large, descriptive rank **0.010**. S includes contemporary potential-resource richness derived from the response and cannot substantiate prospective or independent ecological effects.

All outputs remain **same-source WCVP/HOSTS potential-consumer-resource models**; the matched null tests source-label structural specificity, *not* field feeding, butterfly occupancy, recruitment, ships or colonization. Independent GBIF observed community H-given-O improvement remains +0.091% in a different 224-species/169-region cohort. Both positive and adverse source results are retained without numerically subtracting disparate endpoints.

[Full qualified null report](../island104_host_family_degree_null_20261119/README_JA.md) · [Fixed exploratory protocol](../island104_host_family_degree_null_20261119/PROTOCOL.md) · [Numeric main receipt](../island104_host_family_degree_null_20261119/VERIFIED_FAMILY_DEGREE_NULL_20261119.json) · [Source-success workflow #44448741452](https://github.com/zuizui0223/chocho/actions/runs/44448741452).

## 19 Nov 2026: independent 224-species observed GBIF geographic transfer uncertainty

Source-checked the original **169/194/133** independent GBIF region-held-out model scores (M2=other introduced plants, M3=M2+known larval-host plants). In the primary original **169-region** cohort, H after O improves region-equal MSE **+0.0913%**, but regional bootstrap 4,999-draw descriptive central 95% range is **−0.0915 to +0.2621%** (world-stratified −0.0764 to +0.2474%; 96/169 regions favor H). Relaxed 194-region cohort: **+0.3325%** (descriptive +0.0712 to +0.5770); stringent 133-region cohort: **−0.0764%** (−0.1554 to −0.0007). The 3 eligibility cohorts overlap and differ in region composition; these are NOT independent experiments or detection-corrected estimates. Leave-world-out summary is from original regional scores, not retraining. Do not convert reused CV units into iid significance. [Protocol/report](../gbif224_region_uncertainty_20261119/README_JA.md) · [numeric receipt](../gbif224_region_uncertainty_20261119/VERIFIED_GBIF224_REGION_UNCERTAINTY_20261119.json) · [source-exact CI](https://github.com/zuizui0223/chocho/actions/runs/44452349265).

**New source location:** Source-locked annual GBIF species×WGSRPD3×year counts and region-year effort tables, missing from main, are present at frozen branch commit `dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0`. They permit the next explicitly approximate temporal dataset-day effort-thinning sensitivity without re-downloading the 14.6M occurrence archive.

## 19 Nov 2026: original GBIF annual butterfly community recording-intensity controls

Recovered source-original GBIF 224-butterfly species×WGSRPD3×year positive records and dataset×year effort from archived immutable analysis commit `dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0`. Reverified full source SHA, 91,882 species×region×year cells, 5,390 year-region effort cells, 224 indexed species, 169 original primary observational regions, and 14,196 regional dyads. Direct raw GBIF Jaccard exactly reproduces earlier source prior **0.114213** and late **0.167880** (difference **+0.053667**).

In 499 simulated observation-intensity-matched repetitions, after thinning each region's indexed records to the lower total **datasetKey×calendar-day** budget, late-minus-prior Jaccard difference median falls to **+0.015064** (conditional Monte Carlo range +0.013041 to +0.017291, **71.9% smaller raw contrast**). If thinning to each region's smaller **indexed occurrence count** budget instead, median difference falls to **+0.002914** (+0.000217 to +0.005503, **94.6% smaller**). Earlier/later original observation intensities differ in 168/169 regions under both proxies. This is a source-conditional **assumption-sensitive recording intensity diagnostic**, not an empirical detection-corrected biodiversity signal, independent replication, or 95% population interval. The two period lengths, human observer distribution, and species-specific detection are not standardized; no host-flora transfer or colonization model is refitted here.

**Interpretation:** Raw observed GBIF butterfly-assemblage similarity increase is highly sensitive to uneven recording effort under explicit random-thinning assumptions. Do not label the original raw delta direct biological homogenisation or link it causally to shipping/host plants. The separate +54.9% GLOBAL potential HOSTS×WCVP resource gain and positive port–alien flora association remain unchanged.

[Successful source-matched Actions run #44453480218](https://github.com/zuizui0223/chocho/actions/runs/44453480218) · [Full numeric JSON artifact #14496313483](https://github.com/zuizui0223/chocho/actions/runs/44453480218/artifacts/14496313483) · [Qualified main report](../gbif224_annual_effort_thinning_20261119/README_JA.md) · [Summary receipt](../gbif224_annual_effort_thinning_20261119/VERIFIED_ANNUAL_EFFORT_THINNING_20261119.json).

## 19 Nov 2026: independent observed butterfly H-after-O re-fitted on effort-thinned populations

**Original four GloNAF/HOSTS/WCVP+GBIF models and all 169 region held-out MSEs reproduced exactly from pinned original source commit dfd49cfcf9fb01a3db47cc611aae8cdd9c53cda0**. Unchanged 224 indexed butterflies, original 169 botanical WGSRPD3 regions and 14,196 nonindependent region dyads, original early and late observed GBIF Jaccards. Freeze original geography, contemporary botanical predictors, environment, colonial history, observation effort covariates and training-only ridge lambda 10; only substitute source-intensity-thinned **prior butterfly Jaccard feature** and **late observed butterfly Jaccard outcome** together on each draw.

Original independent **H after O** H-vs-unlisted-alien prediction gain is +0.0912908% held-region mean MSE. Under 29 deterministic draws matched on within-region prior/later **datasetKey×calendar-day** proxy, gain median is **+0.214956%**, conditional Monte Carlo central 95% **+0.065888 to +0.470555%**; 29/29 draws positive, 102/169 heldout regions improve at the median. Under 29 **indexed butterfly record-count** proxy draws, gain median **+0.437833%**, conditional range **+0.154318 to +0.736163%**, 29/29 positive, 104/169 regions improve at median. No pair union was undefined in the source reconstructions or simulations. The thinned late−prior Jaccard mean differences remain small (+0.015407/+0.003008 in these 29 draws), consistent with earlier 499-draw effort-intensity diagnostics.

**Scientific interpretation:** A small positive *independently recorded community prediction* increment from known larval-host-alien geography survives both explicit intensity-matching assumptions. This improves on a purely WCVP/HOSTS construction-only response, but **does not identify true occupancy, population size, causal colonisation, historical plant introduction or feeding quality**. The source original late-period effort/record covariates are held fixed while the simulated community outcome changes; the source is **retrospective**, with no standardized visits or species-specific detection model. Neither 29 Monte Carlo draws nor their percentile ranges constitute independent biological replications or population confidence intervals. Do not infer from 0.2–0.4% MSE gains that the much larger +54.9% source *potential* resource increase is realized demographically. Preserve original ports, island native/alien floras and endemic direct larval-feeding positive evidence.

[Reproducible complete original-source CI #44456355685](https://github.com/zuizui0223/chocho/actions/runs/44456355685) · [All 58 model-score repeats in original artifact](https://github.com/zuizui0223/chocho/actions/runs/44456355685/artifacts/14498146982) · [Human-readable interpretation](../gbif224_thinned_HO_transfer_20261119/README_JA.md) · [Numerical receipt](../gbif224_thinned_HO_transfer_20261119/VERIFIED_THINNED_HO_TRANSFER_20261119.json).

## 19 Nov 2026: stricter across-biogeographic-world butterfly host transfer

The 224-species independently GBIF-recorded community H-after-O prediction test now also uses **entire TDWG Level1 groups** as unseen geographic units rather than only one held-out region. The **original eight world-holdout M2/M3 source scores exactly reproduced**, as did the four original source models and their 169 region-holdout fits. Using the **same 29×2 effort-thinning simulations, same 169 regions and same 14,196 pairs**, equal-weighting the eight world MSEs gives an H-after-O improvement of **+0.6945%** (datasetKey×calendar-day proxy; conditional 29-draw range +0.1470 to +1.1570) and **+0.7793%** (record-count proxy; +0.4355 to +1.6183), positive in 29/29 simulations for each proxy. Original unthinned strict world-equal benefit was +0.0971%.

**Important heterogeneity retained:** each simulation improves only **5–7 of the 8** world groups (median 6/8), including a **negative South America (TDWG Level1=8)** result. Its H-after-O MSE increment is **−2.2517%** under dataset-day matching (0/29 positive; conditional range −4.7596 to −1.0282) and **−2.0700%** under record-count matching (1/29 positive; −4.1587 to +0.0465); the original raw source was −2.9806%. Level1=6 comprises only two original regions and is unstable. Do not claim a universal region-invariant ecological effect because the overall average is positive.

**Interpretation:** accepted larval-host plant geography provides a *small*, directionally positive aggregate predictor of independently recorded butterfly community geography, sometimes surviving transfer to unseen biogeographic regions under strong, explicit **observation-intensity simulation** assumptions. This is **not** a causal or detection-corrected effect on true butterfly occupancy, colonization, populations or larval host quality. All 29 draws share one observational source, and cross-world dyads appear in more than one world test. Source positive +54.9% potential resource supply, port–alien plant flora associations and independent endemic larval feeding remain separately preserved.

[Source exact CI #44461719514](https://github.com/zuizui0223/chocho/actions/runs/44461719514) · [all 58×8 signed held-world scores](https://github.com/zuizui0223/chocho/actions/runs/44461719514/artifacts/14500267175) · [main geographic-transfer report](../gbif224_thinned_world_transfer_20261119/README_JA.md) · [numeric summary](../gbif224_thinned_world_transfer_20261119/VERIFIED_THINNED_WORLD_TRANSFER_20261119.json).
