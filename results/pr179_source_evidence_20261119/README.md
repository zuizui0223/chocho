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
