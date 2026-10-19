# Original-source butterfly life-stage readout — 19 October 2026

**Status: ORIGINAL SOURCE VERIFIED; PRIOR-ART REANALYSIS, NOT NEW FIELD ECOLOGY; NO CROSS-SPECIES RECRUITMENT LAW DEMONSTRATED.** Authors' original CC0 data were obtained directly from Zenodo's published Dryad mirror, not simulated, all five original MD5 and file sizes verified in successful [Actions run 39202355562](https://github.com/zuizui0223/chocho/actions/runs/39202355562), which created the [full machine-readable JSON receipt artifact 12164785204](https://github.com/zuizui0223/chocho/actions/runs/39202355562/artifacts/12164785204). This document is a scientific interpretation of a reproducible original-data audit on separate `analysis/cross-butterfly-lifestage-raw-20261019` branch. It does not change PR #38 (frozen GEB paper) or original experimental publications.

## Pre-specified question

A plant that provides an *apparently* geographically new host does not necessarily generate viable butterflies. Can genuinely measured performance on host plants be linked across butterfly species to environmental context, rather than treating potential host-geography unions as a recruitment outcome? **What actual component of the sign-reversal hypothesis is identified** by existing public butterfly cohorts?

A total of **3,458 original individual-larva records in three previously published experimental datasets representing four butterfly species** were audited:

| Original study | Butterfly | Initial larval records | Endpoint | Status |
| --- | --- | ---: | --- | --- |
| Saldivar & Wilson-Rankin, doi:10.5061/dryad.wdbrv15v5 | *Vanessa cardui* | 1,154 on five plants; 530 diet-switch (=1,684) | Reported time to larval death; ambiguous pupation codes | Endpoint quality NOT IDENTIFIABLE for adult/pupation from CSV |
| Ferreira & Rodrigues, doi:10.5061/dryad.3bk3j9kk2 | *Danaus erippus*, *Danaus gilippus* | 40 laboratory; 54 field (=94) | Lab survival and development time; field first/second-instar presence | Available within lab; field host and site confounded |
| Verspagen et al., doi:10.5061/dryad.kwh70rz2m | *Melitaea cinxia* | 1,680 | Death vs reaching diapause (not adult) | 14-family matched host×temperature contrasts computable |

**Never pool all 3,458**: they are neither independent temperature-treatment clusters nor equivalent life-stage endpoints.

## 1. Danaus: survival versus conditional development time

Exact original `lab_performance.txt`, original n=10 assigned per species×host; *Asclepias curassavica* is native in the focal field geography, *Calotropis procera* exotic. Source-coded survival is 1/0. Development time averages condition on recorded success; the source explicitly does NOT include died-early individuals in developmental mean.

| Butterfly | Native Asclepias survival | Exotic Calotropis survival | Native survivors mean development time | Exotic survivors mean | Exotic minus native |
| --- | --- | --- | --- | --- | --- |
| *Danaus erippus* | 9/10 | 7/10 | 13.00 d (n=9) | 17.43 d (n=7) | **+4.43 d** |
| *Danaus gilippus* | 10/10 | 9/10 | 14.20 d (n=10) | 16.33 d (n=9) | **+2.13 d** |

Both species show delayed successful development on the exotic in this small lab experiment. This direction was **already in the original paper**. The survival differences are only 2 and 1 larvae of ten; **no evidence for a universal exotic survival disadvantage or geographic population decline**. Survivor-only developmental times must not be treated as unbiased developmental speed across all assigned larvae.

Field original n=54 *D. erippus* entries are distributed: Asclepias unbagged 3/12 observed survive (1 of 13 missing), bagged 8/10 observed survive (3 of 13 missing); Calotropis unbagged 1/14, bagged 7/14. Four NAs were reported as cattle-grazing missing, not deaths. **Different host plants at two different sites**: within-site exposure contrast is descriptive; host identity, vegetation patch, abiotic conditions and site are inseparable. Bagging may alter environmental conditions. The difference between field host values is NOT a randomized native-vs-exotic host plant effect.

## 2. Melitaea: stage-valid thermal/family reanalysis

The source `Survival_data.csv` contains 14 mothers (18 egg clutches), 4 day temperatures (28,30,32,34°C, nighttime 15°C), 2 host plants (*Plantago lanceolata* / *Veronica spicata*), 15 larvae per family×temperature×plant group. Exactly **112 groups ×15 =1,680**; all rows have binary `Dead=0` (arrived alive at diapause) or `Dead=3` (died). Adults/reproduction are not measured.

| Temperature (day) | `Plantago` arrival alive | `Veronica` arrival alive | Mother-paired `Veronica–Plantago` |
| --- | --- | --- | --- |
| 28°C | 192/210 = **91.43%** | 199/210 = **94.76%** | +3.33 pp |
| 30°C | 194/210 = **92.38%** | 203/210 = **96.67%** | +4.29 pp |
| 32°C | 186/210 = **88.57%** | 201/210 = **95.71%** | +7.14 pp |
| 34°C | 185/210 = **88.10%** | 199/210 = **94.76%** | +6.67 pp |

Across four temperatures equal-weighted within the 14 mothers, `Veronica–Plantago` **+5.357 pp**, and the **exploratory 9,999-draw maternal-family bootstrap 95% percentile interval: +1.071 to +10.000 pp** (source-verified script seed 20261019). Between-family signs: 9 positive, 4 negative, 1 exactly zero. Thus **a host effect on reaching diapause is visible**, but the estimate is from a single source experiment.

The specific extreme-temperature interaction `[(Vs–Pl) at 34°C] – [(Vs–Pl) at 28°C]` is **+3.333 pp**, 9,999-family-bootstrap 95% **−4.286 to +13.810 pp**, includes zero. Therefore **an increasing or sign-reversing temperature effect is not confirmed by this selected survival-only diagnostic**. Developmental physiology and original time-to-event analyses may show other host×temperature responses, as already published; this one contrast does not refute them. Seven maternal families show raw positive and negative effects across the four temperatures, but 15 larvae per cell and rare death events mean this is **not proof of genotype-specific biological sign reversals**.

Data integrity: four `Dead=0` rows have numeric `Time2=4` though README says alive `Time2=NA`. Event-coded `Dead` was retained, and all four exceptions are in the JSON receipt. A sensitivity coding those four as censored/unknown rather than survived is necessary before stronger survival inference; current observations **do not justify silently treating them as deaths**.

## 3. Vanessa: individual time events are not a success denominator

The original `larv.plants.csv` contains 1,154 rows and `larv.switched.csv` 530; all **1,684** have `status=2`. Original paper separately reports four reached pupation phase, two reached a chrysalis. Those cannot be identified as successful cases from the supplied two tables' `status`. Consequently no observational **pupation success rate, true zero-success rate, adult production or cross-host recruitment** is independently reconstructed here. Species/trial/source and large batches of same event times must not masquerade as individual stage-confirmed survivorship.

## 4. Which studies actually overlap global 239-species panel?

Frozen `data/frozen/figure_sources/anthropogenic_species_metrics.csv` contains:
- *D. gilippus*: 17 host species, native potential host-resource geography 154 WGSRPD3 cells, contemporary 170, **16** introduced-added
- *V. cardui*: 104 host species, native 282, contemporary 360, **78** introduced-added

*D. erippus* and *M. cinxia* are not entries in this frozen 239-species S1 panel. **Zero** of these occurrence-level global geographic counts are a measured local cohort host-area link. The fact that a species appears in both data systems **does not permit** multiplying 16 or 78 resource cells by one laboratory larval survival fraction. Seasonal host choice, locality, native fallback, actual introduced plant presence and egg-to-next-generation endpoints are all missing.

## 5. Interpretation and clear next target

**Supported within original studied systems (prior art, now source independently recalculated):**
- Two Danaus species' successful development on exotic versus native hosts is delayed (small n=10 per condition).
- Melitaea diapause-arrival probability depends on which of two host plants was assigned, with 14-family replication; a large high-vs-low-temperature interaction has *not* been estimated precisely.
- Vanessa life-stage recoding is not recoverable by interpreting raw `status=2` as a clean failure of pupation.

**Not identified:**
- Within-patch oviposition diversion `q_t` from a viable native host to a worse introduced one;
- Same egg-cohort comparisons of native vs introduced hosts followed through flight-capable adults or next-generation offspring;
- Plant-availability seasonality and resource quantity, movement/source–sink contribution or geographically matched butterfly population growth;
- Repeated independent butterfly-species and site/season units proving nonmonotone total recruitment under host redistribution.

**Decision**: The original global butterfly paper can still claim *potential larval-resource geography* and source-dependent novelty, not fitness. The **independent ecological follow-up** should target repeated same-species × patch × season egg allocation and *tagged* larval cohort-to-adult outcomes, paired native and introduced hosts accessible **within** sites, before any inference that adding exotic host plants lowers butterfly population growth. Primary benchmark is host-presence-only; evaluate whether quality×choice×timing improves held-out cohort production. Prior published effects in these three source papers cannot be sold as a newly discovered pan-butterfly effect.

## Reproduction

Run from root:
```sh
python scripts/download_exact_cross_butterfly_sources_20261019.py
python scripts/analyze_cross_butterfly_lifestage_raw_20261019.py --output artifacts/original_butterfly_lifestage_receipt.json
```
Fail-closed if Zenodo source MD5 or source length differs. If unable to download exact source, report NOT EVALUABLE, not a fictitious value. See `data/exploratory/realized_host_20261019/README_SOURCE.md`. The independent Github Actions receipt is a machine-readable output, not a new dataset sampled in October 2026. Frozen main GEB manuscript PR #38 unchanged; PR #43 alternate manuscript remains draft.
