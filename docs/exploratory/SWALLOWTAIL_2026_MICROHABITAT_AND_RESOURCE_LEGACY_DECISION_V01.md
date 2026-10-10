# Same larval host, different microhabitats: new 2026 butterfly evidence changes the Kyoto mechanism test

**Date checked: 2026-10-10**  
**Status: independent peer-reviewed evidence + published-count audit. No unpublished biological result, no inferred causal competition effect.**  
**Repository scope:** exploratory branch `exploration/realized-host-window-v01` only. GEB manuscript PR #38 unmodified.

## Independent 2026 publication located

**Jang JY, Choi H, Jeong SJ & Kim JG (2026),** *Niche partitioning of Sericinus montela and Atrophaneura alcinous in small scale riverside habitats*. *Journal of Ecology and Environment* 50:15, published 6 August 2026, DOI [10.5141/jee.26.017](https://doi.org/10.5141/jee.26.017), [full text](https://www.e-jecoenv.org/journal/view.html?uid=1281&vmd=Full).

Four river-bank habitats in South Korea with a **common native host plant *Aristolochia contorta***. Observations of butterfly use were recorded at two nested levels: host ramets and quadrats; larval preference and performance were also measured independently with sun and shade grown leaves.

**Published ecological evidence (not discovered by us):**

- The paper's **ramet-level** sample has 542 neither, 204 *Sericinus* only, 124 *Atrophaneura* only and 5 both, summing to **875 sampled ramet observations**. The 5 both occupy **0.57%** of ramet records, not a percentage of all *Sericinus* foodplant opportunities across a year. This does not imply exclusive coexistence failure, competition, or complete absence of earlier host use.
- *S. montela* occurrence is positively associated with relative light intensity (OR per +1 RLI percentage point **1.009**, 95% CI **1.002–1.016**), whereas *A. alcinous* occurrence is negative (OR **0.979**, CI **0.973–0.986**). The analyses are observational associations adjusted for survey region/month and partly clustered by quadrat, not randomized sunlight experiments.
- *A. alcinous* consumed a lower proportion of sun-exposed leaves (published leaf-choice experiment); in the independently manipulated diets, larval **growth differences between sun and shade leaves were not clearly detected** for either species, while there were species differences in RGR, assimilation and consumption. Thus local host use can differ even without an experimentally demonstrated shade-vs-sun larval-performance tradeoff. This preference/performance distinction is **already the paper's result**.
- The species were not proven to interact or alter one another's occupancy by manipulation. Measured preference, intrinsic habitat distribution and plant architecture are equally plausible sorting causes.
- Korea's host is ***A. contorta*** and *S. montela* is part of the local habitat assemblage, whereas the 2023 Japan butterfly density experiment used ***A. debilis*** and an introduced-to-Japan *S. montela*. Do not transfer numerical effects or botanical traits across these different plant populations, land-use histories and country contexts.

## Unexpected, reproducible denominator discrepancy in the published 2026 article

A NEW, publication-only numerical consistency audit was performed, **NOT a source-data reanalysis**. Exact paper table fragments transcribed and cross-checked in `scripts/audit_2026_swallowtail_microhabitat_published_counts.py`. Workflow [2026 swallowtail count audit](https://github.com/zuizui0223/chocho/actions/runs/38045210348) passed **4/4 software tests** and stores a JSON receipt.

| Published study count | Printed value | Independent arithmetic using the SAME publication |
| --- | ---: | ---: |
| Methods: total quadrats | **255** | Site counts 50 + 51 + 57 + 57 = **215** |
| Table 1: *S. montela* only quadrats | **155** | Month counts 12+21+21+31+30 = **115** |
| Same Table 1, same category | **155** | Site counts 26+28+32+29 = **115** |
| Other quadrat categories | 17 neither / 72 *A. alcinous* / 11 both | Both site and month breakdowns agree |
| Quadrats accounting from the internally consistent category margins | **255 as printed header total** | 17+115+72+11 = **215** |
| Ramet categories | 542 / 204 / 124 / 5 | **875**, internally consistent |

**Best-supported suspicion:** a transcription error in the original quadrat sample size and header `S. montela only N=155` versus the consistent internal 215 and 115. This is **NOT author-confirmed**: the original raw quadrat-level data are not publicly deposited, only available on request from the corresponding author. The author should confirm before corrected denominators are used as a formal reanalysis. The published coefficient estimates cannot be validated using the article aggregates alone. We must NOT invent 215 raw records or claim to have reproduced the models.

The **11 both-species quadrats** is separately consistent across the site/month table and corresponds to 11 / 215 = 5.12% **if** 215 is indeed the real quadrat denominator. The underlying 255/215 discrepancy means this percentage must be qualified; the five both-species ramets in 875 are independently internally consistent.

## Independent source availability: a critical gate, not missing evidence of a mechanism

| Study | Original biological contrasts | Data accessibility verified | What it can and cannot establish |
| --- | --- | --- | --- |
| Hashimoto & Ohgushi **2017**, DOI 10.1007/s10144-016-0568-8 | Host growth and compensation following past *Sericinus* vs *Atrophaneura* herbivory; new larvae on regrown plants | 2017 paper and supplementary **PDF** publicly listed; original leaf-level raw CSV with rearing fates **not verified accessible** | Authors' regrowth finding is prior art; no data-identified short post-damage metabolite trajectory accessible here |
| Hashimoto & Ohgushi **2023**, DOI 10.1002/ece3.10164 | Competition density, repeated larval survival stages, **final** host defoliation scores, one-molecule AAI short bioassay | Author's original **CSV files independently retrieved and MD5 verified** in the chocho repository workflows | Food access or quality THROUGH the larval stages is absent; no newly estimable mediation effect |
| Jeong et al **2024**, DOI [10.3390/plants13111456](https://doi.org/10.3390/plants13111456) | *A. contorta* **simulated herbivory × CO2**, multiple plant tissue primary metabolite responses, including systemic responses | Source article and supplementary ZIP publicly present; **raw plant-level metabolite data not public per author Data Availability Statement** | Published positive evidence that the plant chemistry/metabolome CAN respond after simulated injury; no recipient *A. alcinous* fitness or real interspecific competition, no direct transport to *A. debilis* |
| Park et al **2025**, DOI [10.1002/ece3.71175](https://doi.org/10.1002/ece3.71175) | *A. contorta* damaged-plant volatiles attract the parasitoid *Ooencyrtus* in an experiment involving specialist *S. montela* | Published results, no jointly identified field co-use and cross-species adult recruitment source verified | VOC/egg parasitoid pathway is already published, not a new inference about the Kyoto butterfly larvae |
| Jang et al **2026**, DOI 10.5141/jee.26.017 | Ramet/quadrat co-use, light-driven habitat associations and diet choice/performance | Full text and tables publicly accessible; **underlying original measurements on author request only**, quadrat count discrepancy unresolved | Important local **encounter-filter** prior art, but not a causal proof that microhabitat partition prevents competition |

This comparison explicitly rules out claiming as new: butterflies partitioning sun/shade, plant metabolite changes after herbivory, plant VOC-attraction of parasitoids, or asymmetric shared-host butterfly competition.

## The remaining distinct ecological question

The **global host-plant network** in the GEB paper estimates *regional potential* use of host taxa. The **2026 independent butterfly study** demonstrates a strong observational example of potential host sharing that need not become simultaneous micro-site use. The 2023 Kyoto experiment on the other hand can produce an asymmetric negative interaction when co-use is imposed under rearing conditions.

The explanatory framework is **two consecutive ecological filters**, not one invented transmission mechanism:

1. **Encounter filter:** Does host plant distribution translate into actual same-ramet or same-day larval co-use, given microhabitat (light, plant size, canopy, phenology, survey effort)? Host taxon presence alone cannot identify this.
2. **Payoff filter after encounter:** Conditional on a real or experimentally assigned shared host encounter, does a competitor depress the native species through stage-specific accessible food, changed host-plant quality, or direct interference? Original 2023 data do not separate those effects.

Critically **observationally low co-use is not evidence of negative competition**: microhabitat choice could create segregation even in isolation, rather than being caused by past or current rivals. Nor does the 2026 paper establish that simultaneous host co-use is necessarily rare in all *A. debilis* habitats in Japan. Such generalization would be unjustified.

### Concrete next prospective test that avoids posthoc source manipulation

**Phase A field-feasibility (Kyoto *A. debilis*):** record independently identified *A. debilis* ramets; site and plant IDs; survey dates, repeat visits and negative/detection effort; two butterfly species' egg and larval instars on every surveyed plant; exact light availability, host accessible area and plant stage; previous presence/herbivory and plant-persistence across visits. Distinguish same plant at different dates from true simultaneous co-use. Compare overlap to conditional expectations matched within site/date, plant leaf supply and light; prefer a preregistered heldout season/site. **No inference of competition** based solely on deficit of co-occurrences.

**Phase B causal study:** retain the previous *whole-enclosure* 2×2 random assignment of *Sericinus* present/absent × natural/supportive accessible *A. debilis* food supply, with matched-age, sham-supplied leaves. **Light conditions must be measured and balanced/blocked before randomization**, not tuned afterward based on outcomes. The adult fates, native neonate denominators, loss reasons and real 3rd–5th instar resource accessibility remain primary data. A follow-up independent prior-herbivore conditioning versus mechanical wounding experiment is justified only if support does not remove competitive suppression.

The novel deliverable would be a *causal* connection between naturally measured encounter probability and competitor payoff across replicated source populations/plant resource conditions, rather than a second publication describing the already-known sun–shade split.

**STOP / access rule:** No new model fitting with imaginary 2026 raw quadrat rows, no inferred corrected counts without author clarification, no reused known metabolite effects presented as a novel butterfly discovery, and no edits to GEB PR #38. This is a literature and printed-table audit plus a new evidence-informed design, **not** a newly observed biological mechanism.
