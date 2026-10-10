# Experimental butterfly ecology after rejecting static resource and parasitoid shortcuts

**Date:** 2026-10-10  
**Decision:** SOURCE-VERIFIED AUDITS COMPLETE; STOP UNSUPPORTED NEW ECOLOGICAL CLAIMS. Prioritize a genuinely randomized temporal resource supply × heterospecific competitor experiment.  
**Branch:** `exploration/realized-host-window-v01` — GEB PR #38 NOT altered.

## A. What was actually validated during this continuation

### Swedish *Urtica* host and parasitoid paper: existing direct prior art

Audusseau et al. (2021), *Oikos*, DOI https://doi.org/10.1111/oik.07953 and public original raw data https://doi.org/10.17045/STHLMUNI.14260211 (Figshare article 14260211).

Authors already tested: same-week co-occurrence of butterfly species, abundance, larval season, county, year, and time since establishment of the newcomer *Araschnia levana*, against parasitism rates. They documented **6,777 butterfly larvae**, **1,508 parasitised**, **11 identified parasitoid species**, and a 19-site design. This already covers substantial parts of any generic "shared resource → parasitoid interaction" claim.

Original source checksum gate succeeded in [GitHub Action 38016704850](https://github.com/zuizui0223/chocho/actions/runs/38016704850): Figshare `Audusseau_etal2021_batch_monitoring.csv`, exact MD5 **5e861203477ce8d2ae3efdfacc7436fc**, **1,080 original batch rows**; original first-observation `Araschnia` table has 18 rows. The original source metadata and CSV checksums match. **No new parasite ecology effect was inferred by verifying the archive.**

The post-source lagged parasitoid *feasibility* audit was performed at [GitHub Action 38017138169](https://github.com/zuizui0223/chocho/actions/runs/38017138169); 3/3 unit tests passed. The original sheet uses **au, aio, alev, va** rather than full binomials. This correction was made before any model fitting. The original table has 249 *Ag. urticae* batches, 204 *Ag. io*, 168 *Ar. levana* and 459 *V. atalanta* batches. **28** `Aglais` batch rows contain zero larvae and are not individual-fate failures; the **425** positive-denominator `Aglais` rows encompass 4,510 butterflies with **144 `Sturmia bella` parasitisms**, matching the paper's 39+105 adult parasitoids attributed to the two resident hosts.

The site-week lag screen (only previous sampling occasion 1–3 weeks earlier) produced:

| Prior *Araschnia* | Current *Araschnia* | Eligible resident butterfly batches | Resident larvae | *Sturmia bella* parasitism | Different source BMS codes |
| --- | --- | ---: | ---: | ---: | ---: |
| absent | absent | 233 | 2,234 | 4 | 17 |
| absent | present | 48 | 565 | 45 | 9 |
| present | present | 34 | 403 | 27 | 6 |
| present | absent | 15 | 298 | 36 | **4** |

**Do not interpret these raw ratios as delayed indirect competition or parasitoid reservoir effects.** *Araschnia levana* is absent in Stockholm, and *S. bella* is also absent from that northern region in the original publication, so spatial geography can largely determine the raw contrast. The critical *previous present/current absent* category is only **four BMS site codes**. Moreover the source contains **21 unique BMS site codes in 2017, 19 in 2018**, while the published study describes a **19-site design**; this requires independent reconciliation, not invention of 21 independent sites. The valid outcome and source gate intentionally FAILS independent temporal-contrast support. **No predictive effect was fitted.**

**Field infection confounding:** In the original study, larvae were fed leaves collected daily from their source locations during laboratory rearing. Microtype-egg parasitoids including *Sturmia bella* can be ingested *after field collection*. Therefore the observed parasitoid emergence cannot be automatically dated to the larval field collection and used to establish a multi-week in-field reservoir mechanism. Even with more sample sites, an untreated shared-leaf exposure control would be needed.

### Kyoto swallowtail two-species competition: exact prior art with original time data

Hashimoto & Ohgushi (2017), *Population Ecology*, DOI https://doi.org/10.1007/s10144-016-0568-8, previously tested repeated herbivory and plant regrowth. Hashimoto & Ohgushi (2023), *Ecology and Evolution*, DOI https://doi.org/10.1002/ece3.10164, found asymmetric competitive effects on the shared **native** host *Aristolochia debilis* between Japan-native *Atrophaneura alcinous* and introduced-to-Japan *Sericinus montela*. The 2023 lab experiment estimated about **2.6×** greater food demand in *A. alcinous*. The presence of *S. montela* worsened *A. alcinous* survival and development time, whereas the reverse effect was not observed. **Food-demand asymmetry, host biomass depletion, prior herbivory and plant regrowth are already published results.**

The exact original Figshare article https://doi.org/10.6084/m9.figshare.23170898 was independently source-audited and confirmed in [GitHub Action 38017296067](https://github.com/zuizui0223/chocho/actions/runs/38017296067), with 3/3 software tests passing. The original code and data are CC BY 4.0. The files include **690** cage × time rows of larval stages/survival, **30** experimental cages with four `defoliation.1–4` FINAL plant-individual scores (four different potted plants per cage, **not four dates**), and a separate **120-row** one-day food-demand assay. Original CSV byte checksums are included in the artifact.

These are real repeated measurements, **not** independent randomization of the *timing of new edible biomass*. The four `defoliation.1–4` columns are **four different plant individuals measured at the final endpoint**, not repeated observations over time. Each is a post-treatment herbivory outcome influenced by caterpillar density and feeding, not an independently assigned plant-regrowth trajectory. Fitting a model of mortality on these same post-treatment defoliation measurements cannot establish a causal temporal-regrowth effect, and a same-source refit would not be a separate ecological discovery.

**Correction (2026-10-10, original codebook cross-check):** The original 2023 plant CSV `defoliation.1–4` labels refer to the four potted plants **within** each cage at its final assessment; they are not sequential plant measurements. The 690 repeated rows concern butterfly larval stages. Previous characterization of four temporal plant assessments was mistaken. The 2012 ESJ59 poster already proposed temporal/spatial variability in *A. debilis* regrowth as a competition-relaxation hypothesis. A novel result must address the explicitly unresolved food-quantity-independent suppression (S. montela density reduces A. alcinous survival despite no clear S. montela-induced food shortage), separating plant conditioning from direct interference with matched edible biomass. See `KYOTO_NONRESOURCE_COMPETITION_MECHANISM_V01.md`.\n\n## B. The new question, and why it is not fixed by the experiment design

> At **identical total season-long edible plant supply**, does the **timing of available fresh host tissues** switch which butterfly suffers from competing with the other species?

Potential **qualitatively different** biological responses: (1) high-demand *A. alcinous* benefits most from fresh tissue before its peak instar demand; (2) its benefit reverses if later instars require fresh tissue most; (3) the smaller-demand *S. montela* can compensate under asynchronous supply; (4) competitor-presence effects remain asymmetric regardless of supply timing; (5) plant tissue quality or behavior, not mass availability, determines the result.

The **target is the interaction between an independently manipulated food-supply schedule and an independently manipulated heterospecific neighbor**, not the tautology that withholding food during a critical period reduces survival.

**Keep host plant identity (*Aristolochia debilis*) fixed.** Here the expanding alien organism is a **butterfly**, not an introduced plant. No botanical globalization mechanism or 239-species phylogenetic generalization follows from this experiment alone.

## C. Causal experiment requiring genuinely new data

Randomize independent potted host plants/enclosures in a factorial design:

- Focal species: *A. alcinous* or *S. montela*.
- Other butterfly species: absent or present, using comparable and predefined original heterospecific densities and independently observing plant biomass depletion.
- Edible-host supply timing: early-pulse / distributed / late-pulse fresh tissue replenishment. **Match total season-long edible biomass across timing arms** and measure actual accessible biomass daily or at instar-specific intervals; apply matched handling/sham operations, standardized tissue developmental stage and plant source.
- When relevant, make separate adjustments for focal/competitor larval age, as their food demand and phenologies differ. Otherwise nominally equal total biomass can mask accessible-stage differences; report measured access.
- Outcome per **originally assigned focal cohort**: successful flight-capable adult emergence (including deaths and complete zeros), larval developmental time, competitor survival, instar-specific tissue consumption, host regrowth and phenological stage. No multiplication of sample size by counting individual larvae on the same plant as independent cages.
- Block/randomize by source population, parental family and experimental date. Inference must treat plant/enclosure, not larva, as the treatment unit; ideally replicate across populations and separate experimental seasons. Perform a pilot to estimate independent-unit variance before calculating final power.
- Require applicable collection/transport/containment permissions and biosecurity, especially for the introduced butterfly. Do not release or translocate individuals into natural habitats as part of this study.

**Main difference-in-differences estimand for each focal species**: let (Y(s,c,r)) be viable focal adults per initial cohort under competitor (c\in\{0,1\}) and resource timing (r\in\{early,late,distributed\}). Then compare ([Y(s,1,early)-Y(s,0,early)]-[Y(s,1,late)-Y(s,0,late)]), with distribution arm as a mechanistic reference. Estimate the contrast across independently assigned plants and source populations.

**Positive gate:** a reproducible interaction effect on *viable adults*, not only egg-laying, herbivory or plant cover, with error intervals and species-specific sign; do not call a negative competition effect "adaptive" merely from its sign.

## D. What NOT to do next

1. Do not claim *Sturmia bella* lag from Swedish raw occurrence; 4 site codes, geographic confounding, and possible post-collection parasitism do not identify it.
2. Do not interpret four FINAL plant individuals per cage as four timepoints, or use post-treatment defoliation outcomes to retrofit a causal plant-regeneration timing effect.
3. Do not reinstate the earlier rejected static spatial null or time-history GBIF detection lag as a novel butterfly mechanism.
4. Do not relabel the existing 2017/2023 research as a novel plant-mediated competition result.
5. Do not merge these notes into the GEB global potential-resource manuscript: it has a different estimand and PR #38 remains untouched.

**Current conclusion:** the sources are audited and the novelty floor is clear. A new empirical paper on temporally modulated competitive asymmetry is **potentially testable**, but **not established by existing public data**. Further progress should produce an independently randomized experimental outcome or a genuinely equivalent exogenous temporal supply contrast, not a selection of favorable post-hoc models.
