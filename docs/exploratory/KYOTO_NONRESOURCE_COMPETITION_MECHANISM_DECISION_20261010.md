# A native butterfly loses despite no apparent food shortage: identify the missing causal pathway

**Date:** 2026-10-10  
**Status:** SOURCE-AUDITED ORIGINAL CASE / DISCRIMINATING EXPERIMENT PROPOSED / NO NEW BIOLOGICAL MECHANISM ESTIMATED  
**Working branch:** `exploration/realized-host-window-v01`; GEB source paper PR #38 **unchanged**.

## 1. Why this question survives repeated novelty tests

In a 2023 randomized-density common-garden cage experiment, *Sericinus montela* (introduced into Japan) significantly depressed survival and delayed pupation of native *Atrophaneura alcinous*, although the presence of *S. montela* did not straightforwardly deplete the shared native food plant *Aristolochia debilis*. The 2023 authors explicitly considered (i) direct interference and (ii) plant-quality changes as explanations of the non-depletion suppression; this is **a published open problem, not a new observation**. Their short AAI addition assay did not show feeding/growth differences (Hashimoto & Ohgushi 2023, DOI [10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164)).

An earlier study found that host plants re-grew after herbivory, and that previously attacked plants did not automatically depress larval growth after regrowth (Hashimoto & Ohgushi 2017, DOI [10.1007/s10144-016-0568-8](https://doi.org/10.1007/s10144-016-0568-8)). The authors also **already proposed** temporal variability in regrowth as a potential competition-relaxation mechanism in [ESJ59 P1-196A (2012)](https://esj.ne.jp/meeting/abst/59/P1-196A.html). Arrival order is a known plant–herbivore effect in other systems (Erb et al. 2011, DOI [10.1111/j.1365-2745.2010.01757.x](https://doi.org/10.1111/j.1365-2745.2010.01757.x)). Later work shows VOC–parasitoid interactions involving *Sericinus* and **another** *Aristolochia* host in Korea (Park et al. 2025, DOI [10.1002/ece3.71175](https://doi.org/10.1002/ece3.71175)); this does not prove the pathway in Kyoto.

**Honest novelty floor:** causally resolve **which part of the introduced butterfly's presence suppresses the native butterfly at a controlled amount of edible plant tissue**, and whether it persists in the absence of the introduced butterfly. New mechanistic identification could be a valid focused butterfly ecology paper. The broad idea of indirect competition or induced defense is NOT new, and a single-population species pair would not justify an unqualified global invasion law.

## 2. The source data were not what earlier exploratory notes assumed

Primary Figshare dataset: [10.6084/m9.figshare.23170898](https://doi.org/10.6084/m9.figshare.23170898).

[Original-MD5 source integrity and design audit](https://github.com/zuizui0223/chocho/actions/runs/38021529338), with three software tests passing:

| Raw source item | Confirmed structure | Scientific meaning |
| --- | --- | --- |
| Larval time series | **690 rows = 30 cages × 23 observation dates each** | Butterfly instar/pupation, not plant regeneration |
| Original cage factorial | **15 initial density combinations × 2 independent cages each** | Randomized butterfly density, NOT randomized food supply |
| `plantdata.csv` | **30 cage rows × 4 FINAL potted plant observations = 120 plant endpoints** | `defoliation.1–4` are different plant individuals, NOT four dates |
| Plant remaining-food scores | Categories **0 (24 plants), <25 (21), 25–50 (8), 50–75 (9), >75 (58)** | Coarse observer ordinal percentages, NOT directly measured continuous edible grams |
| Pupation endpoint | `cumul.sp`, `cumul.ap` | Cumulative pupae; no individually documented viable adult emergence |
| Direct interference or plant quality measures | **None** | Not estimable from existing source |

The initial two-species larval design contained 120 *S. montela* and 120 *A. alcinous* larvae across cages; these are **not 240 independent cage replicates**. No new field/experimental mechanism has been estimated by auditing this source.

**Correction logged:** An earlier note referred to four repeated plant defoliation measurements. This was **wrong**. Only larvae were measured repeatedly; four plant scores were taken at the END for four plants within each cage. The old route was corrected in `BUTTERFLY_SHARED_HOST_TEMPORAL_PAYOFF_GO_NO_GO_20261010.md` and `KYOTO_TEMPORAL_REGROWTH_COMPETITION_PROTOCOL_V01.json`. The source checks uncovered this before any new causal fit was claimed.

## 3. Causal design: plant conditioning versus direct interference

The primary contrast holds plant species, target quantity and developmental stage of the food constant, and randomizes **who previously fed on the plant** rather than merely whether a second butterfly is present.

### Phase I: does suppression persist in the competitor's absence?

Randomize independent potted *A. debilis* plants to the following history groups and remove any donor before introducing a fresh, independently assigned recipient larva.

| Donor plant history | Mechanism separated |
| --- | --- |
| No-donor sham | Handling/environment baseline |
| Mechanical damage at matched leaf-area loss | Quantity loss and generic wound response |
| Previous *S. montela* larvae, matched feeding damage | Species-specific plant/residue legacy beyond tissue removal |
| Previous *A. alcinous* larvae, matched feeding damage | Conspecific legacy versus heterospecific legacy |

Measure how much **edible, accessible biomass and leaf stage actually remain** when the recipient arrives and during development. Standardize plant age, moisture, nutrient status, date and enclosure handling. The immediate recipient is *A. alcinous*; repeat with *S. montela* if both independently sourced populations and host material are feasible.

**Primary estimand:** difference in **flight-capable adult emergence per initially assigned recipient cohort** between Sm-conditioned and equally damaged mechanically conditioned plants, tested with donor absent. Also report pupation, development time and the duration-specific response. Do not count all larvae on the same plant as independent randomizations.

Interpretation: a negative effect that persists without a live competitor argues against direct contemporaneous larval interference, **but not for a specific chemical**. Leaf-age/quality changes, residues, induced responses and microbiota remain competing paths. A mechanical-damage effect equal to Sm-conditioned treatment argues against donor-species-specific induction.

### Phase II (trigger only if Phase I warrants): isolate plant-tissue and direct-contact paths

- **Plant-mediated/leaf-borne pathway:** compare fresh tissue from conditioned versus sham/mechanically damaged donor plants transferred into identical clean recipient units at equal edible amounts, with donors absent and transfer handling matched. This distinguishes a tissue-associated effect from persistent social contact, but does not identify a molecule without further manipulation.
- **Contemporaneous interference:** manipulate a live *S. montela* neighbor with physical contact allowed versus physically separated while maintaining equivalent food accessibility, with randomized independent enclosures and sham partitioning. An effect limited to contact-allowed conditions supports a direct-interference route.
- **Post-treatment timing:** a SHORT no-donor post-feeding window versus a predefined later regrown state can reveal a transient mechanism, with *actual plant developmental state* measured. Do not switch windows after seeing survival results.
- **Mechanism-specific measurements:** objective leaf mass, tissue stage, nitrogen/proxies, other measured metabolites or microbes with separate controls. A detected metabolite that correlates with outcomes does not by itself prove mediation.

### Counterfactual outcomes

| Observed result in newly randomized trial | Narrow valid conclusion |
| --- | --- |
| *A. alcinous* performs worse on Sm-conditioned than mechanical-damaged plant, Sm absent | Non-quantity donor-history effect exists; specific quality mechanism still unknown |
| Sm-conditioned ≈ mechanical damage < sham | Generic damage or amount/handling more likely; not Sm-specific |
| All no-donor histories similar, but contact-allowed Sm suppresses native larvae | Direct coexistence/interference candidate |
| Equal food supply eliminates all deficits | New mechanism not supported; reinterpret 2023 density effect as context/amount-mediated |
| Strong opposite effect (Sm conditioning improves native survival) | Plant-mediated facilitation; investigate without forcing competitive narrative |

### Replication and claims

Randomize plants/enclosures (not individual eggs on one plant), block by experimental batch, plant maternal genotype and butterfly source family/population, and retain failed larvae in denominators. Use an *a priori* pilot of independent plants to estimate variance and failure rate before numerical power calculations. Do not invent an exact required sample size from 30 cages with only two replicates per density treatment. Absent cross-population replication, state the result applies to the tested plants and butterflies only.

Butterfly transport/captive maintenance and nonnative species containment are subject to local permits and facilities. No open-environment releases or unauthorized deliberate spread.

## 4. Relationship to chocho global paper

The 239-butterfly GEB project reconstructs worldwide geographic opportunity added by **human-moving host PLANTS**. In this Kyoto mechanism, the expanding alien is a **butterfly** and the shared host plant is **native**. The two should be independent manuscripts and should not be presented as direct causal validation of WCVP/HOSTS range projections.

**Go/no-go as of now:** Source-valid original factorial data; unresolved mechanistic alternative is real but already noted in the published paper; experimental design written, **not implemented**. Stop turning ordinal final food-remaining scores into temporal exposure or hand-selected mediation evidence. A genuine result requires independently randomized plant history versus mechanical damage and adult offspring fate.
