# Kyoto nondepletion paradox — evidence, causal ambiguity, and executable next experiment

**Date:** 2026-10-10. **Status:** SOURCE-AUDITED FACT + EXPERIMENTAL AND ANALYSIS CONTRACT; NO NEW BUTTERFLY BIOLOGICAL EFFECT.  
**Scientific isolation:** exploratory branch only; GEB submission PR #38 unchanged.

## 1. What is known from the real 2023 experiment

Hashimoto & Ohgushi (2023), *Ecology and Evolution*, DOI https://doi.org/10.1002/ece3.10164, studied two specialist butterflies on the native *Aristolochia debilis* host in Kyoto: invasive-to-Japan *Sericinus montela* and native *Atrophaneura alcinous*. *A. alcinous* requires ~2.6 times as much plant biomass for larval development as *S. montela*. Experimentally increasing *S. montela* initial density **decreased** *A. alcinous* survival and prolonged development, but did **not** deplete more measured food at the **final** cage visit; the authors' fitted final-food model showed an *increase* in remaining foliage with increasing `S. montela` density in some density contexts. The ecological puzzle is published and not a new finding of chocho.

The first experimental source audit was corrected and reproduced with original Figshare file checksums: https://github.com/zuizui0223/chocho/actions/runs/38021529338 (3/3 tests). It confirms **30 cages**, **15 randomized initial-density treatment groups**, **two cages per group**, **690 cage×23-day butterfly stage records**, and **120 independent host plants nested as four plants in each cage**. The four `defoliation.1–4` measures are **four distinct terminal ordinal leaf-remaining observations (NOT four different days)**. Their observed bins were 0:24, `<25`:21, 25–50:8, 50–75:9, `>75`:58. No plant biomass time course or experimental early-versus-late replenishment is available.

## 2. Why final leaves cannot diagnose food competition

Both of these can produce many leaves at the end:

**Biological alternative A — no shortage:** *S. montela* affects *A. alcinous* by a plant-quality or behavior mechanism. The surviving native caterpillars have abundant edible leaves throughout.

**Biological alternative B — early bottleneck with survival–consumption feedback:** during a brief instar-sensitive period, the food is temporarily inaccessible/insufficient; native larvae die, consume fewer leaves later, and the plant regrows. Thus the terminal plot looks rich in remaining leaf area even if feeding conditions were poor when mortality happened. This is *a hypothesis*, not evidence that the original cages had such an event.

A schematic directed causal graph is:

```text
Randomized Sericinus initial density
   |                    \
   v                     v
Current food availability  Plant state, chemistry, interference
   |                     |
   v                     v
Atrophaneura survival and growth
   |
   v
Subsequent food consumption by surviving larvae
   |
   v
Final ordinal host-food-remaining score
```

It is methodologically invalid to match on, adjust for, or stratify by `final_food_remaining` as if it were a *baseline* cause of survival: it can itself be affected by native butterfly survival and is a post-treatment variable. This may produce post-treatment selection/collider bias. The original experiment's exact 15 density groups and only two cages per group do not fix this by adding more regressions.

## 3. Executable new randomized comparison

Freeze a **2×2** factorial assignment before any recipient larvae enter cages:

- Competitor *S. montela*: absent / present.
- Plant food-access strategy: natural regrowth / actively monitored fresh edible-host-area floor with sham handling.

The experimental unit is the independently randomized *Aristolochia* plant/enclosure. Each arm has the same **predeclared original native larval cohort size**. Measure accessible edible leaf area throughout the larval instars, not just the start and end. Measure food supplementation, leaf age, microclimate and direct larval contact. End with tracked successful native **flight-capable adult emergence**, including deaths and complete cage-level failures.

**Primary biological estimand**:
```text
(clamped supply: competitor − no competitor)
    −
(natural supply: competitor − no competitor)
```

A positive interaction with a negative natural competitor effect is compatible with short-term food restriction, *if* the clamp actually removes it. But supplementation may change nutrition, plant chemistry or handling, so one cannot claim a pure food-mass mechanism from this intervention alone. A negative effect that persists under an achieved clamp calls for a **second** donor-history/standardized-mechanical-wounding or direct-contact design rather than asserting a plant defense.

## 4. Newly implemented reproducible workflow

The following work is executable **once genuine experimental measurements exist**:

- `KYOTO_HIDDEN_RESOURCE_BOTTLENECK_FACTORIAL_V01.json` — study question and 2×2 manipulation.
- `KYOTO_FOOD_CLAMP_ANALYSIS_CONTRACT_V01.json` — prespecified observational unit, fates, attrition bounds, interaction and limitations.
- `KYOTO_HIDDEN_BOTTLENECK_FIELD_SCHEMA_V01.md` — exact allocation, original cohort, scheduled resource-access, and adult-fate CSV field definitions.
- `scripts/randomize_kyoto_bottleneck_factorial.py` — auditable, balanced, fixed-seed block randomization; **4/4 unit tests passed** in [run 38021938106](https://github.com/zuizui0223/chocho/actions/runs/38021938106).
- `scripts/analyze_kyoto_food_clamp_outcomes.py` — source-blind randomized cage-level interaction, natural/clamped competitor effects separately, block bootstrap, objective repeated food-floor status, missing outcomes retained; **6/6 unit tests passed** in [run 38026863764](https://github.com/zuizui0223/chocho/actions/runs/38026863764).
- The latter reports correctly constructed worst-case missing-outcome bounds for an interaction with *both positive and negative coefficients*. Simply treating every unknown adult outcome as success is **not** necessarily the upper bound of such an interaction, and has been explicitly corrected.
- These test fixtures are **synthetic for code integrity only**. They are NOT claimed to be observed butterfly outcomes or evidence of biological effects. No experimental cohort or statistical sample size has been invented.

## 5. Decision

**Research question worth testing:** Does temporally adequate *edible* host supply remove a native butterfly's sensitivity to a newly arrived competitor, or does nondepletion suppression persist through a plant-quality/behavior pathway?

**Still unidentified:** plant-quality mediator, food shortage in the 2023 cages, viability beyond newly tracked adult emergence, and generality beyond the Kyoto two-species system.

**Immediate future real-world prerequisites:** permits/containment for introduced butterflies; genuine accessible host plant inventory; source-family matched plant/cage assignment; resource target and independent plant-level variance pilot; only then execute the randomized field/indoor experiment and analyze individual fates.

**Manuscript boundary:** Butterfly invader × native host plant is not a direct causal validation of the unrelated chocho global *introduced HOST PLANT* resource expansion paper. Never merge this exploratory note or new pre-experiment protocols into the existing GEB manuscript without independent editorial/scientific justification.
