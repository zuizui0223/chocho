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

## 2a. Source-driven advance: measure the high-demand instars, not only remaining leaves

The published article's own **Table 1** (original openly accessible [authors' PDF](https://www.ecology.kyoto-u.ac.jp/~ohgushi/ja/achievements/PDF/Ohgushi271.pdf)) gives original leaf consumption (dry mass per individual):

| Species | instars 1–2 (mg) | instars 3–5 (mg) | total (mg) | fraction consumed during instars 3–5 |
| --- | ---: | ---: | ---: | ---: |
| *A. alcinous* | 21.66 | 1,003.23 | 1,024.89 | **97.9%** |
| *S. montela* | 1.88 | 394.25 | 396.13 | **99.5%** |

These are **published consumption measurements**, not standing leaf biomass or known starvation thresholds. The larvae in instars 1–2 were originally reared in groups of 10; the later-stage consumption estimates have reported SEs (approximately 39.56 mg for A and 9.10 mg for S).

To locate the recipient's high-demand stage without inventing a vegetation time series, the actual MD5-identical `Hashimoto_and_Ohgushi_2023_larvaldata.csv` was independently checked through [GitHub Action 38035260979](https://github.com/zuizui0223/chocho/actions/runs/38035260979) (**2 software tests passed**). Of **30** cages, **24** initially included *A. alcinous*; late instars (3–5) appeared in original census records for **23**. The first observed late-instar day was **day 3 in 3 cages, day 5 in 14, day 7 in 6**. This is an observed census interval, not exact molt dates; the remaining one cage cannot be classified by mechanism from this count alone.

**Practical change:** Start recording accessible *Aristolochia* tissue **before** day 3, throughout third to fifth instars and before pupation, adjusting the actual calendar after a source-independent pilot under the new conditions. The original raw `a1–a5` count columns describe larval stages, but no synchronous accessible leaf mass, age/chemistry or stage-specific starvation exposure is recorded.

**Not new:** The 2023 original authors themselves explicitly suggested in the Discussion that added *S. montela* might leave more foliage at the end **because native larval survival fell and subsequent consumption declined**. This feedback is an already published interpretation, not a hypothesis discovered by our reanalysis. What remains empirically unresolved is whether an undetected stage-specific shortage ever occurred, or whether tissue quality/contact drove the negative interspecific effect.

### Source-independent stage-aware experimental ledger validation

The updated `KYOTO_FOOD_CLAMP_ANALYSIS_CONTRACT_V01.json` and `scripts/analyze_kyoto_food_clamp_outcomes.py` require each scheduled food-access observation to contain **living native larvae** and **living third–fifth instars**, accessible area, separately added accessible leaf area and both current/added leaf-age categories. The original ITT whole-cage estimand does **not** condition on survival to third instar; a cage in which all natives died earlier cannot be reclassified as evidence of adequate food. Source-independent synthetic tests of the ledger and randomized cohort integrity: [GitHub Action 38035377637](https://github.com/zuizui0223/chocho/actions/runs/38035377637), **10/10 successful**.

A useful sham where feasible is to add *the same age and provenance of leaf material* in all four arms, making the additional leaves inaccessible behind matched perforated barriers in the natural-food arms and accessible in the clamp arms. This partially equalizes fresh-leaf arrival, odors and handling, but could still modify plant microclimate or larval movement; monitor those directly. The true experimental unit is one **entire cage** even when it contains multiple potted plants. The CSV column `plant_id` is a legacy name for that randomized cage ID.

> **2026-10-10 independent-study update:** Jang et al. (2026), DOI [10.5141/jee.26.017](https://doi.org/10.5141/jee.26.017), studied these **same two butterfly species** on a DIFFERENT host plant (*Aristolochia contorta*) in South Korea. Strong opposite light/microhabitat associations and only **5 both-species out of 875 original ramet observations** are reported. This is **published observational prior art**, not a novel effect of chocho and not causal evidence that competition caused spatial partition. Their printed quadrat denominator has an internal **255 versus 215** conflict requiring author confirmation; see [source-audited published-count decision](SWALLOWTAIL_2026_MICROHABITAT_AND_RESOURCE_LEGACY_DECISION_V01.md) and [reproducible arithmetic audit](https://github.com/zuizui0223/chocho/actions/runs/38045210348). A Kyoto test should now measure **actual co-use opportunity before manipulating payoff** and balance light within experimental blocks. **No edits to GEB PR #38.**

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
- `scripts/analyze_kyoto_food_clamp_outcomes.py` — source-blind randomized **whole-cage** interaction, natural/clamped competitor effects separately, block bootstrap, *stage-specific* repeated food-floor and leaf-age diagnostics, missing outcomes retained; **10/10 current tests passed** in [run 38035377637](https://github.com/zuizui0223/chocho/actions/runs/38035377637).
- The latter reports correctly constructed worst-case missing-outcome bounds for an interaction with *both positive and negative coefficients*. Simply treating every unknown adult outcome as success is **not** necessarily the upper bound of such an interaction, and has been explicitly corrected.
- These test fixtures are **synthetic for code integrity only**. They are NOT claimed to be observed butterfly outcomes or evidence of biological effects. No experimental cohort or statistical sample size has been invented.

## 5. Decision

**Research question worth testing:** Does temporally adequate *edible* host supply remove a native butterfly's sensitivity to a newly arrived competitor, or does nondepletion suppression persist through a plant-quality/behavior pathway?

**Still unidentified:** plant-quality mediator, food shortage in the 2023 cages, viability beyond newly tracked adult emergence, and generality beyond the Kyoto two-species system.

**Immediate future real-world prerequisites:** permits/containment for introduced butterflies; genuine accessible host plant inventory; source-family matched plant/cage assignment; resource target and independent plant-level variance pilot; only then execute the randomized field/indoor experiment and analyze individual fates.

**Manuscript boundary:** Butterfly invader × native host plant is not a direct causal validation of the unrelated chocho global *introduced HOST PLANT* resource expansion paper. Never merge this exploratory note or new pre-experiment protocols into the existing GEB manuscript without independent editorial/scientific justification.
