# Stage-specific food limitation versus plant-mediated interference in Kyoto swallowtails

**2026-10-10 | EVIDENCE-BASED EXPERIMENT PLAN, NOT A NEW BIOLOGICAL FINDING**  
Worktree: `exploration/realized-host-window-v01` only. GEB submission PR #38 remains unmodified.

## Why this is the next distinct question

In [Hashimoto & Ohgushi (2023), *Ecology and Evolution*, DOI 10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164), raising introduced *Sericinus montela* density depressed native *Atrophaneura alcinous* survival and prolonged development. Yet the study did not detect a corresponding decrease in food remaining **at the end of the experiment**. The authors highlighted unmeasured plant quality as a plausible mechanism. Direct offensive/interference behavior was not seen and naturally synchronous occupation of the same plant was relatively infrequent, making a plant-mediated mechanism plausible but **not demonstrated**.

### Critical inferential point

**No corresponding reduction in FINAL ordinal food-remaining categories is not proof that fresh edible foliage remained continuously sufficient at every vulnerable instar.** Plants can be depleted during a critical developmental window and later regrow; leaf age and quality can change independently of final percentage cover. The 2012 [ESJ59 P1-196A](https://esj.ne.jp/meeting/abst/59/P1-196A.html) had already proposed temporary/spatially variable plant regrowth as a competition-relaxation process, and [Hashimoto & Ohgushi (2017)](https://doi.org/10.1007/s10144-016-0568-8) measured plant recovery after damage. **The concept of regrowth is not novel.** The currently untested causal discriminator is whether equalizing *time-resolved edible access* eliminates the 2023 negative density effect, or whether a residual quality/direct-interference path persists.

## Original-data audit—fully reproducible and audited, no mechanism claimed

- Original Figshare dataset: https://doi.org/10.6084/m9.figshare.23170898.
- Latest passing exact-MD5 source audit: https://github.com/zuizui0223/chocho/actions/runs/38021529338, 3/3 software tests, SHA256 and source MD5 retained in run artifacts.
- Exactly **30 cages × 23 larval observation occasions = 690 rows**, **15 randomized initial competitor/recipient density combinations × 2 independent cage replicates each**.
- Within each cage **four separate potted host plants**, each scored at a **single endpoint** for remaining food. The 120 ratings are categories: **0 (24), <25 (21), 25–50 (8), 50–75 (9), >75 (58)**. They are **NOT** a plant time series and NOT 120 independent cages.
- Original cumulative pupation and instar stage are recorded; independently documented flight-capable adult production, food-access trajectories, plant induced chemistry and direct contacts are not.
- Therefore **NO new plant-mediated, aggressive-interference or hidden-bottleneck mechanism is statistically identified by the old data**.

## Focused randomized intervention before more chemistry or diet modelling

The frozen machine-readable causal protocol is `KYOTO_HIDDEN_RESOURCE_BOTTLENECK_FACTORIAL_V01.json`.

Randomize each **independent potted plant/enclosure** (the primary experimental unit) into this 2×2:

| | Unaltered natural edible supply | Independently monitored and supplemented edible supply |
| --- | --- | --- |
| *S. montela* absent | NATURAL_NO_COMPETITOR | CLAMP_NO_COMPETITOR |
| *S. montela* present | NATURAL_COMPETITOR | CLAMP_COMPETITOR |

Fixed initial *A. alcinous* neonate count per cage. Choose positive *S. montela* larval density and resource-support threshold **before outcomes**, using a separate blinded pilot and available biology, not arbitrary 2023 cells. All treatments get comparable enclosure checks/sham handling. Measure actual edible leaf surface/mass and plant tissue developmental state repeatedly THROUGHOUT native larval instars. Maintain an adequate, comparable accessible food target in both supported arms, not just the same initial or endpoint biomass.

**Primary endpoint:** number of **flight-capable native butterfly adults per originally assigned native larval cohort**, counting all failures/complete zeros. Report pupation and development separately.

**Key interaction:** `[Y(clamped, competitor)−Y(clamped, no competitor)]−[Y(natural, competitor)−Y(natural, no competitor)]`, with absolute success rates and uncertainty across randomized cages and biological source blocks. The two within-resource competitor contrasts must be reported separately.

- If the natural-supply competitor effect is negative and clamp neutralizes it **AND** repeated monitoring verifies it removed a real transient shortage, that is evidence compatible with a hidden timing-specific supply constraint. However, fresh supplement may also change nutrition/plant chemistry; a pure food-mass path still needs additional matched-tissue controls.
- If suppression remains under achieved adequate supply, food limitation alone is insufficient. Move to the distinct **donor removed × plant prior herbivory / damage-matched mechanical wounding / direct-contact separated** design `KYOTO_NONRESOURCE_COMPETITION_MECHANISM_V01.json` rather than ascribing the effect automatically to aristolochic acid.
- If suppression is absent in both natural and clamped conditions, a putative general novel mechanism is not supported; verify whether the original density contrast was successfully reproduced and stop outcome-chasing.
- If clamp never achieves a real difference in access through larval development, the experiment is **inconclusive** rather than proof against transient competition.
- If all assigned focal larvae are not followed to adult eclosion, restrict inference to pupation, not actual viable recruitment.

**Predeclared unit:** cage/potted plant. Four final host plants in the old paper are not four replicates. Future larval offspring within one plant remain nested. Collect daily leaf accessibility and treatment manipulation checks, adult eclosion, parental source and handling controls using `KYOTO_HIDDEN_BOTTLENECK_FIELD_SCHEMA_V01.md`.

## Operational artifact created, not a biological experiment

`scripts/randomize_kyoto_bottleneck_factorial.py` takes an independently prepared genuine `plant_id,block_id` eligible-plant inventory and a recorded random seed, rejects duplicate/blank plants, requires complete 4-treatment blocks and assigns arms in balanced random order. A source-blind workflow checked **4/4 test cases**: [GitHub Action 38021938106](https://github.com/zuizui0223/chocho/actions/runs/38021938106). No user plants, experimental larvae or outcomes were generated. Do not mistake a randomization-test pass for ecological hypothesis support.

## Scientific originality and limits

1. A randomized **instar-specific resource-access × competitor interaction** is experimentally independent of the 2023 finding that crowding depresses survival while final food scores did not change, and of the 2017 plant-regrowth work.
2. Positive experimental resolution of a hidden bottleneck is a contribution to mechanism in the *specific* system. It is NOT the discovery that temporal regrowth exists or that butterfly competitors can asymmetrically affect one another.
3. If the effect is a tissue-mediated legacy, a second randomized experiment is necessary to isolate plant induction, donor residues and direct contact. A metabolite correlation cannot stand in for mechanism-specific manipulation.
4. This is competition with a nonnative BUTTERFLY on a native *Aristolochia* plant. It does not validate the separate GEB hypothesis about **human redistribution of plant host ranges**. General claims need additional independent species/population replication.
5. Permit, biological containment and plant-material availability checks precede any live work; no uncontrolled introduced-butterfly release or reintroduction.

**Decision:** Prioritize the small causal bottleneck experiment, with measurable objective and explicit failure criteria, over another retrospective regression on the same 30 cages. Distinguish deliverable research DESIGN from any actually measured ecological RESULT. The original GEB paper remains unchanged.
