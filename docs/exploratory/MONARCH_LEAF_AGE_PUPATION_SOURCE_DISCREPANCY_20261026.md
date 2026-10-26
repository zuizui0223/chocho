# Named-butterfly host-leaf resource quality and conversion to actual pupation: source discrepancy verdict

**2026-10-26 — source-verified posthoc original-experiment audit.** This is a new original-data QA analysis of Fisher et al. (2026), DOI [10.3390/insects17020215](https://doi.org/10.3390/insects17020215), **not a new causal biological experiment** and not a butterfly population increase after introduced host relocation. No change to frozen GEB main or submissions PR #38/#70.

## Why the prior butterfly macroecology cannot prove larval production

The 239-butterfly resource geography reconstruction observes **potential presence** of known host plants across mapped regions, not whether juveniles feed, develop, pupate or become reproductive adults. The 2010–2024 Caterpillars Count field dataset validated actual juvenile use of native vs alien woody plant branches but only **15 photographed source sightings from a single identified butterfly species** matched the original 239 species. The field source largely includes moth caterpillars and cannot provide a general butterfly-specific success rate.

The Fisher et al. individually recorded **Danaus plexippus** larval feeding experiment on upper/younger (Top) versus lower/older (Bottom) **Asclepias syriaca** leaves is a controlled butterfly life-stage calibration. Original author GitHub [MilkweedLeafAgeMonarchSuccess](https://github.com/kelseyefisher/MilkweedLeafAgeMonarchSuccess) source commit **acbb93c75e1f4810355c6ce03e6076e2303b2a98** was checked out and six input CSV Git blobs were independently verified. Eight source treatment summaries represent **480 originally assigned larvae, 4 YearTrial blocks (2 in 2023, 2 in 2024)**. Primary fate = reaching *pupation*, not adult eclosion.

## Critical original-data discrepancy, independently replicated

**2024 trial 1 (YearTrial=20241)** is the sole source conflict:

| YearTrial20241 source | Top pupated of 60 assigned | Bottom pupated of 60 assigned |
|---|---:|---:|
| Original `SurviveToPupation.csv` aggregated count | **30** | **36** |
| Original `2324_GD_PupalDays.csv` individually measured pupae | **36** | **30** |
| Original `2324_GD_PupalMass.csv` individually measured pupae | **36** | **30** |

The 3 other experimental blocks match their treatment-specific survival and individually recorded pupal counts exactly. Total pupae **262** match across all three files. **This is not proof the treatment effect is real or a treatment switch did occur**: the original data contain conflicting assignment labels, either in the summary or the underlying individual pupal files. We never silently relabel original records. The first CI failed closed upon detecting it. A separate frozen [two-scenario amendment](../../provenance/reviewer_defenses/protocols/MONARCH_LEAF_AGE_PUPATION_STAGE_CONVERSION_20261026.md) was committed before re-estimation, and corrected CI succeeded. Until original study authors resolve this discrepancy, it is **undetermined which treatment-specific pupation estimate is justified**.

## Exact two-scenario original-data audit

Mantel–Haenszel common odds ratio (Top/Bottom), stratified by 4 real YearTrial blocks. The conditional exact p-value is a 4-hypergeometric convolution fixing treatment group size and total pupae per block. Its calculation does not make four trial environments independent biological replicate years. A **10,000-resample YearTrial-cluster bootstrap** carries the four whole trial blocks and shows between-trial uncertainty.

| Source-trust scenario (both preserved) | Top pupated | Bottom pupated | Trial-stratified OR | 95% four-block bootstrap OR | Fixed-margin exact two-sided p |
|---|---:|---:|---:|---:|---:|
| **S0: use original survival summary exactly** | 140/240 (58.3%) | 122/240 (50.8%) | **1.4024** | **0.868–2.213** | **0.0986** |
| **S1: hypothetical trust the 2 individual pupal files only for YearTrial20241** | 146/240 (60.8%) | 116/240 (48.3%) | **1.7745** | **1.449–2.213** | **0.00478** |

S0 source itself includes a negative Top−Bottom pupation difference in 20241 (−10 percentage points); S1 reverses that specific block to +10 percentage points while the other 3 remain unchanged. Leave-one-experimental-block-out ORs: S0 **1.210–1.898**, S1 **1.624–1.913**. Four-cluster uncertainty S0 Top−Bottom pupation absolute risk-difference bootstrap **−3.3 to +16.7 percentage points**; S1 **+8.3 to +16.7 points**.

These are sensitivity interpretations of an unresolved source discrepancy, **not** two independent experiments, and the smaller p-value is NOT a basis for selecting S1. Under the original summary as recorded there is no convincing treatment-wide success advantage given four real trial blocks. Even if original authors correct YearTrial20241 in S1's direction, the finding belongs to the original Fisher et al. experiment and may not apply across other plants, landscapes or animal taxa.

## Performance of measured juveniles: a different life-history layer, not emergence

| Original source measured stage | Rows | Equal-trial Top minus Bottom mean among measured animals | Identification |
|---|---:|---:|---|
| Six-day larval mass | 270 | **+0.0160 source mass units** | 3/4 blocks only; 2023 Trial1 absent |
| Pupal mass | 262 | **+0.04025 source mass units** | Only caterpillars that reached pupation; treatment ambiguity in 20241 |
| Days to pupation | 262 | **−0.7961 days** | Only successfully pupated larvae; not time to adult eclosion |

From separately sampled original leaves, mean nitrogen content **4.824% Top vs 3.068% Bottom** (10 samples each). Top leaf force mean **0.730**, Bottom mean **0.893** (108 original leaf observations each). These are *leaf phenotype descriptions* with no verified individual caterpillar mediator linkage. There is no evidence in this source that the tested treatment changed **adult emergence**, adult population size, host plant arrival or interspecific competition.

Thus a butterfly can show an advantage in a measured **growth** endpoint among sampled larvae without necessarily showing a well-determined **survival-to-pupation** advantage across independent trial blocks. Source treatment-label contradiction further prevents an unqualified pupation conclusion.

## Scientific conclusion for chocho

The general biological mechanism of interest is a **conversion filter**, not "introduced hosts create more butterflies":

`Mapped host opportunity → observed larval feeding → juvenile growth/condition → pupation → adult emergence/recruitment`.

Our geographic macro paper measures the first stage. Caterpillars Count documents broad Lepidoptera branch-level larval occupancy but not butterfly species-specific recruitment; Fisher measures named-butterfly growth and pupation in an experimentally controlled leaf-age context with an unresolved YearTrial source-label discrepancy; the Illinois source measures adult count fluctuations but not host biomass or newly introduced host arrival. **These are distinct cohorts, species assemblages, sampling grains and mechanisms, not one joined causal sequence.**

Do not multiply observed potential resource additions by the Fisher laboratory treatment survival rate to infer global butterfly population changes. No unsupported synthetic product/experiment results in chocho GEB. Future original ecological test needs independently recorded local milkweed leaf quantity/age, herbivore egg/larval fate and adult emergence across multiple host types and sites, plus control for year/season and competitors.

## Reproducibility and branch boundary

- [Fixed posthoc protocol and discrepancy amendment](../../provenance/reviewer_defenses/protocols/MONARCH_LEAF_AGE_PUPATION_STAGE_CONVERSION_20261026.md)
- [Source-locked Python code](../../scripts/analyze_monarch_leaf_age_stage_conversion_20261026.py)
- [Machine-readable per-block/2-scenario source receipt](../../provenance/reviewer_defenses/results/monarch_leaf_age_pupation_stage_conversion_20261026.json)
- [Exact original-source GH Actions success](https://github.com/zuizui0223/chocho/actions/runs/40321446363)
- Source six original CSV Git blob SHAs included in machine receipt. No source raw individual IDs committed.
