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


## 2026-10-26 follow-up: published survival percentages and test method source-checked

The article's publicly accessible [Section 2.6 and Section 3.2](https://www.mdpi.com/2075-4450/17/2/215) report pupation rates **60.8% Top and 48.8% Bottom**, χ²=1.3139 and p=0.2517, and no statistically detectable survival difference. In the exact original author R script (Git blob de410164236730994c559bcf115a1e64e187be18), the corresponding function call is literally `chisq.test(c(48.8,60.8))`, and the Day-6 survival call is `chisq.test(c(74.2,75.4))`.

An [independent fixed-source statistical-method audit](../../provenance/reviewer_defenses/results/fisher2026_published_percent_chisq_source_audit_20261026.json) reproduced the reported pupation χ² **1.3138686** and p **0.2516954** exactly as a **chi-square goodness-of-fit test whose two input frequency values are percentages**, not treatment-specific numbers surviving/failing out of **240 assigned Top and 240 assigned Bottom larvae**. The original R code thus does **not implement the two-treatment sample-denominator-aware pupation survival test described in the paper's Methods**. This is a specific methodological source mismatch, not a judgement of the authors' intent.

The actual original counts do not uniquely reproduce the paper's published rates:

| Original data source | Top pupated / 240 | Bottom pupated / 240 | Correct count-based pooled Fisher exact two-sided p | Four-block conditional exact p |
|---|---:|---:|---:|---:|
| S0, literal SurviveToPupation.csv | 140 (58.3%) | 122 (50.8%) | **0.1190** | **0.0986** |
| S1, hypothetical correction supported by both individual pupal files | 146 (60.8%) | 116 (48.3%) | **0.00778** | **0.00478** |
| Published prose/code (`chisq.test(c(48.8,60.8))`) | implies 146 when Top denominator =240 | **48.8% is neither source count; closest 117/240** | **not computed from original individuals** | **not computed from original four blocks** |

The original published p=.2517 therefore cannot resolve the original source contradiction; it reproduces from percentage numbers even if the true treatment numerator labels are different. Treat p-values in the second and third columns as **different statistical objects**, not as proof S1 is the correct sample. Only original authors can establish whether `SurviveToPupation.csv` 20241 had inverted labels, or the two raw pupal outcome files did, and why the published Bottom survival 48.8% does not exactly match either.

The original 2026 article's **growth and earlier pupation** results remain independently measured among relevant source animals; the source-integrity finding applies to treatment-specific survival-to-pupation inference, and does not overturn unrelated original experiments. The scientific scope remains leaf age of *A. syriaca* within a controlled *D. plexippus* study, **not human-introduced larval hosts causing new butterfly populations**.

**Action-ready neutral questions for study authors (not sent):**
1. Which YearTrial=20241 Top/Bottom pupation assignment is correct: 30/36 in original survival summary, or 36/30 in both pupal individual-stage files?
2. What original successes and denominators give Bottom 48.8% quoted in the article, given pupal stage records show 116/240=48.3%?
3. Would a reproducible survival-to-pupation comparison based on the original binary successes/failures and four trials, rather than the two percentages passed to `chisq.test`, better reflect the design?

[Published-method audit code](../../scripts/audit_fisher2026_published_percent_chisq_20261026.py) · [fixed audit protocol](../../provenance/reviewer_defenses/protocols/MONARCH_FISHER_PUBLISHED_PROPORTION_TEST_SOURCE_AUDIT_20261026.md) · [successful Actions audit](https://github.com/zuizui0223/chocho/actions/runs/40326629096).

**Stop rule:** Preserve all source files and both statistical scenarios. Do not silently "correct" treatment labels or tell readers a firm survival effect exists pending original author/source verification. No additional parameter sweep; identifying plant quality effects in field populations still needs actual local host availability and tracked juveniles.


## Distinct individual-stage finding: day-six larval mass predicts recorded later pupation within original treatment strata

**Posthoc within-original-source prognostic test, 2026-10-26.** The original article quantified mean six-day mass, mean time to pupation and source treatment survival separately. We located an actual *individual larva ID* bridge in original `2324_GD_AllMass.csv` (Git blob `3905f64bfe2667ce566a09d291e7c848378a11bf`) and validated each mass against the two independent **file copies** `2324_GD_6DayMass.csv` and `2324_GD_PupalMass.csv` from the same original experiment. These are not independent biological experimental replications.

The cross-file agreement is exact: **532 measured stage records**, of which **270 larvae were weighed on Day 6** and **262 pupae were weighed**. The oldest trial `20231` lacks Day-6 weights entirely, but contains 72 recorded pupae. In the other three trials, **all 190 subsequent pupal IDs appear among the 270 Day-6 measured larvae**, leaving 80 original Day-6 weighed larvae with no later pupal record. This gives a direct within-cohort stage **measurement link**, not actual adult eclosion. Individual IDs were linked by `YearTrial×Treatment×ID`, so the same numeric ID in another trial is never misjoined.

The [fixed source-fate protocol](../../provenance/reviewer_defenses/protocols/FISHER_MONARCH_DAY6_MASS_TO_PUPATION_PROGNOSIS_20261026.md) set the primary prognostic statistic as an **equal-weight average of within-six-(trial×recorded leaf)-stratum Mann–Whitney AUCs**. It asks whether a Day-6 larva that is *later seen pupating* had higher Day-6 mass than a same-trial, same-recorded-leaf larva without a pupal record. AUC=0.5 is no ranking information. Group-specific measurement is essential because 2023/2024 larval masses differ greatly.

| Source YearTrial × recorded leaf label | Day-6 weighed | Pupa ID matched | Within-stratum Day-6 mass prognostic AUC |
|---|---:|---:|---:|
| 2023 Trial 2, Top | 34 | 19 | 0.7965 |
| 2023 Trial 2, Bottom | 35 | 15 | 0.7800 |
| 2024 Trial 1*, Top | 42 | 36 | 0.6806 |
| 2024 Trial 1*, Bottom | 46 | 30 | 0.7313 |
| 2024 Trial 2, Top | 59 | 50 | 0.5911 |
| 2024 Trial 2, Bottom | 54 | 40 | 0.5179 |

Equal-six-strata primary AUC = **0.68288**. A fixed-null 20,000 permutations (seed 20261026) shuffling recorded pupation labels **within each original trial×leaf group** while preserving its original count produced null 95% `[0.4184,0.5844]` and two-sided exploratory conditional tail **0.00005** (none of 20,000 random shuffles was at least as extreme). Mean AUC by true YearTrial varies **0.7882, 0.7059, 0.5545**: the third experimental trial shows only a weak rank relationship. Thus a within-sample association is strong against the *individual-level* conditional independence null, but generalization beyond **three real Day-6 trials** is unestablished.

**Crucial source-identity robustness:** Swapping the **names** Top and Bottom in the disputed `20241` source group, without moving an animal's ID, Day-6 mass or pupation observation, changed the primary equal-six-strata AUC by only machine rounding (~2.2×10⁻¹⁶). Thus the **individual-state→subsequent-recorded-stage association** is identifiable *within the recorded source* even though the **experimental treatment-specific pupation advantage** remains unresolved. This is a *prognostic* result: Day-6 mass is measured **after** the leaf-age treatment and may reflect condition, genetics, handling and unmeasured larval vulnerabilities. We **cannot infer that making a larva heavier causes survival**, that plant nitrogen mediates the link, or that the finding predicts regional butterfly population increases.

[Source-checked numeric result](../../provenance/reviewer_defenses/results/fisher_monarch_Day6_mass_pupation_ID_AUC_20261026.json) · [Exact individual linkage/permutation implementation](../../scripts/analyze_fisher_monarch_day6_mass_to_pupation_AUC_20261026.py) · [Successful original-source workflow](https://github.com/zuizui0223/chocho/actions/runs/40328362969).

**Biological conclusion:** Candidate host presence, observed larval feeding, Day-6 juvenile condition, later pupation, and successful adult eclosion are separate evidence levels. Only the Day-6→pupation association within the subset can be assessed directly with individual IDs here. Neither source label contradictions nor strong rank concordance justify multiplying a global potential host-region gain by a fitted laboratory survival probability.
