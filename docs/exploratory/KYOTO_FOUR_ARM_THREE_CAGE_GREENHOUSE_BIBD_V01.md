# 京都蝶：3ケージ／温室で4処理を実験する釣り合い型不完備ブロック設計

**2026-10-14｜SOURCE-GROUNDED DESIGN REPAIR, NO NEW ANIMAL EXPERIMENT.**

**Exploratory branch** `experiment/kyoto-4arm-three-cage-bibd-v01`, derived from `experiment/kyoto-six-arm-identity-food-v01`. The main GEB submission paper and PR #38 are unchanged. This is a prospective candidate design, **not a registered, powered or authorized experiment**, and the 2012 physical facility is NOT known to be available now.

## The exact biological question — not a statistical methods paper

> **Can making fresh edible leaves continuously accessible rescue the reduction in native *Atrophaneura alcinous* adult production caused by newly introduced *Sericinus montela* larvae on the shared native host *Aristolochia debilis*?**

Hashimoto & Ohgushi (2023) [DOI 10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164) established suppression of the native butterfly in a density experiment; it did not time-resolve accessible young/mature host biomass and had no randomized leaf-access intervention. Its original **2012 setup held 30 whole experimental cages in 10 greenhouses, each with only 3 cages**, with multiple host plants nested in each cage. The present research **does not have a confirmed current greenhouse inventory**.

An earlier experimental protocol allowed 4 complete-treatment blocks, requiring *four cages in every independent block*. That layout is impossible **within the original three-cage greenhouse**. It does not follow that the four-arm causal question is nonidentifiable under the historical size restriction. The answer is a connected **balanced incomplete block design (BIBD)**, standard statistical experimental-design theory (Purdue STAT514 lecture notes on connected incomplete blocks; Penn State STAT503 Section 4).

## 1. Exact BIBD combinatorics

- 4 treatments `v=4`.
- 8 greenhouses `b=8` (each preselected as legally and physically eligible for 3 independent whole cages).
- Each greenhouse `k=3` different treatment cages; **one of the four is omitted**.
- Each arm used exactly `r=6` cages.
- Each distinct pair of arms shares `lambda=4` greenhouses.

Identities: `bk=8×3=24=4×6=vr`; `lambda(v−1)=4×3=12=r(k−1)=6×2`.

| Greenhouse schedule (schematic, not literal assignment) | Three arms included | Arm deliberately omitted |
|---|---|---|
| G1 and G5 | NATURAL_NONE, NATURAL_ADD_S, SUPPORTED_NONE | SUPPORTED_ADD_S |
| G2 and G6 | NATURAL_NONE, NATURAL_ADD_S, SUPPORTED_ADD_S | SUPPORTED_NONE |
| G3 and G7 | NATURAL_NONE, SUPPORTED_NONE, SUPPORTED_ADD_S | NATURAL_ADD_S |
| G4 and G8 | NATURAL_ADD_S, SUPPORTED_NONE, SUPPORTED_ADD_S | NATURAL_NONE |

**Actual allocation is randomized twice**: assign the eight balanced arm-omission labels to the eight prescreened greenhouses using one recorded seed; then permute the remaining three treatment labels across the three genuinely eligible whole-cage IDs in each greenhouse. No selection of omitted arms after observing caterpillar survival, condition of leaves or temperature. Where ten greenhouses exist, the choice of the eight main-study greenhouses must be fixed **before randomization**, blinded to outcomes; two can be prospectively designated for a separate pilot but cannot be casually pooled as main-study replicates.

A single greenhouse is not independently assigned to one treatment. Its three cages are the treatment units, with **greenhouse adjustment as blocking**. Cages are not interchangeable replicates of a greenhouse-level intervention. This design needs actual independent cage containment, foliage handling and no treatment leakage.

## 2. Statistical estimability and efficiency: verified algebra

Original response `Y`: the count of **flight-capable A adults** observed by the prespecified census date divided by the entire **initially assigned focal A cohort** in each cage, including zeros. Living diapause or pending pupae are not called dead and remain as possible eventual successes in upper-bound diagnostics.

Primary four-arm contrast:

```text
D = (Y_SUPPORTED_ADD_S − Y_SUPPORTED_NONE)
  − (Y_NATURAL_ADD_S − Y_NATURAL_NONE)
```

The estimator fits **greenhouse fixed effects** and three treatment indicators to **whole-cage risks** (one row per cage), never 24*larval_n independent rows. The contrast is estimable because all arm-pairs co-occur and the treatment incidence matrix has rank 3 after greenhouse adjustment.

The adjusted 4×4 treatment information matrix is

```text
C = diag(6,6,6,6) − (N' N)/3
  = (16/3) I − (4/3) J.
N' N has diagonal 6 and off-diagonal 4.
```

For the signed `c=(+1,−1,−1,+1)`, **assuming homoscedastic independent residuals of variance sigma² at the cage level after greenhouse adjustment**:

```text
Var(c' treatment_estimate) = sigma² * c' pinv(C) c = 0.75 sigma²
Residual df = 24 − 8 greenhouse indicators − 3 treatment indicators = 13.
```

Comparison with **hypothetical** 6 complete greenhouse blocks of four cages each (also 24 cages) gives variance `(4/6)sigma² ≈ 0.6667sigma²`; thus the 3-cage BIBD incurs **12.5% higher variance**, or approximately **88.9% relative efficiency**, under this *specific homoscedastic fixed-block model*. This is a conditional design calculation and not an observed source outcome.

**Avoid a false “24 cages are enough” claim.** With purely assumed within-greenhouse residual SD = 0.20 of an adult-success **proportion**, the BIBD standard error is `sqrt(0.75)*0.20≈0.173`. A nominal Student-t `df=13` half-width is `t(13,0.975) × 0.173≈0.374`, very broad. This **does not provide 95% coverage without normal/independent/homoscedastic cage residual assumptions**, and no empirical residual SD/power are available. Binomial larvae in one cage cannot be counted as independent cages.

If greenhouse-specific treatment responses exist, greenhouse effects are nonadditive or a common original maternal family links blocks, the nominal model interval is unreliable. Source-family/donor batches should be blocked or balanced prospectively and a genuinely independent variance/recruitment pilot is required. The code provides leave-one-greenhouse-out estimator sensitivity **as descriptive diagnostics**, not 8 independent statistically powered replication proofs.

## 3. What the implemented analysis actually enforces

- Unique original `cage_id` and `greenhouse_id`; exact **8×3** randomization, 6 per treatment, omission balanced 2 per arm and all 6 treatment pairs sharing exactly 4 greenhouses. An actual 10×3 inventory is **rejected** rather than silently dropping 2 greenhouses.
- Explicit original `focal_a_initial_n` unchanged over the trial, separate original added `S. montela` individuals and fates in relevant arms.
- A living `diapause_alive` or pending pupa is **not** coded dead, and all unknown/possible eventual adult fates are retained for **signed, worst-case adult contrast bounds** using exact OLS contrast weights (not a naïve all-unknown-as-success difference-in-differences bound).
- Fixed greenhouse effect estimates on 24 whole-cage adult risks; 13 residual degrees of freedom; analytic BIBD contrast variance factor 0.75; treatment-specific observed adult percentages separately reported.
- Repeated actual edible host leaf area and stage through 3rd–5th instars, current versus supplemented leaves, leaf donor source and sham handling, with missing area different from genuinely zero accessible area. **Never select cages based on successful supplementation**, high residual food, or native larvae surviving late instars.
- Source hash/seed and precise input identity, no archive data copied or false 2012 recruitment estimates.

Even a positive primary `D` is the effect of an assigned **food-access-support strategy** and does not prove leaf mass alone mediated the rescue: supplemented young leaves can change chemistry/quality, humidity, encounter opportunity or plant age exposure. A secondary manipulation of leaf stage or prior damage is needed to discriminate those mechanisms.

## 4. Future real-data input contract

Run the randomizer only after local species/rearing/handling permissions, full containment of introduced S and **current actual greenhouse capacity** are independently verified. A possible source list of eight approved greenhouses with three usable whole cages apiece uses precisely:

```text
eligible_whole_cages.csv
cage_id,greenhouse_id
```

```bash
python scripts/randomize_kyoto_four_arm_greenhouse_bibd.py \
    --inventory eligible_whole_cages.csv \
    --allocation assigned_24_cages.csv \
    --receipt original_randomization.json \
    --seed <integer_recorded_before_outcomes>
```

The prospective analysis needs the **following four original CSVs**, not constructed data:

```text
allocation.csv
cage_id,greenhouse_id,treatment,competitor_present,food_support,random_seed,design_version

assignments.csv
cage_id,greenhouse_id,treatment,competitor_present,food_support,focal_a_initial_n,added_s_initial_n,expected_n_food_visits,food_floor_cm2,donor_batch_id

fates.csv
cage_id,larva_id,species,cohort_role,fate

resources.csv
cage_id,observation_datetime,focal_a_alive_n,focal_a_late_instar_alive_n,added_s_alive_n,accessible_edible_leaf_area_cm2,accessible_leaf_added_area_cm2,leaf_development_stage,added_leaf_stage,leaf_source_batch_id,sham_handling_done
```

`cohort_role` must be `focal` (`Atrophaneura alcinous`) or `added` (`Sericinus montela`). Every originally assigned focal/additional larva has a fate row. Fates include flight-capable adult, confirmed immature/pupal death, nonflight adult, living diapause/pending pupa and genuine unresolved follow-up. Timestamps have timezone offsets; measured zero foliage is `0`, missing observation is an empty field, not zero.

```bash
python scripts/analyze_kyoto_four_arm_greenhouse_bibd.py \
    --allocation allocation.csv \
    --assignments assignments.csv \
    --fates fates.csv \
    --resources resources.csv \
    --receipt results/kyoto_BIBD_original_fates.json
```

**No initial cage inventory, measured leaf supply, butterfly fates or intervention permission exists in this branch. All code tests use explicit artificial cage and larval records.**

## 5. Ecological decision rule

This fixes one **real feasibility blocker** from the original 3-cage greenhouse arrangement. It does not validate the hypothesis. 24 cages are unlikely to give a narrow main-effect interval unless genuine between-cage residual variability is small; a real variance pilot is critical.

- **GO for new feasibility work**: preverified set of 8 suitably isolated 3-cage greenhouses, independent source families/donor plant materials, reliable repeated edible leaf access measurement, originally assigned native larvae tracked to adult or alive diapause, introduced species contained and lawful handling.
- **HOLD**: 10×3 can be confirmed only from published historical facilities, not currently eligible; no independent prior SD; mother/donor history or greenhouse conditions cannot be balanced.
- **STOP/revise**: source blocks are not genuine (shared units or carryover), or the feasible allocation collapses all treatment groups onto different greenhouse environments. Then use a different predeclared blocked design or independent 4-cage greenhouse/scheduling, not an outcome-based redesign.
- If a leaf-stage pilot confirms resource use depends on young/mature tissue composition, treat it as a **follow-up ecological mechanism experiment**, not a post-hoc explanation forced by the same support manipulation.

**Main finding supported by publications:** nontrivial asymmetric native-invader butterfly competition, and prior neonate leaf-stage differences. **New contribution now:** physically feasible conditional allocation and original-outcome analysis *method*, not an empirically supported new fitness result. Keep the separately frozen global GEB potential host-resource paper untouched.
