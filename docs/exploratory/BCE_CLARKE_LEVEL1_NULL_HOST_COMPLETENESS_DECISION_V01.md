# BCE/Clarke accepted host additions under a continent-preserving resource-homogenization null

**Date:** 2026-10-08. **Status:** COMPLETE VERIFIED GITHUB OUTPUT AFTER NONSTATISTICAL METADATA REPAIR; ADDITIONAL DETERMINISTIC-ORDER REPRODUCIBILITY REGRESSION IN PROGRESS. **Post-hoc resource-network sensitivity, not ecological causality.** Frozen GEB paper PR #38 unchanged.

## Scientific question

The GEB reconstruction from 239 butterflies, fixed HOSTS and plant native/introduced ranges showed higher regional potential butterfly-resource similarity than a fixed-margin null. A separate BCE/Clarke evidence-ranked (1–3) foodplant checklist identifies **1,027 candidate accepted-plant-ID host edges** missing from the pinned HOSTS for **81 European butterfly species**. The one-direction addition alters the global 239-butterfly expansion ratio **54.85% → 50.14%** and the native resource baseline **26,530 → 29,852** butterfly × WGSRPD3 units. A conventional fixed-row/column-margin null on the fixed original 355-region domain retained small nonrandom homogenization excess **+0.008604 → +0.007906**, both one-sided Monte Carlo p=0.002 (run [37781612222](https://github.com/zuizui0223/chocho/actions/runs/37781612222)).

This follow-up asks whether that residual is explained by the broader *continent-level* pattern of butterfly-specific introduced resource additions.

## Frozen design and inputs

The source-blind protocol `BCE_CLARKE_LEVEL1_NULL_HOST_COMPLETENESS_PROTOCOL_V01.json` was committed **before** the strict null outcome. It fixes:
- **239 exact butterfly species and the same original 355 native-active WGSRPD3 regions in both scenarios**, excluding the one newly native-active BCE region from the paired comparison;
- original HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`, WCVP commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`, WGSRPD Level3/Level1 mapping commit `52da7828aba9d461dd133c27b3bd7a4407161f54`;
- no cherry-picked hosts: all same 1,027 candidate accepted-ID additions for their originally documented butterflies; no new links for remaining species;
- **499** conditional double-edge-swap samples with seed `20261007`, burn-in `max(10000,15×number of added links)`, interval `max(2000,3×number of added links)` and eligible swaps drawn **within the same WGSRPD Level1 unit**;
- invariants: preserve each butterfly's added-region count, added butterfly incidence per region, and each butterfly's added counts **within each Level1 region**, keeping native cells unavailable for introduced additions.

The new `scripts/analyze_bce_clarke_homogenization_level1_null.py` mirrors the original `scripts/analyze_butterfly_resource_homogenization_level1_null.py` method and explicitly tests row, column, species × Level1 and native-exclusion invariants. All **three synthetic unit tests passed** in the original job [37784230410](https://github.com/zuizui0223/chocho/actions/runs/37784230410). Original 239 resource sums and four frozen regional/species overlap observables were validated against reference before sampling.

## Final output and reproducibility resolution

The first workflow [37784230410](https://github.com/zuizui0223/chocho/actions/runs/37784230410) completed all 499 permutations per scenario but failed during final JSON creation because of an incorrect provenance-only key; calculations and output were printed. The corrected [GitHub Actions rerun 37786325837](https://github.com/zuizui0223/chocho/actions/runs/37786325837) **completed successfully**, uploaded the source JSON artifact `chocho-bce-clarke-continent-preserving-host-source-null-v01` (artifact ID **11554734608**) and passed all 3 unit tests. Its exact regional excesses are **+0.0036249015** (original HOSTS) and **+0.0034393111** (BCE-augmented), with one-sided Monte Carlo `p=0.002` for both. They differ by approximately 1×10^-5 and 1×10^-6 from the first calculation, which suggests hash/insertion order affected the exact pseudorandom null sampling despite identical seed and scenario. The algorithmic conclusion is unchanged. Follow-up deterministic edge/group sorting was implemented in [commit 858c347](https://github.com/zuizui0223/chocho/commit/858c347aa387e2c8cc9dafa5c71883351ae9eca9), with reverse-insertion order regression test in [commit 0398e53](https://github.com/zuizui0223/chocho/commit/0398e5372220d2c75573f6ea37e69e7db589e1ee). The stable-order rerun is separate; **the statistics below are from the verified successful 37786325837 receipt**, not the first failed workflow.
| Same 355 regions, each scenario's continent-constrained null | Original frozen HOSTS | BCE/Clarke accepted-host augmentation |
| --- | ---: | ---: |
| Native mean regional Jaccard | 0.2769818977 | 0.2993793713 |
| Contemporary mean regional Jaccard | 0.4620834005 | 0.4850609467 |
| Mean contemporary Jaccard under within-Level1 null (median) | **0.4584584990** | **0.4816216356** |
| **Regional observed−null excess** | **+0.0036249015** | **+0.0034393111** |
| One-sided Monte Carlo p (499 samples) | 0.002 | 0.002 |
| Butterfly resource geography observed−null excess | +0.0027159214 | +0.0028002081 |
| Butterfly resource geography one-sided Monte Carlo p | 0.002 | 0.002 |
| Accepted swaps across burn/samples | 1,023,872 | 985,706 |
| New resource butterfly × region incidences within fixed 355 regions | 14,121 | 14,444 |

The raw regional Jaccard gains are approximately **+0.18510** and **+0.18568**, respectively. Only about **1.95%** and **1.85%** of these raw changes are represented by the observed-minus-strict-null residual in this specific comparison; this is a descriptive effect-size ratio, not a variance explained decomposition. The stricter geographic null absorbs much of the nonrandom structure remaining beyond the ordinary row/column margin control:

| Conditional null | Original observed−null regional excess | BCE-augmented observed−null regional excess |
| --- | ---: | ---: |
| Butterfly × region margins only | +0.008604 | +0.007906 |
| Add butterfly × Level1 regional margins | +0.003625 | +0.003439 |

## Scientific inference

The stricter source-completion test still yields a positive, Monte Carlo-resolved structural residual under the frozen model, without requiring the original HOSTS-only foodplant inventory. **This supports limited robustness of potential butterfly-resource geographic homogenization beyond continent-preserving introduction opportunity margins**. It does **not** show that butterfly populations have actually homogenized, colonized introduced-host areas, improved fitness, undergone host switching or entered greater ecological competition.

The small residual is important: **most** reconstructed homogenization is accounted for by the amount, geographic concentration and continental allocation of potential resource opportunity. Only a small conditional structure remains. The BCE/Clarke source is Europe-biased and literature-derived; its 1,027 added accepted taxon edges are source-conditional candidate host-use links, not 1,027 newly witnessed developmental successes. The scenario-specific MCMC nulls condition on different margins, and their identical p=0.002 values do not prove equivalence of effects.

## Decision

Keep the GEB paper's primary **butterfly** question focused on anthropogenic larval-host redistribution expanding and homogenizing *potential resource geography*. The new result increases robustness of the **structure conditional on bibliography** but cannot promote it to a causally realized ecological mechanism. Avoid post-hoc new geography threshold searches or presenting the source-augmented results as corrected global estimates. Preserve original PR #38 and all original frozen figures. Subsequent independent ecology must measure local realized host access, development and colonization or demography.
