# Candidate GEB Supplementary Methods and Table S9 — independently compiled host-inventory sensitivity

**Status (2026-10-08): EDITORIAL CANDIDATE ONLY — NOT YET INCLUDED IN THE GEB SUBMISSION BRANCH.** This is optional text for scientific/author review; it does not revise the frozen 239-species manuscript or claim new functional host-use measurements. Its stricter continent-constrained null results were printed in the finished calculation for [run 37784230410](https://github.com/zuizui0223/chocho/actions/runs/37784230410), but that run failed during a non-statistical output-metadata lookup; the successful final JSON artifact is pending [rerun 37786325837](https://github.com/zuizui0223/chocho/actions/runs/37786325837).

## Supplementary Methods — host-interaction knowledge sensitivity (candidate)

Because butterfly–host compilations are incomplete and differ in their inclusion criteria, we assessed the dependence of our resource-geography reconstruction on independently **compiled** larval foodplant knowledge, rather than assuming the original HOSTS database is exhaustive. A European butterfly host checklist (*Clarke* 2024; DOI [10.1002/ece3.10834](https://doi.org/10.1002/ece3.10834)) is presented by Butterfly Conservation Europe as species foodplant lists retaining evidence ranks 1–3 (wild feeding and rearing, wild larvae, or wild oviposition). We retrieved these **derived** species pages, not the original relational Dryad tables, whose automated source download was unavailable. Among the original 239 butterfly species, 92 names matched BCE exactly and 83 yielded plant species-level foodplant lists. We reconciled the additional names against fixed WCVP accepted taxon identities, recovering **1,027 accepted-ID butterfly–host pairs absent from the original HOSTS compilation** for **81 focal butterflies**, involving **679 accepted plant species**. These are published *candidate host pairs*, not independent feeding observations. We added those pairs to a *separate* sensitivity reconstruction, preserving the other 158 butterfly host lists and leaving the main dataset unchanged.

We recalculated plant-native and contemporary butterfly × WGSRPD3 resource opportunities, verifying all original 239 butterflies' frozen resource counts first. For regional similarity, we retained the *same 355 original native-resource regions* in both conditions to exclude changes in the regional analysis domain. We calculated mean pairwise Jaccard similarity of resource assemblages and sampled separately scenario-conditioned fixed-margin nulls (499 permutations; seed 20261007) that preserve each butterfly's and each region's added opportunity counts, and additionally used within-WGSRPD Level1 swaps to preserve each butterfly's introduced-added counts by continental group. The original 239-species/355-region baseline and all four original Jaccard quantities were reproduced exactly.

## Supplementary Table S9 — knowledge sensitivity and regional resource similarity (candidate)

| Measure | Original HOSTS | BCE/Clarke-derived additional host pairs |
| --- | ---: | ---: |
| Evaluated focal butterfly species | 239 | 239 (81 host lists augmented) |
| Total native butterfly × WGSRPD3 resource units | 26,530 | 29,852 |
| Total contemporary resource units | 41,083 | 44,819 |
| Introduced-added resource units | 14,553 | 14,967 |
| Relative expansion over native baseline | 54.85% | 50.14% |
| Butterflies with any resource expansion | 206 | 213 |
| Regional mean Jaccard — native, fixed 355 regions | 0.27698 | 0.29938 |
| Regional mean Jaccard — contemporary, same 355 regions | 0.46208 | 0.48506 |
| Unconstrained added-row-and-column fixed-margin null median | 0.45348 | 0.47716 |
| Regional observed-minus-unconstrained-null median | +0.00860 | +0.00791 |
| One-sided Monte Carlo p (499 draws) | 0.002 | 0.002 |
| WGSRPD Level1-constrained null median | 0.45847 | 0.48162 |
| Regional observed-minus-Level1-null median | **+0.00362** | **+0.00344** |
| One-sided Monte Carlo p (499 draws) | 0.002 | 0.002 |

*Important denominators:* all-regions aggregate expansion counts use all 239 resource footprints. Regional Jaccard statistics condition on the 355 original native-active WGSRPD3 regions for comparability; the new 356th native-active region is excluded from the paired Jaccard analysis. Accordingly, these are **not** two independent estimates of actual local host utilization. The small positive null excess measures host-resource **opportunity**, not observed butterfly population homogenization or competition.

## Candidate restrained Results wording

> A secondary, evidence-filtered European host checklist exposed substantial knowledge dependence of the original butterfly–host inventory: 1,027 additional accepted species-level candidate host associations among 81 butterflies. Including these relationships increased native resource opportunity more than introduced-added opportunity, reducing the reconstructed proportional expansion from 54.9% to 50.1%. Nevertheless, regional resource-assemblage convergence remained above fixed-margin expectations under a 355-region matched comparison, including a stricter continental-configuration-preserving null (observed-minus-null median +0.0034 versus +0.0036 in the original; both permutation p=0.002). This is a source-specific, European taxonomic-completeness sensitivity and cannot be interpreted as a corrected worldwide estimate or local evidence of larval development.

## Editorial gate

Do **not** insert this candidate into the blinded GEB submission until (i) the complete successful [rerun 37786325837](https://github.com/zuizui0223/chocho/actions/runs/37786325837) verifies the Level1 artifact and (ii) the authors decide whether to include substantial post-hoc, Europe-only host-source information in Supplementary Information. If promoted, carry all limitations into the supplement and Data/Code Availability, add an original publication bibliographic reference, and rerun paper CI and anonymous-bundle preflight. Do not increase the main text above the GEB word limit.
