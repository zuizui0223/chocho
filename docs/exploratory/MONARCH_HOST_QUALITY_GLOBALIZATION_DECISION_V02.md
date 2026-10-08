# Monarch larval host quality versus anthropogenic plant geography: decision v0.2

**Date:** 2026-10-08. **Scope:** exploratory `analysis/resource-homogenization-v01` only. Main GEB manuscript unchanged.

## Biological question and data provenance

Are host-plant species independently classified as **high** or **low** performance for monarch larvae systematically redistributed into more non-native botanical regions, after controlling the number of regions in their native range?

- Original performance evidence: Greenstein, Steele & Taylor (2022), PLOS ONE DOI [10.1371/journal.pone.0269701](https://doi.org/10.1371/journal.pone.0269701), **S1 Appendix Table 2** (127 classified plants: H=34, L=42, N=33, U=18). Original supplement downloaded in GitHub Actions [run 37736225347](https://github.com/zuizui0223/chocho/actions/runs/37736225347).
- Original botany: rWCVPdata pinned commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`. Intro and native WGSRPD3 codes, no inferred local host use.
- Frozen analysis protocol `MONARCH_HOST_QUALITY_GLOBALIZATION_PROTOCOL_V01.json`; before final WCVP crosswalk outcome. Correction note `MONARCH_HOST_QUALITY_WCVP_ACCEPTED_NAME_CORRECTION_V02.md`, written **before rerunning** with corrected matches.
- Reproducible corrected source + 9,999-permutation run: [37745236295](https://github.com/zuizui0223/chocho/actions/runs/37745236295), artifact `butterfly-monarch-host-quality-globalization-v01`. Source-matched output includes all 127 names and reasons for any excluded taxon.

## Critical matching repair

Initial successful [run 37744576643](https://github.com/zuizui0223/chocho/actions/runs/37744576643) matched only 108/127 because it declared 16 names ambiguous wherever exact accepted and synonym spelling rows had multiple destination IDs, including `Asclepias curassavica`. Its original primary p=0.5561 is an **incorrectly truncated target panel** and is retained only for audit.

A deterministic correction prioritizes a **unique exact accepted species record**, and uses a synonym only when the accepted ID is unambiguous. Corrected crosswalk matches **124/127** original plants, with 3 unmatched (Brassica oleracea var. capitata, Gossypium arboretum, Pachycarpus grandifloras). The key positive control `Asclepias curassavica` matches WCVP accepted ID 500848, 41 native and 97 introduced regions. Original H3 evidence includes adult emergence on this plant.

## Results

| Test (intro count controlled for native region count) | Distinct resolved taxa | H / L taxa | Permutation p, two-sided |
| --- | ---: | --- | ---: |
| **Asclepias high vs low, primary** | 55 | 26 / 29 | **0.5644** |
| Asclepias H3 vs L3, strong evidence sensitivity | 35 | 9 / 26 | 0.5972 |
| Asclepias high vs low excluding `A. curassavica` | 54 | 25 / 29 | 0.5376 |
| Apocynaceae high vs low, descriptive | 72 | 34 / 38 | 0.6303 |

Only four of the 55 Asclepias species had any recorded introduced region: `A. curassavica` (H3, 97), `A. syriaca` (H3, 28), `A. incarnata` (H3, 1) and `A. speciosa` (H3, 1). All four are high classified, **but both high and low groups have median zero introduced regions**, the observed class contrast is not unusual under the native-breadth-constrained randomization, and no predictive quality effect is established. Interpreting four introduced taxa as 55 independent introduction events would be misleading.

Source quality categories are graded **evidence**, not standardized experimentally measured survival probabilities (H1/H2 may have less direct evidence than H3). This is not proof that any WCVP region contains a viable local monarch population. Taxonomic/phylogenetic and botanical-introduction biases remain.

## HOSTS knowledge-selection diagnostic

Among the original high-quality Greenstein species, **23/34** have an exact accepted-name match in the frozen monarch HOSTS accepted sidecar; among low-quality species only **8/42** do. The cross-class record representation odds ratio is **8.89** (two-sided Fisher exact p≈2.22×10^-5). This diagnoses strong **ascertainment/host-list representation bias**, not positive selection on butterfly fitness. The HOSTS list nonetheless includes some low-performance and nonhost records, and omits the widely naturalized H3 `Asclepias curassavica`.

## Decision

1. No reliable evidence that higher-quality larval hosts are preferentially spread globally **conditional on native plant range** in this data.
2. Do **not** substitute this source-quality association for a causal ecological mechanism, or claim globalization improves larvae' fitness.
3. Historic botanical occurrence and citizen-science garden records established a real cultivated-host mapping gap, but a 915-larva independent source-field probe found no explicit `Phoebis×S.polyphylla` feeding link; see `CULTIVATED_SENNA_ORIGINAL_SOURCE_DECISION_V04.md`. These negative **coverage** results are not proof of biological absence.
4. Return to genuine **demographic response** instead of more static geography permutations: a separate, two-year common-garden `Asclepias speciosa`/ `A. fascicularis` adult-emergence interaction test is frozen in `MONARCH_H3_HOST_ENEMY_SURVIVAL_PROTOCOL_V01.json`; its [source and analysis workflow](https://github.com/zuizui0223/chocho/actions/workflows/butterfly-monarch-H3-host-enemy-survival.yml) was launched. The host-dependent growth effects are **already published by Diethelm et al. (2026)**, DOI [10.5061/dryad.51c59zwgx](https://doi.org/10.5061/dryad.51c59zwgx). Do not claim a new discovery without an independently supported, distinct outcome.

**No modifications to the main GEB manuscript.**
