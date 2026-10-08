# Butterfly host-identity geographic selectivity: robustness decision v0.2

**Date:** 2026-10-08  
**Status:** EXPLORATORY / POSITIVE CONDITIONAL STRUCTURE, NOT YET A DISTINCT ECOLOGICAL DISCOVERY  
**Scope:** `analysis/resource-homogenization-v01` only; do not alter `main` or the GEB submission manuscript.

## Scientific claim under examination

Does human redistribution of known larval hosts preferentially connect distant regions with similar climates because of the **particular butterfly–host species identities**, beyond floristic redistribution and coarse butterfly/host phylogenetic association?

This is a potential-resource geography reconstruction, not observed butterfly colonization, larval use, competition or fitness.

## Verified strict conditional link null

- Frozen panel: 239 butterflies, 355 WGSRPD3 regions, 1,706 documented host plants, 2,596 butterfly–host links.
- Distant region pairs: different WGSRPD Level1 codes, ≥3,000 km apart.
- Four climate variables: z(BIO1), z(BIO7), z(log1p(BIO12)), z(BIO15).
- Climate analogues: bottom quartile of climatic distances; discordant: top quartile. Each group: 12,944 region pairs.
- Contrast: (contemporary − native butterfly-resource Jaccard gain in climate analogues) minus the same gain in discordant pairs.
- Strict link randomization: exact botanical host distributions fixed, double-edge switches **only within the same butterfly family × botanical host family**, preserving each butterfly's number of hosts in each plant family and each plant's consumer count within each butterfly family.
- Both independently seeded Markov chains used 499 permutations and accepted >1.13 million swaps.

| Run | Observed contrast | Strict null median | Excess | One-sided Monte Carlo p |
|---|---:|---:|---:|---:|
| Seed 20261008 | 0.098437 | 0.087899 | +0.010538 | **0.030** |
| Seed 20261009 | 0.098437 | 0.088901 | +0.009536 | **0.026** |

**Interpretation:** exact butterfly–host identity has a modest conditional association with climate-selective resource-geographic convergence even after preserving butterfly and botanical family affiliations. The effect is **borderline** and post-hoc; do not claim causal consumer response or novel realized trophic links.

Reproducible GitHub Actions runs:
- strict original: 37707320944 (artifact `butterfly-host-identity-lineage-null-v02`)
- six-family leave-one-out and independent seed: 37711152248 (`strict-host-identity-family-robustness-20261008`, `strict-host-identity-family-robustness-20261009`).

## Leave-one-butterfly-family-out audit

After excluding each complete butterfly family, all observed-minus-null effects retained positive sign under **both** seeds. This is a sensitivity audit, not six separate independent hypothesis confirmations.

| Omitted family | Butterfly species omitted | Excess seed 20261008 | p seed 20261008 | Excess seed 20261009 | p seed 20261009 |
|---|---:|---:|---:|---:|---:|
| Hesperiidae | 51 | +0.009422 | 0.068 | +0.008405 | 0.050 |
| Lycaenidae | 25 | +0.012603 | 0.016 | +0.011873 | 0.012 |
| Nymphalidae | 103 | +0.013175 | 0.030 | +0.012322 | 0.038 |
| Papilionidae | 16 | +0.011123 | 0.026 | +0.010229 | 0.026 |
| Pieridae | 43 | +0.006826 | 0.116 | +0.005934 | 0.156 |
| Riodinidae | 1 | +0.010577 | 0.030 | +0.009568 | 0.026 |

Conclusion: not a Nymphalidae-only artifact, but removal of **Pieridae** substantially weakens the signal. Since the response is a nonadditive, normalized Jaccard contrast, differences between leave-one-out estimates are *not* additive family contributions.

## Distinguish pre-existing host structure from anthropogenic increment

Recomputed the same strict null separately for native and contemporary resource geography, holding the same seed and the exact edge-swap schedule; the primary original strict p=0.030 was reproduced.

| Conditional contrast | Actual-minus-strict-null median |
|---|---:|
| Native host-resource climate selectivity | **+0.040536** |
| Contemporary host-resource climate selectivity | **+0.050996** |
| Additional selectivity in contemporary minus native | **+0.010538**, p=0.030 |

The exact native and contemporary median residuals are not algebraically additive because they are medians of different null distributions, but their difference is close to the directly tested incremental residual.

**Critical interpretation:** much of the host-identity-associated climate selectivity **predates** anthropogenic redistribution. Human introductions add a smaller positive contrast to an already structured butterfly host network. Do NOT present all contemporary network selectivity as newly caused by plant globalization.

Decomposition reproducibility:
- code: `scripts/audit_butterfly_host_identity_native_current.py`
- workflow: `.github/workflows/butterfly-host-identity-decomposition.yml`
- run: 37711546754 (GitHub replication running when this note was written).
- local exact-input result: same strict original p=0.030 and 1,134,859 accepted swaps with seed 20261008.

## The strongest previous results are not themselves novel consumer biology

1. Cross-continental climate-analog regional butterfly-resource assemblage convergence exceeded fixed-margin swaps by +0.03742 in the analogue-minus-discordant gain contrast (p=0.002, 499 swaps). However **Yang et al. (2021)** had already reported that naturalized plants homogenize geographically distant floras especially in climate-similar regions.
2. The host-plant-only assemblage descriptive analogue-minus-discordant gain contrast was about 0.0924 versus about 0.0984 in the butterfly-resource projection. These are different assemblage units, so their arithmetic difference is not a formal consumer-amplification effect.
3. The stricter within-butterfly-family and within-host-family network null finds a *smaller* residual +0.0105 (p≈0.03). This is the best positive consumer-resource-specific structural evidence here. It is not equivalent to realized food-web homogenization.

## Published prior art and boundaries

- Yang, Q. et al. (2021). *The global loss of floristic uniqueness*. Nature Communications 12:7290. https://doi.org/10.1038/s41467-021-27603-y — globalization/plant naturalization promotes floristic similarity over long distances especially for climate-similar regions. This **overlaps** the coarse geography claim.
- Daru, B.H. et al. (2021). *Widespread homogenization of plant communities in the Anthropocene*. Nature Communications 12:6983. https://doi.org/10.1038/s41467-021-27186-8 — large-scale taxonomic and phylogenetic plant homogenization.
- Fricke, E.C. & Svenning, J.-C. (2020). *Accelerating homogenization of the global plant–frugivore meta-network*. Nature 585:74–78. https://doi.org/10.1038/s41586-020-2640-y — global interaction-network homogenization established using realized local frugivore interactions, so network homogenization per se is **not** a new general phenomenon.
- *Unifying host-associated diversification processes using butterfly–plant networks* (2019). Nature Communications. https://doi.org/10.1038/s41467-018-07677-x — host associations are phylogenetically structured and not exchangeable between butterflies. This qualifies interpretation of link-swap nulls.

## Decision and next genuine ecological hurdle

**Retain as a promising exploratory structural pattern, not an independently publishable ecological mechanism yet.** The result supports the statement:

> The exact identity of butterfly–host associations channels the existing biogeography of globalized plants into a modest, nonrandom amplification of climatically selective *potential larval-resource geography*.

But it still cannot establish whether butterfly populations use those introduced resources locally, experience colonization, alter abundance or reproduction, or are affected by competition. The HOSTS interaction database and WCVP contemporary introduced distributions do not measure those endpoints.

Before an independent ecology-focused manuscript, prioritize **an independent, non-circular biological endpoint** (e.g. independently sourced local larval host-use events or verified butterfly occupancy, ideally with dates and effort controls). Do not add more static network metrics or post-hoc bin choices to maximize p-values.

The GEB main manuscript is unchanged.
