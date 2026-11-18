# 239-butterfly paired island/mainland host-gap filling: source-verified trait heterogeneity

**15 November 2026 | Executed original data re-analysis | Source target: potential resource geography, NOT butterfly colonization**

## Decision

The large **+118.9% island** relative gain compared with **+40.9% mainland** is real under the original frozen native-versus-contemporary HOSTS×WCVP resource definition, but **is not caused by more efficient native-host-gap filling on islands**. Island butterflies started with fewer native-host-covered botanical regions, hence a large relative increment can coexist with weaker conversion of uncovered regions into newly usable regional resources.

This is a property of **butterfly × WGSRPD3 region potential host coverage**, not observed larval fitness, abundance, sustained populations, actual feeding, ship transport, or range expansion.

## Original source panel, reproducibility

- **239 named butterfly species**, each evaluated against the **same fixed panel** of 126 strictly insular and 237 mainland WGSRPD3 Level-3 botanical regions.
- Original source-locked [by_species.csv](../island_resource_gap_filling/by_species.csv), [pooled.csv](../island_resource_gap_filling/pooled.csv), frozen [butterfly trait table](../../data/frozen/figure_sources/anthropogenic_species_metrics.csv), and [original gap-filling summary](../island_resource_gap_filling/summary.json).
- Original resource distributions and known HOSTS identities were kept **unchanged**. No new ecological records, origins, or reconstructed butterfly presences were generated.
- **GitHub Actions source test passed**, with exactly 239 paired rows, sums reproducing native 4,829 and introduced-added 5,741 insular combinations versus native 21,215 and added 8,686 mainland combinations, and the previous 65/141/33 split.
- [Successful source-locked workflow](https://github.com/zuizui0223/chocho/actions/runs/43668672593), [full 239-row paired CSV, summary statistics and SHA-256 artifact](https://github.com/zuizui0223/chocho/actions/runs/43668672593/artifacts/14142524691).

## What changed: finer butterfly ecology, not rewriting positive findings

| Paired estimand (same butterflies) | Original eligible taxa | Result |
|---|---:|---|
| Island native-host coverage within all 126 insular ×239 possibilities | 239 | 4,829/30,114 (**16.0%**) |
| Mainland native-host coverage within all 237 mainland ×239 possibilities | 239 | 21,215/56,643 (**37.5%**) |
| Added hosts filling initially uncovered island opportunities (pooled across butterfly × region) | 239 | 5,741/25,285 (**22.7%**) |
| Added hosts filling initially uncovered mainland opportunities (pooled across butterfly × region) | 239 | 8,686/35,428 (**24.5%**) |
| Island versus mainland **species-equal** fill fraction | 239 | **−5.440 percentage points**, descriptive species bootstrap interval **−7.168 to −3.753** |
| Species with greater island / mainland / identical gap-fill fraction | 239 | **65 / 141 / 33** |
| Species with taxonomically adequate HOSTS coverage, same island-minus-mainland mean | **215** | **−5.589 pp**, descriptive bootstrap **−7.411 to −3.703** |
| Original 15 matched area × broad biogeography strata, NOT species-equal | 239 | island 24.12% vs mainland 22.85% (**+1.265 pp**, descriptive permutation p=0.5622) |

The different pooled, equal-species and area-stratified comparisons use **different weights** and can have different signs without contradiction. They do **not** imply a higher rate of host-gap repair on islands. Non-indexed or regional plant absence is not a proven field-level ecological absence.

### Does the direction depend on dietary specialization?

| Original frozen HOSTS host-family breadth | Butterflies | Island greater | Mainland greater | Tied | Species-equal island-minus-mainland fill difference |
|---|---:|---:|---:|---:|---:|
| 1 host family | 97 | 29 | 53 | 15 | **−4.36 pp** |
| 2 host families | 59 | 8 | 39 | 12 | **−6.15 pp** |
| 3–5 host families | 50 | 16 | 28 | 6 | **−7.83 pp** |
| ≥6 host families | 33 | 12 | 21 | 0 | **−3.74 pp** |

The direction is **negative in all four original breadth categories**. Species with positive island-specific opportunity nevertheless exist, and that variation should not be treated as uniform ecological success.

Illustrative **island-favored potential resource coverage**, NOT butterfly colonization:
- *Junonia lemonias*: island-minus-mainland gap-filling fraction **+41.6 pp**;
- *Elymnias hypermnestra*: **+36.7 pp**;
- *Calpodes ethlius*: **+33.5 pp**.

The original global result of **206 of 239 butterflies** with at least one added regional opportunity and **+54.9% globally** survives; likewise the original no-clear-host-family-breadth prediction for the *overall relative global expansion* (rho approximately 0.008) concerns a different estimand and is not contradicted by breadth-specific gap filling.

### Counterfactual accounting under equal species-specific fill rates

Take each species' **observed mainland gap-filling fraction** and apply it to that same species' **observed number of native-host-uncovered island regions**, holding native island baseline fixed:

- Actual island introduced-host-added butterfly×region combinations: **5,741**.
- Hypothetical additions at the same species' mainland fill fraction: **7,135.26** (a fractional expectation, **not actually observed regions**).
- The hypothetical − actual difference is **+1,394.26** combinations, illustrating that mainland-rate gap filling would have filled **more** island opportunities than the actually reconstructed host plant geography did.
- In that accounting scenario, new island opportunities/native island opportunities would be **+147.76%**, compared with the observed **+118.89%**. The baseline multiplier persists even when island gap filling is not enhanced.

This exercise is algebraic standardization, **not an experiment or model of plant dispersal, establishment or butterfly population rescue**. It cannot support statements about which plant or which harbor imported the food source.

## Statistical and ecological boundary

- The 239 butterflies are related phylogenetically and share hosts; **species bootstrap intervals and an unclustered binomial sign diagnostic are exploratory descriptions**, not valid phylogenetically independent confirmatory tests. The audit output marks this caveat explicitly.
- The island and mainland regions differ in size, climate, botanical record quality and native floristic history. The original 15-block area×Level-1 comparison gives no robust island advantage for the gap-filling *fraction*.
- Habitat suitability, introduced plant density, food quality, seasonal overlap, larval development, island breeding, shipborne vectors and dispersal trajectories are **not** measured by these tables.
- The original WCVP origin classes apply to geographic regions and a known-host list; a butterfly species' newly gained potential region is not a verified new butterfly locality.
- The *ports* positive correlation operates through **whole introduced flora**, not a verified host-specific harbor enrichment after floristic controls (original 98-region beta=−0.0415, descriptive p=0.093). The *marine* circuits are retrospective community predictors under a separate 51-island cohort, not the source of the 239×region opportunity gains.

**GEB interpretation:** Human movement of host plant ranges reshapes butterfly **resource opportunity** most visibly where native host coverage is sparse. The novelty is a **resource-base asymmetry** across island and continental botanical regions, with butterfly-species heterogeneity, rather than proof that ports specifically increase host plants or that islands preferentially attract suitable food. This distinction keeps the butterfly's dietary opportunities as the object of study.
