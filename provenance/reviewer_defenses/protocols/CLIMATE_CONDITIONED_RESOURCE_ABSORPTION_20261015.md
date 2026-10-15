# Climate-conditioned host-range null — exploratory protocol

Date: 2026-10-15. **Post-hoc protocol written after the unconstrained plant-range null and its absorption results were known, but before the climate-conditioned numerical result was computed.** This protocol is not a prospective/confirmatory biological study. No change to frozen GEB submission PR #38.

## Biological issue and estimand

In the previously frozen 239-butterfly panel, the plant-level introduced-range shuffle held 75,346 recorded host × butterfly × region incidences fixed. It found 46,387 incidences falling within butterfly regional resource cells already covered by native host plants, versus about 45,300 under plant-family × continent constrained swaps. Resulting unique new resource opportunity (14,553 butterfly × WGSRPD3) was below the null median (~15,385).

**Question:** Does this excess *potential native-resource absorption* persist when destination climatic context is also approximately held constant? It would be premature to claim a plant-specific or butterfly-specific ecological sorting mechanism without this check.

Primary endpoint: `introduced_incidences_in_native_supported_cells` (absolute count **A**). Preserved total H=75,346 means A/H is a meaningful normalization. Secondary directional consequence: `added_units` (U) under exact H=A+N and U=N-R accounting. Other network overlap, single-host dependency, and co-user outcomes are descriptive and must not be selectively promoted after inspection.

## Climate data and fixed preparation

Use **existing independently produced climate region values**, not species outcome responses:
- TTF GitHub Actions run **36270581743**, successful artifact **10919985332**, `butterfly-climate-release-postgate-independent-result-v01`.
- `unit_climate_crossfit.csv` has CHELSA V2.1 1981–2010 `bio1`, `bio7`, `bio12`, `bio15` median values across deterministic WGSRPD3 within-region sample points (up to 16). Although rows are repeated per butterfly, the climate vector must match *exactly* across rows for each region; otherwise stop.
- A different host-reconstruction path is not allowed. Pinned HOSTS commit `808e0b869f9ec1adf8efff87cf6a395adda103e0`; pinned WCVP v13 commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`; original S1 descriptors run **36222306369**; fixed WGSRPD3/Level1 lookup.
- Four variables are retained for provenance and one strict sensitivity. Missing climates are **not imputed**, and their plant×region incidences stay fixed during randomization. Report missingness and coverage.

## Null designs (defined before opening climate result)

All designs keep recorded butterfly–host links fixed; they keep every plant's introduced-region count within WGSRPD Level1 and each region's number of introduced host plants per plant family. Native plant geography is held fixed and introduced assignments to native cells are prohibited.

Under these constraints, shuffle introduced *plant × WGSRPD3* incidences with degree-preserving double-edge swaps **within exactly matched climate-region categories**:

1. Primary: `bio1` (annual mean temperature) and `bio12` (annual precipitation), each divided at its median **within each Level1**; 2×2 climate bins.
2. Sensitivity: the same two variables in within-Level1 tertiles (3×3 bins).
3. Stricter sensitivity: `bio1`, `bio7` (annual temperature range), `bio12`, `bio15` (precipitation seasonality), each split at the within-Level1 median; 2×2×2×2 bins.

Every member of a swap has the same plant family, WGSRPD Level1 region and climate bin. Quantiles are region-unweighted and defined solely from the earlier CHELSA region summaries. This controls only broad co-occurring climate **categories**, not continuous climate distance, regional topography or native/introduced plant physiological niches. Avoid causal language.

## Quality-gate and execution

- At least **85%** of all focal introduced plant×region edges must have a climate signature.
- At least **40%** of introduced plant×region edges must fall in viable (>=2 plant species and >=2 regions) double-swap groups.
- Burn-in proposes **5× E** swaps (E=total introduced plant×region edges). Must accept >=100 and move >=**5%** of original incidences compared with initial observed graph. Report accepted swaps and moved fraction, including when below cutoff.
- For each design, run two RNG seeds **20261015** and **20261016**, each with **299** draws separated by **2× E** swap proposals.
- Every chain **restarts from original data**. Do not pool serial chain draws as independent permutations; record lag-one draw correlation as a diagnostic and report null median, central 95% randomization interval, nominal one-sided probabilities and differences, per chain.
- If a gate fails, return `NOT_EVALUABLE` and do not treat absent randomization as evidence that the ecological pattern disappeared.

## Decision framework

- If two evaluable climate-matched chains retain a substantial positive native-absorption excess, describe **persistent *conditional spatial overlap*** and acknowledge coarse climate bins; this is neither plant choice nor larval feeding evidence.
- If the excess shrinks toward zero or changes sign, demote the unconstrained-null absorption result as potentially climate-composition-driven.
- If the chains fail coverage or mixing criteria, state **not evaluable** and do not infer support/non-support.

The analysis is an exploratory sensitivity selected in response to the observed result, so a nominal p<0.05 is **not** confirmatory evidence. The amount of environmental correction is limited by the resolution and coverage of the geographic climate artifact.

## Reproducibility
- Code: `scripts/analyze_butterfly_climate_matched_plant_range_null.py`.
- Synthetic invariance tests: `tests/test_butterfly_climate_matched_plant_range_null.py`.
- Isolated workflow: `.github/workflows/butterfly-climate-matched-structural-null.yml`.
- The workflow must emit a full JSON artifact with source SHA, null coverage, movement/acceptance, and all results. Do not modify original v0.1/v0.2 receipts.
