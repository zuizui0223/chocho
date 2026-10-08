# Next ecological test: when do anthropogenic host plants produce butterflies?

**Date:** 2026-10-08. **Status:** independent-source review, NOT a new experimental result. Analysis branch only; main GEB submission unchanged.

## Current verified result and barrier

The frozen global reconstruction shows potential resource geography. For the two Florida species *Phoebis philea* and *P. sennae*, the Florida WGSRPD3 cell **already has respectively 1 and 4 known native Fabaceae hosts** in frozen HOSTS/WCVP, plus respectively 5 and 4 mapped introduced known hosts. Therefore newly sampled cultivated *Senna polyphylla* adds local exposure but **zero regional envelope presence** at this coarse scale. This is a concrete scale mismatch, not an inferred colonization event.

The source-audited Florida pilot found 59 photo-supported cultivated *S. polyphylla* reports in 46 approximate 5-km cells; WCVP has no Florida distribution row for that host. A published rearing experiment records adult emergence from cultivated plants. These independently establish a *known mapping omission*, but **not garden-level larval occurrence in the 59 plant records**.

The held-out 50-butterfly GloBI panel produced 128 high-criterion larval records but **zero** in introduced-only reconstructed butterfly×region cells. This is absence of strict evidence, not evidence of ecological absence.

Koptur et al. (2024) rearing source has 176 unresolved juvenile fates. Worst-case native/exotic adult emergence contrast ranges from -34.6 to +15.6 percentage points, so sign remains unidentified over full collected cohorts. Both origins were cultivated; it was NOT a planted-versus-wild study.

## Three types of independently sourced evidence

1. **Manipulated survival:** Diethelm et al. (2026), Ecological Entomology, DOI [10.1111/een.70094](https://doi.org/10.1111/een.70094), original files on [Dryad](https://doi.org/10.5061/dryad.51c59zwgx). Their 2018/2019 two-plant common-garden predator exclusion experiment followed larvae to adult emergence (493 larvae across the two years), with plants that are both *native*. This is the best publicly identified causal test of food-plant identity × predator exposure, but it does NOT directly test alien-host origin. The paper already reports host-dependent non-lethal predator effects and butterfly morphology.

2. **Regional seasonal ecology:** Erickson, Schultz & Crone (2025), Ecosphere, DOI [10.1002/ecs2.70259](https://doi.org/10.1002/ecs2.70259), source dataset/code [Figshare v1](https://doi.org/10.6084/m9.figshare.25648644.v1). Monthly East Bay surveys covered 15 urban routes in 2022; four routes were followed in 2023. Native and nonnative milkweeds show distinct seasonal profiles. The published result finds separate origin-composition model comparisons for egg counts (p=.010), larvae (p=.122), adults (p=.006), OE infection (p=.251). Eggs peak in summer; larvae peak in winter. **The p-values are not a test of differences between stages.** The article attributes seasonal effects to several possible pathways, including predation, local density, phenology and OE. It does NOT measure stage-specific mortality directly.

3. **Alien-host exposure vs local conditions:** Ferreira & Rodrigues (2021), Ecology & Evolution DOI [10.1002/ece3.7821](https://doi.org/10.1002/ece3.7821), dataset [Dryad](https://doi.org/10.5061/dryad.3bk3j9kk2). Two Danaini species in lab native-vs-exotic feeding trials, plus in-field bagging/exposure on hosts. While valuable because it follows real larval responses, the field study placed each plant species in **one different location**, perfectly confounding host and site. The protection treatment also partly changes abiotic exposure. It cannot isolate a causal natural-enemy×host-origin term.

Detailed machine-readable provenance and limits: `docs/exploratory/BUTTERFLY_RESOURCE_REALIZATION_SOURCE_MANIFEST_V01.json`.

## One focused falsifiable question

**Does the fraction of nonnative year-round larval host plants modify how strongly egg abundance corresponds to larval abundance, conditional on resource density and season?**

Do not call this larval survival: repeated, independently detected egg and larval stage counts are not cohorts, and adults move.

Planned *before fitting any new raw data*:
- source: Figshare v1 with original route×month egg/larva, native/exotic milkweed and search effort counts;
- source gate: actual original data accessible, all included routes/years audited, stage-specific denominator verified;
- primary model comparison: stage (egg vs larva) × nonnative host fraction interaction, with total visible host abundance and surveyed length, season/day-of-year and repeated route in the model; model family chosen from dispersion checks before looking at interaction;
- robust alternate: analyze only 2022's 15 routes (main), then four repeatedly surveyed routes in 2023 (temporal check), with native/non-native plant phenology described openly;
- require a single interpretable effect-size contrast, uncertainty interval and held-out route/season robustness; simple sign differences in individual stage models are not enough;
- **stop rule:** If the dataset lacks comparable plant-accessibility or stage-specific effort data, or if stage×host-type effect is non-estimable due to season collinearity, do not call the model an ecological mechanism.

This can become a useful empirical validation *of the scale/phenology gap*, but **not a universal statement about all butterflies from one monarch system**. A broad inference needs additional butterfly lineages, paired sites and measured adult emergence.

## Runs and open blockers — updated 2026-10-08

- **Original source gates resolved, not evidence of local fitness:** 59 cultivated *Senna* records and the *Phoebis* direct-larva/observer audit completed in [run 37728156494](https://github.com/zuizui0223/chocho/actions/runs/37728156494); a separate 915-larva independent field probe did not establish an explicit *Phoebis × S. polyphylla* feeding link. The CRG rearing endpoint feasibility audit completed in [run 37730024907](https://github.com/zuizui0223/chocho/actions/runs/37730024907). Feasibility/coverage is **not** an independent observed fitness effect.
- **Figshare original route-level stage counts are now accessed:** source gate [37735892534](https://github.com/zuizui0223/chocho/actions/runs/37735892534); source-audited 2022 15-route/3,001-patch analysis and host-identity triad/zero-inclusive controls succeeded in [run 37763626497](https://github.com/zuizui0223/chocho/actions/runs/37763626497). The 91 three-host route-month presence cells and 40 fully active strata show strong native-vs-native stage composition differences, thus botanical origin is **not a sufficient explanation**. Full numerical results and uncertainty are in `URBAN_MILKWEED_STAGE_IDENTITY_ECOLOGICAL_DECISION_V01.md`. These are post-hoc independently detected stage counts, **not** tracked offspring or survival.
- **True demographic experiment still lacks accessible individual data:** original Diethelm et al. Dryad DOI `10.5061/dryad.51c59zwgx` publicly lists source CSVs and adult-eclosion/survival variables, but [run 37745988292](https://github.com/zuizui0223/chocho/actions/runs/37745988292) failed to download originals: Dryad v2 file and DOI-export endpoints HTTP 401; public file_stream HTTP 403. No experiment outcome was fitted or invented. Its published host-dependent growth/morphology findings are **prior art**, not our discovery.

**Updated decision:** Global host redistribution changes reconstructed butterfly resource supply, while its ecological value is not inferable from botanical presence, binary plant origin or cross-sectional egg/larva ratios alone. Continue only with **directly accessible, independently informative demographic contrasts**; no extra static geographic cutoffs or recycled post-hoc stage contrasts to manufacture novelty. Main GEB manuscript and PR #38 remain scientifically separate.
