# Species-level coupling test: introduced-host resource gains vs added port-connectivity information

**15 November 2026 | Executed actual original per-species prediction losses | Exploratory observational analysis.**

## Main decision

**New, source-verified result:** Butterflies with a larger global increase in mapped known-host resources do **not** show a clear trend toward greater *additional* GBIF recording-prediction benefit from adding a port-connectivity feature to an already host-informed model.

This is a different and more relevant-to-ports estimand than the **M3 climate+landcover → M5 host-informed gain**. That original per-species M3/M5 contrast is **not** present in the publicly committed source snapshot and is **NOT reconstructed here**. The successfully recovered per-species records are **M5 already including host identities/resources**, **T0** adding port location and an unknown-port-support indicator, and **T1** adding baseline historical port connectivity. The source extension's M5 predictions were explicitly checked against the original 155-species baseline by its original evaluator. Those three per-species errors had been frozen in an independent *transport_forecast_temporal_20261030* results directory before this secondary resource-gap analysis.

The analysis is a retrospective prediction diagnostic, not a test of host introduction chronology, larval performance or actual inter-island establishment.

## Source-verified scope and raw positive/negative contrasts

Exact source crosswalk: **155 originally eligible butterflies**, the same 2015–2018 predictor → 2019–2022 GBIF target scores (**208,482** source species×cell evaluation rows), plus the **95 butterflies** with original strict insular spatial support. We joined each species to the original 239-species fixed HOSTS×WCVP global resource expansion and source-recording support table by exact accepted species names.

| Retrospective source comparison | Species with per-species losses | M5 species-equal MSE | T1 species-equal MSE | Relative MSE improvement T1 vs M5 | Butterfly taxa with T1 better |
|---|---:|---:|---:|---:|---:|
| All originally modeled taxa | **155** | 0.030742795 | 0.030673697 | **+0.2248%** | 96/155 |
| Source-matched strict island subrows | **95** | 0.022127825 | 0.022259568 | **−0.5954%** | 57/95 |
| Strict island subrows excluding Britain and Ireland | **95** | 0.028672672 | 0.028817657 | **−0.5057%** | 47/95 |

For all 155, including port **location** alone (T0) slightly worsened MSE (0.030809744, **−0.2178%** compared with M5), while including connectivity in T1 improved MSE relative to T0 by about **+0.442%**. This original positive T1-vs-T0 contrast and the tiny T1-vs-M5 positive global score are both preserved. The **island T1-vs-M5 score is negative** despite some individual taxa improving; a count of positive species is not the species-equal error-size difference.

**Critically:** only **5** strict-island source landmasses had indexed port-feature support, with **38,623 of 38,648** supported island evaluation rows from Britain/Ireland. The remaining original non-British/Irish island group comprised **25 directly port-supported rows across 3 cells on 3 islands**; the wider 95-species non-British/Irish evaluation mostly had **missing port feature support**, not measured zero connectivity. Consequently the island-only contrast is NOT a well-supported transport null or a negative biological port effect.

## Does known-host expansion predict extra port-connectivity model gain?

Primary preregistered-here descriptive exposure: **original number of globally added butterfly×WGSRPD3 resource regions** due to introduced ranges of already known larval hosts. Primary response: **species-level absolute M5 squared-error − T1 squared-error**, positive for benefit. The null/covariate tests below are newly exploratory and postdate the original transport forecast.

| Species-level source cohort | Spearman rho: added resource units × T1-vs-M5 absolute gain | 4,999-draw taxon bootstrap 95% interval | Source family × native-resource-quartile conditional permutation p |
|---|---:|---:|---:|
| **155 source-matched butterflies** | **−0.0526** | **−0.200 to +0.100** | **0.499** |
| **95 strict-island-observed butterflies** | **+0.0340** | **−0.164 to +0.233** | **0.696** |
| **95 butterflies, non-British/Irish island source subrows** | **−0.0155** | **−0.204 to +0.168** | **0.703** |

Alternative full-155 rank correlations are similarly weak: **log relative known-host resource expansion vs absolute port predictor gain rho=−0.111**; added host-resource units versus **T1 relative species-MSE improvement rho=−0.098**; added host-resource units versus the *connectivity-only* **T0→T1** increment rho=**−0.040**.

A rank-residualized full-155 correlation adjusting original *global native-host resource abundance, number of known hosts, GBIF record-source count, original M5 model squared error* and three host-diet family breadth indicators is **−0.0665**. This is an **exploratory partial correlation**; no phylogeny, climatic region, individual GBIF observation process, or species-specific plant transport trajectory is controlled by this model.

### Source-fixed introduced-resource quartiles

The **155** butterfly species were ordered by global introduced-host-added potential WGSRPD3 region units, then split into four similarly sized groups. The incremental MSE benefit T1 vs M5 was:

| Introduced-host-added resource quartile | Approximate added resource span | Source species | MSE improvement |
|---|---|---:|---:|
| Q1 | 0–25 | 38 | **+0.548%** |
| Q2 | 25–65 | 39 | **+0.403%** |
| Q3 | 66–91 | 39 | **−0.019%** |
| Q4 | 93–259 | 39 | **+0.128%** |

No monotonic dose–response increase is apparent. The quantiles are sorted empirical groups; they are not prespecified biological effect thresholds and contain overlapping cutoff ties.

### Source-scope sensitivity and limits

Restricting the 95 strict-island-included butterflies to species with at least 10, 30 and 100 observations **on the original target-51 strict island subset** changes the resource-added × MSE-improvement Spearman rho to approximately **+0.091** (66 species), **+0.065** (60), and **+0.020** (58). It does not produce a strong positive association. These thresholds refer to the *earlier 51-target island attribution*, not exact 25-island coverage for the transport model, so they are only source-quality proxies.

The original port extension used 2019 port coordinates and modern source-name matching alongside baseline-quarter connectivity aggregates, and a static contemporary host catalogue. Therefore **not** a genuine historical forecast of shipping-causally-mediated butterfly establishment. A positive added resource is only a mapped potential food plant distribution and not an observed larval resource. Port connectivity is not ocean currents, adult wind trajectories or shipborne insects.

The bootstrap treats butterflies as independently resampled units even though host networks and ancestry are correlated; the confidence intervals and stratified permutation p-values are **descriptive, not phylogenetically independent or multiplicity-adjusted**. Species-held-out model transfer and geographically strict ocean-transfer checks have separate limitations.

## Scientific synthesis

1. **Retain all positive results for their actual endpoints:** 239 species mapped global resource opportunity +54.9%; island native-baseline amplification +118.9% versus mainland +40.9%; known-butterfly M3→M5 recording prediction +4.12% later period and +3.00% held space; island Circuitscape geometry-only retrospective MSE gain +5.08%.
2. This new, actual per-species **M5→T1** join does **not** support the stronger *species with most human-redistributed hosts benefit most from harbor connectivity* monotonic mechanism. The distinct positive global T1-vs-M5 +0.225% remains but is tiny and not transferred to island-specific test rows.
3. Without direct host-feeding, source-dated local plant arrivals, directional marine exposure and independently confirmed butterfly breeding on a newly reached island, the mechanistic arrival × resources × settlement chain remains an interesting hypothesis, **not a newly validated result**.

## Exact source and re-execution

The audit reads **eight pre-existing GitHub source files**, verifies exact 155/95 per-species MSE agreement with the original transport forecast published means to 1e-12, checks original modeled-species selection against the 225 candidate ledger, and preserves their SHA-256 hashes. It emits **345 source-qualified species×analysis-scope rows**, 4,999-draw species bootstrap intervals, 4,999 conditional permutations, alternative rank and 7-source-covariate partial correlations, species resource quartiles and island data-support sensitivities.

- [Successful source-locked execution and 155/95/95 original-MSE validation](https://github.com/zuizui0223/chocho/actions/runs/43674229891).
- [Complete immutable source crosswalk, full result JSON and input SHA-256 receipt](https://github.com/zuizui0223/chocho/actions/runs/43674229891/artifacts/14145049953).
- Executable source: `scripts/audit_port_connectivity_gain_vs_host_expansion_20261115.py`.
- Source archived species errors: `provenance/reviewer_defenses/results/transport_forecast_temporal_20261030/transport_temporal_all155_species_errors.csv`, companion 95-species island and island-excluding-Britain-Ireland panels; transport source code and comparator tables also committed.

**Editorial placement:** one clear paragraph in Supplementary Results S.C and a precise statement in Limitations, not a new causal claim in the headline abstract. The original GEB single integrated draft and all earlier positive source measures remain retained.
