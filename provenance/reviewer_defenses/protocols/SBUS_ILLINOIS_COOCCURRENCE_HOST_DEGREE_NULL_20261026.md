# Outcome-blind follow-up: observed adult co-occurrence and documented larval-host overlap (26 October 2026)

**Scope:** a NEW descriptive field co-occurrence question after the Illinois lagged coefficient failed robustness/degree-null checks. This source follow-up **does not re-test competition** and does not change main GEB paper, PR #38 or PR #70. The endpoint (species-pair field co-detection) and all thresholds below were fixed before computing co-occurrence contrasts.

## Original fixed source, without source additions
Original Henry et al. (2025), *Scientific Data*, US butterfly monitoring DOI `10.6084/m9.figshare.27934602.v1`, exact 55,066,096-byte CC-BY ZIP checksum MD5 `3161b1c8802a0d67771b2044e18da928`. Original focal butterfly panel: 239 exact species names, blob SHA `9f83d2d9ac8e46f1cfd4ded572c7a023273503dd`. Original accepted HOSTS `808e0b869f9ec1adf8efff87cf6a395adda103e0` and accepted WCVP `65bed76bae9d644ccb6ad200c05f9f5071d89e05`. Exact LepTraits butterfly family source, no fuzzy taxa aliases.

**Geographic/time population** = Illinois Butterfly Monitoring Network, same original source-validated transect method, same accepted survey/site-year inclusion in the preceding pilot: source positive and finite visit duration, ≥5 distinct survey events in a site-year, ≥3 distinct 31-day DOY bins, sites with ≥3 consecutive eligible calendar years. Require reconstruction of **121 eligible sites and 851 eligible site×years** and **47 original panel butterflies positive in that panel**. Fail on drift; other programmes are excluded from this first test because they failed the locked ≥50 eligible sites gate for lagged analysis.

## Primary co-detection estimand
For each original butterfly species in Illinois, create a binary vector over the exact 851 qualified site×years: 1 if the species was positively **recorded** that site-year, else 0 **only as an observational non-detection proxy**. This is not a larval host-plant observation and it is not independently verified biological absence. Restrict to species positive in ≥10 distinct qualified site-years and absent/not recorded in ≥10 site-years; exact threshold fixed before new test.

For each unordered species pair among these retained taxa:
- observed site-year co-detection similarity `J_ij = #both positive / #either positive` (site-year Jaccard), a descriptive response.
- **primary prevalence-standardized binary association**: phi (Pearson correlation of two species' site-year detection vectors). This subtracts the effect of marginal detection frequency algebraically, but does not control species habitat, site effort or spatial autocorrelation. If either variance zero, exclude with count.
- `H_ij = 1` when at least one accepted larval-host *plant ID* is shared; `H_ij = 0` for none in currently recorded incomplete HOSTS catalogue.

Predeclared primary contrast: **equal eligible species-pair weighted** `Δphi = mean(phi | H=1)−mean(phi | H=0)`. Positive means documented host-sharing pairs have more synchronous *geographic adult detection patterns* conditional only on individual detection prevalence; it does **not** prove shared host filtering, local host presence or facilitation. Report sample sizes, baseline shared fraction, mean Jaccard contrast and pair prevalence distribution (aggregate only).

## Strict host-network structural null (computed only AFTER fixed contrast)
199 independently restarted host-graph randomizations, seed `20261026+7919×replicate_index`, 5×E proposed swaps per replicate. Use the audited 239 butterfly × accepted plant bipartite network and exact botanical family / butterfly family assignments: double-edge swaps **only within butterfly-family × plant-family groups**. Preserve each consumer's accepted-host species count within family, each host's consumer count within butterfly family, every butterfly-family × plant-family margin, all 47 Illinois adult site-year vectors and eligible pair identities. No source butterfly/site/year resampling in null.

Require ≥190/199 accepted replicates with ≥100 accepted swaps, ≥10% originally viable edges changed and exact degree checks. The **same species-pair set** must be used in every null state (no eligibility re-selection). For each graph, classify each original pair as sharing or nonsharing and recompute Δphi from the original **fixed** vector correlations. Require both categories contain ≥20 pairs in all qualifying graphs, else mark that draw NOT_EVALUABLE. Archive all 199 per-graph contrasts but no site or event IDs.
- Compare original observed Δphi with null 2.5% / median / 97.5%, as a *posthoc source-conditioned structural diagnostic*, not a causal randomized field test or a species-level confidence interval.
- The reported quantity is observational adult **co-detection**, which may reflect abiotic niche filtering, geographic ranges, species' shared phylogeny, adult mobility and source effort more than larval resources. This null does not preserve climate niche, larval host abundance or within-site larval encounter. The stronger inference would require localized plant abundance and surveyed host plants.

## Secondary: within-site adult rate co-fluctuation (only if identifiable)
For sites where two butterflies were both observed positive at least once across eligible years, compute each species' annual log(1+adult-count/survey hour) residualized by its *site-specific mean*. For each pair require ≥5 jointly eligible sites, ≥30 paired site×years, and strictly positive summed squared residual for both species. Aggregate pair-specific residual Pearson correlations, then the same shared-minus-nonshared contrast and degree-preserving null without changing retained pair set. Report `NOT_EVALUABLE` if fewer than 20 host-sharing pairs meet fixed support, rather than relaxing thresholds or selecting programs posthoc.

## Stop rules
- Do not interpret even a strong positive Δphi as evidence of facilitation, actual population growth or competitive release. Conversely negative correlation does not prove competition or exclusion.
- Do not rerun multiple thresholds/season windows, choose taxonomic families or select only apparent positive species based on outcomes.
- Any new field co-detection pattern is a **separate exploratory source-conditional result**; original GEB paper's host globalization contrast remains conceptual and cannot be assigned a population effect from this undated botanical source.
