# Global butterfly habitat-affinity filtering under host redistribution
**Frozen before model fitting: 19 October 2026.** Exploratory macroecology, not literature reanalysis. Inputs are the existing pinned global `chocho` species panel, with independent broad LepTraits habitat descriptors. No changes to the GEB initial-submission PR #38 or `main`.

## Biological question
**Does anthropogenic redistribution of already-known host plants systematically privilege disturbance-affiliated butterflies over disturbance-avoidant butterflies in *potential larval resource geography*, beyond consumer diet breadth and their original native resource range?**

The naive union of plant occurrences necessarily grows when introduced occurrences are added. By contrast, a **species-trait gradient** in the magnitude of that growth is NOT entailed by the union definition. A positive association would be evidence of *ecological sorting of potential resource opportunity*, not proof of butterfly population, recruitment or abundance change.

## Source population and selection risks
- Frozen `data/frozen/figure_sources/anthropogenic_species_metrics.csv`: 239 butterfly species with native host resources, 365 WGSRPD3 regions, host richness, native and introduced-added cell counts, wing size proxy, voltinism. No new butterfly measurements.
- Independently link `data/frozen/figure_sources/s1_resource_descriptors.csv`: 339 butterfly trait candidates, 239 resource-eligible in this version. Label classes are qualitative natural-history descriptors, not an experimental measurement of disturbance tolerance. Exact species matching only, no imputed trait.
- Taxonomic family from pinned `data/external/leptraits_consensus_v1.0.csv`, as in the earlier family-adjusted analysis. Group genera by first token of accepted study name.
- **Primary exposure** compares only predeclared strong classes: `Disturbance-associated (strong)` versus `Disturbance-avoidant (strong)` (expect n=150 vs n=33). Exclude weak, variable, generalist descriptions from primary; broad strong+weak exposure is sensitivity only.
- Audit 339→239 resource-panel inclusion fraction by disturbance categories. No imputation of missing resource geography, no pretending selection weights recover missing response.

## Primary outcome and contrast
For butterfly i, native region count N, contemporary C and added A=C−N, 365 possible regions:
- `Y = log(C/N)=log_resource_expansion`. Confirm arithmetic against source.
- Estimate the **difference in Y between the two strong disturbance categories**, first unadjusted then controlling: nonlinear native range (`log1p(N)` and its square), `log1p(resolved_host_species)`, `log1p(host_family_count)`, `log(wing_size_proxy)`, voltinism B/M/U, a taxonomic-family fixed effect, and host inventory adequacy indicator. Report actual n, coefficients, finite-sample uncertainty.
- Fit OLS using all primary eligible species (do not select significant species/regions); genus-cluster sandwich CI and a genus-resampling bootstrap 95% interval. Disturbance affinity is an *association*, not a randomly assigned cause.
- Declare primary positive directional biological prediction. **No result promoted solely because nominal p<0.05**.
- Report major-family-specific descriptive sign and commonality; no pretending Nymphalidae/Pieridae family estimates are all independent.
- Also fit a nested model including a 3-level edge-association score (edge-associated/edge-avoidant/other), because edge ecology may covary with disturbance; this is a robustness model, not a new fishing outcome.

## Predeclared sensitivity responses
1. Native-area ceiling: `fraction of previously unsupported regions potentially opened = A/(365−N)`, with continuity-corrected `logit((A+0.5)/(365−N+1))` as alternative outcome; retain native range covariate. These are **resource regions**, not actual range occupancy.
2. HOSTS completeness: repeat strict analysis only for `host_taxonomy_lower_bound_adequate=True`; if direction collapses, demote.
3. Trait ambiguity: broad associated strong+weak vs avoidant strong+weak; keep middle groups excluded from binary analysis and explicitly report n.
4. Extreme ranges: predeclared leave-one-genus-out and leave-largest-genus-out; no cherry-picked species.
5. Selection pressure: all-339 trait category inclusion rates, and the host-species richness distribution by category. Treat source ascertainment as a limit, not as a causal finding.

## Distinguish effect size and predictive utility
- Compare baseline covariate model vs baseline+disturbance **grouped out-of-genus 5-fold cross-validation**, repeated 10 deterministic group assignments; report held-out mean squared error and whether the exposure adds useful predictive information. A fitted significant coefficient *without* out-of-genus gain does not establish practical prediction or universality.
- As an auxiliary null, use 4,999 genus-wise Rademacher wild score sign flips under the reduced model, preserving genera as blocks, to assess the residualized exposure association. This is **model-dependent approximate null**, not exact exchangeability or an introduction randomization test.
- A planned cross-family direction check (not significance optimization) reports which major families contain both strong categories and whether the signs align.

## Stop/interpretation gates
- If the unadjusted disparity is explained by native host ranges, resolved HOSTS richness, taxonomic Family, source adequacy or the coupled ceiling response, **report no independent ecological habitat-filtering evidence**.
- If out-of-genus prediction does not improve, write "description of this panel; limited predictive value" instead of "predicts butterfly globalization".
- If the 339→239 eligibility screen is strongly skewed against disturbance-avoidant species, do not generalize to butterflies worldwide.
- Without region-level native-versus-introduced plant abundance, historical introduction dates, and observed butterfly population change, NEVER claim actual butterfly range expansion, local competition, invasive-plant benefits or causal extinction threat.
- No new literature collection, no separate literature-derived case studies, no manuscript overwrite. Place all outputs on dedicated exploratory analysis branch and optionally draft a PR to the **alternative synthesis branch only** once computations and CI pass.

### Claim permitted if predeclared gates survive
**Potential butterfly larval-resource gains from redistributed plants are ecologically uneven, conditional on disturbance affinity, even after controlling taxonomic host breadth and native-range opportunity in this frozen sample.** Qualify the sampling, family proxy and observational nature; all claims must be conditional on actual estimates.
