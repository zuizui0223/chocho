# Planned BCE/Clarke inventory-sensitivity diagnostic for the resource-gain association
**24 October 2026, exploratory extension defined after observing the fixed HOSTS-only gain–single-host regression.** This is NOT an independent discovery or registered test; an explicit posthoc source robustness analysis.

The original 206 expanded species showed a negative source-conditional association between log resource expansion and proportion of new regions supported by only one recorded introduced host. That fraction cannot be reconstructed solely from BCE result counts H/A/R/U. We must not claim a replicate of the same response.

## Pinned independent alternative host catalogue
Read the exact pre-existing BCE/Clarke host-inventory source artifact from GitHub Actions run `38988757959`, named `butterfly-bce-resource-novelty-sensitivity-v01` (GitHub artifact 12084636746). Require its original source cohort (239), 81 BCE candidate-augmented species, 1027 extra accepted plant–butterfly links, and the invariant totals original `H=75346,A=46387,U=14553,R=14406`; augmented `H=103349,A=67740,U=14967,R=20642`. Reuse the exact 239 frozen species descriptors for baseline log region sizes, host breadth, family breadth and genus. No updated external citation list or observed response used.

## What is mathematically identifiable from both source scenarios
Per species/scenario:
- `E = log((native + U) / native)` (positive native and U required)
- `support_intensity = (H−A)/U = 1 + R/U` (number of candidate introduced host incidences *per* newly opened butterfly-region cell)
- `log_support_intensity = log(support_intensity)`.

This **mean support intensity is NOT the fraction of cells with ≥2 hosts** and cannot demonstrate actual population insurance. Do not replace one metric with the other when reporting.

## Fixed paired comparison
Use exactly the **81 known BCE-augmented butterflies**, and for a paired before/after regression subset require U>0 in BOTH scenarios, reporting exclusions. Estimate species-equal-weight OLS separately for original and augmented source `log_support_intensity ~ 1 + E + log(native_regions) + log1p(resolved_host_species_for_that_scenario) + host_family_count`; adjusted host count adds exact number of new accepted candidate host IDs, if source consistency checks pass. Report coefficient on E and 5,000 genus-cluster percentile CI (fixed xorshift32 seed 20261024) for each, descriptive percent differences, and within-species paired change in E versus paired change in support intensity. Require explicit warning that both are derived from same source catalogue and can be mathematically coupled. An opposite sign or zero-crossing CI falsifies any claim of source robustness.

## Not eligible to become biology on its own
The BCE plant-host catalogue is European-biased, not an independent world sample, and merely changes the definition of known hosts. An increased support count is not independent larval survival or local resource availability. Source-based genus bootstraps do not correct phylogeny, geographic dependence, introduction timing or the shared matrix. Keep frozen `chocho/main` and GEB PR38 unchanged. 
