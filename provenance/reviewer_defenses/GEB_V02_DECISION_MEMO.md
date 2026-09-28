# GEB v0.2 reviewer-defense decision memo

**Date:** 2026-09-28  
**Submission surface:** `manuscript/butterfly_specialization_ecology_v0.2.md`

This memo records how the v0.1/v0.2 scientific concerns changed the allowed claims. It is provenance, not manuscript text.

## 1. Matched-host null was confounded by plant prominence

### Concern
Uniform sampling from all alternative HOSTS plants in the same plant family over-represents minor, sparsely recorded plants. Actual butterfly hosts could appear unusually anthropogenically expanded simply because frequently used/recorded host plants are geographically broad and commonly introduced.

### Resolution
The strict post-hoc null retained exact butterfly host-species richness and plant-family composition, then added:
- probabilistic matching on host native WGSRPD3 breadth;
- weighting by the number of **other** Lepidoptera species using each candidate plant in HOSTS;
- a combined native-range + other-Lepidoptera-use null.

The strict-null panel contained 207 butterflies.

### Result
- Observed total introduced-added units: **12,853**.
- Native-range-only null median: **7,347**; observed remained above all 999 randomizations.
- Other-Lepidoptera-use-weighted null median: **15,236**.
- Combined null median: **12,971** (95% interval 12,084–13,843), **p = 0.584**.
- Mean log expansion: observed **0.439**, combined-null median **0.409**, **p = 0.073**.
- Median log expansion: observed **0.333**, combined-null median **0.309**, **p = 0.079**.

### Decision
**The original “actual butterfly host identities are unusually globalized” claim is not allowed.**  
The simple same-family host-identity excess is absorbed by plant native range plus network-wide host prominence/recording intensity.

### Positive ecological replacement
Across **8,909 host plants in 278 families**, plant network degree was positively associated with anthropogenic geographic expansion:
- Spearman log degree vs log expansion: **0.313**;
- partial rank association controlling native breadth: **0.267**;
- within plant-family × native-breadth-quintile strata: **r = 0.295**; **0/4,999** permuted correlations were as extreme (Monte Carlo **p < 0.001**).

Allowed interpretation: **network-prominent host plants are disproportionately redistributed**, but HOSTS degree is a joint measure of ecological prominence/commonness and study/recording intensity.

## 2. Occurrence validation was potentially dominated by one species

### Concern
`Pyrgus communis` contributed 22/66 recovered species × region units, while several complete-recovery cases had only 1–4 outside-native units. Unit-level pooling also ignored species heterogeneity.

### Resolution
Species was made the replication unit:
- region-matched null for mean and median species recovery;
- 50,000 species-cluster bootstrap resamples;
- leave-one-species-out analysis for all informative species;
- explicit `Pyrgus communis` and `Pieris brassicae` exclusions;
- unrecovered cases reported.

### Result
- 23 informative species; pooled recovery **66/115 = 57.4%**.
- Mean species recovery: **0.583** vs region-matched null median **0.385**; **0/99,999** null draws as large.
- Species-cluster bootstrap pooled recovery 95% interval: **0.349–0.774**.
- Excluding `Pyrgus communis`: **44/93 = 47.3%**, null median **30** recovered units (24–36); **0/99,999** exceedances.
- Excluding `Pieris brassicae`: **62/111 = 55.9%**, null median **46** (40–52).
- Heterogeneous failures retained: `Erynnis tristis` **0/11** recovered; `Historis acheronta` **5/13** recovered.

### Decision
Occurrence alignment remains a supported **secondary** result, but it is not evidence of local larval use or causal host-facilitated range expansion.

## 3. Occurrence panel was not designed as a validation panel

The 32-species panel was originally stratified for the climate test. The manuscript now states this explicitly and describes occurrence validation as a **secondary reuse** of that panel.

## 4. Phylogenetic non-independence

### Coarse model
Butterfly-Family random-intercept model:
- standardized rank host-breadth coefficient **0.021**;
- 95% CI **-0.111 to 0.154**;
- **p = 0.752**.

### Species-level PGLS
Kawahara et al. (2023) time-calibrated global butterfly tree:
- exact species matches: **124/239**;
- Brownian rank-PGLS beta **0.140**, 95% CI **-0.036 to 0.317**, **p = 0.119**;
- estimated Pagel lambda **0.033**;
- Pagel rank-PGLS beta **-0.051**, 95% CI **-0.234 to 0.131**, **p = 0.581**.

### Decision
There is no supported broad-generalist advantage after species-level phylogenetic correction. Unmatched taxa are not phylogenetically imputed.

## 5. Finite geographic ceiling and regional bias

The host-family-breadth vs proportional-expansion relationship remains near zero after:
- controlling starting native resource breadth;
- progressively excluding species with broad native envelopes;
- butterfly-Family adjustment;
- species-level PGLS;
- dominant native-resource-region stratification.

Europe shows a modest positive regional association, but the global near-zero result is not produced solely by the Northern American component.

## 6. Climate analysis

The n = 24 climate test is now secondary and its figure is Supplementary Figure S1.

Distance-matched filtering scores remain above 0.5 for:
- **21/24** species at 250 km;
- **23/24** at 500 km;
- **22/24** at 1,000 km.

The pre-specified host-breadth release prediction remains unsupported:
- partial rho **-0.166**;
- one-sided **p = 0.2237**;
- bootstrap 95% interval **-0.583 to 0.261**;
- approximate 80%-power detectable magnitude **|rho| ≈ 0.505**.

The manuscript does not claim equivalence or a zero effect.

## 7. Portfolio architecture

The raw family-breadth architecture gradient is retained only as a structural diagnostic:
- raw rho with effective contributors **0.492**;
- exact host-species-count-adjusted correlation **0.059**, permutation **p = 0.502**;
- raw rho with maximum single-host share **-0.486**;
- exact-count-adjusted correlation **-0.071**, **p = 0.415**.

The original specialist-versus-generalist portfolio mechanism is not an allowed headline claim.

## Current paper-level conclusion

The defensible positive result is now:

> Anthropogenic redistribution has substantially expanded butterfly larval-resource geography, and this expansion is concentrated in plants that are prominent across the Lepidoptera–host network. Broad butterfly taxonomic diet breadth does not predict proportional gain. Contemporary introduced-host geography also aligns with butterfly occurrences beyond structural-overlap expectations, although that validation is secondary and non-causal.

## Current main display logic

1. **Figure 1:** resource expansion + strict host-bias nulls + direct plant network-prominence association.
2. **Figure 2:** secondary occurrence recovery + species-cluster/leave-one-out robustness.
3. **Figure 3:** ceiling and regional robustness of the absent generalist advantage.
4. **Supplementary Figure S1:** climate filtering, distance sensitivity and imprecise negative test.
