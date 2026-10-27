# Frozen exploratory documentation-bias model for HOSTS/WCVP plant-introduced host-use reports

**Frozen 2026-10-27 before examining any regression coefficient.** This is a posthoc statistical audit of the already established 239-butterfly Location/WCVP source classification; not a prospective confirmatory experiment or a host-switch/establishment proof. Stack on PR #76, independently of GEB main manuscript.

## Outcomes and sampling units

Full original-239 butterfly species (unique binomials) is the numerical cohort. Source WCVP host-report geography originally uses a 95% whole-country introduced/native/unlisted area rule, safe >=95% geographic mapping, plus exact-vetted Hawaii and West Malaysia botanical units. All 257 original HOSTS Location labels retained in the source denominator. DO NOT silently discard source-unmatched, unsafe or mixed cases. An apparent absent record is always 'no positive report in this indexed source', not a butterfly-host biological absence.

Primary response: binary 1 when that species has >=1 triple of the form butterfly×accepted host plant×original HOSTS Location with source WCVP plant introduced-area status; binary 0 otherwise. 62 source-positive species observed, not 62 invasive butterflies. Secondary responses: any plant-native-area documented host report (positive-control reporting) and any WCVP-unlisted-location host report (negative-control/alternative source reporting), and explicitly record their overlapping species support. No response is 'adopted a novel host'.

## Predictor hierarchy (fixed)

Primary predictor: log(1 + number of globally documented WCVP-accepted HOST PLANT families) per original frozen 239 species. Original GEB result near-zero association with potential resource geography remains a separate estimand; do not claim equivalence.

Adjusted model nuisance predictors:
- log(1 + original accepted species×host×Location unique source triple count): source reporting volume / study effort proxy.
- log(1 + contemporary minus native reconstructed resource WGSRPD3 cells from the original study): potential introduced-plant resource opportunity proxy.
- intercept.

All three continuous predictors standardized once over all 239 source butterflies to original cohort mean and sample SD; coefficients report odds ratio per **one original-cohort standard deviation of log(1 + documented host-family breadth)**, not per raw host family. Fit deterministic weak-ridge logistic maximum penalized likelihood with intercept unpenalized and penalty 1e-3 on all other terms; do not fit a complex model if source events or rank fail. Report model convergence, correlation and VIF diagnostics; a sign cannot be interpreted when identification is unstable. Genus-cluster bootstrap 499 resamples, same genus labels, fixed seed 20261027, use percentile 95% interval. No naive per-row inferential p values.

Predictor-swap sensitivity uses log(1 + accepted host PLANT species richness) rather than family breadth with identical nuisance predictors. Species-level analysis with source-native-area report presence (>0) only is a recording adequacy sensitivity, NOT a claim that the butterflies are native or that native report completeness has been proven.

## Reporting-bias falsification

Fit identical adjusted family-breadth predictor with outcome 'native plant-area HOSTS positive source report' and with outcome 'WCVP-unlisted location positive source report'. For all three responses compute genus-cluster bootstrap differences in breadth coefficient introduced minus native, and introduced minus unlisted. If the introduced slope is similar to source-native or source-unlisted slopes (or CIs include zero), do not interpret introduced reporting as specific evidence of a new ecological mechanism; the result could reflect general host-recording intensity and plant geography.

## Biological inference firewall

Predictors are traits and *catalogue size proxies*. WCVP introduced status is for host plant only, not butterfly native/introduced range. The original butterfly-native historical climatic niche, native HOSTS original citation effort and actual larval survival/fitness are not measured. Crop status and plant phylogenetic isolation are NOT justified from plant family counts or ad hoc crop-name screening.

Result status remains SOURCE_REPORTING_ONLY; statistical associations of source recorded status cannot establish host acceptance mechanism, butterfly colonization, human-mediated invasion, niche adaptation or competition. Source data is source-frozen, no new original HOSTS/WCVP raw interaction pairs published, aggregate per species/source only.
