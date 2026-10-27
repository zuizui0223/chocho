# Breadth predicts where host-use REPORTS exist, not uniquely an introduced-area butterfly host-adoption mechanism

**2026-10-27. Complete retrospective source-index macro analysis, deliberately descriptive.** Draft stacked PR #77 on PR #76. Main GEB v0.2 manuscript untouched. [Original source-pinned real analysis](https://github.com/zuizui0223/chocho/actions/runs/40440717141) SUCCESS: **11 tests** passed, all 239 accepted butterfly species and exact original WCVP/HOSTS source 11-member hashes verified, fixed 969 disjoint original country area atoms. The full derived species-aggregate and model source-only files are archived in [Actions artifact](https://github.com/zuizui0223/chocho/actions/runs/40440717141/artifacts/12705598677), not original HOSTS interaction matrices.

## Why original question must distinguish source observation from ecological process

Previous chocho GEB main (frozen) found a near-zero relationship between documented HOST-PLANT family count and newly reconstructed *potential resource geography* (Spearman rho ~0.008). Our new source WCVP Location work recovered a different phenomenon: **179 geographically explicit host-use reports in introduced-source-plant areas**, involving **156 unique butterfly–accepted plant pairs / 62 source original butterfly species**. Of those pairs, 145 have **no separately geographically classifiable source-native-area HOSTS record**. This is an indexed-evidence gap, NOT 145 proven new plant-host adoptions.

The test is whether this documented introduced-area signal is more specialized to introduced hosts than to source reports in native or WCVP-unlisted plant areas, after accounting for source recording effort and potential plant redistribution.

## Frozen adjusted model

Species-level binary indexed source outcomes over **all 239 original butterflies**:
- any documented host-use in a WCVP introduced plant area (62 species);
- any documented host-use in a WCVP native plant area (122 species);
- any documented host-use in an original WCVP-unlisted location (107 species).

Main exposure log(1+host-family breadth), standardized across the original 239. Log(1+number of exact accepted HOSTS species×plant×Location reports) controls citation/report volume; log(1+newly added host-resource WGSRPD3 cells) controls opportunity. Deterministic logistic ridge (lambda 0.001), fixed 499 genus bootstrap draws over **134 genera**; sample-identity and design checks pass. Adjusted design condition number **2.11**, moderate source effort vs family-breadth predictor correlation **0.496**. The models predict *whether a report exists*, **not how frequently a butterfly uses that host in nature or whether the butterfly established**.

### Main adjusted association

| Indexed HOSTS Location reporting outcome | Positive species | Odds ratio for 1 SD wider source family breadth | Genus-bootstrap 95% percentile interval (OR) |
|---|---:|---:|---:|
| WCVP introduced plant-area source report | 62 / 239 | **1.49** | **1.04–2.27** |
| WCVP native plant-area source report | 122 / 239 | **1.79** | **1.31–2.63** |
| WCVP-unlisted plant region source report | 107 / 239 | **1.46** | **1.05–2.17** |

A positive source-breadth slope is **not introduced-area-specific**. The introduced minus native breadth log-odds slope difference is **−0.184**, bootstrap 95% **−0.677 to +0.305**; introduced minus unlisted difference **+0.020**, bootstrap 95% **−0.453 to +0.632**. Neither comparative estimate is identified as a distinct introduced-area ecological adoption response; all outcomes come from the same HOSTS source-reporting system.

### Sensitivities and falsification

1. Restricting to the **122 original butterflies with ≥1 classifiable native-plant-area HOSTS report** reduces the introduced-area source report breadth effect to OR **1.34**, log-coefficient 95% bootstrap **−0.251 to +0.845** (OR approx **0.78–2.33**), compatible with no association. Native-recorded presence is fixed by selection here and intentionally NOT refitted as a negative-control dependent outcome.
2. Replacing host-FAMILY count with total source documented accepted HOST-PLANT SPECIES count yields high predictor–source-effort correlation **0.981** and condition number **11.18**, producing unstable negative coefficients (introduced-area OR point 0.20, CI spans null broadly). This is multicollinearity / statistical suppression, not evidence that specialization causes introduced-host reporting.
3. The original 95% whole-country WCVP origin gate is stringent: of 3586 unique accepted butterfly×host×raw-Location triples, **only 881 (24.6%)** are strictly native/introduced/unlisted, while broad biogeographic, unsafe country and mixed categories are still retained in the source denominator. At 80% geographic dominance the introduced triples change **179→197**, and the source-native co-reported triples **12→14**. Country geography uncertainty must remain visible alongside fitted source models.

## Ecological inference (the manuscript-worthy distinction)

These results support a **new, reproducible observation-level fact**: host-use catalogue reports from introduced plant regions are documented for 62 of the original chocho 239 butterflies, and there is large geographical asymmetry in native versus introduced **documentation** of the identical pair. They also reveal that diet breadth cannot currently be promoted as the mechanism that converts potential geographic host opportunities into genuinely adopted new larval hosts: similar breadth effects occur for native and unlisted host reports, and source effort/sampling are entangled with the same underlying HOSTS trait catalogue.

This is a **source-aware falsification** of an overstrong intuitive ecological narrative, not proof that host switching fails. It redirects true invasion niche work toward independently verified invasive *butterfly* native/alien range boundaries, well-dated colonization, local larval feeding or experimental fitness, and effort-controlled source-native reference coverage. In the PR #75 eight-species invasion candidate list, only five belong to the original 239 and geographic invasion-area HOSTS feeding records are very sparse.

**Defensible current conclusion:** anthropogenic plants reshape butterfly *potential host-resource geography*, while verified plant-introduced-region source host-use reports exist but do not yet identify an introduced-host-specific shift in diet breadth or butterfly realized niches. Do not state that 145 pairs are 'new' interactions.

## Reproducibility

- [Recording-bias source-blind protocol](../../provenance/reviewer_defenses/protocols/HOSTS_LOCATION_BREADTH_REPORTING_BIAS_20261027.md).
- [Original 239 country-WCVP interaction reconstruction](../../scripts/analyze_hosts_location_wcvp_origin_20261027.py).
- [Deterministic 499 genus-bootstrap logistic analysis](../../scripts/analyze_hosts_location_breadth_reporting_bias_20261027.py) and [synthetic controls](../../tests/test_hosts_location_breadth_reporting_bias_20261027.py).
- [Original source only CI and numeric results](https://github.com/zuizui0223/chocho/actions/runs/40440717141).
