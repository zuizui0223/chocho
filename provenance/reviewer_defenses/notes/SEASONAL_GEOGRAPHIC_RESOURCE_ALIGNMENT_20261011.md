# Seasonal and subcontinental conditioning leaves a very small resource-alignment residual

Date: 2026-10-11. Status: completed exploratory post-hoc sensitivity, not preregistration, causal inference or a demographic result.

## Decision

Do not promote the residual native-resource alignment into the main novelty claim. It remains positive under the tested constraints, but becomes much smaller when seasonal proxies and subcontinental introduction geography are preserved. This does not change the existing +54.9% descriptive resource-expansion reconstruction or establish that biological effects are absent.

## Same observations, progressively restricted comparisons

Primary panel: 191 butterflies with positive observed resource expansion and a minimum host-taxonomy consistency filter. This filter does not establish complete host records. There are 73,748 introduced host-by-butterfly-by-WGSRPD3 incidences; 46,096 are in regions already covered by at least one native known host (B).

All nulls conserve known butterfly-host identities, native plant distributions, and introduction counts per host and per plant family-region inside the permitted blocks. Mean climate controls use within-continent BIO1 and BIO12 tertiles. Seasonal controls add within-continent median splits of BIO7 (annual temperature range) and BIO15 (precipitation seasonality). Geographic controls add the official WGSRPD level2 subcontinental unit. These are not exact physical-distance matching or plant phenology.

| Null definition | Expected B, mean | Observed excess B | Excess / all 73,748 incidences |
|---|---:|---:|---:|
| Mean climate, matched input guard | 45,720.76 | +375.24 | +0.509 percentage points |
| Plus seasonal climate proxies | 45,935.07 | +160.93 | +0.218 percentage points |
| Plus subcontinental geography | 46,020.62 | +75.38 | +0.102 percentage points |
| Plus both: primary | 46,045.29 | +50.71 | +0.069 percentage points |
| Both, independent seed and doubled burn-in/spacing | 46,045.32 | +50.68 | +0.069 percentage points |

Each model sampled 199 constrained-chain states. Upper-tail Monte Carlo p is 0.005 in each row, the attainable minimum with the +1 correction. The primary combined null's central 95% range is 46,007.90–46,077.05; this is a null-distribution interval, not a causal-effect confidence interval.

Relative to the matched-input mean-climate baseline, the B excess is 86.49% smaller after adding both constraints. This is a change between conditional null models, not an estimate that seasonality/geography causally explains 86.49% of an ecological effect. The numerical reduction is greater with the geographic refinement alone than with the seasonal refinement alone, but that contrast is not a causal variance decomposition.

## New resource geography

Observed new butterfly-by-region resource cells U remain 13,529. The mean-climate null expects 13,823.10 (2.13% shortfall); the combined null expects 13,588.87 (0.441% shortfall). The longer chain gives 13,589.59 (0.446% shortfall). Lower-tail p is 0.005 in both combined runs. The data support a small conditional shortfall, not a large extra geographic constraint.

The primary combined accounting differences are B +50.71, D +9.17 and U -59.87, which sum to zero up to rounding because total introduced incidences T=B+D+U are fixed. D, duplicate introduced-host incidences outside the native envelope, is not unusually high (upper-tail p=0.305). The pooled single-host share is 57.16% versus 57.27% expected (upper-tail p=0.835); the extraordinary single-host-vulnerability narrative remains unsupported.

Direction in overlapping sensitivity panels: the U shortfall is 0.436% for all 239 butterflies and 0.567% for the 80-species high-host-richness subset. These are not independent biological replications.

## Constraints reduce the scope of inference

The combined model has 246 regional blocks, including 175 with one region. Of 74,055 introduced plant-region edges, 33,901 lie in blocks eligible for swaps. It accepted 1,160,067 swaps; the longer chain accepted 2,319,023. At the final checkpoint only 1,768 of 30,373 focal-host edges across the full 239-species panel differ from their original positions (5.82%); the longer chain gives 1,737 (5.72%).

Fixed margins and native-range exclusions passed. The B lag-1 autocorrelation in the primary run is 0.118, with an approximate AR1 effective sample size of 157 of 199. Agreement across seeds supports numerical stability of this constrained statistic. It does not establish complete state-space mixing, full climatic equivalence, or applicability beyond the exchangeable subset. A shrinking residual under tighter constraints must not be equated with a complete biological explanation.

## Conservative input guard

The frozen 5,300-row climate file contains 365 unique regional vectors, versus 369 regions in the fixed lookup. Four missing-climate regions remain individually frozen. The stored BIO12 value for CPI is exactly 65,535, a possible storage ceiling that has not been verified as nodata against original raster metadata. CPI is therefore individually frozen and excluded from cutpoint calibration in EVERY current comparison, without removing or imputing any of its records. Its 30 introduction edges, including 12 focal-host edges, remain present.

This common input guard and independent simulation seed explain why the current mean-climate reference (+375.24) is not numerically identical to the earlier reference (+366.49). Both are about half a percentage point. Never contrast the new constrained outcome with an unmatched-support baseline as a causal effect.

## Verification and reproducibility

Execution: GitHub Actions run 38123136476, exact commit d9387319a0c174616f4c2eb83538f69adc051890, five matrix jobs. All source digests matched. Full repository logs show 95 passed for each job; local red-green testing also passed all 95 tests after the new implementation.

All five result ZIP digests were verified against GitHub metadata. A separate SciPy integer sparse-matrix reconstruction, without the null-engine imports, reproduced observed T/B/D/U and singleton counts in all three panels. An independent NumPy audit checked every stored mean, quantile and tail probability and all 2,985 panel-by-replicate accounting identities.

The machine-readable receipt is ../results/butterfly_seasonal_geographic_null_20261011.json. Full traces and model-specific results are retained in the five Actions artifacts referenced there and in the downloadable verified source/results bundle. No manuscript text or scientific estimate on GEB PR #38 was changed.

## Ecological interpretation boundary

These results concern the spatial configuration of documented larval resources. They do not show larval host use, reproductive success, local abundance, seasonal synchrony, competition, persistence or extinction risk. The next biologically informative evidence is realized breeding or larval-use records in regions that have introduced known hosts but no mapped native known hosts; further subdivision of the same spatial null should not substitute for that evidence.
