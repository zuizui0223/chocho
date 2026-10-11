# Climate conditioning strongly attenuates resource-geography alignment

Date: 2026-10-11. Status: **completed exploratory, post-hoc sensitivity; not a preregistered or demographic result**.

## Main finding

Introduced known hosts overlap butterflies' already resource-covered regions more than expected under the tested climate-conditioned spatial null. **However, the excess becomes about 71% smaller when coarse temperature and precipitation constraints are added.** The remaining excess is reproducible numerically but small: about **0.50 percentage points of introduced host-by-butterfly-by-region incidences** in the primary panel.

The defensible result is therefore not that climate is irrelevant, nor that a new ecological mechanism has been identified. It is that a large part of the original model-relative excess is sensitive to climatic constraints, while a small residual remains after the two regional climate variables used here.

## Primary panel and quantitative comparison

The primary panel contains 191 butterflies with introduced-resource expansion and at least as many resolved host species as recorded host families. That filter is a minimum taxonomic consistency check; it does not establish complete host records. The panel contains 73,748 introduced host × butterfly × WGSRPD3 incidences. Of those, 46,096 lie in a region already covered by at least one native known host.

All rows below compare the **same observed 46,096 incidences**, and all nulls preserve butterfly-host identities, native plant distributions, host introduction counts, and regional introduced plant-family richness.

| Spatial null | Expected overlap, mean | Observed minus expected | Excess as share of all 73,748 incidences | Upper-tail Monte Carlo p |
|---|---:|---:|---:|---:|
| Plant family × continent, no climate bins | 44,818.50 | +1,277.50 | +1.732 percentage points | 0.005 |
| Additionally temperature × precipitation, 2 × 2 bins | 45,409.47 | +686.53 | +0.931 percentage points | 0.005 |
| Additionally temperature × precipitation, 3 × 3 bins: primary | 45,729.51 | +366.49 | +0.497 percentage points | 0.005 |
| Independent 3 × 3 seed, doubled burn-in and spacing | 45,725.02 | +370.98 | +0.503 percentage points | 0.005 |

Each row uses 199 sampled null matrices; 0.005 is the minimum attainable +1 Monte Carlo p value. These are constrained Markov-chain samples, not 199 independent field observations. For the primary 3 × 3 run, the central 95% null range for the overlap is 45,652.95–45,805.15. This is a **null-distribution interval**, not a confidence interval on a causal effect.

The 71.31% attenuation is calculated as `1 - 366.4874 / 1277.4975`. It is a change between null definitions, **not the percentage of a biological effect causally mediated by climate**.

## Added resource geography also remains slightly below expectation

The 191 butterflies gain 13,529 unique butterfly × region resource cells. The no-climate null expects 14,447.41: a 6.36% shortfall. The primary climate-conditioned null expects 13,823.19: a **2.13% shortfall**, lower-tail p = 0.005. The independent longer-mixing run gives 13,823.95, essentially the same scale of shortfall.

Accounting uses the exact identity `T = B + D + U`, where T is all introduced known-host incidences, B is overlap with the native-host resource envelope, D is extra introduced-host incidences beyond the first outside that envelope, and U is unique new resource geography. Under the primary climate null, observed-minus-null-mean differences are **B +366.49; D −72.30; U −294.19**. They sum to zero because T is held fixed. The data do not show an excess of introduced hosts piled into the same previously resource-free region.

The direction also remains in the full 239-species panel (14,553 new cells versus 14,855.99 expected, 2.04% fewer) and the 80-species high-host-richness subset (6,224 versus 6,371.25, 2.31% fewer). These overlapping panels are sensitivity analyses, not independent biological replications.

## The single-host vulnerability narrative is not supported

The primary panel's pooled single-host share is 57.16%, versus 57.32% expected under the primary climate null (upper-tail p = 0.81). In the 80 species with at least 10 recorded hosts, it is 45.36%, versus 45.21% expected (upper-tail p = 0.365). A high descriptive single-host fraction must not be promoted to unusually poor redundancy, population dependence, or extinction risk.

## Climate coverage and mixing checks

The frozen climate CSV contains 5,300 butterfly-region rows but only 365 unique regional climate vectors; vectors were deduplicated before binning. The WGSRPD lookup has 369 regions. Four missing-climate regions (ANT, BOU, HMD, SSA) were frozen separately in every model. **None of the 74,055 introduced plant-region edges lacks climate data**, so no introduced edge was discarded for missing climate. All available regions had 16 valid points in the frozen extraction.

BIO1 and BIO12 regional medians were placed into within-continent quantile bins. Butterfly occurrence, recording effort, and fitted niche-mismatch scores did not enter the bins. This is coarse climate conditioning, not exact matching of local environments.

The primary 3 × 3 chain accepted 4,085,853 switches; the longer run accepted 8,175,284. Approximately 6,000 of the 30,373 focal-host introduction edges differed from their observed locations at the recorded checkpoints. Fixed margins and native-range exclusions passed. Primary-overlap lag-1 correlations were 0.156 and 0.226 in the two 3 × 3 runs; approximate AR1 effective sample sizes were 145 and 126. These diagnostics and agreement across seeds support numerical stability at this resolution, not a proof of complete state-space mixing.

## Independent verification and correction of a misleading CI status

The exact input and output ZIP digests were checked. A separate SciPy integer sparse-matrix reconstruction, not the geographic-null engine, reproduced T/B/D/U and single-host counts in all three panels. It also checked all **2,388 panel-replicate accounting identities** and independently recomputed stored means, null quantiles, and Monte Carlo tail probabilities.

The original scientific workflow displayed success, but its captured full-suite log actually contained **86 passed, 1 failed, 1 skipped**. The failing test was `tests/test_repository_layout.py::test_active_protocol_and_result_surfaces_are_whitelisted`: three earlier post-hoc notes had been placed on the active protocol surface. The success label was misleading because pytest output was piped to tee without explicit failure propagation.

The three notes were moved to `provenance/reviewer_defenses/notes/`; the whitelist was not weakened. Explicit Bash/pipefail verification and the full optional test environment then produced **89 passed in 4.56s**, run 38115148981 at commit `32891e207081c9c9e61144af24765afcfa4d5d0a`. The climate workflow itself was also changed to fail fast and preserve pipeline exit codes. These repairs did **not** change scientific inputs, the null engine, or the reported numerical results.

## Reproducibility

Numerical execution: run **38114331740**, job **114396102659**, exact commit `d14eedb4239fa269e888be5701984712a59b72d9`. Full results, four raw trace files and reconstructed native/contemporary source CSVs: artifact **11692453601**, `chocho-climate-conditioned-geonull-20261011`.

Artifact SHA256: `2ec1749445d382246324c74a661ecbd7065eb75cb6f1024cf340cac42b765eb5`.

Full result JSON SHA256: `30ea62cfc8bad716b21c214937e12dbfcfe08cacaa0b1da544d979543d2eacb1`.

Permanent compact receipt: `provenance/reviewer_defenses/results/butterfly_resource_climate_geonull_20261011.json`.

Independent checker: `scripts/verify_butterfly_climate_result.py`.

The submission PR #38 remains at `88209b8a215ac1bc9c4c5e67622fa104179c8895`; it has not been merged, amended with these post-hoc findings, or submitted by this workflow.

## Remaining scientific boundary

Climate seasonality, fine spatial proximity, within-region climatic variation, introduction pathways and recording effort are not controlled. No local host abundance, host quality, oviposition, larval success, population persistence or historical introduction sequence was measured. The next discriminating test is finer climatic and geographic matching with adequate swap freedom; a residual that vanishes under those controls should not be redescribed as a novel host-specific ecological mechanism.
