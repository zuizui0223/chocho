# Does latitude absorb AIS and wind prediction gains?

Exploratory original-51-island macroanalysis, completed 2026-11-07. Stacked on [PR #133](https://github.com/zuizui0223/chocho/pull/133). The existing candidate-host, port, pre-2011 botanical recording effort, Circuitscape and source-target data are unchanged.

**Main result:** Across all original 1,275 island pairs with both endpoints unseen at training time, adding native AIS source resistance retains a small positive contribution after the four frozen latitude controls. Adding mean wind-speed resistance does not retain its earlier positive predictive improvement. The result depends on latitude/geographic specification, not proven transport mechanisms.

## Later butterfly community Jaccard (first indexed 2011–2022)

| Added marine predictor | MSE gain before latitude controls | MSE gain after controls | MAE gain after controls | Between-longitude-block MSE gain after controls |
|---|---:|---:|---:|---:|
| AIS only +2.2226% +1.0020% +0.8485% +1.3347% |
| Wind only +1.3593% -4.6618% -1.6892% -6.9606% |
| Wind after AIS +2.3030% -5.1453% -1.7995% -7.8725% |
| AIS after wind +3.1581% +0.5446% +0.7410% +0.4936% |
| AIS + wind +4.4744% -4.0918% -0.9357% -6.4327% |

Positive percentages indicate lower held-out prediction error. Original-model denominators and latitude-controlled denominators differ; percentage figures should not be subtracted and interpreted as biological effect sizes.

After latitude controls, AIS alone retains **+1.0020% MSE** and **+0.8485% MAE** improvement. AIS beyond wind retains **+0.5446% MSE** and **+0.7410% MAE**. In contrast, wind alone worsens MSE **-4.6618%**, and adding wind conditional on AIS worsens MSE **-5.1453%**.

## Why this is more than a story about latitude

The four extra response-free pair controls were the endpoint sum of signed latitudes, endpoint sum of absolute latitudes, absolute latitude difference and strict opposite-hemisphere indicator. They were taken from the source-audited original WPI port representative coordinates, not inferred from butterfly responses. They enter BOTH reference and augmented models for every latitude-controlled circuit contrast.

The controls themselves affect forecast difficulty unevenly. HRE→HRE_LAT improves all-pair MSE **+2.2786%**, but within-longitude pairs **+14.2008%** and between-longitude pairs **-5.5861%** have different signs. In the full AIS+wind model, adding latitude controls produces an all-pair **-6.4845% MSE change**, with within-block **+10.8110%** and between-block **-17.9750%**. These opposing values suggest geographical transfer and model stability problems, not conclusive proof that latitude causally explains wind.

Indeed, the original unadjusted HCRE MSE over all pairs is 0.007107995, which is smaller than latitude-adjusted HRE_AIS MSE 0.007198528. Thus the latitude-adjusted AIS model is NOT automatically the best standalone predictive model, and the extra controls should not be selected by looking for a desired positive result.

## Earlier butterflies (2001–2010): negative results retained

| Marine addition, with latitude controls | MSE gain | MAE gain |
|---|---:|---:|
| AIS only -8.1571% -3.5753% |
| Wind only -1.7476% -1.0998% |
| Wind after AIS +0.3498% -0.1114% |
| AIS after wind -5.9276% -2.5627% |
| AIS + wind -7.7788% -3.6907% |

These older records were not observed under the AIS 2015–2021 surface. They are not a historical null that must have zero correspondence with modern geography, and differences between response windows are not evidence of changes in transport mechanisms.

## Design, data and verification

The protocol was [frozen on this branch](PROTOCOL.md) before performing the new latitude-controlled fits, but all relevant parent PR133 outcomes were already known. This sensitivity is **post-hoc, observational and exploratory**, not independently preregistered or independent geographical replication.

Identical source population: **51 original islands, 194 source ports, 182 marine-qualified ports and all 1,275 unordered source pairs**. The previous source-only longitude quintiles of 11+10+10+10+10 islands form 15 training/test block pair splits. **Both endpoints** of each test pair are absent from its training fold. Five within-block tests total 235 pairs and ten between-block tests total 1,040 pairs. Train-only z-scoring, OLS intercept and prediction clipping to the [0,1] Jaccard range were held constant.

Four previous models (HRE, HRE_AIS, HRE_WIND and HCRE) and four latitude-controlled variants were fitted for both butterfly target periods. The parent 30 source-factor MSE/MAE scope comparisons were source-locked and independently reproduced before interpreting the new results. All original and latitude-controlled source-factor comparisons, plus two direct latitude additions, are retained in the [complete 72-row numerical summary](contrast_scope_summary.csv).

The independent [auditor](verify_latitude_controls.py) reconstructs **20,400 source-pair predictions, 240 model folds, 360 model-comparison folds and 72 summary rows**, validates all 1,275 original pair identities and 2,550 target observations, and independently confirms no source endpoint leakage. Original source data and all prior positives and adverse results remain unchanged. The GitHub Actions artifact stores the full unrounded predictions and source/fit receipts.

## Scientific interpretation boundaries

**Supported:** AIS-derived spatial connectivity conveys some conditional predictive information even after these coarse latitude controls. The apparent wind-speed connectivity improvement is not robust to this latitude parameterization; it becomes adverse on held-out MSE and MAE. **Not supported:** direct proof of vessels transporting larval hosts, butterfly arrivals or colonization, actual directional wind dispersal, or a unique latitude-mediated mechanism.

The original host predictor still combines 2001–2024 candidate host records, potentially POSTDATING the 2011–2022 butterfly target. The 2001–2010 plant effort covariates are general all-plant recording counts, NOT dated host-species identities. The study does not resolve that temporal information issue. The longitude blocks and test pairs are mutually dependent observations from the same selected 51-island panel; there are no naive independent-fold p-values or confidence intervals. Any conclusion of generalization beyond these islands requires new geographic and historical data.

Previous GEB results such as the +54.9% resource-opportunity reconstruction and original island–port analysis keep their own populations and estimands and are not superseded by this sensitivity.
