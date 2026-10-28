# Follow-up falsification: is weather-based leaf predictability temperature-driven and transferable across years?

**Date: 2026-10-28.** This is a **posthoc secondary prediction validation**, designed AFTER PR #93 original primary held-out-site weather result, but BEFORE inspecting any temperature-only, precipitation-only or year-holdout model outcomes. Keep frozen original PR #93 result unchanged. Model is a test of observed *Leaves* YES/NO on monitored Asclepias plants, **not an identified physiological or causal warming response**.

## Exact frozen empirical panel
- Reuse original 72 NPN sites selected on monitoring effort, same original six named species, 2018–2025, 27,057 original unique plant×date statuses. Guard original NPN each yearly source SHA-256 and Daymet original 8-year daily source combined SHA256 against [PR #93 frozen result](../results/USANPN_DAYMET_HELDOUT_LEAF_20261028.json). This prevents posthoc cherry-picking, date changes, data increments, taxon substitutions, and site relabeling.
- Require the original complete analysis matched plant statuses **20,770** at 72 stations, 23 states and 5 source taxa. Missing Daymet weather, NPN status or monitored site is unscorable; never assign biological zero.

## Pre-fixed component tests
Maintain exactly the original PR #93 two-harmonic seasonal/site geography/species/year reference, train-fold-only standardization, site equal weights, logistic L2 C=0.1, fixed 5 SHA-hashed held-out-site folds.

Compare baseline against:
1. **Temperature-only** increment = 30-day mean air temperature **plus 30-day frost-days**, no precipitation covariate.
2. **Precipitation-only** increment = log1p total 30-day daily precipitation, no temperature covariates.
3. **Original full** increment = temperature + precipitation + frost, as a mandatory original-results reproducibility check.

Primary component gate: temperature-only must beat baseline on **site-equal log-loss for at least 4 of 5 held-out site folds AND for the 72-site aggregate**. Precipitation-only is a competing explanation, not a posthoc chosen winning model. Never interpret logistic coefficient sizes as physical photosynthetic sensitivity, leaf mortality or a causal gradient.

## Independent time-window control
- **Train** on complete original observations 2018–2023 only; **test** on 2024–2025 observations, at the original fixed sites. All status/model scaler/regularization is trained only on earlier years. Report site-equal log-loss and Brier for baseline, temperature-only and full. Require ≥30 sites with >=3 source species, both explicit YES/NO test statuses, and >=1000 withheld observations; otherwise fail source information gate.
- Temperature-only must also beat baseline in this time holdout to authorize a *separate, conditional* future temperature perturbation simulation. If temporal skill fails even though site folds pass, do **not** claim future leaf phenology predictable.
- This is **not** leave-future-years-out site-independence: the time holdout allows same sites in early/late years, which the five spatial folds separately address. These are complementary source-diagnostic gates, not causal identification.

## Future CMIP6 transfer firewall
Positive predictive skill is not automatically future climate transfer. A separate future experiment must source-lock the **same WorldClim v2.1 baseline/future monthly tmin/tmax archive SHA256s as PR #90**, calculate **monthly future minus baseline** temperature deltas at the frozen source sites, and apply only those deltas to historical daily Daymet Tmin/Tmax in a hypothetical future replay. Seasonal precipitation remains observed/historical **FIXED**, never silently replaced with future precipitation or called a drought forecast. Compute new frost-day counts from shifted daily Tmin, never infer daily frost counts from 20-year monthly Tmin directly.

Temporal response is learned from 2018–2025 observation status at nonrandom monitored sites. Model extrapolation outside historical weather support must be explicitly audited; future data has no matching host site/irrigation changes, nor observed leaf responses, nor egg-to-adult survival. Conditional future leaf probabilities may be presented as **model scenarios only**, and only if both component and held-out-year skill pass. **No monarch abundance/range-generation estimates or causal claims.**

## Data and reporting controls
All original NPN station IDs/coords/plant IDs, raw observation dates, Daymet daily rows and predicted site-specific labels stay temporary in the GitHub Action process. Public git and artifacts contain only source counts, five source-independent folds, time holdout performance, source SHA provenance and scientific interpretation. Retain USA-NPN mandatory participant attribution and Daymet DOI 10.3334/ORNLDAAC/2361. Original GEB manuscript/main and CMIP6 PR #90/#92 unchanged.
