# Post-source physical compatibility amendment — preserve original daily diurnal range

**Date 2026-10-28. Status: explicit posthoc model revision made after the first source-only daily-climate QA failed, but before ANY CMIP6 conditional future leaf-probability results were calculated.**

## Failure and cause
At [source workflow #40621864779](https://github.com/zuizui0223/chocho/actions/runs/40621864779), all 72 original Daymet 2018–2025 source station histories were retrieved and rechecked, but the initial CMIP6 scenario construction stopped with `CMIP6 monthly anomaly creates inverted daily tmin/tmax; no silent repair`. The source original WorldClim Tmin/Tmax monthly future-minus-historical changes, applied **separately** to original historical daily Tmin and Tmax, sometimes imply Tmin of a synthetic future *day* > Tmax of that day. This reflects a fundamental **timescale-transfer failure**, not a bad measured Daymet day or the absence of milkweed leaves. The original source data, 18 WorldClim ZIP source hashes, original 72 stations, all model coefficient choices and leaf statuses have NOT been changed.

## Replacement scientifically distinct conditional climate perturbation
- Define monthly `delta_daily_temperature_mean = 0.5 × (delta_CMIP6_tmin + delta_CMIP6_tmax)` at each frozen station's WorldClim 10-minute pixel.
- Add this **same monthly mean anomaly to BOTH** original Daymet Tmin and Tmax on each original 2018–25 day, so original **daily diurnal range remains exactly fixed**, and a physically valid historical daily Tmin≤Tmax pair remains valid after a common warming shift.
- Compute hypothetical 30-day mean temperature and frost-day count from those physically compatible hypothetical daily temperatures; original date, plant species, observed leaf status, Daymet daily precipitation and daily temperature variability are all fixed.
- This method does **NOT** represent the separately predicted future change in diurnal temperature range. It is a common delta-change sensitivity for **historical daily mean temperature and frost**, not a future daily Tmin/Tmax projection. Any future model output must be labeled `POSTHOC_DTR_PRESERVING_MONTHLY_MEAN_TEMPERATURE_DELTA_SCENARIO`, never retroactively presented as the originally frozen separate-Tmin/Tmax daily scenario.
- Validate with source-independent synthetic test that an extreme difference between monthly ΔTmin and ΔTmax still preserves original daily Tmin≤Tmax and daily DTR, and that identical ΔTmin=ΔTmax exactly reproduces the original method. Zero anomalies must reproduce original Daymet 30-day features and model probabilities byte-exact/floating 1e-12.
- Keep all other precommitted gate rules from `USANPN_CMIP6_DYNAMIC_LEAF_DELTA_CONTRACT_20261028.md`: original full source/fold integrity, common 72-station source panel, eight genuine CMIP6 futures, 95% valid climate-source station gate, ≥10 winter stations/≥100 visits, 1st–99th original temperature/frost support reference and **≤10%** accepted out-of-historical-support visit share, otherwise `UNSUPPORTED_EXTRAPOLATION`.

## Scope
This is still **a model scenario about predicted leaf YES at original monitored visits**, not directly observed future leaf phenology, edible biomass, a demonstrated causal warming response, drought response, monarch breeding, survival, future generations or population/range expansion. Keep original GEB frozen manuscript and `main` unchanged.

NASA/ORNL Daymet DOI:10.3334/ORNLDAAC/2361; WorldClim v2.1 original source. Required USA-NPN attribution retained. No individual stations/lat/lon/observer/raw statuses are exported.
