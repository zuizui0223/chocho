# Spatial-null mixing qualification v0.2 — frozen before extended-chain outcomes

Date: 2026-10-15. **Motivation and disclosure:** The original geographically local climate null run **38699068001** passed basic swapping geometry and reported 199 chain draws per seed at three distances, but the primary 1,500-km null of native resource absorption had lag-1 correlation around **0.78**, estimated effective sample size (ESS) **~9–10**, and early-versus-late quarter drift of **~−34 to −44** native-absorbed incidences. Consequently its nominal Monte Carlo p-values are *not valid independent permutation p-values*. This qualification procedure is chosen **after seeing those mixing diagnostics**, and must not be described as pre-registered inferential planning.

Scientific estimand, fixed data/source, geographic constraints and original-footprint tether remain unchanged: primary outcome **A = introduced-host incidence inside already native-supported butterfly×region resource envelope**, secondary outcome **U = unique newly added butterfly×region cells**, 239 butterflies and 365 WGSRPD3 regions.

## One final, finite recalibration

- Restrict to the previously declared **primary 1,500-km** original-footprint-tethered climate-bin null. Do **not** choose a favorable geographic radius after new results. The 750- and 3,000-km sensitivity results remain explicitly exploratory and unqualified.
- Restart independently from exactly the original observed host-range configuration at **two new seeds: 20261017, 20261018**.
- Increase burn-in from 5×E to **30×E** proposed swaps. Increase between-draw proposals from 1×E to **20×E**. Retain **99 drawn worlds per seed**, meaning total simulation effort increases substantially but the raw number of correlated draws drops.
- Record full draw series for **A and U** for each seed, lag-1 correlation, initial-positive-pair integrated autocorrelation time (IACT), approximate ESS=N/IACT and first-vs-last-quarter mean drift.
- **Reliability gate** for each primary A and secondary U chain:
  1. `lag1 <= 0.50` (and finite),
  2. `ESS >= 30` estimated by nonnegative initial autocorrelation pairs,
  3. absolute difference between first and last quarter means <= `0.5 * full-chain sample SD`.
  4. each chain independently passes original one-step/geographic mobility, burn acceptance and moved-fraction gates.
- Additionally, the **between-chain primary null median discrepancy** must be <= `0.5 * larger within-chain sample SD`.
- If either chain fails any criterion, end with **NOT_EVALUABLE_MCMC_CALIBRATION**. Do not relabel the initial simulation's tiny nominal p-values as statistical significance. Do not launch new parameter sweeps to rescue the effect.
- If qualified, report both chain medians, central 95% simulation intervals, observed-minus-null A and U, and the same inference ceiling. Even qualified, these are *post-hoc MCMC sensitivity estimates* rather than 99 independent exact-permutation samples; the simulation p-values are not confirmatory.

This is primarily a **methodological correction to prevent false confidence** in the already observed possible spatial residual. Neither host range redistribution nor butterfly competition/colonization is experimentally identified.
