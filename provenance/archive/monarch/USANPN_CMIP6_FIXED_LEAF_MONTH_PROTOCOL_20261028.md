# Counterfactual future thermal opportunity at **observed** USA-NPN milkweed leaf monitoring sites

**Date 2026-10-28.** This follow-up is an exploratory source-informed model after inspecting [PR #91] original leaf statuses. It is NOT prospective pre-registration or future leaf phenology prediction. The eligible site-winter panel is held source-fixed **before inspecting CMIP6 outcomes at those NPN sites**.

## Study panel
- Exactly two source-eligible milkweed species according to original prespecified winter criterion: *Asclepias syriaca* **30** and *A. tuberosa* **21** site×winter-years, each with ≥2 distinct leaf observation days in **each** Nov, Dec, Jan and Feb, across 9 source US states per species. Four other original species remain unscorable, never imputed as absent.
- Original USA-NPN source-accepted individual×day statuses 2018–2025: SHA-locked per calendar year as in `provenance/archive/monarch/results/USANPN_MILKWEED_STATUS_2018_2025.json`; original registry `Leaves` phenophase ID 488. Match source station coordinates in memory from the **same** source observation records, reject inconsistent stations exceeding 0.03° coordinate spread. Never export exact site IDs, site coordinates, individual IDs, observers or response bodies.
- The same NPN site may be monitored in multiple winters; the primary result is **equal-weighted by unique site**, averaging per-winter indicators within the site before averaging across sites. Report distinct sites, site×winter-years and climate grid pixels separately.

## CMIP6 climate source and scenario plan
- Use exactly **WorldClim v2.1 10 arcminute original `tmin` and `tmax` archives** and the source SHA values frozen in PR #90's `monarch_cmip6_thermal_window_20261028.json`. Baseline 1970–2000, MIROC6 and IPSL-CM6A-LR, SSP1-2.6/SSP5-8.5, 2041–2060/2081–2100. No invented temperature offset or dynamic range extrapolation.
- Physically invalid source grid pixels are excluded on **common nine-scenario support**, not corrected. Fail closed if >10% of NPN sampled climate pixels are lost (source quality guard for this small panel).
- Biological screening constants fixed from PR #90: lower 11.5 °C, upper 33 °C, winter 323 °C·days over 30+31+31+28 climatological days. This is a coarse thermal screen, **not** verified monarch generation completion.

## Exactly three counterfactuals on every fixed NPN site×winter source
1. `thermal_only`: all Nov–Feb days treated as if leaves were present.
2. `optimistic_month_leaf`: for any month with ≥1 recorded direct **YES** on its ≥2 sampled observation dates, assume the whole climatological month is leaf-available. This is an optimistic *source-recorded monthly switch*, not a measured all-day upper bound.
3. `conservative_month_leaf`: assume whole-month leaves only where **all sampled dates** were YES (zero NOs); else none. This is a conservative source-date switch, not a rigorous physical lower bound.

Require for every sample/month/scenario `conservative_GDD <= optimistic_GDD <= thermal_GDD`; do not fill winter months without direct observations. The source observed leaf switches are **held unchanged under future climates**, explicitly demonstrating the effect of assuming persistent versus nonpersistent leaves, not predicting plant phenology under warming.

## Gates and inference
- Do not combine these source sites with the original 840 iNaturalist photo 100-km hexes: their geographical sampling frames are different and were not co-located.
- Baseline climate 1970–2000 vs observed NPN leaf site winters 2018–2025 has a historical epoch mismatch. This is an **exploratory physical counterfactual**, not validated past winter reproduction.
- Neither leaf masks nor warming effects are observationally assigned, and there are no egg-to-adult survival, sampled host stem absence, dispersal, host plant adaptation or competition data. No projected monarch abundance, population growth, range expansion, actual winter generation number or benefit of cultivated/alien milkweeds is justified.
- Report failure if no scenario can distinguish the counters; zero eligible climate sites means nonidentifiable, not absence of monarchs.
- Keep main, frozen GEB submission and parent PR #90 unmodified.

**Attribution:** “Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature's Notebook program.” USA-NPN raw site/observer data are not redistributed. Also cite original WorldClim/CMIP6 raster sources and authors.
