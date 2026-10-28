# Observed milkweed leaf months constrain conditional monarch winter thermal opportunity under CMIP6

**2026-10-28 — completed source-locked exploratory counterfactual. Not future plant phenology, monarch recruitment, climate-niche adaptation, or an expected change in populations.** Stacked on chocho draft PR #91; the original frozen GEB manuscript and `main` remain unchanged.

## Exact original sources and analysis cohort

- **USA National Phenology Network**, Nature's Notebook, live unfolded `Leaves` status (`phenophase_id=488`), 2018–2025; original source [official service](https://services.usanpn.org/npn_portal/) and [reproducible source code](https://github.com/usa-npn/rnpn). Eight complete source years and exactly 27,057 accepted unique individual×date leaf statuses were re-retrieved and matched to the frozen per-year source response SHA receipts. Missing observations are never treated as leaf absence. Original observer, individual and station identifiers and coordinates were not exported.
- Before CMIP6 outcomes were inspected, the PR #91 **full-winter source gate** retained only *Asclepias syriaca* (30 complete site×winter-years) and *A. tuberosa* (21). A complete winter requires **at least two distinct leaf-status observation dates each in November, December, January and February**. Other four source species lack sufficient complete repeated winter sampling; no status is imputed for those species.
- *A. syriaca* contributed **14 monitored sites, 30 site×winter-years, 12 distinct 10-arcmin climate pixels**. *A. tuberosa* contributed **9 sites, 21 site×winter-years, 9 distinct climate pixels**. These species-level site counts should not automatically be added as a verified union of unique sites. Across the source sample, 19 distinct original climate pixels have physically valid temperatures in all nine panels, with **zero pixel exclusions**.
- Exactly the original 18 WorldClim 2.1 10-arc-minute tmin/tmax archive SHA256s from chocho CMIP6 PR #90 were checked: one historical baseline (1970–2000) plus MIROC6 and IPSL-CM6A-LR × SSP1-2.6 and SSP5-8.5 × 2041–2060 and 2081–2100. No invented warming increments, botanical presence maps, site coordinate uploads, or refitting to future outcomes.
- Egg-to-adult thermal screen fixed **323 °C·days**, developmental threshold 11.5°C, upper cap 33°C. The November–February monthly climate GDD is only an approximate *thermal opportunity*. It is neither daily temperatures nor a demonstrated generation completed in winter.

## Three **fixed observed leaf-month** counterfactuals

At each already monitored source site×winter-year, the same observed leaf records are applied to every climate scenario. The model counts a whole month as thermally available under three explicitly assumed leaf masks:

1. **Thermal-only:** all four months are treated as if edible leaves were present.
2. **Optimistic observed leaf-month mask:** at least one direct live-leaf YES observed on ≥2 observation dates in that month (other visits may be NO). Treat the entire month as available; **not proof of uninterrupted leaves**.
3. **Conservative observed leaf-month mask:** all sampled dates in that month were live-leaf YES. Treat the entire month as available; **not a demonstrated biological minimum**, because unobserved days remain unknown.

Results are first averaged over winter-years **within each station**, then across stations equally. Source years are not 51 independent randomized biological replicates; sites share climate grids and nonrandom observation effort. Hence the shares in a host-month scenario are **site-equal averaged indicators**, not exact counts of sites expected to breed.

## Observed numerical results

Percentage of *source-monitored sites* with winter degree-day opportunity ≥323 °C·days, site-equal (years within site averaged), shown as **thermal-only / optimistic leaf-month / conservative leaf-month**:

| Climate scenario | *A. syriaca* (14 sites) | *A. tuberosa* (9 sites) |
|---|---:|---:|
| Historical 1970–2000 | **0.0 / 0.0 / 0.0%** | **11.1 / 11.1 / 11.1%** |
| MIROC6 2041–2060 SSP1-2.6 | 7.1 / 7.1 / 7.1% | 33.3 / 22.2 / 11.1% |
| MIROC6 2081–2100 SSP1-2.6 | 7.1 / 7.1 / 7.1% | 33.3 / 22.2 / 11.1% |
| MIROC6 2041–2060 SSP5-8.5 | 7.1 / 7.1 / 7.1% | 33.3 / 22.2 / 11.1% |
| **MIROC6 2081–2100 SSP5-8.5** | **7.1 / 7.1 / 7.1%** | **33.3 / 33.3 / 22.2%** |
| IPSL-CM6A-LR 2041–2060 SSP1-2.6 | 7.1 / 7.1 / 7.1% | 22.2 / 11.1 / 11.1% |
| IPSL-CM6A-LR 2081–2100 SSP1-2.6 | 7.1 / 7.1 / 7.1% | 33.3 / 22.2 / 11.1% |
| IPSL-CM6A-LR 2041–2060 SSP5-8.5 | 7.1 / 7.1 / 7.1% | 33.3 / 22.2 / 11.1% |
| **IPSL-CM6A-LR 2081–2100 SSP5-8.5** | **7.1 / 7.1 / 7.1%** | **33.3 / 33.3 / 33.3%** |

Mean thermal-only winter GDD (site-equal) increases in *A. syriaca* from **18.1** historical to **91.9** (MIROC6) or **146.5** (IPSL) in 2081–2100 SSP5-8.5, but the focal 323 threshold is still reached in just **7.1%** of its sample sites. For *A. tuberosa* mean thermal GDD rises from **116.3** to **333.5** or **395.4**, with strong differences between the two site-level winter-leaf masks.

**Ecological inference**: higher temperatures can increase climatic thermal opportunity at some sites, yet *the same* thermal environment yields different apparent host-gated opportunity according to how the observed month leaf records are converted to an availability assumption. This is an empirical-data-informed counterfactual sensitivity, **not** an identified host limitation coefficient: the leaf mask is held constant in the future and the unobserved days are not known. The two geographically different milkweed species are **not a controlled native versus introduced host experiment**. The observed host species were selected for temporal source completeness and need not represent the range of monarch.

## Fundamental remaining dependencies

- NPN leaf records from 2018–2025 are compared to a 1970–2000 WorldClim baseline and 2050/2090 future climatologies. No historical equal-year validation or climate-driven leaf phenology function was fitted. **Future leaf persistence may change with warming, irrigation and drought**, but this analysis holds the observed month masks fixed by design.
- USA-NPN records measure *live unfolded leaves yes/no*, not leaf biomass, leaf quality or host accessibility to eggs/larvae. Do not confuse a month with one YES date for 30 days of leaves; the modeled monthly masks are assumed switches.
- There is no observation of larval egg-to-adult survival, overwintering completion, population growth, host movement or climate range shift. There is also no externally representative sampling-effort or plant stem absence denominator. No probability interval is justified by just two GCMs.
- The NPN site panel was **not co-located with the 840 iNaturalist monarch photo hexes** from PR #90. Do not multiply these 14- and 9-site leaf-mask proportions by that 840-hex climate thermal fraction or promote them to global prevalence.

## Reproducibility and provenance

- [Actual original NPN × CMIP6 source workflow #40615388524](https://github.com/zuizui0223/chocho/actions/runs/40615388524), **SUCCESS**, 4 scientific regression checks, no original site identifiers or coordinates serialized.
- [Aggregate-only source JSON artifact](https://github.com/zuizui0223/chocho/actions/runs/40615388524/artifacts/12784334491); archive ZIP SHA256 `f9f1a6d2c4648bb0eb6137614e8689c5f60c5c6cb61bec9ab56cf033e40ce4b0`; inner JSON SHA256 `69737ce9ff14a5c23daa2b78453eaf1795494f2b7dc436814294b6bb8fd7d0ff`.
- Permanent aggregate-only JSON: `provenance/archive/monarch/results/USANPN_CMIP6_FIXED_LEAF_MONTH_20261028.json`. 18-row summary CSV: `data/exploratory/usanpn_monarch_cmip6_fixed_leaf_site_summary_20261028.csv`.
- Exact original source fingerprint protection and re-run scripts: `scripts/simulate_usanpn_milkweed_cmip6_leaf_window_20261028.py`, `tests/test_usanpn_cmip6_leaf_window_20261028.py`, `tests/test_usanpn_cmip6_frozen_receipt_20261028.py`.
- Mandatory USA-NPN attribution: “Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature's Notebook program.” WorldClim/CMIP6 sources also require original research credits.

**Decision**: The independently observed leaf-month constraint is useful to falsify unconditional winter development opportunity as a sufficient mechanism, but these data cannot identify real winter breeding, future leaf dynamics or a native-vs-alien causal effect. This draft continuation is not merged into the original GEB paper, the parent CMIP6 thermal analysis, or `main`.
