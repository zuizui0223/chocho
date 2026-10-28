# Real climate improves held-out-site prediction of directly observed milkweed leaves

**2026-10-28 — Completed USA-NPN × Daymet source-verified field observation test.** This is predictive evidence about live-leaf YES/NO at already-monitored plants, **not** an experimentally estimated weather effect or a future plant, butterfly population, or winter-generation forecast. The original GEB main/manuscript and CMIP6 #90/#92 are untouched.

## Provenance and cohort

- USA-NPN Nature's Notebook 'Leaves' records: original six Asclepias taxa and 2018–2025 source, 27,057 distinct source plant×date statuses; eight original response SHA256 receipts revalidated. A directly monitored 'NO' is real; unvisited sites/months are **not** NO. [NPN source](https://services.usanpn.org/npn_portal/).
- [NASA/ORNL Daymet 1 km daily weather](https://daymet.ornl.gov/web_services) (DOI 10.3334/ORNLDAAC/2361): 8 years × 365 daily records at **all 72 of 72 selected stations**, with **zero source rejections**. Daily Tmax/Tmin/precipitation was used only during the source workflow. No exact sites, observers, individuals or raw daily climate are published.
- Sites were selected **before weather outcomes**, purely from monitoring frequency: >=20 distinct observation dates, >=2 calendar years, >=4 sampled months, valid stable site coordinates in the Daymet domain. From **117** original eligibility-qualified distinct sites, up to 12 per source milkweed species were prioritized, then filled to a fixed maximum 72 by date effort. The original site-eligibility sample had no *A. fascicularis* meeting this source information criterion (not plant absence).
- Exactly **20,770 leaf YES/NO records** at 72 selected stations in **23 states** and **five Asclepias species** survived a complete past 30-day weather gate. YES = 10,742, NO = 10,028. Species-index-record totals: syriaca 11,733; tuberosa 5,351; incarnata 2,502; speciosa 928; curassavica 256. In leap years Daymet does not provide Dec 31; observations with an incomplete prior 30-day weather window were excluded, not interpolated.

## Model comparison and independent-site testing

L2 logistic C=0.1, exact identical sample for both models, feature standardization learned only on each training fold, equal site weights. Baseline: plant species, two seasonal harmonics, latitude/longitude, year. Extended model adds **mean temperature, log1p precipitation, frost days over the preceding 30 observed calendar days**. Entire NPN monitoring stations were held out in five deterministic folds, with no train/test reuse of a station. Metric: station-equal predictive log-loss (lower is better), with Brier as a secondary proper score.

| Held-out fold | Stations | Plant×date observations | Baseline log-loss | +Real weather log-loss | Improvement |
|---|---:|---:|---:|---:|---:|
| 1 | 12 | 5,387 | 0.3668 | 0.3416 | +0.0253 |
| 2 | 18 | 3,800 | 0.3061 | 0.2823 | +0.0239 |
| 3 | 15 | 3,309 | 0.3170 | 0.2838 | +0.0332 |
| 4 | 14 | 3,935 | 0.3994 | 0.3852 | +0.0142 |
| 5 | 13 | 4,339 | 0.3455 | 0.3147 | +0.0308 |

**5/5 sites-held-out folds improve.** Site-equal baseline **0.343748** → full weather **0.318342**, difference **+0.025406**. Brier score improvement **+0.008281**. This satisfies the original >0 total and >=4 positive folds **prediction gate**, without any causal p-value or host-plant physiological interpretation. Twenty-three distinct states are covered, but site selection is **not representative** of North American milkweed population availability.

## What is supported and what remains unresolved

**Supported:** Actual past 30-day temperature, precipitation and frost information improves prediction of directly observed live leaves beyond a calendar/species/geography/year reference on independently held-out USA-NPN monitoring stations. This provides a data-backed reason to *consider* climate-sensitive host timing rather than keeping host season static by assumption.

**Not supported yet:** a temperature-specific effect, climate-driven date of leaf senescence, causal warming response, future CMIP6 plant phenology, monarch reproduction or population growth. Weather components are correlated, observations come from nonrandom monitored plants, and nearby held-out sites may still share observers/climate. Distinguish daily weather-driven status prediction from 20-year monthly CMIP6 climatological scenarios: daily frost counts cannot be substituted with monthly mean Tmin. A distinct posthoc temperature-versus-precipitation component test, complete-year temporal holdout and feature-domain transfer gate are needed before projecting future leaf status.

## Reproducibility

- [Original source-complete action #40618930634](https://github.com/zuizui0223/chocho/actions/runs/40618930634) SUCCESS; 7 source-independent synthetic regression tests pass.
- Exact only-aggregate source JSON SHA256 `52a0f857202f79db874c512549af54c05b001beea65d68e67b9f2d81b2e1c014` at `provenance/archive/monarch/results/USANPN_DAYMET_HELDOUT_LEAF_20261028.json`.
- Script `scripts/audit_usanpn_daymet_leaf_prediction_20261028.py`; original preweather frozen protocol `provenance/archive/monarch/USANPN_DAYMET_LEAF_RESPONSE_PREDICTION_PROTOCOL_20261028.md` and independent tests.
- Mandatory acknowledgment: “Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature's Notebook program.” Cite NASA/ORNL Daymet Single Pixel Extraction DOI 10.3334/ORNLDAAC/2361.
