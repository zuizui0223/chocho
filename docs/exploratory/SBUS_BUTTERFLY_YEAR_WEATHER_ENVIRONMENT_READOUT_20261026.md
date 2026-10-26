# Host-sharing butterfly adult count synchrony after broad environmental forcing: source-gated analysis

**2026-10-26 | exploratory continuation only.** The 239-butterfly GEB `main` and PR #38/#70 remain unchanged. An independent 2000–2020 U.S. adult butterfly transect source (Henry et al. 2025, `10.6084/m9.figshare.27934602.v1`) was already proven to contain 1,321,882 species-survey count records. Illinois fixed panel: 121 monitored sites, 851 site×years; **640 qualified within-site temporal butterfly species-pair correlation values**, including 62 literature host-sharing pairs, 578 nonsharing. These are **adult counts per reported survey hour**, not breeding population size, larval success, or experimental competition.

## Pre-environment results and their evidence boundaries
- Source-published 10-hour survey duration cap removed five long survey events, changing effort denominators in only two eligible site×years; original adult annual correlation difference shared-minus-nonshared became **+0.143325885**, 199 source family-stratified host network null median **+0.088676607**, 95% range **[+0.045865,+0.143163]**, nominal right-tail **0.030**.
- A second **independent source-constrained seed** reproduced marginal direction (the original uncapped figure +0.143231, second network null upper +0.143418; right-tail .030). The fixed 121-site jackknife stayed positive but deletion of one site reduced the difference to +0.05749; the 38-butterfly leave-one-out remained positive. Pair time-coverage weighting did not strengthen beyond the host null.
- **Three** spatial/temporal related endpoints had already been examined. The nominal .030 tail on secondary within-site synchrony is not confirmatory after considering multiplicity and shared data source. Even under valid independent randomization a simple Bonferroni cutoff `0.05/3=0.0167` would not accept .030 as a confirmatory family-wise test.

## Spatial validity of a separate weather source
The independent NASA/ORNL Daymet daily gridded weather API was confirmed accessible with complete 2015 data (tmax, tmin, precipitation and solar radiation), DOI `10.3334/ORNLDAAC/2129`. But **all 121 Illinois public monitored sites map to only 14 distinct rounded coordinate classes** (`lat-lon precision: Nearest 1`). Therefore original **1km Daymet pixels cannot be paired with unique actual monitoring site microclimates**. Precise site-specific temperature, precipitation, leaf water content and host leaf phenology **NOT IDENTIFIABLE** from this public source. Original source/QC [receipt](../../provenance/reviewer_defenses/results/sbus_illinois_daymet_site_source_gate_20261026.json).

A separate **coarse regional forcing** dataset was obtained only at the 14 public rounded positions, 2000–2020, one model-pixel per coarse source class and **365 days per year**. Daymet returned a complete 14×21-year source and 100% 851 site×year mapping; raw daily weather and even rounded source coordinates are not copied into GitHub. [14-class model-weather source provenance](../../provenance/reviewer_defenses/results/sbus_daymet_14region_daily_source_gate_20261026.json). No attempt is made to call those model pixels directly measured site climate.

## Independent common-year-only environmental sensitivity COMPLETE
A site-leave-one-out, species-specific yearly anomaly adjustment was frozen before its results and applied to the **identical 640 temporal butterfly pairs**, using the original fixed source adult monitoring data but NO independently measured weather. For every accepted species×calendar year, the annual species count-rate anomaly was estimated across **other sites**, excluding the target site, only with ≥4 other contributing sites. Adjusted original source species-site-year entries **15,978/16,295 (98.05%)**; 317 sparse entries retained unadjusted year anomaly.

| Same fixed cohort | Original source | Species-specific other-site year effect removed |
|---|---:|---:|
| Mean within-site adult rate correlation (shared recorded hosts) | +0.19934 | +0.19030 |
| Mean within-site adult rate correlation (nonshared recorded hosts) | +0.05611 | +0.05005 |
| Shared-minus-nonshared difference | **+0.14323** | **+0.14025** |
| Null host-link degree-preserving median | +0.08860 | **+0.09412** |
| Null 95% central range | +0.04584..+0.14304 | **+0.05864..+0.14910** |
| Posthoc empirical right-tail fraction | .030 | **.065** |

The adjusted association is **inside** the conditional fixed-host-degree graph null. Calendar-year regional common adult responses therefore reduce the perceived exceptionality, without establishing whether meteorology or host-plant phenology caused it. This is **not** a test of environmental causation, only a model-falsification of source-wide annual shocks. [Exact source/run results](../../provenance/reviewer_defenses/results/sbus_illinois_year_shock_niche_synchrony_20261026.json), [replay workflow](https://github.com/zuizui0223/chocho/actions/runs/40303050049).

## Distinct independent coarse-seasonal climate test
A separate protocol defined one regional Daymet model pixel per public rounded position class and specifically chose **spring DOY 60–181 accumulated GDD5** and **summer DOY 182–273 total precipitation** before retrieving its covariate values. Original 640 pairs, same 62/578 host-sharing classification, and 199 accepted-host family-degree fixed null graphs. Source method notes:
- seasonal GDD5 is an **abiotic proxy for developmental thermal forcing**, NOT measured botanical leaf-out;
- a single 1km pixel at a broad approximate coordinate is not representative of its entire field bin. Reanalysis is regional-weather-proxy sensitivity only;
- no available source identifies when a particular host plant appeared, flowered, leafed out, became edible or was fed on in each local transect;
- climate, phylogeny, detection error, species traits and local density remain confounded.

[Predeclared-for-this-test coarse regional source protocol](../../provenance/reviewer_defenses/protocols/SBUS_DAYMET_14BIN_COMMON_CLIMATE_PROXY_20261026.md). The independent final numeric weather-adjusted adult cofluctuation readout lives at `provenance/reviewer_defenses/results/sbus_illinois_coarse_daymet_weather_host_niche_20261026.json` when source and exact-cohort workflow qualify. Do NOT consider this claim complete until that receipt has been checked.

## Biological decision before weather adjustment
No demonstrated **increase of butterfly breeding populations from introduced host plants** nor demonstrated **negative density-dependent interspecific competition** has emerged. The earlier lagged abundance contrast is source-structural and not robust to missing-zero or covariate changes; the spatial adult co-detection difference lies inside 199 family-constrained host graph nulls. The only surviving environmental question is whether marginal real annual adult rate cofluctuation is due to shared forcing versus observed local host resources. With 14 coarse public weather classes and no host phenology, causal discrimination remains limited.

The appropriate ultimate question remains: *When does a newly redistributed host plant become an effective breeding resource versus a shared environmental cue or competitive bottleneck?* Actual test of host-mediated synchronization needs verified plant taxon and age/tissue use, timed plant biomass/leaf phenophases, larvae/egg counts, adult follow-up and replicated sites, rather than reusing a region-level static host-list model.

All source counts and code kept on `analysis/sbus-environment-shock-20261026` and prior exploratory source branch; never merge into GEB submission without independent scientific author assessment.
