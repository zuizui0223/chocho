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

## Final qualified Daymet coarse-seasonal weather sensitivity (E1b)

The pre-specified regional-weather-proxy test **COMPLETED** and passed source identity, all 14 daily gridded-weather response SHA256 comparisons, 121-site/851-site×year source cohort, all original **640** temporal butterfly pair membership checks and all **199** accepted-host family×plant-family degree-preserving null graphs.

The Daymet seasonal predictor set was fixed before weather values were inspected: 21 years of daily modeled data (2000–2020) at each of the **14 public rounded source coordinate classes**; spring GDD5 (DOY 60–181) and summer precipitation (DOY 182–273). Unlike E1a's common-year effect, these are **independent meteorological model inputs**. However they are **not** 1km measured weather at 121 actual monitored sites, and GDD5 is a seasonal environmental proxy, not an observation of host plant phenology.

| Same 640 pairs (62 documented-host-sharing, 578 nonsharing) | 10h source QA original | After coarse regional spring-GDD5 & summer-rain adjustment |
|---|---:|---:|
| Shared group mean within-site annual adult encounter-rate correlation | +0.1994716 | +0.1998833 |
| Nonsharing group mean | +0.0561458 | +0.0522420 |
| Shared minus nonsharing difference | **+0.1433259** | **+0.1476413** |
| Host-degree null median | +0.0886766 | **+0.0930785** |
| Host-degree null central 95% interval | [+0.045865,+0.143163] | **[+0.052267,+0.149143]** |
| Exploratory empirical upper-tail fraction | .030 | **.040** |

The corrected contrast increased only **+0.0043154** instead of disappearing. The median proportion of each partner's pooled within-site temporal adult rate variance statistically captured by both weather predictors was **2.01%** (14.14% of partners had in-sample seasonal weather R² >0.10). The weather-adjusted observed difference remains **inside** the host-degree structural 95% null envelope; empirical `.04` is a posthoc-tail diagnostic, not evidence of climate-immune, causally host-mediated synchrony.

[Machine numeric E1b receipt](../../provenance/reviewer_defenses/results/sbus_illinois_coarse_daymet_weather_host_niche_20261026.json). [Successful source-exact workflow](https://github.com/zuizui0223/chocho/actions/runs/40304049483). An independent [Daymet daily source coverage workflow](https://github.com/zuizui0223/chocho/actions/runs/40303635788) had already archived 14 response hashes and confirmed complete 21-year daily meteorology, with no local site coordinates persisted.

### Scientific stopping decision
This independent climate model's inability to remove the positive adult cofluctuation cannot establish that common forcing is absent. The public coordinates collapse 121 survey sites into only 14 region classes; summer and spring weather are modeled pixels at rounded source coordinates, not measured microhabitat; only two prespecified seasonal dimensions were controlled, not local host biomass or species-specific phenology. The original adult-count zero assumption, shared-sites/pairs pseudoreplication, 3-endpoint multiplicity and literature HOSTS incompleteness all persist.

**Main ecology conclusion:** *Butterflies sharing recorded larval hosts have an exploratory tendency toward more synchronized within-site annual adult count indices in one monitoring programme; this is not conclusively exceptional under degree-constrained host network nulls and not explained away by the limited broad-year or coarse seasonal model-weather controls used here. Local larval-host mediated dynamics and direct interspecific competition remain unidentified.* Under the original tested source-null conditions, field co-detection and next-year competition proxy both failed to exceed expected source structure. A new statistically convincing ecological mechanistic result would require **independently measured, localized, annual larval host-plant phenophases/available leaf biomass and actual larval occurrence or survival**, with more than 14 independently resolved local climate histories.

Do not replace current GEB fixed-source paper's Results or infer butterfly colonization. This analysis remains exploratory in its own source branch. No further source-window cherry-picking or weather threshold retuning.
