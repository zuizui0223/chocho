# Source-first independent US butterfly count data: niche filtering vs competition
**25 October 2026 — fixed BEFORE inspecting new count data.** This is an exploratory extension of `chocho` and DOES NOT alter or submit the existing GEB butterfly-resource paper.

## Original independent monitoring source
Henry E.H. et al. (2025), *Twenty years (2000–2020) of butterfly monitoring data across the contiguous United States*, *Scientific Data*, DOI 10.1038/s41597-025-05513-8. Figshare data+code DOI 10.6084/m9.figshare.27934602; original API `https://api.figshare.com/v2/articles/27934602`. According to author release: 19 programs, 64,984 surveys, >12.6 million adults counted, and 1.2–1.4 million species×survey count records. Dataset includes programme, site/event IDs, survey year/date, harmonized taxonomy crosswalk, counts, and effort. Original authors report source-specific effort and absence reporting differences. Our aim is not to reproduce their US decline model but to test **niche-sharing structure in local adult abundance dynamics** by linking documented larval-host relationships to repeated monitored population records.

Data-license, source-file identity and survey restrictions must be verified on original Figshare metadata/manifest before downloading. Do not redistribute raw site data, recorder IDs, anonymous site IDs or precise coordinates in GitHub receipts, even if the source is open; aggregate diagnostics only.

## Primary new test
**Does documented larval-host sharing shift butterfly interactions from positive spatial abundance association (resource filtering) to negative lagged abundance association (potential competition)?**

Predefine original 239 butterfly names using source Git blob SHA 9f83d2d9ac8e46f1cfd4ded572c7a023273503dd (no alternative taxa selected based on responses). The HOSTS exact accepted plant ID network is pinned to globalbioticinteractions/HOSTS commit 808e0b869f9ec1adf8efff87cf6a395adda103e0 and WCVP 65bed76bae9d644ccb6ad200c05f9f5071d89e05. Resolve synonym crosswalk by the author-provided harmonized dictionary; **no fuzzy matching or guesses**. Any focal species missing reliable taxonomy, a resolved host list or effort fields is excluded with explicit reason.

**Prerequisite source gate, before numerical models:**
1. Fetch Figshare article original metadata and exact file list, DOI, modified version, license, per-file checksum and bytes with HTTP status.
2. Identify a verifiable count dataset `complete-share.csv` under `4_res/` or a ZIP containing that file. Inspect code map, taxonomic dictionary, and select a bounded exact-source data member only after metadata reveals available bytes. If access fails/large files require special intervention, stop: do not invent observations.
3. Restrict to **structured repeated transect monitoring** with site and event/date and measured effort; exclude count-circle aggregate and opportunistic field-trip programs from confirmatory count dynamics.
4. Record event/site uniqueness, distinct programmes/species/sites/years, missing zero conventions, effort fields, and exact/mapped overlap with focal 239; requires >=15 focal species, >=100 repeated site-years, >=3 consecutive years and >=100 informative focal site-species transitions before fitting.
5. Dataset counts refer to adults (mainly), not larvae; host-sharing indicates only documented resource niche, not local shared feeding. No WCVP regional plant introduction status is a timestamp for site-level plant arrival.

**If gate permits**, estimate separate (i) shared-host vs matched nonshared species contemporaneous co-fluctuation conditional on site, year, climate/habitat and monitoring effort and (ii) focal `Δlog(1+count/effort)` vs lagged host-sharing competitor abundance, own abundance, same-site non-host-sharing comparator, species/site and year factors. Use lagged competitor measures only at same fixed site and same standardized survey method; a shared species environment can drive correlated and delayed trends (weather, pathogens, management). Species–host network fixed-degree null and placebo time-shuffled competitor test are required before claiming excess niche signal.

**Stop and inference ceiling:** any cross-lag coefficient is observational and does not prove resource competition; ecological mechanism requires local host biomass, stage co-occurrence, egg/larval success or manipulated competitor removals. Separate source feasibility from biological result. Do not repeat the already published 2025 US national trends as a new claim.
