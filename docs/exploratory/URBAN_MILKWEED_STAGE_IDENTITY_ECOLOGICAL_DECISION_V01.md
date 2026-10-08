# Seasonal milkweed host identity and monarch stage composition — decision v0.1

**Date:** 2026-10-08. **Scope:** exploratory analysis branch only; GEB main submission unchanged.

## Research question

Does the **same number of putative larval resources** imply the same ecological response when different milkweed species are present at different times? Specifically, does the relative representation of monarch eggs and later larvae across native versus exotic milkweed reflect botanical status, species composition or seasonal survey mixture?

This is a reanalysis of already published observational data, **not** an experimental estimate of larval survival.

## Provenance

- Original source: Erickson, Schultz & Crone (2025), *Ecosphere* 16: e70259, DOI https://doi.org/10.1002/ecs2.70259.
- Raw public Figshare: https://doi.org/10.6084/m9.figshare.25648644.v1, file \`Milkweed_dat_with_coords_DEC2024.csv\`. Frozen source gate GitHub run **37735892534**.
- Code: \`scripts/analyze_urban_monarch_stage_origin_matched.py\`. Workflow \`urban-monarch-stage-origin-matched.yml\`; successful final **GitHub run 37754705269**, artifact \`butterfly-urban-monarch-stage-origin-matched-v01\`.
- The main matched protocol \`URBAN_MILKWEED_ORIGIN_STAGE_MATCHED_PROTOCOL_V01.json\` was frozen **after** a preliminary matched result was viewed; the species pair addendum \`URBAN_MILKWEED_SPECIES_PAIR_STAGE_PROTOCOL_V01.json\` is likewise post-hoc. **Neither is an independent preregistered discovery.**

## Exact field panel

4,518 milkweed-patch records from 15 urban survey routes during 2022–2024. Complete numerical nonnegative \`Eggs\`, all \`Instar_1 ... Instar_5\`, and \`Num_plants > 0\` leave 4,480 observations. Full 2022 contains **3,001 complete patches**, **131 route×month groups** with both native and exotic milkweeds, and **88 route×month pairs** in which both categories have at least one egg or 4th/5th-instar larval observation.

Analysis statistic is Mantel–Haenszel OR of (instars 4+5) relative to eggs, exotic over native, conditioned on route×month. The events are **cross-sectional reports, not longitudinally followed offspring**. Stage ratios are not larval transition survival, adult recruitment or fitness.

## Temporally induced reversal, not a fixed exotic/native advantage

| Comparison | Relative late-instar/egg OR (exotic vs native) |
|---|---:|
| Full year without adjustment | **1.7624** |
| Route adjustment only | **1.9448** |
| Calendar month adjustment only | **0.5313** |
| Same route and same calendar month | **0.4137** |

The matched annual estimate has **4999-route-bootstrap 95% interval 0.1243–1.0285**, so a nonzero overall native/exotic effect is **not secure**.

2022 4th+5th instar/egg conditional associations by fixed season:

| Season | Active paired route-month groups | MH OR |
|---|---:|---:|
| January–March | 0 | not estimable |
| April–June | 27 | **3.8699** |
| July–September | 32 | **0.1737** |
| October–December | 29 | **0.1332** |

Spring/fall OR ratio **29.05**, route bootstrap 95% **6.16–93.68**, descriptive only. The original 2025 study **already** documented seasonal changes in urban monarch larva:egg counts. This is a new host-origin-stratified reanalysis of published data, **not** new butterfly tracking or an independent mechanism test.

Additional stage-defined MH OR: first+second instar versus eggs **0.395**, all larvae versus eggs **0.467**, fifth instar alone versus eggs **0.369**. Year 2023 has only 15 active paired route-month groups (OR **0.231**); the 2024 carryover has six (OR **6.351**), so neither is a complete comparable 2022 replication.

## The botanical species identity diagnostic is the key result

All nine observed milkweed species occur in exactly **one** native-status category; native/exotic status is perfectly confounded with host plant species. The top three 2022 species by plant counts are exotic *Asclepias curassavica*, native *A. fascicularis* and native *A. speciosa*.

Using exact-route×month contrasts with 4th+5th instars relative to eggs, and 4999 full-route bootstraps:

| Species 1 / species 2 | Active pairs | Matched OR | Route-bootstrap 95% CI |
|---|---:|---:|---:|
| **Native A. fascicularis / native A. speciosa** | 43 | **3.766** | 1.241–28.716 |
| Exotic A. curassavica / native A. fascicularis | 70 | **0.178** | 0.064–0.384 |
| Exotic A. curassavica / native A. speciosa | 55 | **2.279** | 0.371–8.112 |

These contrasts imply that a single native/exotic indicator masks pronounced **host-species differences**. They do NOT prove a physiological cause or compare identically positioned plants. Changing species pairs changes available region/month strata; the effect contrasts are not additive. In spring, *A. speciosa* has many egg counts relative to later larvae; in autumn *A. fascicularis* has disproportionately many later-instar reports. This may reflect asynchronous oviposition/larval cohort timing, field detectability, host phenology or plant architecture.

### Access-field quality trap

In original CSV \`Too_far\`, the string \`NA\` is a literal missing code (428 rows in 2022), separate from blanks and the exact explicit \`No\`. Including blank/NA/No yielded matched OR **0.380**; using only explicit \`No\` produced apparent OR **2.822**, but \`No\` is **absent from the entire autumn block**. Thus the strict-No filter is **not** a comparable full-season effort correction; it simply removes the season where the origin contrast is lowest. This is important nonrandom observer/reporting metadata.

## Independent experimental context, not direct validation

DuBose, Hoogshagen & de Roode (2025), *Journal of Insect Science*, DOI https://doi.org/10.1093/jisesa/ieaf061, compared **actual adult eclosion** on native *A. incarnata* vs exotic *A. curassavica* in summer and fall in Georgia. Their original source is public GitHub \`gabe-dubose/w3mp\` commit \`0c556b075ec3a45bacfb96ead0accaaffa4e4030\`.

The authors' \`analysis/proportions.R\` reports adult-emergence outcomes (summer 29/49 incarnata vs 31/50 curassavica; fall 3/70 vs 11/80). The corresponding aggregated host×season logit interaction LRT has p≈0.128, not strong evidence for a differential seasonal *relative* response, though the paper already reports extension of development into autumn. **Data reconciliation warning:** \`data/eclosions/adults.tsv\` itself has 78 unique adult rows (summer 31 incarnata + 32 curassavica; fall 3 + 12), versus 74 in \`analysis/proportions.R\`. Without a source-confirmed exclusion rule, do not silently replace one with the other or promote an exact new experimental effect. The Atlanta study uses a different native species and geography from the Northern California field data; it cannot validate the field stage ratios numerically.

## Ecological decision and stop

A global botanical host-presence map or an annual egg/larva abundance average misses important **seasonal host-identity turnover**. However these published data **do not** establish that exotic plants increase or reduce cohort survival, nor that a general physiological mechanism was discovered. The most defensible hypothesis is:

> Seasonally complementary host-plant identities can produce different oviposition and advanced larval-stage footprints within the same local resource landscape, so larval-resource opportunity is time- and stage-dependent rather than a static scalar.

The hypothesis is not a discovery until direct within-species host use, repeated individual fate/survival and plant seasonality are observed under comparable conditions. Do not turn this source reanalysis into a new ecology manuscript on the strength of a post-hoc p-value or the 29-fold seasonal descriptive contrast. Any new test must separate plant origin from species identity, account for stage-specific detectability, and follow cohorts or measure adult eclosion/parasitoids.

**Submission main/GEB untouched.**

## Identical three-species route-month support — supplementary diagnostic, 2026-10-08

An additional **post-hoc negative-control sensitivity** fixes the survey-comparison universe to route × calendar-month cells where all three dominant species (*A. curassavica*, *A. fascicularis*, *A. speciosa*) are reported. This avoids comparing different available route-month universes for each botanical pair. Reproducible source-audited GitHub Actions [run 37763626497](https://github.com/zuizui0223/chocho/actions/runs/37763626497), commit `4114edc`.

There are **91 route-month cells** with all three plants present, of which **40** have nonzero egg or late-instar records for every species (spring 14, summer 19, fall 7, winter 0). On those *identical 40 active strata*, the matched fourth+fifth instar to egg odds ratio (first/second species) is:

| Same-route-month pair | MH OR | 4,999 route-resample 95% CI |
| --- | ---: | ---: |
| Native *A. fascicularis* / native *A. speciosa* | 3.389 | 1.030–26.924 |
| Exotic *A. curassavica* / native *A. fascicularis* | 0.300 | 0.093–0.672 |
| Exotic *A. curassavica* / native *A. speciosa* | 2.038 | 0.130–8.530 |

**Selection-on-observation sensitivity:** requiring *all three* to have observed stages can select unusually active patches/routes. Repeating on all 91 *plant-present* cells, allowing zero stage counts and excluding only entirely event-free pairs, gives 79/88/86 stage-informative route-month strata and OR **3.766 / 0.220 / 2.237**, respectively (route-bootstrap 95% CI **1.241–28.716 / 0.075–0.470 / 0.371–7.934**). The direction of all three botanical pair contrasts is unchanged under this additional support definition.

**Ecological interpretation:** an order-of-magnitude late-instar/egg contrast exists **between two native host species** on the same matched calendar and route support, while the two exotic/native pair contrasts point in different directions. Thus **nativeness is not a sufficient botanical explanatory variable** for the observed stage composition. This is a limitation of the categorical indicator, *not* proof that specific host physiology caused survival differences. Host identity remains confounded with location within route, planting phenology, plant architecture and observation/detection; cross-sectional stage counts are not tracked survival. The matched ratios are noncollapsible across groups and should not be multiplied to derive other comparisons. With only seven active fall triad cells and broad route-bootstrap intervals, do not fit a stronger biological mechanism or create a new stand-alone manuscript from this post-hoc reanalysis.

The main GEB manuscript and PR #38 remain scientifically separate and unchanged.
