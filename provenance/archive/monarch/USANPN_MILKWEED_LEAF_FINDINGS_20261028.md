# Six milkweeds: independent observed live leaves, and the limits of CMIP6 projection

**2026-10-28 — REAL USA-NPN source data; original GEB main and previous PR #90 untouched.** This is not a demonstration that warming leads to sustained monarch winter breeding. The observed plant phenophase is **Leaves** (one or more live unfolded leaves, not biomass or host suitability); plant cultivation category is not recorded.

## 1. Source identity and data completeness

- USA National Phenology Network *Nature's Notebook* public observations: [observational data](https://data.usanpn.org/observations/), [API documentation](https://www.usanpn.org/data/code), [milkweed Leaves definition](https://naturesnotebook.usanpn.org/npnapps/species/Asclepias/syriaca).
- Exact six registered *Asclepias* taxa, registry source IDs `199,197,195,224,1329,1326`, eight complete years **2018–2025**, exact `Leaves` phenophase ID **488** in every year. The source registry was frozen [before any observation outcomes](USANPN_LIVE_LEAF_OBSERVATION_PROTOCOL_20261028.md).
- Full original source responses 8/8 **HTTP 200**, **30,476 Leaves-status records**. After site+individual+date exact de-duplication, omission of non-geolocatable state, conflicts and status unknown, **27,057 accepted individual×day observations**; **5,620 observed site×month instances**. Unvisited months are UNKNOWN, not leaf-absent months.
- [Source-verified acquisition](https://github.com/zuizui0223/chocho/actions/runs/40614049701) passed ten independent/offline tests; [independent full-four-month monitoring sensitivity](https://github.com/zuizui0223/chocho/actions/runs/40614463921) passed three additional tests.

## 2. Source observation support, not winter host availability at randomly sampled sites

| Original species | Winter site×years with at least 2 dates (not necessarily 4 months) | Explicit winter YES / NO observation-days | Site×winter-years monitored **all four** Nov–Feb months on ≥2 dates **each month** | At least one observed YES in **each** winter month | At least one winter month with **all** observed visits NO |
| --- | ---: | ---: | ---: | ---: | ---: |
| *A. syriaca* | 148 | 56 / 1,749 | **30** (9 states) | **1/30** | 29/30 |
| *A. tuberosa* | 73 | 184 / 852 | **21** (9 states) | **5/21** | 16/21 |
| *A. incarnata* | 44 | 24 / 326 | 7 (3 states) | 1/7 | 6/7 |
| *A. speciosa* | 16 | 5 / 146 | 6 (2 states) | 0/6 | 6/6 |
| *A. fascicularis* | 2 | 14 / 4 | 0 | not estimated | not estimated |
| *A. curassavica* | 19 | 98 / 43 | **1** (1 state) | 0/1 | 1/1 |

The counts are **site×winter-years**, not distinct independent plant populations. YES in each monitored month is only an optimistic *sampled-date* screen; nothing establishes leaf continuity between days. Repeated NO statuses document real winter leaf absence in the monitored plants. Monitoring locations differ substantially among species and are citizen-science selected.

**Source information threshold (≥10 fully sampled winter site×years across ≥2 states)** passes *A. syriaca* and *A. tuberosa* only. All other species cannot receive a winter-availability parameter for future prediction from this source. The more permissive gate requiring ≥15 site×winter with any two dates across ≥3 states passed five taxa; it was insufficient for four-month leaf continuity.

## 3. Geographic independence check: FAIL

The original *A. curassavica* vs *A. syriaca* same state×winter-year×month union has **12** candidate shared source strata, but **0** contain at least two independently repeated survey site-months in each species. Therefore the published-looking large contrast in raw winter YES frequency is **not** a geography-controlled species-native/invasive/evergreen effect. Do NOT compare their raw winter percentages as a causal coefficient or conclude cultivation causes more winter leaves.

## 4. What this can and cannot feed into the CMIP6 experiment

The existing [PR #90 real CMIP6 thermal scenario](https://github.com/zuizui0223/chocho/pull/90) is a climate-only source localities test: historical 74/840 winter GDD-threshold source hexes vs 2090 SSP5-8.5 **23.45–25.60%**. Its `all_winter` vs `half` leaf persistence classes were **explicitly hypothetical**.

These new NPN observations establish **actual live-leaf YES/NO phenophase sources** and identify a weak but scorable *A. syriaca*/*A. tuberosa* four-month repeated-monitoring panel. They **do not** map to the exact original 840 monarch photo hexes, nor to annual future plant leaf projections, and their observation years 2018–2025 do not align with WorldClim historical 1970–2000. Do not multiply the USA-NPN YES fraction by the global 840-hex fraction or replace missing winter months with NO.

The only defensible immediate coupled next test is to take the **same real NPN monitored site geographies**, apply the real CMIP6 monthly temperature scenario and hold *observed month-status flags* fixed as **counterfactual leaf persistence**, reporting source and geographic matching separately. A genuine future leaf response would require modeling phenology versus actual annual weather with out-of-site validation and uncertainty.

## Attribution

"Data were provided by the USA National Phenology Network and the many participants who contribute to its Nature's Notebook program." Cite USA-NPN data service and exact API request provenance. Observer, individual and raw site IDs/coordinates are not redistributed.
