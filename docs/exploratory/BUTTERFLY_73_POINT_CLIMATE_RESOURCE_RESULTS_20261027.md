# 2026-10-27 butterfly niche opportunity: complete 73-point-source and TDWG3 co-location readout

**Status:** source-locked observational global macroecology, NOT realized butterfly climate-niche shift. Draft PR #74 stacks on #73 and #72. Original GEB PR #38/#70 and main unchanged.

## 1. Cohort reconciliation — source-based cases, not independent individuals

Original validated PR #72 country opportunity classes, exact WCVP/HOSTS/BCE version and 2010–2025 accepted GBIF index:

| | Novel-host-climate dominant source countries | Analogue-host-climate dominant source countries |
|---|---:|---:|
| Species×country GBIF-positive units | **32** | **41** |
| Original indexed butterfly occurrence records | **516** | **5,560** |
| Distinct butterfly species represented | 15 | 15 |
| New full-source accepted species/country queries | 32/32 complete | 41/41 complete |
| Adult lifeStage source labels | 206 | 2,846 |
| Juvenile/immature source labels | **0** | 87 |
| Unspecified lifeStage labels | 299 | 2,520 |
| Point coordinate uncertainty reported ≤1km | 260 | 2,918 |

The complete 73-case live independent retrieval exactly reproduced the original 6,076 indexed GBIF records, with **88/88 successful API calls, zero new country count drift, and 10 independent source-control tests passed** in [GitHub Actions #40432323278](https://github.com/zuizui0223/chocho/actions/runs/40432323278). The original complete novel+analogue roster was separately reconstructed with 5 PASS source independent tests [#40431411694](https://github.com/zuizui0223/chocho/actions/runs/40431411694). [Original source-only derived artifact](https://github.com/zuizui0223/chocho/actions/runs/40432323278/artifacts/12701987908) contains case counts, source fingerprints, source metadata only (no individual coordinates).

**Source dependency:** two observed butterfly×country cases, *Pseudozizeeria maha* / Taiwan (2,220 indices) and *Catopsilia pomona* / Taiwan (2,198 indices) account for **4,418 / 6,076 = 72.7%** of pooled index records, both in the analogue-dominant source country class. Their removal must be displayed rather than treating records as i.i.d. wild butterflies.

## 2. Tested the claimed mapping using original source-verified WCVP/TDWG polygons

[GitHub Actions #40433784232](https://github.com/zuizui0223/chocho/actions/runs/40433784232) **SUCCESS**, 14 regression tests PASS, full exact original author's offline source ZIP, SHA256 `408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323`, same 239 butterflies, 235 native-host climatic supports and 365 WGSRPD3 host-resource regions. **73/73 positive case queries complete with 6,076/6,076 accepted point indices;** 88 source requests. The raw original plant distribution polygons and GBIF point coordinates were never published.

Primary **native-host CHELSA regional-climate 5–95%** rectangle:

| Original country class | In introduced-added regional analogue resource | In introduced-added regional climate-novel resource | No mapped added regional resource | Other unmapped / native / ambiguous |
|---|---:|---:|---:|---:|
| Novel country (516 point indices) | 4 | **450** | 3 | 59 unmapped |
| Analogue country (5,560) | **5,306** | 1 | 0 | 225 unmapped + 27 native-region + 1 ambiguous |

**Both source host country labels are mostly spatially corroborated at WGSRPD3 region resolution**, but this is partly *built into the selection criterion*: to be called novel/analogue the original case must have ≥80% of the source introduced-host **area** inside/outside the source rectangular climate window. It is not an independent causal association. For source-novel cases, the source country-area novel-host fraction averaged **94.68% when weighted by the observed butterfly case record counts**, versus **87.21%** of the observed butterfly indices geographically falling in corresponding novel host regions (including the source unmapped points in the denominator). A uniform-area country null is not a valid wildlife survey expectation.

The point-to-region result also shows nonrandom source-coverage failures: **59/516 = 11.4%** of recorded butterfly points in source-novel countries cannot be matched to a valid modelled TDWG3 polygon; source-analogue countries have 225/5560 = 4.0% unmapped plus one ambiguous. Mapping missingness should never be treated as no butterfly or no host.

**Critical sensitivity**: changing the climate analogue from 5–95% to 10–90% on the SAME original regional CHELSA source changes its biological meaning sharply:
- Novel-dominant source countries: 450 → 453 points classed in climate-novel introduced-host regions.
- Original analogue-dominant source countries: **1 → 4,253** points classed in climate-novel introduced-host regions, because the originally labelled analogous region climate lies outside the narrower rectangle for many source taxa.

The result is not evidence of 4,253 new host introductions or butterfly climatic niche evolution: it is an **analyst's threshold-driven reclassification of unchanged occurrences and unchanged plant distributions**. As a result the binary “novel niche” label must always be accompanied by its climate threshold and full geographic/source uncertainty.

## 3. What is still necessary

The **actual point climate** (rather than a broad botanical-region CHELSA mean) has not been assessed by this result. A strictly separate, same-product WorldClim v2.1 1970–2000 BIO1/BIO12 point and source native-host range raster test is now implemented; it independently recomputes the reference from the same WorldClim grid before comparing butterfly pixel climates. This is an **exploratory sensitivity**, not direct comparability to CHELSA 1981–2010 percentages. The correct biological claim will remain *potential host resource geography* even if butterfly point climatic novelty is found: true native-butterfly climatic baseline, plant introduction dates, local host use, independent wild populations, breeding and standardized effort remain absent.

Source values: [Exact region-source and all 73 point-index action artifact](https://github.com/zuizui0223/chocho/actions/runs/40433784232/artifacts/12702868554), [executable script](../../scripts/audit_butterfly_73_point_tdwg_resource_20261027.py), [frozen protocol](../../provenance/reviewer_defenses/protocols/BUTTERFLY_73_POINT_WGSRPD_RESOURCE_20261027.md).
