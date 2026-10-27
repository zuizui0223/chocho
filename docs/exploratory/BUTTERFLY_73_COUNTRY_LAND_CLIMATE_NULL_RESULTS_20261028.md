# 2026-10-28 — The apparent butterfly/introduced-host climatic contrast is closely matched by country climate opportunity

**Source-complete exploratory environmental null — no inferred climatic niche expansion, colonization, larval feeding, or host mediation.** Dedicated independently executed [GitHub Actions #40534253657](https://github.com/zuizui0223/chocho/actions/runs/40534253657) **SUCCESS**, same 73 source butterfly × country positive cases as PR #74 and PR #87. This is a test of available geography, **not** a survey-effort-matched observation null. `main`, the manuscript and the GEB peer-review materials are unchanged.

## Fixed real-data inputs and outcome-independent null

- **73** original source-positive butterfly species×country cells; **20** butterfly species, **51** source countries, **6,076** original GBIF indexed point records. The observed species-paired comparison from PR #87 was **10 informative species**. Countries with both source classes for *different* butterfly species: CU, DO, GT, LA, LK, VE (only six); they do not yield a country-fixed causal design.
- Frozen WorldClim v2.1 **1970–2000 BIO1/BIO12 10-arcminute**, original SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`. Original source WorldClim butterfly-record and native-*HOST-PLANT* envelope JSON SHA256 `5d3ca4137821e60aaca7e1bf965c509b1cddee76c643fb00ee15df44ae775937` (from [Actions #40435313884](https://github.com/zuizui0223/chocho/actions/runs/40435313884)).
- Original audited Natural Earth 1:50m map-unit source SHA256 `b8d421aca6e9e08e8cdf09cc26af111cc3e0deba4fe915611d58ade71e8a4db0`. Use **valid country land-climate pixel centers** with `cos(latitude)` weights as climate availability; same species' unchanged native-host regional BIO1/BIO12 5–95% rectangle classifies both the observed records and every environmental control pixel. All 73 country cases have ≥10 valid reference pixel centers; no omitted country in primary.
- This is **not** drawing random butterfly localities from the same observation effort. It only estimates what a random valid land pixel's climatological novelty would be within that country. Original CHELSA area-defined *host-resource source class* remains fixed, no class switching or species postselection.

## Empirical result — geographic sorting is sufficient to reproduce nearly all of the unadjusted contrast

| Contrast among paired butterfly species | Equal-species mean percentage points |
|---|---:|
| Observed butterfly point-climate: novel-host-source country minus analogue country | **+42.42** |
| Expected from the country land-area WorldClim climate distribution alone | **+41.07** |
| Residual after subtracting that country-area expectation | **+1.35** |

The ten paired species show **5 positive versus 5 negative adjusted residuals**. Descriptive equal-species bootstrap 95% percentile span (50,000 seeded resamples) for the mean residual is **−24.12 to +29.39 pp**; this is **NOT** a population confidence interval, a causal effect interval or an equivalence test. Case-level and species-level ecological dependence remain unmodelled.

### All prespecified sensitivity checks

| Sensitivity of original cases / geography | Paired species | Observed difference (pp) | Country-environment expected difference (pp) | Residual (pp) |
|---|---:|---:|---:|---:|
| Original 5–95% native-host support, >=1 record/case | **10** | +42.42 | +41.07 | **+1.35** |
| >=3 records/case | **2** | +33.05 | +32.07 | **+0.98** |
| >=5 records/case | **2** | +33.81 | +31.37 | **+2.44** |
| >=10 records/case | **2** | +39.24 | +32.13 | +7.12 |
| Narrower 10–90% native-host climate bounds | 10 | +63.45 | +43.72 | **+19.73** |
| Equal-pixel 5–95% native-host bounds | 10 | +46.39 | +41.54 | +4.85 |
| Butterfly coordinates with known uncertainty <=1 km | **9** | +36.57 | +47.84 | **−11.27** |
| Country polygon all-touched raster boundary sensitivity | 10 | +42.42 | +41.72 | +0.70 |
| Exclude countries with <30 10-arcminute valid land pixels | 10 | +41.44 | +40.33 | +1.11 |

At >=5 original butterfly records/country, only *Mechanitis lysimnia* (**+16.69 pp** adjusted) and *Zizula hylax* (**−11.80 pp** adjusted) retain both climate-source classes. They show **opposite** geographic-adjusted directions. The narrower 10–90% interval changes the residual sharply; the primary point is therefore a **lack of robust excess beyond environmental availability**, not evidence of a exactly-zero biological response.

## What follows biologically — claims retained, rejected, and next test

**Retain:** geographic host redistribution changes potential larval resource space (separate frozen GEB evidence). The source-positive butterfly observations are preferentially in climates unlike mapped native larval HOST-PLANT geography under the source definition, but most of this observed climate contrast is matched by country-level available-climate differences. **The high unadjusted contrast should not be described as interaction-mediated abiotic niche release.** This rules out a simple unadjusted macro signal as compelling evidence of the proposed mechanism, not the mechanism itself.

**Do not claim:** non-native larval HOST ingestion; wild juvenile recruitment; a shift beyond each butterfly's own original breeding-climate niche; net climatic adaptation; increased abundance; new generations; historical colonisation after a known host arrival; local host vs nonhost selective preference; equivalence of observed and background distribution. Climate class stems from host-plant ranges, while a country's full land climate includes unsampled agricultural/urban habitats and inaccessible terrain. The climatological reference is fixed 1970–2000 not annual contemporaneous microclimates. Focal butterfly GBIF records are not spatially effort-matched to country land pixels.

**Next independent primary study:** test preselected species×country×year wild records vs *other butterfly observations from the same data source/time/administrative area*, using focal butterfly native-breeding-climate support and fine-scale introduced-host presence where possible, and verify actual larval recruitment separately. Test explicitly whether introduced-host landscapes are exceptional relative to **non-host plant** landscapes and to survey effort; the current country-area null cannot answer these.

## Reproducibility and source receipts

- Pre-outcome fixed protocol: `provenance/reviewer_defenses/protocols/BUTTERFLY_73_COUNTRY_LAND_CLIMATE_NULL_20261028.md`.
- Analysis: `scripts/analyze_butterfly_73_country_climate_area_null_20261028.py`.
- Frozen 20 host-derived WorldClim bounds: `data/exploratory/butterfly_73_native_host_worldclim_reference_bounds_20261028.json`.
- Source-complete aggregate result: `provenance/reviewer_defenses/results/butterfly_73_country_land_climate_null_20261028.json` (JSON SHA256 `5e4d4c31536ee18ec0b907f6d086e32bb1d8d7bf0d1de05836d7baeebe1597ee`; GitHub Action ZIP artifact `12747279612`, digest SHA256 `e8895fb81a5838d94aa0c8960b4c932248f07b9e02e2917082abcb4942194e11`). No butterfly raw point coordinate / observer, botanical author-only polygon, or unapproved climate grid exported.

Command requires only public frozen sources and the vendored derived reference (the Actions workflow downloads/verifies both external source archives):

```bash
python scripts/analyze_butterfly_73_country_climate_area_null_20261028.py \
  --worldclim-zip <verified-WorldClim-2.1-10m-BIO.zip> \
  --natural-earth <verified-Natural-Earth-50m-ISO-map-units.geojson> \
  --output /tmp/chocho-country-climate-null.json
pytest -q tests/test_butterfly_73_country_climate_area_null_20261028.py
```
