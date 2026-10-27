# Pre-outcome geography-only null: can receiving-country climates explain the 73 butterfly point contrast?

**Frozen before opening new WorldClim country-area outcomes on 2026-10-28.** This is an independent *exploratory follow-up* to PR #87/74, not a prospective pre-registration or a new butterfly occupancy dataset. Leave frozen GEB main/submission materials unchanged.

## Question and selection

The original source-locked 73 butterfly × country positive cases (20 species, 6,076 GBIF point indices) contrast 32 **host-resource-climate-novel country** cases against 41 analogue cases. The existing posthoc species-paired observed contrast is +42.42 pp for 10 species. The same geographical classification is partly defined by climate; therefore an observed climate difference can be environmental availability alone. Can we subtract the climate distribution of the *receiving country itself* without using focal butterfly positions to define the null?

The cohort, country labels, counts, and native-**host-plant**, not butterfly, WorldClim BIO1/BIO12 reference intervals are **identical** to validated source run #40435313884. These were extracted before this null, and their original JSON/SHA, WorldClim raster ZIP SHA and original country vector SHA are fixed; no species/country can be selected based on the new null outcome.

## Reference and estimand

1. Download identical **WorldClim 2.1 BIO1/BIO12 1970–2000 10 arcminute** full raster; SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`. Reject any other product. Use the original 20 species' unchanged **native HOST-PLANT regional** area-weighted 5–95% BIO1/BIO12 rectangular bounds, not a native butterfly tolerance baseline.
2. Download originally audited Natural Earth **1:50m administrative map units** using commit `117488dc884bad03366ff727eca013e434615127`, SHA256 `b8d421aca6e9e08e8cdf09cc26af111cc3e0deba4fe915611d58ade71e8a4db0`; dissolve by original ISO_A2_EH rule (never generic country names). Only valid WorldClim land pixel **centers inside** receiving-country boundary contribute; approximate equal-area weights as `cos(latitude)`. Require ≥10 valid country climate pixels. Zero pixels or insufficient pixels is unscorable, **not zero climate opportunity**.
3. For source case `s,c`, observed fraction `O_sc = outside source native-host rectangle GBIF indices / all original indexed points`. Calculate **country-land geographic opportunity** `G_sc = area-weighted country WorldClim pixels outside the same rectangle / valid country land-climate pixels`. Compute `R_sc = O_sc − G_sc`. Both source classes are scored on exactly the same climate product/reference. `G_sc` is available climate, **not observed butterfly sampling effort**.
4. For each butterfly represented in BOTH original climate classes, average the country case values `O_sc` and `G_sc` equally **within each class**; calculate `D_obs,s = mean(O|novel)−mean(O|analogue)` and `D_geo,s = mean(G|novel)−mean(G|analogue)`, then `D_resid,s = D_obs,s−D_geo,s`. Average these differences equally over paired species. **Do not** treat 6,076 points as 6,076 independent units or compute ecological exchangeability/permutation p-values.

## Locked before result: sensitivities, stop criteria, wording

- Minimum butterfly GBIF point indices ≥3, ≥5 and ≥10 per case, each requiring at least one case in each source class **within the same species**.
- Native-host WorldClim 10–90% rectangle and equal-pixel rather than area-weighted 5–95% native-host reference, without reoptimizing thresholds. Strict original documented butterfly coordinate uncertainty ≤1km as a separate sensitivity; no made-up uncertainty for missing metadata.
- Natural Earth boundary `all_touched=True` sensitivity, and minimum 30 country climate pixels. Report omitted cases and missing geography. Exclude no taxa or countries after looking at outcomes; if the primary panel becomes insufficient, report that explicitly.
- Full environmental area baseline includes urban areas, crops, forests, alpine climates and all valid pixels regardless of site visit, butterfly ability to reach a pixel, actual planted host establishment, phenology, elevation microclimate or season. **It cannot correct observation effort or supply a real nonhost-plant negative control.** Positive residuals may simply mark observer hotspots; zero/negative residuals cannot prove failure of host-mediated climate entry.
- No baseline treated as an independent native butterfly breeding-climate envelope, no adult index as breeding, no plant occurrence as larval ingestion, no temporal colonization claim. The comparison is diagnostic for **geographic opportunity confounding only**, not host causality. Only aggregate source/derived country/case fractions and SHA receipts can leave the run; NO exact GBIF observation coordinates, observers, or author-only plant range polygons exported.

## Next test after the null

If sufficient paired species remain, use an independent **survey-effort-matched reference** from other butterfly sightings within the same country / year / data source, and ultimately actual larval recruitment on a local introduced host vs non-host landscape. The country-area null does not substitute for those biological endpoints.
