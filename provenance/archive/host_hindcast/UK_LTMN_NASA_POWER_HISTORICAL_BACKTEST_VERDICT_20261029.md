# UKBMS historical 1981–95 → 2005–21 NASA POWER out-of-time weather baseline: failed predictive gain

**29 October 2026; exact original-data run SUCCESS but the weather predictive hypothesis FAILS.** This is an exploratory climate-only benchmark based on actually visited UK butterfly monitoring transects. **No historical host-plant mediation was measured; no future CMIP6 prediction is authorized by this result.**

## Actual observational and climate sources

- **19** source-verified fixed UKBMS traditional transects originally visited in both eras, each with precisely matched visit date/time and source butterfly species counts, georeferenced through the Natural England LTMN transect metadata. Training era **1981–1995**, held-out test era **2005–2021**; this replaced the original 1976–95 ERA5 protocol following Open-Meteo HTTP 429; neither the model era nor the climate product is interchangeable with the original protocol.
- **27** fixed chocho butterfly taxa with UKBMS common-name crosswalk, including taxa not positive in a particular period. Original survey visit denominator requires **≥3 other positively recorded non-ambiguous butterfly names**, independent of whether focal butterfly is detected; missing focal count is at most visit-conditional non-detection, NOT site absence or extirpation.
- NASA POWER daily T2M and PRECTOTCORR, 1981–2021, **14,975 dates for each of 19 transects, zero invalid weather values**, approximately coarse 0.5°×0.625° MERRA-2-based field, not a true transect microclimate. All 19 source response SHA256 values retained without original coordinates, station identifiers or person records. 115,249 early and 138,577 late eligible **site×species×visit labels**; 23,060 and 29,456 respectively were detected. There are only 19 geographically independent transects; 253,826 is not a biological sample size for inference.

## Source-tested temporal prediction

Both models are L2 logistic (C=0.1), fixed taxon dummies, flight-day seasonal harmonics, transect latitude/longitude, linear calendar year, training-only feature scaling and transect×species equal total weights. The weather model alone adds actual **preceding 30 days** NASA temperature and precipitation to exactly the same historical source records.

| 2005–2021 true observed test | Without weather | With prior-30-day NASA weather | Weather gain |
|---|---:|---:|---:|
| Site×species-equal log-loss (lower better) | **0.410004** | **0.418668** | **−0.008664** |
| Site×species-equal Brier (lower better) | **0.132748** | **0.135459** | **−0.002710** |

Five independent site-held-out test groups have **1/5** positive weather increments: fold 1 −0.017925, fold 2 −0.001843, fold 3 −0.006077, fold 4 −0.011404, fold 5 +0.000073. The prospectively fixed criterion was ≥4/5 positive folds AND overall improvement. **Fails**. Original real-data run [#40802452296](https://github.com/zuizui0223/chocho/actions/runs/40802452296) SUCCESS, independently repeated initial run [#40802418044](https://github.com/zuizui0223/chocho/actions/runs/40802418044) also SUCCESS. The corrected final-run source cohort uses 19 site panels and exactly the original 3-other-positive-butterflies recording gate.

## Why this is important, without overclaiming

- The outcome is **whether a focal species was detected on an actual monitored visit**, not range-edge tracking or a true biological absence. A failed extra weather predictor can reflect coarse NASA POWER grid, mismatched historic/modern butterfly phenology, changed survey effort/land use, original linear-year extrapolation, confounding, temporal novelty or weak actual weather effect. It **does not establish that temperature was irrelevant to butterfly population or warming**.
- It is scientifically inappropriate to tune new warming features or change the scoring period in response to the negative result, then present an optimized future model as prospectively validated. Maintain this result permanently and separately.
- **Main interaction test is not yet performed**. The actual question is whether independently *pre-test-era* plant host availability helps predict historical butterfly detection better than M0 climate/effort/habitat, and whether uniquely *human-planted/maintained hosts* can be separated from wild flora. Static 2010–2025 chocho host geography is not legitimate for predicting 2005–2021; broad habitat categories are not proof of a host species or management.
- A NEW independent source became available: the **BSBI Plant Atlas 2020** 10km×survey-era atlas (Zenodo DOI 10.5281/zenodo.10998003, original file PA2020_10by10kmSquareData.zip) and original UK/EU butterfly–larval-plant relation tables (Clarke 2024, Dryad DOI 10.5061/dryad.1vhhmgr12). Using pre-2000 plant records only can build a distinct *historical observed host presence* source, with equal-epoch nonhost plant recorder-intensity controls. These are plant presences, **not established absences and not human horticultural maintenance**. Do not add them to the main paper until historical held-out predictions pass.

## Provenance / scientific restrictions

The complete source-verified aggregate JSON is frozen alongside this report. GitHub Actions artifact 12866801448 compressed ZIP SHA256 `08eafca509d689317a97171dd680ba472cb07a502bafbd8a6f16e334ec2a2734`; inner original JSON SHA256 `18daee3d584983e03822a121601607359da7ed6760d37075bfd33227c895803a`. The independently obtained original data were never published as individual site/person/butterfly records. UKBMS © Butterfly Conservation, UKCEH, BTO and JNCC; NASA POWER acknowledgement and terms as cited in protocol.

**Decision:** Future CMIP6 requires a new calibrated and held-out validated explanatory model. This run passes source reproducibility, **fails weather prediction**, and does not identify host benefit or climate tracking.
