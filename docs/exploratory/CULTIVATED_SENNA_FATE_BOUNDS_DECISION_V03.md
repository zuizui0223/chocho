# Horticultural host exposure versus demographic effects — fate audit v0.3

**Date:** 2026-10-08. **Branch:** `analysis/resource-homogenization-v01`. **Status:** exploratory; not an independent butterfly-fitness result.

## Where the ecological claim stands

In frozen WCVP, the cultivated Caribbean species *Senna polyphylla* has no mapped Florida occurrence. Public iNaturalist records include 59 labelled cultivated observations in 46 approximately 5-km reporting cells across ten observation years between 2010 and 2025. Six cells contain cultivated reports in at least two different years, and one cell has reports in six years (2010, 2011, and 2021–2024). These are **reports**, not a planting census, independent plant-persistence proof, or larval-use observations.

This repeated horticultural footprint is compatible with plant-resource exposure outside the *naturalized-plant geographic inventory*. Importantly, the original Koptur et al. (2024) study documents larval rearing on cultivated alien *S. polyphylla* in Florida but **both native and nonnative Senna in that study were sampled in cultivation**.

Sources:
- Frozen plant distribution reconstruction: [GitHub Actions run 37716723493](https://github.com/zuizui0223/chocho/actions/runs/37716723493).
- Frozen botanical source record pilot: [GitHub Actions run 37723666255](https://github.com/zuizui0223/chocho/actions/runs/37723666255).
- Independently sourced rearing article: Koptur et al. (2024), *Insects* 15:123, [DOI 10.3390/insects15020123](https://doi.org/10.3390/insects15020123).

## Resolution of fitness outcomes matters

The rearing paper reported 305 immature Phoebis on native Senna and 389 on introduced Senna. Of these, only 235 and 283 respectively had a known adult or parasitoid emergence. Based on the **published site rows**, not the inconsistent grand total, there were 16 native-host and 49 introduced-host parasitoid emergences. Adult emergence was therefore confirmed for 219 native-host and 234 introduced-host immatures. **176 fates remained unresolved** (70 native, 106 introduced).

| Outcome | Native | Introduced |
|---|---:|---:|
| Immatures found | 305 | 389 |
| Outcomes resolved | 235 | 283 |
| Confirmed adult emergence | 219 | 234 |
| Confirmed parasitoid emergence | 16 | 49 |
| Unresolved fates | 70 | 106 |
| Adult-emergence fraction of all initially found, logical bounds | 71.8–94.8% | 60.2–87.4% |
| Parasitoid emergence fraction of all initially found, logical bounds | 5.2–28.2% | 12.6–39.8% |

The **introduced-minus-native adult-emergence risk difference** could mathematically range from **−34.6 to +15.6 percentage points**, given only these outcome counts. The difference in parasitoid outcomes could range from **−15.6 to +34.6 percentage points**. Both signs are possible, so the ecological effect on whole-cohort fate is **not identified** by these aggregate data.

This does **not** mean all extreme combinations are biologically plausible; these are transparent worst-case bounds on what the incomplete observed fates can establish. Within the subset of individuals whose outcome was resolved, parasitoid emergence is more frequent on introduced Senna at each of the three monitored sites. Because species identities, sites, host quality, survey timing and repeated shrubs are confounded, that association is not an experimental origin effect.

Reproduction: `scripts/audit_senna_rearing_fate_bounds.py` reads the source-verified frozen benchmark `benchmarks/exploratory/pieridae_fabaceae_empirical_quality_decision_v0.1.json`. The original paper's grand-total native parasitoid count (17) disagrees with site-row addition (16); both are explicitly documented rather than silently reconciled.

## The outstanding independent consumer endpoint

A separate, **pre-outcome** original-source audit of the 59 cultivated *S. polyphylla* observations and wild `Phoebis sennae` / `P. philea` larval observation fields is specified in `CULTIVATED_SENNA_OBSERVER_AND_DIRECT_LARVAL_USE_PROTOCOL_V02.json`. The launched [GitHub workflow run 37728156494](https://github.com/zuizui0223/chocho/actions/runs/37728156494) was still **queued at this note's preparation**. Until its artifact has been checked, no claim can be made about observer independence or directly witnessed larval association in these gardens. Even explicit field annotation remains an observer assertion, not verified adult emergence.

A prior GloBI feeding-record panel failed its strict direct-use coverage gate; do not lower that gate now. The new observer-source audit is an opportunity to understand provenance, not a substitute for rearing data.

## Research implication and publication boundary

The biologically meaningful, falsifiable hypothesis is that **cultivated and naturalized larval hosts can provide distinct pathways of exposure, while realized butterfly demographic benefit depends on host species identity, local accessibility, phenology and enemy pressure**. This is a hypothesis, not a demonstrated universal mechanism.

Published garden-placement experiments already show that host accessibility matters independently of plant count (Baker & Potter 2019, DOI 10.3389/fevo.2019.00474), and host-quality / exotic-trap outcomes are already known (Braga 2023, DOI 10.1016/j.cois.2023.101074). The new *measurement discrepancy* here is real; however **a distinct, causally demonstrated consumer-level ecological finding has not yet been produced**. Real progress would require replicated direct host-use and adult-emergence endpoints, plant density and management status, detection effort and predator/parasitoid outcomes across multiple sites and butterflies.

**Stop:** If the direct-source larval-use branch fails, retain the horticultural mapping discrepancy as a limitation/sensitivity of the GEB resource reconstruction, not a second ecology paper. Do not modify the GEB main submission manuscript.
