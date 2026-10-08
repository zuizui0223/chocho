# Revisited urban monarch milkweed patches: prospective stage-detection gate

**Date:** 2026-10-08. **Status:** EXPLORATORY NEGATIVE / STOP AS ECOLOGICAL DISCOVERY. **Scope:** `analysis/resource-homogenization-v01` only; main GEB manuscript unchanged.

## Why this was tested

Earlier analyses found a dramatic native/exotic seasonal Simpson reversal in cross-sectional ratio of late monarch instars to eggs, and large differences between two **native** milkweeds, *Asclepias fascicularis* versus *A. speciosa*. The key ecological concern was whether those descriptive ratios provide additional, transferable information about later larval detection, conditional on earlier egg input, patch size, season, current late instars and timing. Repeated visits to the same address with the same milkweed species allow prospective tests that one-time Jaccard maps do not.

Original field data: Erickson, Schultz & Crone (2025), Figshare DOI https://doi.org/10.6084/m9.figshare.25648644.v1. Source archived from original Figshare in GitHub Action 37735892534; no new observations. Exploratory prospective protocol `URBAN_MILKWEED_FORWARD_HOST_IDENTITY_PROTOCOL_V01.json` was written before fitting the forward prediction, **but after the earlier species-stage descriptive contrast was known**.

## Prespecified prediction

- Valid original field entries: 4,480 (all required eggs, 4th/5th instars and positive plant counts observed).
- Unique address × species × day aggregation: 4,412 observations, including 68 duplicate-within-day combinations summed before forward links.
- Unit: next consecutive recorded visit to identical route + address + milkweed species, within the **same calendar year**, between **10 and 35 days** after index observation.
- Selected from source record frequency only: `A. curassavica` (exotic), `A. fascicularis` (native), `A. speciosa` (native) and `Gomphocarpus physocarpus` (exotic). At least 100 2022 revisit transitions each.
- Outcome: probability of detecting any 4th/5th-instar monarch at next revisit. **Not larval survival or verified offspring from the earlier eggs.**
- Both models: current egg count, late-instar count, index and destination plant abundance, index calendar month (cyclic), revisit days.
- Compare botanical Native/Exotic classification against exact milkweed species identity; identical fixed-L2 logistic fit settings.
- 2022 leave-one-route-out validation; 2023 temporal holdout with 2022-only fit; 2024 descriptive diagnostic.
- Predefined pass: 2022 improvement with route-bootstrap 95% CI above zero **and** positive 2023 improvement.

Workflow [37763731931](https://github.com/zuizui0223/chocho/actions/runs/37763731931) succeeded. Artifact `butterfly-urban-monarch-forward-species-v01`, JSON `forward_species_prediction_v01.json`.

## Results

| Heldout dataset | Revisit pairs | Positive later-stage detections | Native/exotic-only model log loss | Exact-species model log loss | Improvement (origin minus species) |
|---|---:|---:|---:|---:|---:|
| 2022 leave-one-route-out (15 routes) | 1,901 | 226 | 0.29156285 | 0.29031472 | **+0.00124813** |
| 2023 fully held-out year | 622 | 30 | 0.14882431 | 0.14885149 | **-0.00002718** |
| 2024 sparse, exploratory year | 179 | 29 | 0.47006499 | 0.47082651 | **-0.00076153** |

2022 route-cluster (not individual-row) 4,999 bootstrap 95% CI for origin-minus-species log-loss gain: **[-0.00241532, 0.00449401]**, crossing zero. Routes show improvements and degradations, rather than a geographically uniform effect. 2023 temporal validation shows no gain. The **frozen gate FAILED**.

Within 2022, `A. speciosa` had 12 positive next-visit late-instar events over 233 repeat visits (5.2%), compared with `A. fascicularis` 42/289 (14.5%), `A. curassavica` 151/1173 (12.9%) and `G. physocarpus` 21/206 (10.2%). These are **raw detection probabilities, not comparable per-plant cohort survival**, and the two native species' difference should not be presented as a new physiological effect.

## Scientific decision

1. The original cross-sectional species-stage differences are valid **descriptive patterns** within the published survey dataset, but exact species identity did **not** improve independently evaluated next-visit late-instar detection in a reproducible, meaningful amount beyond native/exotic + prior eggs, stages, season, revisit delay and patch size.
2. Do not rescue via choosing favorable routes, tighter lag windows, different species subsets, selecting summer/autumn after seeing outcomes, increasing flexible model complexity, or treating 1,901 repeated rows as independent site replicates.
3. **Longitudinal patch revisits are not individual fate tracks.** Several cohorts and new eggs can exist between visits; transient plant removal/maintenance and stage-specific observer effort also matter. No recruitment, larval physiological development or fitness is identified.
4. Since all static homogenization, time-split hindcast, pseudo-host lag, field-use and dynamic predictions failed stronger ecological validation, this independent butterfly ecological manuscript **should not be promoted on present data**. The existing main GEB descriptive globalization paper is a different claim and stays unchanged.
5. The next genuinely relevant *new* ecology test needs individually followed offspring, a control plant, adult emergence and documented predators/parasitoids, plus plant identity and naturalization/cultivation status at matched sites. Published Diethelm et al. (2026) Dryad DOI 10.5061/dryad.51c59zwgx has such individual fates and a 2×2 host×cage experiment, but the official Dryad downloads returned HTTP 401/403 in GitHub workflow 37745988292 and the host×predator growth interaction was already published, so reanalysis is not inherently novel.

**No changes to the GEB main manuscript.**
