# First source-locked port-host × marine-access pilot (14 November 2026)

**Status:** Executed in [Actions run #43577791613](https://github.com/zuizui0223/chocho/actions/runs/43577791613) (SUCCESS). An independent in-memory JavaScript rebuild of the same frozen CSV tables reproduced the Python point estimates and approximate uncertainty intervals to numerical precision. JSON with exact source SHA-256s is in [the workflow artifact](https://github.com/zuizui0223/chocho/actions/runs/43577791613/artifacts/14106977646).

**No confirmation of arrival, larval feeding, island establishment or shipborne dispersal.** This exploration is intentionally **not** merged into the main GEB Results, which retain the original global host resource and Circuitscape community evidence.

## Biological question

Does having at least one known larval host *indexed as a regionally introduced candidate near a destination island's port in 2001–2010* strengthen the relationship between subsequent butterfly first indexing and modeled marine accessibility **beyond great-circle distance**?

The proposed story distinguishes a historically **indexed resource** at the destination (candidate establishment environment) from the **geographic access** of that species based on earlier GBIF indexing on other islands. The port is **not assumed to be the butterfly's transport route**.

## Actual input crosswalk

- Frozen original 225 butterflies × 51 selected GSHHG landmass groups = 11,475 cells for 2001–2010 / 2011–2017 / 2018–2022.
- Frozen older 2001–2010 botanical records in 5-km WPI port buffers, candidate native/introduced larval-host plant identities, and port-specific plant record effort. The 225-species historic botanical crosswalk is present for **41/51** islands. Unrecorded host taxa do **not** establish host absence.
- Frozen 1,275 undirected pairwise Circuitscape uniform-ocean resistances, based on the same 51 geographical landmass IDs. Original AIS and wind fields were **not** used in this first time-ordered crosswalk because their dated predictors could not establish the history of actual shipping or flight vectors.
- A species-island candidate enters the source-conditioned test only if no previous record on that destination was indexed and a conspecific was indexed on another of the 51 islands in an earlier full period. The *nearest such indexed source* supplies island-to-island geographic distance and uniform-ocean effective resistance.
- Access = the **negative residual** from a quadratic regression of log(1+uniform-ocean resistance) on log(1+distance) over all 1,275 geography pairs. Greater values indicate lower resistance than expected at the given distance. The residualization is butterfly-outcome-blind.

## Result: no independently detected positive interaction

The response is **first GBIF indexing within a multi-year period**. A linear probability model includes destination-island and butterfly-species fixed effects, native-host port-index indicator, candidate introduced-host port-index indicator, geographic distance, sea-access residual, and their **introduced-host × standardized sea-access** interaction. Approximate uncertainty uses a two-way island/species sandwich; it does not establish confirmatory causal inference.

| Period | Analysis population | Islands | Species×island risk cells | First-index records | Host-positive first-index records | Interaction coefficient (percentage points per 1 SD of access), approximate 95% interval |
|---|---|---:|---:|---:|---:|---|
| 2011–2017 | All indexed ports | 41 | 4,148 | 160 | 13 | −1.105 [−3.920, +1.709] |
| 2011–2017 | **≥10 historical plant records** | 36 | 3,625 | 140 | 13 | **−0.964 [−3.845, +1.917]** |
| 2011–2017 | Within-species AND within-island host-status support | 29 | 2,590 | 93 | 13 | −0.942 [−3.708, +1.825] |
| 2018–2022 | All indexed ports | 41 | 4,584 | 79 | 9 | −0.434 [−2.064, +1.197] |
| 2018–2022 | **≥10 historical plant records** | 36 | 4,006 | 62 | 8 | **−0.290 [−1.881, +1.302]** |
| 2018–2022 | Within-species AND within-island host-status support | 29 | 2,731 | 35 | 8 | +0.114 [−1.367, +1.596] |

These are different source-conditioned subsets and **not** six independent replications. The exploratory interaction is not positive in the two principal periods, and the confidence intervals permit negative or positive modest effects. There are only 13 and 8 first-index events with a candidate introduced port host in the primary ≥10-record panels. Thus it would be incorrect to claim the general hypothesis is disproved or that ocean resistance lacks information about later butterfly assemblages.

## Relation to the existing positive results

This pilot **does not overwrite or refute other endpoints**:

1. The original 239-butterfly fixed-diet host reconstruction added +54.9% *potential* resource opportunity globally, +118.9% on islands and +40.9% on mainland botanical regions.
2. More ports were associated with more non-native floristic taxa, including candidate butterfly hosts, across the distinct 98 insular botanical regions, but there was no host-selective increment after whole-introduced-flora adjustment.
3. In the separate 51-island community similarity target, adding uniform-ocean Circuitscape to geographic controls reduced held-island MSE by 5.082%. Adding AIS/wind fields beyond pooled 2001–2024 host inventories reduced MSE by 6.488%, with documented strict geographic holdout failures.
4. Earlier full route-free annual port-host analyses found a positive point-estimate first-index contrast in a different 74-island 2011–2024 population, sensitive to overlap and source inclusion. This is not an exact re-test of those periods, eligibility criteria or treatments.

## Biological and statistical boundaries

- Original 51 islands were selected through prior port/host/butterfly indexing; **not a probability sample of oceanic islands**.
- Plant origin is based on modern, **regional** WCVP status; a 2001–2010 GBIF plant index is **not** documented island-local introduction, abundance, persistence or larval use.
- Plant **non-index** is unknown, not a biological absence. Excluding zero botanical effort does not eliminate heterogeneous recording probability. The ≥10-record threshold is a sensitivity cutoff, not an effort-matched randomized comparison.
- Species-specific earlier indexed source is only the nearest **among these 51 selected islands**; neither a verified exporting source nor the nearest true population.
- First GBIF indexing is **not** actual new arrival, successful breeding, population persistence or onward colonization.
- Sea resistance is undirected modeled geography, not measured wind-assisted flight, ocean current drift or ship-borne cargo. Available wind is undirected speed, not a resolved wind-vector path.
- First-index events within each time window are not dated to a single year, and the regional origin status remains modern.
- Fixed-effect LPM cluster intervals are approximate, especially with sparse introduced-host first-index events. No experiment, multi-outcome correction or externally validated predictive interaction is claimed.

## Reproduce, inspect or extend

From repo root, with NumPy 1.26+:

\`\`\`bash
python scripts/analyze_port_host_ocean_interaction_20261114.py \
  --out-dir results/port_host_ocean_interaction_20261114
\`\`\`

The generated \`source_join_pilot.json\` preserves four input SHA-256 hashes, exact populations, coefficients, interval approximations and interpretation limits. The dedicated CI workflow checks frozen original source dimensions and both periods.

**Next experiment should not be response-driven fishing among AIS, wind, butterfly wing size and multiple access metrics.** A genuine test of the proposed mechanism needs prior-dated *island-local* host abundance/use, independently measured dispersal capability, equal butterfly observation opportunity, and breeding/persistence events. Monthly **directional** winds and currents can then be added to an independent arrival model; the current undirected circuits cannot substitute for such trajectories.
