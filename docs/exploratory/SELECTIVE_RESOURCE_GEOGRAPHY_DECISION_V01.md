# Selective geographic butterfly resource bridging: exploratory decision v0.1

**Date:** 2026-10-08. **Scope:** `analysis/resource-homogenization-v01` only. **Do not change the GEB main manuscript.**

## Decision status

A global increase in reconstructed regional butterfly-resource similarity is largely an arithmetic consequence of added opportunity. However, spatial selectivity in that increase is not determined by the numbers added alone. A fixed-margin geographical null and a distinct within-host-plant-family link-rewiring null reveal nonrandom spatial structure.

The first two tests have been computed against 499 draws locally; the latitude-only test has independently completed in GitHub Actions run 37706116468. The climate and host-identity GitHub Actions workflow replications were queued at the time of this note. **All analyses remain exploratory and post-hoc relative to the GEB submission.**

## Predefined test 1: geographic placement at constant butterfly-level addition margins

Frozen protocol: `butterfly_selective_geographic_bridging_protocol_v0.1.json`.

- Panel: 239 butterflies; 355 native-resource WGSRPD3 regions.
- Null: preserve each butterfly's added WGSRPD3 count within each WGSRPD Level1, each region's added butterfly-resource count and all native resource incidences. 499 swap samples; seed 20261008.
- Distant: >=3000 km and different WGSRPD Level1 codes.
- Latitudinal analogue: absolute-latitude difference <=10 degrees; discordant: >=20 degrees.
- Primary contrast = (analogue mean contemporary-native Jaccard gain) minus (discordant mean gain).

| Endpoint | Observed | Fixed-margin median | Observed minus null |
|---|---:|---:|---:|
| Analogue gain | 0.255845 | 0.227856 | +0.027988 |
| Discordant gain | 0.170663 | 0.186736 | -0.016072 |
| **Difference in gains** | **0.085181** | **0.041111** | **+0.044070** |

Monte Carlo one-sided **p=0.002** for the primary contrast, 499 permutations. GitHub Actions output: `butterfly-selective-geographic-bridging-v01`, run 37706116468. Latitude alone is **not climate**.

## Predefined test 2: independent four-variable climate analogue

Frozen protocol: `butterfly_climate_analog_bridging_protocol_v0.2.json`, before the actual climate group outcomes were inspected. Climate source: frozen TTF run 36270581743 `unit_climate_crossfit.csv`, with the four region-constant BIO summaries (BIO1, BIO7, log1p(BIO12), BIO15), standardized across valid regions. Distant cross-Level1 climate-analogue pairs are the bottom quartile in 4-D climate distance, discordant the top quartile.

Both groups contain **12,944** region pairs. The same butterfly-region fixed-margin null is applied.

| Endpoint | Observed gain | Fixed-margin median | Observed minus null |
|---|---:|---:|---:|
| Climate analogue | 0.261893 | 0.236605 | +0.025288 |
| Climate discordant | 0.163456 | 0.175572 | -0.012116 |
| **Difference in gains** | **0.098437** | **0.061050** | **+0.037387** |

Monte Carlo one-sided **p=0.002** (local 499-run calculation). Climatic selectivity corroborates the latitude-based pattern, but not a realized butterfly response.

## Consumer specificity gate

With the SAME 1,706 plant species that are documented hosts of the 239 butterflies, descriptive plant-assemblage Jaccard gain is 0.145757 in distant climate-analogue region pairs and 0.053349 in distant climate-discordant pairs: difference **0.092408**. The butterfly resource-assemblage contrast is **0.098437**. These are different resolutions and must NOT be interpreted directly as an estimable amplification coefficient. The qualitative pattern already exists in the host-plant flora.

An additional null fixes every host plant's native and contemporary geography and performs bipartite double-edge swaps ONLY among butterfly-host links within the SAME botanical host-plant family, preserving each butterfly's known-host counts within every botanical family and each plant's number of focal butterfly consumers (2,596 links, 107 botanical families).

- Observed climate contrast: **0.098437**.
- Within-botanical-family link-null median: **0.083458** (95% null interval 0.072520–0.092296).
- Excess **+0.014979**, Monte Carlo one-sided **p=0.002** over 499 swaps (local calculation).
- This is conditional evidence that exact butterfly–host identity contributes to the reconstructed resource-geographic pattern **beyond floristic distributions**, but this null DOES NOT preserve butterfly phylogenetic family structure or the geographic amount added for each butterfly.

A more restrictive null that only swaps links within the same **butterfly family AND host-plant family** was frozen at `butterfly_within_consumer_family_host_identity_null_protocol_v0.2.json`. Its GitHub workflow (`butterfly-host-identity-lineage-null.yml`) was queued at this note's creation. If this stricter null fails, the claim must be narrowed to phylogenetically structured host repertoire rather than consumer-specific surplus beyond lineage context.

## Literature and novelty boundary

- Floristic homogenization between distant climate analogues is already established (e.g. https://www.nature.com/articles/s41467-021-27603-y). That pattern alone is NOT a novel discovery.
- Network-mediated biotic homogenization is established in plant–frugivore meta-networks (https://www.nature.com/articles/s41586-020-2640-y).
- Butterfly–host associations carry phylogenetic conservatism and network structure (https://www.nature.com/articles/s41467-018-07677-x).
- Novelty candidate, conditional on stricter tests: **Specific butterfly–host associations channel the floristic effects of plant globalization into uneven potential herbivore resource geography across climatic analogues.**

The project has NOT measured local butterfly occurrence, realized host use, competition, abundance, fitness, or demographic rescue in the newly supplied regions. Continent codes and within-region sampled climates are coarse approximations. Region-pair values are dependent; do not use region-pair counts as independent ecological replicates. Sources and null rules were selected during exploratory work, so minimum Monte Carlo p-values should not be described as independent confirmatory discoveries.

## Stop rule

If the consumer-family-constrained null is weak or its Markov chain does not mix, do **not** promote the original butterfly-specific interpretation; report spatial heterogeneity with the clearer preexisting floristic explanation. Do not rescue results through additional unplanned cutoff sweeps, competition-alpha simulations or temporal-lag variants.
