# Island ports, introduced flora and butterfly resource opportunity — 2026-11-02

**Status:** exploratory source-audited analysis on `analysis/island-mainland-host-opportunity-20261102`. No change to the frozen GEB manuscript. Numerical associations are not proof of transport, introduction timing, or butterfly establishment.

## The ecological comparison

**Are ports associated with the amount of introduced plant resources, or the *fraction* of butterfly opportunities that depend on introduced plants?** These are different estimands. A naturally host-poor small island can have a high introduced-host *share* with few port facilities and few distinct introduced plants, while a large connected island can support many introduced plants but retain many native host opportunities.

The fixed source population is 126 strict WGSRPD3 island regions (including Canary Is., CNY, restored from an erroneous early mainland classification). For butterfly resource ratios and HOSTS-selected plant taxa, 123 islands have contemporary butterfly resource opportunities. The 2023 World Port Index identifies 982 mapped indexed ports in 91 of those 126 regions; absence of a cataloged port is **not** proof of no shipping or no plant movement.

## Newly audited results

| Metric | Source population and port comparison | Area + Level1 association | Geographic permutation or limit |
|---|---|---|---|
| Relative **introduced-host butterfly resource share** | 123 host-active islands. Mean 54.4% when WPI-indexed port present (91) versus 74.9% when none indexed (32). Unadjusted Spearman ports/share ρ = −0.516. | β = **+0.0102 fractional share points** per unit of log2(1+port count), controlling log area + Level1. Approx +1.02 percentage points. | Within Level1×within-Level1 area half, two-sided p = **0.3646** (4,999 port-label permutations). The adjusted port-share signal is weak, not an established effect. |
| **Distinct introduced larval-host taxa**, deduplicated plant × region | 564 WCVP accepted introduced plant taxa represented among the 239-butterfly HOSTS links; 5,611 plant×island combinations across 123 regions. Mean 53.9 taxa with indexed port vs 21.9 without; unadjusted ρ = +0.393. | β = **+0.1675** for ln(1+taxa) vs log2(1+port count), adjusting log area + Level1 + log native butterfly resource baseline. Approx **18% higher (1+taxa)** per doubling of (1+indexed ports). | Conditional Level1×area-half port-label permutation p = **0.0002** (4,999); sign remains positive excluding Pacific/Antarctic regions, zero native-resource baselines, biggest port region, and restricting to regions with ≥75 contemporary resource cells. |
| **All introduced WCVP accepted plant species**, independent of HOSTS selection | WCVP accepted names census 360,064 species; on 126 islands, 39,024 introduced species×region and 184,997 native species×region records (non-extinct, non-doubtful). Mean 397 introduced species with indexed port (91) versus 83 without (35); unadjusted ρ = +0.628. | β = **+0.2214** for ln(1+all introduced species) vs log2(1+indexed ports), adjusting log area + Level1 + log native flora species richness. Approx **25% higher (1+introduced species)** per doubling of (1+indexed ports). | Conditional Level1×area-half port-label permutation p = **0.0002** (4,999); sign remains positive excluding Pacific/Antarctic regions and highest-port-count island; highest-port exclusion coefficient β = +0.1772, p = 0.0004. |

The all-flora result **tests a different, HOSTS-independent botanical outcome** with the same WCVP geography and the same islands. It is not an independent replication in an entirely new system.

### Why a negative raw ratio can coexist with positive absolute plant numbers

- **Port-poor, naturally host-poor regions** can have extremely high shares of newly enabled butterfly resources simply because the native baseline is near zero.
- **Port-rich larger regions** can contain many more introduced plant species but also retain a broad native host portfolio. Their absolute introduction richness and proportional introduced-resource reliance need not track one another.
- The adjusted port-share effect is small, while the introduced *taxa counts* are positively associated with indexed ports after simple area and broad-biogeographic controls.

This is a potentially interesting **scale/denominator effect** in ecological opportunity, not demonstrated port-mediated ecosystem causation.

## Current three-gate state

| Gate | Source-backed positive observation | Not yet identified |
|---|---|---|
| **Host/resource availability** | 239 butterflies: 26,530 native resource cells vs 41,083 contemporary (global +54.9%). In corrected 126 island regions, 4,829→10,570 (+118.9%) vs mainland 21,215→29,901 (+40.9%). 5,741 island additions supported by 11,573 butterfly×introduced-host×island edges, 564 host taxa. Port-rich regions also have more distinct introduced host taxa (above). | Dated *actual introduction and establishment* of particular host on specific physical islands; plant quality and realized larval feeding. |
| **Climate** | NASA POWER 1981–2025 annual T2M obtained for all 15 literature-case WGSRPD3 regions. Ten pre-first-record annual T2M histories are available for 21/47 literature leads. | Butterfly-specific thermal/degree-day/overwinter limits and timing of actual climatic suitability crossing. |
| **Access / self-dispersal** | 2023-WPI catalog mapped 982 indexed ports to 126 island regions. Frozen LepTraits wingspan information for all 239 GEB species. The case-literature seed includes 47 records (13 species), but only 11 lead rows / 2 species / 7 species×WGSRPD3 cells have host-before-butterfly documentary index plus available climate, port, wing proxies. | Historical ship/plant shipment intensity, butterfly transport mode, native-to-new-island introduction status, reliable arrival and breeding dates. |

**Realized butterfly observation is not established.** In the independent 31-species DOI subsample, only 13/920 newly enabled butterfly×island opportunities have indexed GBIF butterfly records. This **cannot** be interpreted as 98.6% unused. The port-index-only score under within species×Level1 permutations was p(two-sided)=0.2458, wing×port p=0.2181. There is no supported butterfly colonization or vessel replacement-of-flight effect from these sparse observations.

The 47 literature leads also contain nonindependent repeated physical islands within WGSRPD3, and the first indexed plant year is not historical arrival. A causal last-satisfied gate cannot be assigned from these observations.

## Priority next source and model

1. Secure the archived **224/225-butterfly × physical-island × year** observation-risk-set source (original approximately 14.4 million records). This file is not currently available in the analysis runtime/repository and must not be invented. If unavailable, construct a smaller **source-dated documented establishment cohort with explicit absence/recording-effort controls**; do not expand the positive-only 47 reports into artificial failures.
2. Add independent island/year butterfly reporting and plant collection effort, provider/digitization changes and native/butterfly-origin status. Treat nonrecords as *not indexed* until detection is modeled.
3. Calibrate butterfly climate envelopes independently, with cold-season and growing-degree-day variables; use NASA POWER 1981 onward for prospective chronological tests rather than historical first records before 1981.
4. Model species×island×year *indexing or independently verified establishment*, with hosts, climate, port accessibility, wing size and their interactions. Compare port×wing and isolation×wing after sampling adjustment and leave-island/leave-species-out checks. Do not label a last gate without a documented threshold sequence and adequate uncertainty.

## Reproducible evidence

- Ratio response: `results/three-gate-port-resource-faucet/island_port_host_resource_association_audit.json` and `region_port_resource_join.csv`; source code `scripts/test_three_gate_port_resource_island_association.py`.
- Butterfly HOSTS-selected distinct introduced plant taxa: `results/three-gate-port-introduced-host-taxa/introduced_host_taxon_port_association_audit.json` and `distinct_introduced_host_taxa_and_port_by_island.csv`; code `scripts/test_three_gate_port_introduced_host_taxa.py`.
- WCVP *all introduced flora*: `results/three-gate-port-all-flora/all_flora_port_association_audit.json` and `all_flora_intro_native_wpi_islands.csv`; source extraction `scripts/extract_wcvp_all_introduced_flora_islands.R`, model `scripts/test_three_gate_port_all_wcvp_flora.py`; pinned WCVP `65bed76bae9d644ccb6ad200c05f9f5071d89e05`.
- Botanical source, climate, wing data and GBIF records: `docs/three_gates_island_butterflies_20261102.md`.

**Statistical boundary:** Geographic label randomization is a descriptive calibration under restricted exchangeability, not a randomized intervention. WPI counts correlate with population, botanical inventory quality, trade, and land use; neither conditional p nor a log-scale coefficient establishes causal transport or shipping.


## Botanical survey-effort robustness (source-audited)

New outcome-blind sensitivity executed on the same branch with `scripts/test_three_gate_port_flora_botanical_effort.py`. All 126 islands have WCVP flora and WPI port counts, but only **98 islands** have successfully queried, taxonomically consistent, spatially indexed GBIF total botanical record counts (97 from original Q4 response-blind source + the independently indexed Canary Is.). GBIF query failures, complex geometries and no-indexed-record regions are excluded, **not interpreted as zero botanical sampling**.

All regressions use log WGSRPD3 area and fixed Level1 botanical regions. Botanical outcome models also control the appropriate WCVP native-richness baseline, and reporting-effort models include log(1+GBIF indexed plant records) or log pre-1900 records. Port predictor: log2(1+2023 WPI port count). Tail calibration: 4,999 port-label permutations within Level1 × within-Level1 area half. This null is an observational descriptive exchangeability benchmark, not causal inference.

| 98-island measured subset | Port association coefficient | Outcome-scale equivalent per +1 log2(1+ports) | Conditional 2-sided p |
| --- | ---: | ---: | ---: |
| All introduced WCVP accepted species, native baseline / area / region controls but **no GBIF effort** | +0.3337 log(1+species) | +39.6% | 0.0002 |
| Same + **all-years indexed botanical reporting** | **+0.2445** | **+27.7%** | **0.0012** |
| Same + **pre-1900 indexed botanical reporting** instead of all years | +0.2225 | +24.9% | 0.0014 |
| Log introduced/native WCVP species contrast + all-years reporting | +0.2429 | Positive introduced-relative association | 0.0036 |
| HOSTS-selected introduced larval host taxa + native butterfly-resource baseline + reporting | **+0.1639** | **+17.8%** | **0.0110** |
| WPI ports vs **GBIF botanical reporting itself** (native flora / area / region controlled) | +0.5239 | +68.9% indexed botanical records | 0.0002 |

Sensitivity of all-flora effort-adjusted model: ports >0 only, 69 regions (β=+0.2354, p=0.0222); indexed botanical record count ≥1,000, 71 regions (β=+0.2193, p=0.0466); excluding highest-port region, 97 regions (β=+0.2303, p=0.0010); excluding Pacific/Antarctic Level1 regions, 69 regions (β=+0.2428, p=0.0136).

**Interpretation:** WPI-rich islands have more GBIF botanical records, yet the *positive conditional association* between indexed ports and both all introduced plant richness and selected larval-host richness persists when this coarse observation-effort proxy is included. This is a stronger positive observational **resource supply** association than port × observed butterfly occupation. It does **not** demonstrate a mechanism or estimate the arrival route of any plant or butterfly.

**Selection warning:** The excluded 28/126 islands have median area ~25,120 km² and median 6.5 ports vs the analyzed 98 islands' ~956 km² and 1 port. Excluded sites include Great Britain, Japan, Philippines, Borneo and Cuba. GBIF failures and polygon complexity are strongly nonrandom; complete-case results cannot be generalized without these missing locations. Moreover botanical reporting may be downstream of ports and human settlement, so adjusting for records can also remove real shared effects. The exploratory P-values are not spatial-independence corrected, are not corrected for selection/multiple models, and do not identify a causal shipping effect.

Reproducibility:
- `results/three-gate-port-botanical-effort/botanical_effort_adjusted_port_flora_audit.json` (source population, every model and permutation)
- `results/three-gate-port-botanical-effort/botanical_effort_flora_island_inputs.csv` (98 complete-case island observations and WPI/GBIF controls)
- `scripts/test_three_gate_port_flora_botanical_effort.py`, workflow `.github/workflows/chocho-three-gate-port-flora-effort-control.yml` (completed GitHub Actions success).
