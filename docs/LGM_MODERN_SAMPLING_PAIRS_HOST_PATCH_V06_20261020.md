# Original modern insect sampling pairs across hypothetical 0 BP versus 21 ka host-climate patches

**Source-verified exploratory ecology diagnosis — 2026-10-20.** This is **not** observed ancient plant occupancy, ancient insect ranges, realized resource corridors or genetic population breaks. TTF's prior known-host-vs-matched-other-plant **negative** v0.5 climate specificity result remains negative and closed.

## Why this additional study was needed

The 416 native-WCVP modern plant climate models and source-authenticated PMIP4 ICE-6G_C land+ice terrains produced 513 insect potential host-climate landscapes (349 distinct recorded host combinations), often showing reduced suitable area and lower share of suitable area in the largest patch at 21ka. However, those are **global potential graph** summaries, not connectivity of the particular original modern insect sampling sites.

A preliminary response-blind experiment projected all fixed modern sample locations into the two climates using PMIP4 icefree-land constraints and 0.5° climatic host masks, with physical-landmass-constrained nearest suitable patch at radius 50/100/150 km. It showed very different numbers of identifiable sampled edge endpoints in each era. Counting unmatched ancient anchors as 'disconnected' would manufacture a glacial barrier.

The present v0.6 analysis therefore **holds exactly the original same edge identity constant**. An original locality pair enters its component-transition contrast only when both endpoint modern geographic coordinates are contactable under BOTH the 0BP and the 21ka PMIP4 and host climate rules.

## Frozen exact sources and eligibility

- Original TTF 641 insects, 642 source geographic-species entities, 17,202 modern geographic locality rows, 71,728 geographic pair rows, 1,037 accepted insect-host associations. No sequence or genetic-distance field opened.
- For the **513 insects** with all their *recorded* accepted hosts modelled, there are **13,569 original modern localities and 57,252 original geographic pairs**.
- All eligible insect host climate landscapes are formed by OR of the exact original 416 source-verified plant masks on 360×720 0.5° source climate-potential PMIP4 ice-free terrain; source main climatic nonanalog censoring remains fixed. Separate previously prespecified *nonanalog extrapolation sensitivity* is compared, but neither has verified ancient occupancy.
- Geographic anchor must itself fall on epoch-specific PMIP4 4-connected ice-free physical land, and nearest potential host climate patch must lie on the **same terrestrial component** within exactly the prespecified 50, 100 or 150 km great-circle distance. No ocean/ice crossing permitted by unconstrained nearest-neighbor assignment.
- Source script `scripts/audit_common_modern_sampling_edges_v06.py` first recomputed and checked **all 6,156 original per-species/per-era/per-scenario/per-radius** anchor counts from v0.5 *exactly*, then did the within-edge both-epoch intersection. The earlier output was not rewritten.
- The 513 species are not 513 independent host resource landscapes: 349 distinct source host combinations; original geographic edges within insects share endpoints. No confirmatory p-values.

## Result: 100 km, primary climate-nonanalog-censored model

Exactly **19,657 of the 57,252 original pairs** (34.3%) could be assigned to potential host climatic components in BOTH modeled ages; most remaining pairs were unavailable at one or both ages, not proven genetically split.

| Same original geographic pair status | Count |
|---|---:|
| Same potential host patch at 0BP and 21ka | **14,697** |
| Same at 0BP, different at 21ka (*hypothetical split*) | **2,198** |
| Different at 0BP, same at 21ka (*hypothetical merge*) | **1,083** |
| Different in both | **1,679** |
| Total exactly paired | **19,657** |

At the **species** grain, at least one paired edge exists in **385/513 species**; their within-species same-component fraction decreases in **178**, increases in **69**, unchanged in **138**, giving an overall median 0 due to many zero/tiny paired-edge denominators.

For **317 species with >=5 common original pairs** (post-hoc descriptive stability reporting, not qualification rescue), the directional counts become **166 decreases / 60 increases / 91 unchanged**, median within-species difference **−0.0396879** (−3.97 percentage points). These map to **220 shared host combinations**, whose mean within-hostset source-insect changes are negative in 115, positive in 44 and equal in 61 (only a dependence diagnostic, not new independent genetic landscapes).

## Matching radius and nonanalog sensitivity

| Host-climate reference | Radius | Species with ≥5 common edges | Decrease / increase / equal | Median within-species same-component change |
|---|---:|---:|---:|---:|
| Main nonanalog excluded | 50 km | 226 | 75 / 42 / 109 | 0 |
| Main nonanalog excluded | **100 km** | **317** | **166 / 60 / 91** | **−3.97 pp** |
| Main nonanalog excluded | 150 km | 352 | 191 / 73 / 88 | −4.77 pp |
| Novel-climate extrapolation included | 50 km | 234 | 70 / 53 / 111 | 0 |
| Novel-climate extrapolation included | **100 km** | **321** | **166 / 63 / 92** | **−3.23 pp** |
| Novel-climate extrapolation included | 150 km | 354 | 192 / 75 / 87 | −4.77 pp |

The signal is **resolution/radius dependent**, notably weak at 50 km, and selected by whether an original modern edge has suitable-climate contact under both ages. This remains a conditioned geographic experiment, not a valid test of ancestral larval resources or independent genetic breaks.

### Why the primary historic genetic inference remains prohibited

The 513 accepted-host lists are only 1–3 recorded identities per insect and may be incomplete; WCVP current native L3 is broad. The 0BP-fit plant model can extrapolate poorly to 21ka and 0.5° PMIP4 land/ice cells can erase narrow routes. A modern locality itself might be under modeled 21ka ice (5,300/13,569 versus 1,536 in the present 100 km source pairing), but that **does not tell where the ancestral population was during the LGM**. Even a source-host climate patch contact within 100 km on the same physical landmass is **not a continuously host-occupied path from the sampling point**.

The original TTF known-host-vs-climate-matched-alternative-real-plant negative test v0.5 found **no host-specific historical-climate advantage**. These new climatic resource geometry patterns do not overturn that conclusion, nor supply a genetic association. The next genuine genetic test still requires truly independent prespecified panel data and historical resource occupancy/dispersal assumptions.

## Source receipts and final status

[GitHub Actions #39404319303](https://github.com/zuizui0223/chocho/actions/runs/39404319303) **SUCCESS**: 513 insects × 2 source assumptions × 3 distance radii = 3,078 per-species same-edge records. Exact original 6,156 previously audited age-specific records replayed and every edge-status count reconciled to its original denominator.

Source result is stored at `results/lgm_513_original_modern_paired_edges_v06/` in this isolated draft branch; the source SHA-bound archival receipt accompanies the 3,078-row CSV and summary JSON. The response contains **no genotype**, ancient insect range or realized glacial host occupancy.
