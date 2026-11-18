# Port × native/introduced island-flora geography: leave-region-out and port-only source audit

**18 November 2026 | Source-locked 126 WGSRPD3 botanical island regions; one integrated GEB PR #179**

## Source and exact estimand

Plant-richness rows: `results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv` (Git blob `53b618fa51a1a418ccc6e4363f4c866c389c90cc`), with 184,997 native and 39,024 introduced accepted species×region records. Controls shared across native and introduced responses: `log2(1+modern WPI port count)`, log island botanical-region area and full TDWG Level1 dummy set. Separate responses are `log1p(native richness)`, `log1p(introduced richness)` and their **same-island difference**. HC3 normal-approximation 95% intervals describe heteroscedastic regional noise; they are not spatial/phylogenetic cluster intervals and not causal effects. 

Source-locked code: [matched origin and sensitivity reproducer](../../scripts/verify_island_port_native_alien_matched_20261118.py). [Successful source-matched CI on the exact code](https://github.com/zuizui0223/chocho/actions/runs/44288043705) (its JSON artifact includes all coefficient/CI rows). **No model in this block uses measured historical ship traffic.**

## Baseline and stronger geographical checks

| Source-defined test | Native β | Introduced β | Introduced − native β |
|---|---:|---:|---:|
| All 126 islands | **+0.0954** (HC3 [−0.0774, +0.2683]) | **+0.2854** ([+0.0622, +0.5087]) | **+0.1900** ([−0.0187, +0.3987]) |
| Only islands with at least one source WPI port (91 regions) | **+0.0844** ([−0.1475, +0.3163]) | **+0.3110** ([+0.0416, +0.5804]) | **+0.2266** ([−0.0117, +0.4649]) |
| Delete each one TDWG Level1 (nine overlapping refits): **minimum to maximum point β** | **+0.0601 to +0.1257** | **+0.2016 to +0.3797** | **+0.1300 to +0.2676** |

### Leave-one-TDWG-Level1-out: nine actual coefficient re-fits

| Omitted Level1 | Remaining regions | Native β | Introduced β | Introduced − native β |
|---|---:|---:|---:|---:|
| 1 | 116 | +0.0716 | +0.2016 | +0.1300 |
| 2 | 108 | +0.1257 | +0.2620 | +0.1363 |
| 3 | 116 | +0.1155 | +0.3106 | +0.1951 |
| 4 | 108 | +0.1121 | +0.3797 | +0.2676 |
| 5 | 119 | +0.0885 | +0.2634 | +0.1749 |
| 6 | 98 | +0.1234 | +0.2892 | +0.1657 |
| 7 | 121 | +0.0933 | +0.3033 | +0.2100 |
| 8 | 107 | +0.0601 | +0.2920 | +0.2319 |
| 9 | 115 | +0.0763 | +0.2864 | +0.2101 |

All nine geographic leaveouts retain a **positive port–alien coefficient**, and the introduced-minus-native coefficient remains *numerically* positive. However these nine overlapping tests are **not nine independent replications**: most retain 98–121 of the same 126 regions. Moreover the original and port-positive **paired cross-class HC3 intervals include zero**; a strict native-zero or alien-only-port mechanism remains unproven. Deleting Level1=4 makes the paired difference interval positive ([+.0471,+.4881]); this is a post-hoc selected source deletion and cannot replace the source-wide inference.

## Scientific role within one paper

**Supported:** conditional geographical association between contemporary indexed ports and whole introduced plant richness is not driven by one botanical continent/region group or solely the distinction between port present and absent islands. Observational native-flora slopes are weaker and poorly distinguished from zero, while differences between the origin classes remain uncertain under the full crossed model.

**Not supported:** that ports transported the specific plants, that ports have no relationship to native plant richness, that different transport routes are established, that ocean-current ancestry explains all island native floras, or that port presence led to successful butterfly establishment. For paired island-flora *similarity* the marine **shortest-water distance** positively predicts both origin classes; the distinct 104-island **electrical circuit** contrast cannot be substituted for an observed vessel or current flow.

**Future statistical gate:** a spatial-block, species-restricted and genuinely historical port-network test of the native/alien coefficient *difference*, which must not be labelled identified from the present WCVP/WPI aggregate data.
