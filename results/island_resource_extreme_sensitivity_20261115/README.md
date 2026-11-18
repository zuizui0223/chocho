# Are disproportionate island butterfly host-resource gains driven by zero-baseline and tiny islands?

**15 November 2026 | Original 239-butterfly × WGSRPD3 regional source; successful source-locked recomputation.**

## Result

**No: the positive relative island resource-opportunity expansion persists after excluding all regions with zero native-host-covered butterfly combinations, the Antarctic biogeographic Level-1 group, or small WGSRPD3 polygons.** However, its magnitude falls as native-host-richer island regions are selected. This identifies **baseline native-host scarcity** as a structural contributor, rather than claiming unusual island host plants were preferentially introduced.

This study reconstructs **potential regional occurrence of known larval host plants**, not confirmed butterfly island occupancy, climatic suitability, feeding, breeding or colonization.

| Island geographic sensitivity | Island botanical regions | Native host-covered butterfly×region combinations | Newly added from introduced known hosts | Relative increase over native base |
|---|---:|---:|---:|---:|
| **Original strict island set** | **126** | **4,829** | **5,741** | **+118.89%** |
| Exclude all zero-native-resource island regions | **113** | 4,829 | 5,309 | **+109.94%** |
| Require ≥5 native-covered butterfly combinations in an island region | **102** | 4,808 | 4,743 | **+98.65%** |
| Require ≥10 native-covered butterfly combinations | **90** | 4,729 | 3,840 | **+81.20%** |
| Require ≥20 native-covered butterfly combinations | **71** | 4,468 | 2,860 | **+64.01%** |
| Require ≥30 native-covered butterfly combinations | **59** | 4,164 | 2,192 | **+52.64%** |
| Exclude Antarctic broad biogeographic Level 1 (level1=9) | **115** | 4,825 | 5,305 | **+109.95%** |
| WGSRPD3 whole-region area ≥1,000 km² | **67** | 4,049 | 3,275 | **+80.88%** |
| WGSRPD3 whole-region area ≥10,000 km² | **38** | 2,263 | 1,988 | **+87.85%** |

Original continental reference, unchanged: **237 mainland botanical regions**, 21,215 native host-covered combinations and 8,686 added, **+40.94%** relative gain. **This is reference context, not a separately area-, latitude- and native-flora-matched mainland control for each filtered island subset.**

The **13** original zero-native-coverage islands contribute **432/5,741 (7.52%)** of the observed island added-resource combinations. Excluding them still leaves **5,309** new resource combinations above a 4,829 native-resource baseline. The Antarctic Level1=9 island regions contribute 436 added combinations, but in **11** regions only four native species×region opportunities, and are extremely denominator-sensitive; the exclusion still preserves the high island relative gain.

These are **overlapping descriptive filters** of the same original 126 island regions, not independent inferential replicates or treatments.

## Physical interpretation and cautions

- The question is whether introduced known larval hosts add *candidate resources* at the geographic scale of WGSRPD3 regions. Such regions can aggregate several islands or mixed habitats and are not exact ecological island areas or sampled populations.
- A climate-inappropriate or genuinely butterfly-free subantarctic island can still have a mapped host plant under this binary definition; its presence does not create a documented butterfly breeding site. The all-source gain should **never** be called an established increase in butterfly range.
- The all-zero-native group is a valid extreme of the original resource definition, but its astronomical region-level percentage gains are undefined with a zero baseline; the pooled 118.9% estimate remains defined across the original 126 islands.
- After selecting native-resource-rich islands, the island gap-fill **fraction** remains roughly 22–25%, whereas the **relative gain** falls sharply because the baseline is larger. This is a genuine denominator mechanism, not an external causal plant-spread experiment.
- The independent paired-239-species analysis shows that **141/239** butterflies have greater native-gap filling on the mainland, **65** on islands and **33** ties; the island 118.9% large gain is therefore not evidence that island introduction preferentially repairs a higher share of native-host deficits.
- No variable here represents an actual port route, adult wind-assisted movement, ocean drift, an immature butterfly on a named host, or successful island colonization.

## Source and reproducibility

[Source-locked analysis and full island region CSV](https://github.com/zuizui0223/chocho/actions/runs/43669063109/artifacts/14142799905) from [successful CI](https://github.com/zuizui0223/chocho/actions/runs/43669063109). The script `scripts/audit_island_resource_extreme_regions_20261115.py` consumes the original exact WGSRPD3 [region metrics](../island-mainland-resource-q1/region_metrics.csv), [geodesic polygon areas](../island_resource_gap_filling/areas.csv) and frozen [pooled resource audit](../island_resource_gap_filling/pooled.csv). It independently asserts 369 original geographic region rows, strict **126 island / 237 mainland / 6 mixed**, full source-native and introduced-resource totals, and the 13 zero-native regions. The output contains original input SHA256s and a complete 126-island status CSV.

**Manuscript inference:** Keep the +118.9% source positive, but explicitly qualify it as potential regional host-resource geography and note robustness to zero-native, Antarctic and small-area filters. Keep port-host geography and separately estimated Circuitscape community association in their own scientific sections; do not convert a host-resource result into butterfly arrivals.
