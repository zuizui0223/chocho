# Geographic host-resource overlap decomposition (post-result)

**Date:** 2026-10-11  
**Status:** explicitly post-hoc extension after seeing the initial 199-swap-null result, not preregistered.

## Motivation

The plant-family × WGSRPD1 constrained introduction null (Actions run 38112016742) preserves host identity, each host's introduced-region count within WGSRPD1, and the number of same-family introduced plants per region. It found that actual introduced-host additions generate fewer novel butterfly × WGSRPD3 resource units than expected: 13,529 observed vs 14,445 null median in the adequate 191-species panel (199 simulations; lower-tail Monte Carlo p=0.005). It did **not** support an unusual single-host fraction: 57.16% observed vs 57.59% null median; 80 host-rich butterflies 45.36% vs 45.38%.

A shortfall of novel geography does **not** identify its cause. This extension splits a fixed total of documented introduced-host-by-butterfly-region incidences, `T`, into three mutually exclusive terms:

```
T = B + D + U

B = introductions whose region already lies in the butterfly's native host union
D = excess introduced-host incidences in regions outside the native union
    (when two introduced known hosts occupy the same new region, one counts as
     a novel geography unit and the other contributes +1 excess overlap)
U = unique butterfly × WGSRPD3 regions newly covered beyond the native host union
```

Because the null preserves each host's introduced-region total and the butterfly's host identities, `T` is **fixed exactly** for each butterfly and each panel. Any deficit of `U` in observed geography must be an excess in `B + D`.

## Distinguishing two mechanisms

- **Native-resource saturation:** observed `B` above null; newly introduced known hosts disproportionately arrive in regions already supported by at least one known native host of that butterfly.
- **Introduced-host coincidence:** observed `D` above null; introduced known hosts overlap one another in the same region even where no native host is mapped.
- Both may matter; one may offset the other.

The comparison uses the **same 199 sampled null matrices and seed** as the original analysis (not an independent replication). Monte Carlo p values are exploratory and may be sensitive to null-chain dependence. The existing tests assert the accounting identity for every randomized matrix.

## Claims explicitly excluded

Neither `B` nor `D` measures local butterfly population use, quality, herbivory, host biomass, introduction timing, community coexistence, resource competition, or selection. Geographic coincidence can arise from human introduction pathways, climatic sorting, historical trade, correlated recording coverage or spatial dependence. WGSRPD3 regions are coarse.

The GEB PR #38 scientific estimates and submission contents remain unchanged; this remains on isolated `analysis/resource-redundancy-20261011`.
