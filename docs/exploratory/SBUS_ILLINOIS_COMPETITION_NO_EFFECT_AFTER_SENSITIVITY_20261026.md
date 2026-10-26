# Illinois actual butterfly adult-count competition claim: CLOSED after three checks (2026-10-26)

**Scientific status: no source-independent evidence that species sharing recorded larval hosts competitively suppress one another's next-year adult abundance.** This is *failure to identify competition*, not evidence that competition is absent. No change to GEB `main`, PR #38 or author-submission PR #70.

## Primary source
The verified 2025 U.S. monitoring data by Henry et al., DOI 10.1038/s41597-025-05513-8, Figshare v1 DOI 10.6084/m9.figshare.27934602.v1 (MD5 3161b1c8802a0d67771b2044e18da928, CC BY 4.0), provide **reported adult butterfly counts**. We used 2000–2020 Illinois standardized transects, not GBIF reporting frequency. The host overlap network was pinned to original accepted HOSTS (commit 808e0b869f9ec1adf8efff87cf6a395adda103e0), WCVP (65bed76bae9d644ccb6ad200c05f9f5071d89e05), original `chocho` 239 butterfly names. Host sharing is *documented at species level, not at Illinois sample sites*.

Eligible fixed original Illinois panel: **121 sites, 851 site×years, 47 butterflies, 91 host-sharing butterfly pairs**, **9,178 focal butterfly×site annual transitions**, 1,718 species×site groups. Focal response: year-to-year change in log(1 + adult count per monitored hour). Site×species and calendar-year FE, own previous count, other host-sharing butterfly prior counts, comparator non-host-sharing prior counts, effort difference, original within-site-year community mean. Site clustered nominal normal confidence intervals.

Original **peer-sharing** coefficient: **+0.02194** (SE 0.01717; 95% interval **−0.0117 to +0.0556**). It was **not negative**. Non-sharing counterpart **+0.30545** (SE 0.12678). Shared-minus-nonshared contrast **−0.28351** (95% normal site-cluster interval **−0.5132 to −0.0538**). This is a *relative* coefficient contrast between ecological trait proxies, not an estimate of competitive mortality.

## Post-result fragility (different assumptions, no parameter rescue)

| Fixed source/response sensitivity | N focal annual transitions | Shared peer lag β | Shared-minus-nonshared β | Site-cluster 95% interval for contrast |
|---|---:|---:|---:|---|
| Original all eligible, with site-year community mean | 9,178 | +0.02194 | −0.28351 | [−0.5132, −0.0538] |
| Same panel without community mean (mechanically related to peer exposure) | 9,178 | **−0.00148** | **−0.08250** | [−0.15049, −0.01450] |
| Only focal species positively recorded in BOTH consecutive years; original controls | 3,404 | +0.02133 | **−0.66321** | [−1.56977, +0.24335] |
| Same recorded-positive focal subset, no community mean | 3,404 | **−0.01997** | −0.25462 | [−0.41097, −0.09827] |

This is NOT a confirmatory multi-model search: the last three checks were explicitly fixed **after seeing** the first negative coefficient contrast to diagnose the zero-coding and compositional-control vulnerabilities. Importantly, restriction to positive in both years removes attrition/non-detection and creates survivorship selection; it cannot prove the effect in actual populations.

Primary input has **7,501 implied site-year butterfly absence zeros out of 16,454 nominal species×site-year records (45.6%)**. Source field programmes mostly report positives rather than comprehensive explicit absences, so missing as zero is a major unresolved measurement assumption. Of 9,178 focal changes, **3,404 (37.1%)** remain if both years require explicit positive focal observation. Other-species terms in this sensitivity still use implied zeros. Michigan (43 qualifying sites) and Colorado (34) do **not** pass the independently defined ≥50-site criterion, so no external programme replicate is reported as passing.

## Stronger network negative control: 199 independent original-host degree-constrained swaps

Original prior observed contrast: **−0.283507**. We randomized recorded butterfly–accepted plant links within *butterfly family × plant family* strata while retaining each butterfly's botanical-family diet/bipartite degree and each host's consumer degree within butterfly family; every original adult site/year record was unchanged. Across 199 independently started networks (5×E swap proposals; E=2,596 links), all 199 numerical and swap gates passed, with ~59.6–64.2% of viable edges moved. Observed contrast vs randomized:
- **2.5% null quantile:** **−0.374856**
- **null median:** **−0.067694**
- **97.5% null quantile:** **+0.185301**
- **descriptive empirical left tail:** **0.110**.

The observed negative point estimate is **inside the 95% structural null interval**. Consequently it does not demonstrate an exceptional ecological signal after accounting for HOSTS network degrees and broad family strata.

**Caveat on this null itself:** the number of included focal species×site transitions varied from **7,070 to 10,423** across networks (the original had 9,178). The predeclared maximum ±25% criterion passed, but the null **did not hold the regression's pair inclusion panel perfectly fixed**. It is a useful falsification of the original claim, not an ideal exact-composition randomization. No claim is accepted based on this null.

## Biological conclusion

1. **Direct competition:** no. Sharing a documented host does not yield a robust negative coefficient of co-consumer abundance; observational source and network structure remain confounded.
2. **Adult population rise caused by introduced plants:** not identified. There is no year of host-plant introduction or local plant-biomass time series in this monitoring source; WCVP global region codes are not causal local exposures.
3. **Co-occurrence/resource filtering:** still a distinct question, not answered by this lagged density test. The source has site/year adults, so a separate fixed-pair co-detection analysis has been **predefined prospectively for that analysis** in [protocol](../protocols/SBUS_ILLINOIS_COOCCURRENCE_HOST_DEGREE_NULL_20261026.md), with a *fixed survey pair sample* under all host-network nulls. Do not use it to 'rescue' the failed competition coefficient or retroactively rewrite this negative finding.

## Saved original-data receipts and code
- [Original lagged adult abundance estimate](../results/sbus_illinois_host_niche_abundance_pilot_20261025.json)
- [Community-control/explicit-positive sensitivity](../results/sbus_illinois_zero_control_fragility_20261025.json)
- [199 host-network structural null](../results/sbus_illinois_host_bipartite_degree_null_20261025.json)
- [Successful original lagged data replay](https://github.com/zuizui0223/chocho/actions/runs/40112937718)
- [Successful post-result fragility replay](https://github.com/zuizui0223/chocho/actions/runs/40113985621)
- [Successful bipartite host degree null](https://github.com/zuizui0223/chocho/actions/runs/40114268631)

No public site IDs, coordinates or individual adult survey rows were copied into this repository. Attribution: Henry E.H. et al. (2025) data © authors, CC BY 4.0, original contributing US monitoring programmes.
