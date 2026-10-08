# Butterfly-resource climate selectivity: botanical-family audit and ecological pivot

**8 October 2026** · **Exploratory, post-hoc** · `analysis/resource-homogenization-v01` only.

## Summary

The strict butterfly-family × botanical-family conditional null originally produced a modest positive excess in climate-selective *potential* larval-resource geography (observed 0.098437; null median 0.087899; excess +0.010538; one-sided p=0.030, 499 draws; independent seed p=0.026).

A prespecified-for-this-audit leave-one-host-plant-family-out analysis excluded each of the eight highest-link botanical families, which were selected from the frozen HOSTS link counts **before botanical LOO outcomes**. Each excluded family was removed from both the observed and the permuted networks. Species (239), regions (355), butterfly family/host family degree constraints and the same fixed climate pair sets (12,944 analogue and 12,944 discordant pairs) were retained. GitHub run: [37714884267](https://github.com/zuizui0223/chocho/actions/runs/37714884267). No cutoff searches were performed.

| Botanical host family omitted | Original distinct host links | Excess, seed 20261008 | One-sided p | Excess, seed 20261009 | One-sided p |
|---|---:|---:|---:|---:|---:|
| Fabaceae | 573 | +0.005236 | 0.180 | +0.004641 | 0.196 |
| Poaceae | 258 | +0.013005 | 0.010 | +0.011953 | 0.004 |
| Brassicaceae | 200 | +0.011537 | 0.024 | +0.010325 | 0.014 |
| Apocynaceae | 139 | +0.009482 | 0.044 | +0.008414 | 0.042 |
| Asteraceae | 134 | +0.008762 | 0.050 | +0.008080 | 0.050 |
| Apiaceae | 115 | +0.010882 | 0.022 | +0.010131 | 0.020 |
| Passifloraceae | 103 | +0.010113 | 0.040 | +0.009212 | 0.032 |
| Rutaceae | 72 | +0.010550 | 0.028 | +0.009703 | 0.022 |

**Decision:** all botanical-family LOO signs remain positive, but **Fabaceae deletion cuts the excess approximately in half and makes the conditional result non-significant under both seeds**. From the independent butterfly-family LOO audit, deletion of **Pieridae** also reduced the excess and made it non-significant (p=0.116 and p=0.156). These two audits do NOT by themselves identify a Pieridae × Fabaceae causal interaction. They indicate that the broad butterfly-wide statement is too strong.

## Ecology-focused scope

- The frozen 239-species panel contains 43 Pieridae, of which 24 have documented Fabaceae host links (238 distinct links).
- Within these 43 butterflies, 2,798 reconstructed species × region units are new relative to the native resource union. The Fabaceae hosts alone are the unique remaining contemporary known-host support for 1,377 of these introduced-added resource units.
- This is **structural opportunity**, not observed larval survival, local occupancy, competition, population growth, or a benefit of retaining invasive plants.
- Several other Pieridae are Brassicaceae specialists; the fact that Brassicaceae deletion does not eliminate the statistical excess is an interpretable but **post-hoc** contrast, not a preregistered genus-level mechanism.
- The native-versus-contemporary strict decomposition shows a pre-existing native-network selectivity residual around +0.040536 and contemporary residual around +0.050996; the redistribution-related increment is much smaller, roughly +0.0105.

The proper new question is no longer simply “does introduced plant geography homogenize butterfly potential resources?” That phenomenon partly mirrors known plant-flora results (Yang et al., 2021, DOI 10.1038/s41467-021-27603-y). The biological question is **which documented imported Fabaceae hosts are actually used as larval food by Fabaceae-feeding pierids in new regions, and whether feeding corresponds to development and fitness rather than only host listing/oviposition**.

Yoon and Read (2016; DOI 10.1007/s00442-016-3560-2) meta-analyzed 76 studies and reported lower larval performance/survival on exotic hosts relative to native hosts, so the distinction between a documented larval host link and a beneficial habitat is already important and well established. Braga (2023; DOI 10.1016/j.cois.2023.101074) reviews both beneficial and trapping outcomes and the genetic/multitrophic contingencies.

## New independent biology test (launched)

The 24-species Fabaceae Pieridae panel and source criteria are frozen at `PIERIDAE_FABACEAE_LARVAL_USE_PILOT_PROTOCOL_V01.json`. `probe_pieridae_fabaceae_larval_use.py` and GitHub workflow `butterfly-pieridae-fabaceae-larval-use.yml` independently query GloBI for geographically and temporally grounded larval eating observations. Native alternative hosts are retained when classifying introduced-only cells.

**Pilot gate**: >=10 source-grounded introduced-Fabaceae larval feeding observations across >=3 butterfly species and >=3 WGSRPD3 regions, with original occurrence provenance to audit. This is feasibility only: a pass cannot establish colonization, host fitness, or causal range filling. A failure terminates the field-observation shortcut; do not relax stage or provenance thresholds.

## Stop rule and publication boundary

1. Do **not** sell p=0.030 as a generic new butterfly-wide ecological law. The evidence is modest, lineage- and botanical-family-sensitive, and entirely reconstructed from potential host links plus plant distributions.
2. Do **not** claim that host presence entails demographic benefit. The previously audited monarch Puerto Rico case was misclassified as introduced-only because a known native host link was absent from HOSTS.
3. Do **not** merge exploratory results into the submitted GEB manuscript unless a separate, defensible independent biological endpoint is developed.
