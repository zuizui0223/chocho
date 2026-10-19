# Cross-species realized-host evidence gate — 2026-10-19

Status: EXPLORATORY PROSPECTIVE ELIGIBILITY PROTOCOL; no new larval outcome or cross-species effect estimated. This file belongs only to the alternative draft PR #43 branch, not the frozen GEB manuscript PR #38.

## Target claim and estimands
Main biological question: **When does anthropogenic larval-host redistribution increase butterfly recruitment, rather than merely mapped host opportunity?**

Separate three estimands, with distinct denominators and evidence:
1. Geographic opportunity: butterfly × WGSRPD3 native/current host-envelope difference, as previously frozen for 239 species. This is a plant-host distribution construct, not observed butterflies.
2. Host performance: individually assigned eggs/larvae reaching a recorded prespecified life stage under a specified host, site, season and experimental design. At minimum, survival success counts with initial N; adult emergence is preferred, with censoring/diapause and plant/plot-family clustering respected.
3. Choice-mediated total recruitment: butterflies' actual allocation of eggs among coavailable hosts **combined within the same compatible units/cohorts** with stage-specific performance, and a measured/defensible counterfactual. Neither a host survival contrast alone nor an oviposition contrast alone identifies this quantity.

Mathematical boundary for fixed eggs and native-host substitution, per season t:
delta B = sum_t E_t q_t (w_I,t - w_N,t).
The sign-reversal hypothesis requires independently observed egg diversion q>0 and poorer performance w_I<w_N over sufficient demand-weighted windows. With q=0 and no habitat/parasite/source-sink feedback, merely adding a bad host does not reduce recruitment. Treat this equation as a falsifiable mechanism *definition*, not a fitted result.

## Response-blind public dataset discovery (not evidence that the tests succeeded)
Metadata-verified candidate sources, including diversity and limitations:

| ID | Public source | Cohorts / response advertised | Possible valid estimate | Key limit |
|---|---|---|---|---|
| D1 | Ferreira & Rodrigues (2022), Dryad doi:10.5061/dryad.3bk3j9kk2 | Two Neotropical Danaini butterfly species; field_performance.txt, lab_performance.txt; native/exotic Apocynaceae, exposed/bagged | Within designed lab host comparison, and field host × exposure descriptive association | Field host species each occurs at a different site; site and host confounded. Do not call field host effect causal, and don't claim egg-diversion |
| D2 | Saldivar & Wilson-Rankin (2024), Dryad doi:10.5061/dryad.wdbrv15v5 | Vanessa cardui; larv.plants.csv, larv.switched.csv; larval survival and pupation on five potential hosts | Assigned-larva host performance within relevant rearing provenance strata | Wild-vs-commercial origin differs; do not infer native/exotic contrast or recruitment without verified host range and cohort |
| D3 | Glanville fritillary host×temperature×family study, Dryad doi:10.5061/dryad.kwh70rz2m | Larval survival/growth/fat; 2 hosts × temperature with family IDs | Within-family host × temperature reaction norm, if raw family/experimental units exist | No documented egg allocation; not direct proof of introduced-host geography |
| D4 | Jones & Agrawal (2019), Dryad doi:10.5061/dryad.8hd6764 | Danaus plexippus; caterpillar_performance.xlsx and competition_oviposition.xlsx | Parallel measures of host-specific choice vs growth with plant/treatment units | Different experiments cannot be multiplied as if individual cohorts were linked |
| D5 | Narango & Straley (2025), Dryad doi:10.5061/dryad.x95x69pw7 | One larval Lepidopteran on local/nonlocal/non-native congeneric trees; survival_format_long.csv | Individual survival contrast within valid randomized/blocked experimental design | Verify focal taxon is a butterfly (not assume all Lepidoptera are butterflies) before including in main panel |

Sources are discovery leads only; actual file rows, identifiers, licensing, and analysis-ready comparability are NOT yet audited. No pooled meta-analytic estimate is claimed.

## Fail-closed extraction for each dataset
Record source DOI, version, file name, checksum, taxonomic identity (butterfly vs moth), geography, host name and local nativity, design and randomization, original initial N, outcome stage, censoring/diapause, monitoring dates, repeated plant/plot/family IDs, and source-defined covariates. Use statuses:
- ELIGIBLE_HOST_PERFORMANCE — independent butterfly host × life-stage counts with valid denominators and clustering.
- ELIGIBLE_CHOICE_AND_PERFORMANCE — additionally, compatible choice/exposure and cohort performance on coavailable native and introduced hosts.
- GEOGRAPHY_LINK_ELIGIBLE — additionally, species×locality host use and local native fallback validation allow comparison with mapped A/R/U.
- NOT_IDENTIFIABLE — absent or confounded host comparison, missing initial N, only same-day egg/larva snapshot, or source inaccessible.

Do not convert field disappearance into deaths; do not classify live diapause as failed adult production; do not use larvae as independent temperature or cage treatment replicates. Do not assume each host is locally native/exotic from genus name. Check if host range is defined regionally and historically; WCVP current status is not a historical introduction date.

## Sequential tests and stop rules
T0 source audit: independently verify actual rows and units against source paper. If less than two independent butterfly taxa have comparable response-stage data, report descriptive cases only; no across-butterfly law.
T1 within-study performance: compute host contrast on **consistent, explicitly declared life-stage endpoint** and report species/host/site/season-level estimates and uncertainty. Do not combine pupation, survival at day 7, and adult emergence as one endpoint. If host and site are perfectly confounded, label association only.
T2 choice–performance: require matching host-availability and egg-allocation measurements, and prospective follow-up of assigned or tagged cohorts; if absent, sign reversal in *total production* is NOT_IDENTIFIABLE even when host survival differs.
T3 cross-study heterogeneity (only when T0–T2 permit): compare within-study effects at the study/site and butterfly-species levels; cluster by study; disallow many host-pair rows from pretending to be many independent butterfly species. Include negative controls for sampling effort, life stage, insect origin, and plant/cage experiment design. Publish failures.
T4 map linkage (only with geographically verified local-host access): distinguish originally native supported A from genuinely novel U using both source catalogues, and recheck source identity turnover. Do not transport globally measured performance from one plant population to all 365 macro regions.

Primary falsification: if the integrated host-quality + seasonal-diversion model does not improve **held-out, independently observed juvenile/adult cohort production** relative to host-presence-only baseline, do not claim the proposed general principle. Any apparent strong effect driven only by host/site confounding fails.
  
## Immediate authorized next computational step
Acquire permitted public files for D1 and D2 first (a two-butterfly field/lab design plus an independent butterfly developmental study). Audit field names, original larvae denominators and clustering; output a provenance/eligibility JSON with one row per study×experiment×species×host×stage, plus an explicit impossibility matrix for causal native/exotic host attribution. Expand to D3–D5 only after the first gate; do **not** retrofit hypothesis or select datasets by statistical sign. This protocol contains zero newly estimated ecological effects.

## Integrity
The existing 239-species HOSTS/WCVP global manuscript and PR #38 remain unchanged. Inference must distinguish source-proven observations from mathematical examples and literature prior art; no new paper or conservation policy is implied by metadata verification alone.
