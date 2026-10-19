# Public-source eligibility audit — 2026-10-19

Status: METADATA_VERIFIED; RAW_ROWS_NOT_RETRIEVED; NO_NEW_EFFECT_ESTIMATED.
Scope: multi-butterfly public cohort evidence, separate from frozen main GEB PR #38 and from the prior-art reanalyses in PRs #39–44.
Protocol: provenance/reviewer_defenses/protocols/CROSS_SPECIES_REALIZED_HOST_EVIDENCE_GATE_20261019.md.

## Source 1 — Ferreira & Rodrigues / Danaus larvae (two species)
- Source: https://datadryad.org/dataset/doi%3A10.5061/dryad.3bk3j9kk2 ; dataset published 2022-12-20; field_performance.txt, lab_performance.txt and README.txt.
- Publisher abstract explicitly reports both Danaini taxa developed more slowly on exotic host under laboratory conditions, with no detectable larval survival difference in the original study.
- Field first-instar experiment uses bagged versus unbagged exposures on native and exotic Apocynaceae. Reported exposure×exotic-host association with reduced survival is **not** a clean causal host-origin interaction because **each host occurs at a different field site**. Site, leaf microenvironment and host identity co-vary; bagging may itself change local microclimate. Do not treat butterflies or host plants as independent replicate sites.
- Current status: METADATA_VERIFIED / RAW_ROWS_NOT_RETRIEVED / FIELD_HOST_EFFECT_CONFOUNDED.
- Testable from raw data if retrieved: laboratory stage-matched time-to-development and survival within species, field exposure effect within each site/host; NOT cross-site field host-origin causality, not oviposition diversion or long-term recruitment.
- Missing checks: original assigned N, fate codes, censoring, laboratory treatment allocation, plant/plot and maternal family IDs.

## Source 2 — Saldivar & Wilson-Rankin / Vanessa cardui
- Source: https://datadryad.org/dataset/doi%3A10.5061/dryad.wdbrv15v5 ; updated 2024-09-20; larv.plants.csv, larv.switched.csv, README.md.
- Published abstract reports only Malacothamnus fasciculatus and Sphaeralcea ambigua among five candidate plant taxa supported successful pupation of commercially sourced larvae in the tested conditions; and 63% higher survival of wild-collected larvae versus commercial larvae. **These are authors' prior-art published claims, not independently recalculated by chocho.**
- Explicit raw variables include status, source (oviposition vs suppliers), trial, plant species, pot, mean temp/humidity, emergence-in-mail, time to death. Hence provenance/source and complete follow-up are mandatory before pooling.
- Current status: METADATA_VERIFIED / RAW_ROWS_NOT_RETRIEVED / CANDIDATE_HOST_PERFORMANCE_TEST / NOT EVIDENCE FOR INTRODUCED-HOST GEOGRAPHY OR CHOICE.
- Critical counterfactual: if a reported host is unable to complete pupation for commercially reared larvae, this does not prove field populations cannot use it; host-origin/source and acclimation can confound.

## Source 3 — Verspagen et al. / Melitaea cinxia
- Source: https://datadryad.org/dataset/doi%3A10.5061/dryad.kwh70rz2m ; published 2020-11-18; Clutch_data.csv, Individual_data.csv, Survival_data.csv plus three README files.
- Original published abstract reports temperature×host and family-specific larval survival/performance, showing within-species reaction norms are plausible. **Previously published result only**.
- Current status: METADATA_VERIFIED / RAW_ROWS_NOT_RETRIEVED / HOST×TEMP×FAMILY_CANDIDATE / NOT AN INTRODUCTION-ORIGIN EXPERIMENT.
- Testable with raw data: factorial survival response with maternal family hierarchy and temperature unit structure. Requires checking if larval rows correspond to independent treated units and if actual egg allocation measured.

## Across studies: identification-first conclusions
1. All 3 sources contain some documented host-performance related outcomes, but endpoint comparability is not verified without raw rows.
2. None of these metadata listings documents a within-site natural/native versus introduced host *egg diversion linked to full adult recruitment*, so the sign-reversal hypothesis currently remains NOT_IDENTIFIABLE.
3. A *published* biological pattern is already relevant: host origins alone are poor proxies for life-stage performance; lab vs field, temperature, larval provenance and family matter. **This is a prior-art contextual synthesis, not novel replicated independent discovery.**
4. Do not claim a three-species or four-species pooled estimate: Source 1 contains two species but its field native/exotic contrast is confounded, and sources 2 and 3 estimate different endpoints.
5. Data download through this session's external file channel was not successful. The exact CSV rows, sample sizes and distributions were not independently inspected; SOURCE DATA AUDIT INCOMPLETE.

## Next PASS/FAIL rule
To upgrade this audit, first access raw Dryad files lawfully and inspect source README, exact file SHA256, unit identifiers, assigned N and biological endpoint. Extract *one source-specific* table per species×host×treatment×site×season×life-stage, preserving actual source/protocol denominators and alive/diapause censoring. Only thereafter estimate uncertainty; no invented outcomes.

## Decision
No promotion of realization nonmonotonicity to supported global ecology; no merge into PR #38/main. Maintain branch-isolated exploratory status. A false positive from mixing site-with-host or field survival with laboratory development would be more misleading than reporting NOT_IDENTIFIABLE.
