# Monarch host identity vs larval quality: biological evidence checkpoint v0.1

**Date:** 2026-10-08. **Status:** EXPLORATORY SOURCE-CHECKED, NOT A NEW GENERAL ECOLOGICAL RESULT. **Scope:** `analysis/resource-homogenization-v01` only. Main GEB manuscript unchanged.

## Scientific question

Does anthropogenic host-plant redistribution generate **developmentally viable larval-resource geography**, or merely extend a map of documented putative host associations? The difference between plants that are listed as hosts and plants that support actual emergence is the ecological issue. Merely removing links and obtaining fewer mapped regions is algebraically guaranteed and is *not* a discovery.

## Frozen source inventory and reproduced existing data

Pinned HOSTS/WCVP accepted sidecars: GitHub run **37628559686**, artifact `temporal-hindcast-sidecars`, sourced to pinned HOSTS/WCVP commits. `Danaus plexippus` has **42 accepted interaction CSV rows, only 40 distinct accepted plant IDs**. The original deduplicated native resource union covers **181 WGSRPD3 regions**, contemporary union **263**, introduced-only union **82**. Duplicated interactions do not imply duplicated area and were deduplicated for this audit.

All-in-one `monarch_quality_host` source code remains in exploratory work; the same fixed sidecar yields:

- **Missing verified, highly viable known host** `Asclepias curassavica`: the exact accepted host species link for monarchs is missing from the frozen HOSTS accepted interaction sidecar. Greenstein et al. (2022) documented **17/20 (85%)** neonates completing development to adult on this species in a no-choice laboratory experiment. The source paper classed it high performance. This laboratory result does NOT prove success in every region.
- **Low-performance listed host** `Araujia sericifera`: present in the accepted monarch HOSTS list (accepted ID **500657**); in the independent experiment only **2/127 (1.6%)** larvae completed adulthood. The pinned plant map lists it in 5 native and 29 contemporary WGSRPD3 regions. Removing *only* this host link changes the nominal introduced-only union from **82 to 81**, entirely through **CLC (Caroline Islands)**. It does not establish actual field mortality in CLC.
- **Disputed host link** `Gossypium arboreum` (accepted ID **2830992**): source paper warns that some `Gossypium` monarch host claims may arise from mistranslation or vernacular misinterpretation; it does NOT independently prove this exact species unsuitable. Excluding the dubious link as a *sensitivity* reduces introduced-only regions **82→79** (current **263→260**), unique lost codes **GAB, TAI, TOG**. This is a counterfactual, not a reclassification.
- **Native host omitted from reconstructed Puerto Rico butterfly resources**: the frozen model has **PUE native=0, contemporary=1**, so flags Puerto Rico as `introduced-only`. Independent Kew POWO treatment says `Asclepias curassavica` is **native to Puerto Rico** and introduced to Florida. Whether the **identical pinned WCVP version** moves the PUE flag to native is being computed by the dedicated source-pinned R audit. Do not pre-report its result.

**Literature:** Greenstein, Steele & Taylor (2022), *PLOS ONE* 17 e0269701, DOI https://doi.org/10.1371/journal.pone.0269701 . Source species-specific independent quality categories in **S1 Appendix B Table**, which must be retrieved and independently crosswalked before any larger quality-stratified geography claim. Detailed study includes 127 potential hosts (34 H1–H3, 42 L1–L3, 33 N, 18 U). Note that H1/H2 may not mean quantitative survival; require H3 for the strongest evidence and do not impute unmeasured species. Kew POWO `A. curassavica`: https://powo.science.kew.org/taxon/94213-1 .

## Three frozen source-dependent follow-ups

1. **Exact four-scenario spatial reclassification**: protocol `MONARCH_HOST_QUALITY_GEOGRAPHIC_SENSITIVITY_PROTOCOL_V01.json`, script `audit_monarch_missing_high_quality_host_geography.R`, workflow `butterfly-monarch-missing-high-host-geography.yml`, run **37736645743**. Fixed WCVP adds native/current accepted `A. curassavica`, compares baseline, independent high-quality-host addition, disputed cotton exclusion, and both. Need result artifact before numeric statements. If pinned data classify PUE differently from POWO, flag mismatch instead of silently substituting sources.

2. **Direct published host-quality classification table**: protocol `MONARCH_EXTERNAL_HOST_QUALITY_SOURCE_GATE_V01.json`; `probe_monarch_host_quality_original_source.py`, workflow `butterfly-monarch-quality-source-gate.yml`, run **37736225347**. Retrieves original PLOS S1 DOCX and Dryad metadata; no unverified genus matches.

3. **Independent urban stage dataset**: protocol `URBAN_MILKWEED_STAGE_SOURCE_GATE_V01.json`, `probe_urban_milkweed_figshare_source.py`, workflow `butterfly-urban-milkweed-figshare-source.yml`, run **37735892534**. Source *Erickson, Schultz & Crone* (2025) DOI https://doi.org/10.1002/ecs2.70259, Figshare https://doi.org/10.6084/m9.figshare.25648644.v1; checks whether matched 2022 route/month native+introduced milkweed counts and distinct egg/larva counts/survey effort support a stage×alien composition test. Even with data, stage counts are **not** direct larval survival because stages are not a longitudinal cohort.

All three were **queued** at the latest audit when this note was written. They are not claimed completed. Pushing documentation alone should not cancel them; do not rerun on every commit.

## Current defensible interpretation

- The original global result measures **putative known-host presence**, not performance-weighted realized resource use.
- Taxonomic host inventory can miss a verified high-development species while including low-development or disputed species. This supports a substantial *biological measurement warning*; no global bias direction is known without complete linkage and performance coverage.
- In particular, `PUE` is a priority consistency audit: if the high-quality native host restores native resources there, a region previously labelled entirely alien-dependent was an inventory artifact.
- A single butterfly and three debated plant links cannot establish a generalized ecological mechanism or independent high-impact paper.
- The future non-tautological test is whether **independent host performance** is associated with the *geographic identity and distribution* of human-added host resources **conditional on host geography, family, knowledge coverage, and introduction status**. It is **not** whether deleting low-quality rows mechanically makes the resource envelope smaller.
- The main GEB manuscript has NOT been modified.
