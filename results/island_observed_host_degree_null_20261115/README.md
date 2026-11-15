# Independent DOI island butterfly observations vs introduced-only larval-resource opportunities

**15 November 2026 | New source-locked fixed-margins null, n=31 GBIF-recordable butterflies across 126 WGSRPD3 island regions.**

## Decision

The reconstructed global HOSTS×WCVP host-range redistribution generates **many potential introduced-only larval resources**, but this selected independently archived butterfly observation sample **does not show a positive excess of butterfly GBIF indices specifically in introduced-only resource regions after conditioning on species, island recording opportunities and broad biogeography**.

**This is a test of georeferenced GBIF indexing, not real butterfly occupancy or successful plant-mediated island establishment.** The result does not erase the frozen +54.9% global and +118.9% island resource-opportunity expansions or the separate positive 155-species host-informed known-butterfly recording forecasts.

## Frozen source crosswalk and denominators

- Original **32** selected butterfly names and **126** strict insular WGSRPD3 Level-3 botanical regions (32×126 = **4,032** resource combinations), from source `results/island_resource_gap_filling/resource_cells.csv.gz`.
- An **independent original GBIF DOI** 10.15468/dl.pp5nc9 with 53,143 available matched records. Only **31/32** selected butterflies have any records in this archived DOI subset. The unrecorded `Apodemia mormo` contributes **13 potential introduced-only host-region cells** but no DOI records, so it is excluded from the observational **denominator**.
- Recordable cohort = **31×126 = 3,906** source-species-region opportunities. Of these, **159** positive species×island-region GBIF archive indices across **26** butterflies and **49** regions were source-attested, with exact matches of **13 introduced-only** indices to the previous frozen Q3 audit. A non-index is **unknown**, never a real biological absence.
- Original full 32-resource cohort includes **933 introduced-only potential cells**; observed-index-evaluable 31-species subset has **920**. Treat these as separate denominators.

### Island butterfly indices classified by source larval-host range

| Plant resource in WGSRPD3 botanical region | All index-evaluable resource cells (31×126) | Positive butterfly GBIF DOI indices | GBIF archival index fraction (NOT ecological occupancy) |
|---|---:|---:|---:|
| At least one **native known host** (may also contain introduced known hosts) | **785** | **138** | **17.58%** |
| No native known host, **at least one regionally introduced known host** | **920** | **13** | **1.41%** |
| No known host catalogued from contemporary source flora | **2,201** | **8** | **0.36%** |
| **Total** | **3,906** | **159** | — |

## Two nulls that distinguish coarse record geography from host-specific excess

The nulls operate on the **same 159 observed-positive cells**, performing disjoint species×region swaps while retaining exact source record totals. Original butterfly-host identity and each resource category are held fixed throughout.

| What is preserved when relocating positive DOI source indices | Observed introduced-only indices | Simulated null mean, 95% range | One-sided excess p | Biological meaning |
|---|---:|---|---:|---|
| Each butterfly's and each island region's number of positive records | **13** | **46.82** (39–55) | **1.000** | Strong global mismatch; host and record distributions differ across biogeographic realms |
| Same margins **plus species-specific numbers of positive indices within each botanical Level-1 biogeographic realm** | **13** | **12.87** (10–15) | **0.632** | **No introduced-only excess** once broad plant/record geography is preserved |

For **native-host-covered** regions under the stricter Level-1-preserving null, **138** positive indices were observed versus **132.16** null expectation (95% null count range **129–136**; one-sided excess p≈**0.0023**). The complementary **no-known-host** group observed **8**, null mean **13.97** (11–17); descriptive two-sided p≈**0.0027**. These are post-hoc **exploratory** MCMC permutation probabilities, **not** multiplicity-controlled, effort-fully-matched or phylogenetically independent claims.

**Why the two nulls differ:** imposing butterfly × realm margins conditions on broad differences in which islands have historically observed butterflies. A globally unconstrained shuffle can place a butterfly in a botanically rich but irrelevant realm, generating a misleading expectation of introduced-host opportunity. Local realm preservation is much closer to the ecological question, although it does not replace species-specific climate suitability or survey-visit data.

### Cross-seed reliability of the strict fixed-margin null

The original strict null used **2,999** samples (80,000 proposed burn-in and 1,600 proposed edge swaps between sampled states). The edge-chain accepted **130,154** swaps (**2.67%** of proposals); its mean introduced-only count was **12.868**.

Two **independent random-seed / longer-thinning** strict null chains (999 states, 120,000 proposed burn-in and 3,200 attempts between samples each) agreed:

- Independent chain 1: introduced-only expectation **12.815** (95% 10–15), one-sided excess p=**0.602**; native expected **132.155**.
- Independent chain 2: introduced-only expectation **12.839** (95% 11–15), one-sided excess p=**0.605**; native expected **132.181**.

These chains are not independent ecological replicate datasets, and convergence diagnostics are limited to agreement of expected cell counts. They reduce concern that the null result is a single seed or thinning artefact.

## What the new result actually adds

The global 239-species resource reconstruction is an original **plant-defined potential**. The 155-species RF model established positive predictive information from **contemporary host resources and composition** for known butterfly recording distributions, even in held geographic blocks. That model does **not isolate introduced-only hosts** from the native host signal.

This new selected 31-species DOI panel separately asks whether archival butterfly source records occupy **introduced-only resource cells more than coarse resource geography and observation-source margins predict**. The answer under the better-controlled null is **not detectably**. Native resource cells retain record overlap, consistent with the expected value of native host geography for recorded butterfly distributions. It is not evidence that introduced plants are ecologically unusable or unable to support butterflies. Species unrecorded in a DOI subset cannot be treated as absent.

Neither the presence of a new host plant in a WGSRPD3 region nor an indexed butterfly observation implies the plant was used by a larva. Neither identifies ports, actual transport trajectories, colonization, host quality or population growth.

## Exact reproduction and original sources

- [Successful fixed-host-source panel and realm-constrained permutation audit](https://github.com/zuizui0223/chocho/actions/runs/43672529430), [full source SHA256 / 31×126 crosswalk / global and realm null / two independent long-thinning chains](https://github.com/zuizui0223/chocho/actions/runs/43672529430/artifacts/14145275581).
- Executable: `scripts/audit_independent_island_doi_host_overlap_null_20261115.py` consumes only committed source-frozen original `results/island_resource_gap_filling/resource_cells.csv.gz`, strict island registry, frozen 32-species list and original 2026 GBIF DOI selected positive cell index.
- Original DOI sample and independent source limitations: [GBIF archival occurrence support](../q3-observed-butterfly-island-opportunity-v02/q3_archival_gbif_overlap_audit.json).
- All 32 species source counts are in the archive artifact, including the single no-DOI-index species; no positive-only sample is silently equated with 31×126 species-region presence–absence data.

**Editorial placement:** Supplementary Results with one concise Limitations note. Keep GEB main Results focused on the supported 239-species resource transformation and the distinct 155-species known-butterfly recording forecast; present this source-limited test explicitly as the boundary on the stronger claim that *introduced-only host resources already yield demonstrable island butterfly distributions*.
