# Hidden phenological bottlenecks in butterfly host portfolios — feasibility and discriminating test

**Date:** 2026-10-10  
**Status:** EXPLORATORY, source-landing-page verified; no new raw-data analysis or ecological effect estimate yet  
**Branch:** `exploration/realized-host-window-v01`  
**GEB submission PR #38:** intentionally unchanged.

## The scientifically interesting question

**If two local butterfly populations have the same number, identity and total abundance of host-plant species, can they nevertheless differ in egg-placement continuity and later recruitment because the hosts' suitable developmental stages occur at different times?**

This asks about *realized functional complementarity*, not the tautology that the geographic union of hosts expands when introduced ranges are added. A butterfly cannot use a botanical checklist outside its oviposition window. However, the geometrical union of phenologically suitable plant dates is also mechanical, so **the primary biological endpoint must be independently observed oviposition or survival, not a calculated length of the union**.

Contrast three ecological hypotheses:
1. **Numerical redundancy:** once host abundance, host identities and sampling effort are controlled, adding a temporal complementarity descriptor does not improve observed egg placement or recruitment.
2. **Phenological insurance:** at comparable total host availability, greater between-host staggering produces a more continuous supply of used hosts and predicts the **distribution of actual eggs** through seasonal bottlenecks. The effect should persist in held-out site-years.
3. **Behavioral or physiological constraint:** staggered available hosts do not provide insurance when ovipositing females avoid those hosts, larvae cannot mature or enemies make them poor choices. The sign may be zero or negative.

A later *colonization or population-resilience* hypothesis would require independently dated introduced-host establishment, cohort-followed egg-to-adult survival and local population time series. This pilot does not have those variables.

## Why this is not already resolved by the chocho expansion figures

- Chocho's frozen 239-species resource reconstruction records whether known larval hosts occur in WGSRPD3 regions, not the dates and developmental stages at which they can be used.
- The butterfly *Anthocharis cardamines* is in the exact 239-species figure-source panel: 11 accepted known hosts, 162 native host-resource WGSRPD3 cells, 272 contemporary cells and 110 introduced-added cells. These are **coarse species-wide modeled opportunities**, not 110 observed colonizations or site-year supports.
- Independent observations of phenology and oviposition can in principle falsify the assumption that additional mapped hosts imply usable opportunities throughout the butterfly's flight period.
- The independent 2018 study itself already demonstrated that oviposition favors hosts nearer the preferred developmental stage. **Showing this again is not novel.** A new result would require a *conditional portfolio buffering effect* or a transferable, held-out ecological prediction beyond the published per-host synchrony pattern.

## First empirical source (verified only at landing-page / metadata level)

**Toftegaard et al. (2018), Oikos** DOI https://doi.org/10.1111/oik.05720  
Original Dryad DOI https://doi.org/10.5061/dryad.jp328r6

- Study organism: *Anthocharis cardamines*.
- Plant species: *Arabis hirsuta*, *Cardamine pratensis*, *Arabis glabra*, *Arabidopsis thaliana*, *Thlaspi caerulescens*, *Capsella bursa-pastoris*.
- Field sampling: three regions, four seasons (2010–2013). Source archive lists **Plant phenology in the field.txt** (phenology and egg/larva occurrences) and **Sampling site areas.txt**.
- Already-published finding: host timing predicts egg placement, both among host species and among individuals within a host species. This is prior art and cannot be repackaged as our discovery.
- Source access gate: original landing page/README checked 2026-10-10, but original downloadable file stream at `https://datadryad.org/downloads/file_stream/79537` returned HTTP 403 in an independent retrieval attempt. **No row-level schema, sample size, exposure denominators or fitted effect are validated yet.** Do not substitute search snippets for original rows.

The data are a **phenological mechanism anchor**, not proof that a *newly introduced* host contributes to local abundance or fitness. Botanical native/non-native statuses must be reconciled per *historical study region and year* before any globalization-specific inference.

## Related independent endpoints and limitations

| Source | Exact observed endpoint | Scope | Cannot demonstrate |
| --- | --- | --- | --- |
| Barton & Bogner (2024), https://doi.org/10.5061/dryad.866t1g1xb | 224 *Vanessa tameamea* larvae across three populations and five host diets; egg-hatch to adult-eclosion success (incl./excl. flight-inhibiting deformities) | One exotic host (*Cecropia obtusifolia*); a population × host contrast already published | Oviposition choice, real local recruitment or generality across butterflies |
| Singer & Parmesan (2018, 2021), https://doi.org/10.1038/s41586-018-0074-6 and https://doi.org/10.1111/gcb.15656 | Local *Euphydryas editha* loss after adopted host/microhabitat change, multi-population host-use histories | Lock-in and reversal possibility **already discovered** | A general causal estimate across species or a pure plant-removal effect |
| Colom et al. (2023), https://doi.org/10.5061/dryad.8w9ghx3td | Site-year abundance/population growth, host *Rhamnus* cover, larval parasitism for two *Gonepteryx* spp. | Explicit enemy/resource/temperature test; study did **not** find biogeographic support for shared-enemy apparent competition | Alien-host-induced competition or host-portfolio phenology |
| Larsson Åberg et al. (2026), https://doi.org/10.5061/dryad.rxwdbrvrv | Individual female oviposition, larval growth and host distributions for *Aricia artaxerxes* | 7 egg-laying females; preference for common but somewhat poorer-growth host already published | Adult recruitment; global pattern; repeated local populations |

**Do not pool** egg-placement odds ratios, larval growth and adult-emergence risks as if they share a denominator or represent independent experimental replications of one butterfly-globalization effect.

## Auditable first-pass analysis (only after raw source gate)

1. Confirm original files are byte-authentic and parse with explicit delimiters/encoding. Freeze a raw column dictionary, species/taxonomy map, calendar-week and stage codebook *before* looking at any performance contrast.
2. Require egg-negative as well as egg-positive host individuals and **known plant examination effort**. Verify sampling dates, site × year clusters, full relevant host community within each site, and independently measurable adult-flight periods (or label the result egg-sampling-window-only).
3. **Primary response:** eggs per surveyed plant, at site × week (or per-host examined plant with effort offset). **Main contrast:** a botanical count/abundance/season baseline vs an additional predeclared between-host phenology complementarity term. Include total available plant abundance, phenological stage distribution, site/season and host identity; do not make a mechanically defined 'suitable-day union' the outcome.
4. Report site-year-blocked predictive log-score difference, region-wise signs, egg-count scale predictions and uncertainty. Do not treat thousands of plants as thousands of independent site-year replicates. Never optimize seasonal windows after seeing oviposition.
5. A genuine buffering signal requires improvements in independent held-out site-years **after controlling contemporaneous phenologically suitable plant abundance**; any residual complementarity must be interpreted as a behavioral/site-year association, not proof of adaptive demographic insurance.
6. If the source lacks zero-egg plants, stage descriptors, repeated dates, independent flight timing or enough clusters, **STOP**: record data feasibility failure. Do not reframe per-host phenology correlations as new population resilience.
7. No changes to the GEB manuscript or source claims on failure.

## What a publishable second experiment must add

Replicate **at least two butterfly species, each in several independent populations**, ideally with native and introduced established hosts at matched sites. Manipulate host-stage synchronization while keeping host species identities, individual abundance and biomass constant (synchronous vs staggered), and cross with focal-host removal to test fallback. Record: individual maternal plant choice; eggs per accessible plant; monitored offspring to **flight-capable adult**; predators/parasitoids, plant microclimate and season; population change or adult recruitment over ≥2 generations if claiming resilience. Use site/population rather than eggs as the independent replication unit; decide sample size via a pilot-derived variance estimate, not an invented power guarantee.

**Decisive contrast:** does staggered host availability preserve actual viable offspring or post-removal recruitment *at equal total host supply*? Positive, null and negative outcomes all discriminate mechanisms. A pure descriptive calendar overlap result does not.

## Immediate decision

**Priority:** source-gated mechanistic pilot on temporal host complementarity, followed by an explicitly separate experimental recruitment test if justified.  
**Not currently established:** phenological insurance, beneficial alien hosts, fitness buffering, resource-tracking debt or consumer population homogenization.  
**Submission boundary:** PR #38 remains the independent global *potential-resource* paper.
