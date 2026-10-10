# Butterfly host-use hysteresis: can adding a host reduce future ability to fall back?

**Date:** 2026-10-10  
**Status:** INDEPENDENT MECHANISTIC HYPOTHESIS / NOT AN EMPIRICAL RESULT  
**Branch:** `exploration/realized-host-window-v01`  
**Existing GEB manuscript:** unmodified.

## Primary question (beyond geographic expansion)

> Holding **contemporary host identities, biomass, phenological access and enemy exposure constant**, do butterfly populations with different histories of using a redistributed host differ in their ability to return to ancestral hosts after the redistributed host declines or disappears?

The distinctive prediction is **ecological hysteresis** (path dependence): populations with the same *present* resource landscape can have different host-choice behavior, developmental outcomes and responses to removal because their *past* resource landscape differed. Increased mapped host richness is not sufficient to predict the direction of the effect.

This is not a claim that all exotic hosts are ecological traps. A new host may benefit consumers, harm them or have neutral consequences. The novel biological target is the **difference in reversible fallback ability at a fixed present host supply**, with demographic consequences experimentally tested.

## Explicit competing hypotheses

- **H0, instantaneous resource tracking:** after matching present host identity, density, phenological stage, maternal experience, temperature and enemies, historical introduced-host exposure predicts neither ancestral-host oviposition nor egg-to-adult production after introduced-host removal.
- **H1, behavioral hysteresis:** historically introduced-host-exposed populations retain disproportionate oviposition on the introduced host or delay switching back to ancestral hosts, even when native hosts are present; this reduces early reproduction after standardized host loss.
- **H2, evolved/physiological hysteresis:** in a reciprocal common garden (after standardized maternal experience), historical exposure still predicts performance on ancestral hosts, potentially after multiple generations; controls are required before calling it evolutionary.
- **H3, context-only apparent hysteresis:** population differences disappear when phenology, host nutritional stage, local temperature, natural enemies and current plant abundance are standardized. The historical association was ecological confounding, not lock-in.
- **H4, maintained fallback:** populations that have used an introduced host *also retain* ancestral-host preference and larval competence. History does not imply fragility; the effect can be absent or opposite.

## Prior art that removes false novelty

1. **Singer & Parmesan (2018), Nature** DOI https://doi.org/10.1038/s41586-018-0074-6: a Nevada *Euphydryas editha* population became dependent on introduced *Plantago*, then was lost after land management changed the local microenvironment while ancestral *Collinsia* was still present. **Host trapping/lock-in is not a new concept and the management/thermal change must not be falsely described as a simple host-removal manipulation.**
2. **Haan et al. (2021), Scientific Reports** DOI https://doi.org/10.1038/s41598-020-80413-y: another *E. editha* population maintained ability to oviposit and develop on ancestral *Castilleja* after adopting *Plantago*. This provides a *different-population contrast*, not a balanced causal comparison of historic exposure.
3. **Singer & Parmesan (2021), Global Change Biology** DOI https://doi.org/10.1111/gcb.15656: repeated host-use records from 15 *E. editha* populations already documented local dietary changes and patterns of narrowing. Do not describe these published trends as a new chocho discovery.
4. **Barton & Bogner (2024), Biotropica** original data DOI https://doi.org/10.5061/dryad.866t1g1xb: *Vanessa tameamea* 224-larva no-choice experiment, five plants including exotic *Cecropia*, three source populations. Different population responses to the exotic host were already published. These data assess **population × diet performance**, not historical exposure, maternal choice or local post-removal population growth.
5. **Toftegaard et al. (2018), Oikos** DOI https://doi.org/10.1111/oik.05720: preference depends on host developmental stage. Thus a valid test of hysteresis must also control *within-season phenological access*. The Dryad source-audit workflow independently confirmed that original files currently return HTTP 403; no data were reanalyzed.

**Generalization beyond these separate prior examples is presently unproven.** Different published populations, experiments and host contexts cannot be treated as replicated conditions of a controlled experiment without checking their comparability.

## Distinguishing empirical design

### Phase A — find matched local contrasts without selecting outcomes

For each of **at least two butterfly species**, locate populations occupying matched contemporary host assemblages where both (i) ancestral and (ii) introduced recorded food plants are locally present. Measure historical exposure using independently dated introduction/planting/use records, not the same present butterfly outcome. Pre-select high/low historical exposure *before* examining response. Populations with no reliable historical chronology are excluded from the historical-causality comparison.

Match/pair by species, climatic regime, host species and phenological stage, host biomass, vegetation management, enemies, maternal condition and time of season. Replicate **populations and sites**; dozens of eggs from one female do not become independent populations. Stratifying on 239-species chocho mapped opportunity is allowed for candidate discovery, but the global added-region number cannot be the history measure.

### Phase B — reciprocal preference and quality experiments

In each population, measure at the same time:
- Female choice among locally available native and introduced hosts under standardized availability; repeated trials can identify learning/experience, but treat female as repeated-measure unit.
- **No-choice** larval growth and **flight-capable adult eclosion** on ancestral and introduced diets; keep all failed larvae in the denominator; track egg family and source female.
- Plant stage, tissues actually accessible, nitrogen/chemistry proxies, enemies, microclimate and seasonal synchrony.
- Follow host removal in a randomized treatment at the **replicate site or enclosed population**, paired with a sham-retention control. Do not remove invasive plants in natural habitats merely to run this test without site permissions and management safeguards.

### Phase C — randomized host removal, with exposure history as an observational modifier

Primary outcome: population-level **change in viable adult production per female/initial cohort after experimental loss of the adopted host**, comparing historically exposed and less-exposed populations that have the same available native hosts. The estimand is the **exposure-history × randomized-removal interaction**, estimated with species and population structure; interval estimates and predicted absolute adult output are mandatory. **Removal can be randomized, but past host exposure cannot be retroactively randomized**: the interaction estimates differential responses to removal conditional on measured history, not a causal effect of past exposure. A stronger mechanistic test would randomize replicated lineages to prolonged exposure/no exposure in a common garden, then cross that treatment with subsequent removal.

A positive-history association with oviposition preferences but no adult-output loss supports *behavioral history dependence*, not a demographic trap. Persisting quality deficits after standardized common-garden generations provide evidence **consistent with** a heritable component, but genomic/adoption-history inference needs additional tests. A negative post-removal interaction is not automatically adaptation: management, microhabitats, parasites and dispersal alternatives must be assessed.

### Key falsification gates

- If historically exposed populations do not differ from controls at fixed present resource supply and context, **reject hysteresis**.
- If any apparent deficit vanishes on matching phenology and microclimate, classify it as **contemporary ecological filtering**, not host-use memory.
- If choice and larval performance differ but viability after host loss does not, claim **preference/performance decoupling**, not demographic fragility.
- If the only support comes from the already-published *E. editha* case, do **not** claim a broad butterfly law.
- If historic host-use exposure cannot be dated independently, do not estimate a historical-treatment effect; restrict to current-population variation.

## Role of the current source-audited pilot

The exploratory `REALIZED_HOST_WINDOW_HYPOTHESIS_V01.md` asks whether host phenology creates *functional* bottlenecks at matched botanical supply. Its job is to identify an alternative explanation for apparent lock-in and to develop plant-stage measures; it is not itself a demonstration of species-wide resilience.

The existing 239-butterfly static network can support **source-blind stratification** and define a baseline of theoretical host choices. It cannot, even with more permutation tests, estimate behavioral hysteresis, cohort survival or population lambda.

## Scientific decision

**Primary new ecology hypothesis:** historical exposure to redistributed hosts generates *path dependence* in realized diet and capacity to fall back to ancestral hosts.  
**Mechanism control:** local phenology, host quality and enemies.  
**Current evidence:** prior examples and study-design feasibility only; **no new effect yet**.  
**Research implication:** move away from adding static spatial-overlap statistics and toward a matched-history × host-removal experimental design with independently followed offspring. Treat any observational archive as pilot/hypothesis refinement, not substitute for the causal test.
