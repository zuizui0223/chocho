# Non-obvious macroecological findings: butterfly resources under plant globalization

**Date:** 2026-10-15  
**Scope:** Evidence-ranked synthesis of previously computed results in [chocho PR #38](https://github.com/zuizui0223/chocho/pull/38), base head 88209b8a215ac1bc9c4c5e67622fa104179c8895.  
**Status:** DISCUSSION ONLY. No new statistical test, data collection, butterfly experiment, alteration of a frozen result, or change to the GEB submission manuscript. Do not merge into the submission branch automatically.

## Central ecological question

**When existing host plants move across the world, how does the *architecture* of potential larval resources change for butterflies, beyond a mechanically larger resource map?**

The estimand is the geography of **potential host-resource opportunity**, reconstructed with butterfly–host identities held fixed. The dataset does not directly estimate larval survival, butterfly colonization, interspecific competition, temporal host phenology, or population growth.

A defensible summary is:

> Plant redistribution can increase resource opportunity for individual butterflies while making potential resources more geographically shared, often supported by only one redistributed host. Most raw geographic convergence is explainable by how much opportunity is added; only a small, reproducible component exceeds constrained geometric expectations. In a limited occurrence panel, added resources coincide with observed butterflies more often than matched geography alone predicts.

This is the extension from known **plant-flora homogenization** to the **resource layer experienced by butterfly consumers**, not a claim that realized butterfly communities became homogenized.

## Evidence hierarchy: results worth highlighting

### A. The usual specialist/generalist axis does not encode geographic resource specialization

- Host-family breadth versus native host-resource geography: Spearman rho = **0.276** in all 239 resource-eligible butterflies, **0.398** in the 215-species host-taxonomy adequacy subset.
- In the latter subset, **33/215** species have taxonomic/geographic percentile-rank discordance >= 0.5, and **18** butterflies recorded from just one host-plant family lie in the **top quartile** of native resource geography.
- Proportional anthropogenic resource increase is only weakly associated with host-family breadth: rho = **0.008**, species-bootstrap 95% CI **[-0.111, 0.128]**. The post-hoc +/-0.10 TOST gave p = **0.078**, so equivalence to zero is **not** established.
- This pattern is not simply due to Poaceae-associated butterflies: excluding all 58 butterflies with Poaceae hosts leaves rho = **0.026** (n=181).

**Ecological reading:** A one-family specialist can still have access to a geographically extensive native host-resource pool, and the plant-side spatial response to globalization is not well predicted by host-family count. This is a decoupling of ecological axes, **not** proof that diet breadth never matters.

**Receipts:**
- [specialization dimensionality](../../benchmarks/exploratory/butterfly_specialization_dimensionality_result_v0.1.json)
- [anthropogenic expansion](../../benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json)
- [equivalence diagnostic](../../provenance/reviewer_defenses/results/butterfly_expansion_equivalence_v0.1.json)
- [grass-guild sensitivity](../../provenance/reviewer_defenses/results/butterfly_poaceae_sensitivity_v0.1.json)

### B. Global opportunity is concentrated on the PLANT side, but a special butterfly-host preference is not established

- **670** introduced-range host plants contribute some of **14,553** added butterfly × WGSRPD3 units; **38** plants receive half the fractionally allocated credit (top 10 = 25.1%; top 50 = 57.0%). The pattern persists at genus level and after crop-host exclusions.
- Across **8,909** plants in **278** families, host network prominence (the count of Lepidoptera species recorded using a plant) is positively associated with plant geographic expansion: rho = **0.313**; stratified, native-range- and family-aware rank association = **0.295**, permutation p < **0.001**.
- Crucial falsification: a same-family host-identity null suggests unusual focal-host expansion, but the observed excess **disappears** when alternative hosts are also aligned by native-range breadth and other-Lepidoptera use (207-species observed 12,853 added units; strict null median **12,971**, p = **0.584**). The importance of that weighting is itself parameter-sensitive.

**Ecological reading:** The *plant-side network prominence axis*, rather than butterfly diet breadth alone, is a promising organizing description. However, recorded host degree conflates ecological host generality/commonness with sampling and literature intensity. Do **not** claim selective evolution toward especially globalized plants or that prominence fully explains the pattern.

**Receipts:**
- [host contribution concentration](../../provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json)
- [plant network prominence](../../provenance/reviewer_defenses/results/butterfly_host_plant_prominence_expansion_v0.1.json)
- [host prominence-weighted null](../../provenance/reviewer_defenses/results/butterfly_resource_expansion_hostbias_null_v0.1.json)
- [null-weight sensitivity](../../provenance/reviewer_defenses/results/butterfly_hostbias_parameter_grid_v0.1.json)

### C. The important homogenization result is the *decomposition*, not the +66.8% headline

- Among **355** baseline-active WGSRPD3 regions (62,835 region pairs), potential butterfly-resource assemblage mean Jaccard changes **0.27698 -> 0.46208** (raw +0.18510; +66.8%).
- A fixed-margin null retaining every butterfly's count of added regions and every region's count of added butterfly opportunities predicts median **0.45348**. Thus observed-minus-null is only **+0.00860**, about **4.6%** of the raw increase; Monte Carlo p = **0.002** with 499 permutations.
- A stricter continent/Level1-constrained fixed-margin null leaves an excess of **+0.003616**, about **2.0%** of the raw increase.
- A post-hoc, Europe-biased BCE/Clarke host-list augmentation adds 1,027 candidate butterfly–host pairs for 81 focal butterflies and changes aggregate proportional expansion **54.9% -> 50.1%**, while Level1-constrained homogenization excess stays positive but tiny (**+0.003436**; 499 draws, p = 0.002). Not a globally corrected estimate.
- Resource-envelope overlap between 28,441 butterfly pairs rises **0.222 -> 0.328**, against matched-null median **0.321**; the same structural-versus-residual distinction applies.

**Ecological reading:** Redistribution of plants *can*, in principle, produce distinct rather than convergent consumer-resource assemblages, so the measured configuration matters. Yet the **vast majority of the raw similarity increase is structurally expected from opportunity additions**. A very small p-value on a tiny margin-corrected residual is not an ecologically large effect. Avoid implying that the full +66.8% is host-identity-specific.

**Receipts:**
- [resource homogenization](../../provenance/reviewer_defenses/results/butterfly_resource_homogenization_v0.1.json)
- [Supplementary Table S9 and Level1 null](../../manuscript/butterfly_specialization_supplement_v0.2.md)

### D. Identical host links are geographically replicated across butterfly species and families

- Exact shared-host butterfly-pair × region units: **62,473 -> 141,885**, or **+127.1%**; **920 of 986** host-sharing butterfly pairs acquire at least one new shared-resource region.
- Of **79,412** added pair × region units, **41,543 (52.3%)** connect butterflies from different families.
- Among the top 38 plant contributors, **26** serve butterflies in at least two focal butterfly families; concentration in a few hosts does not mean the homogenization signal is attributable solely to those hosts.

**Ecological reading:** Existing trophic associations can be *spatially replicated* by movement of plants, increasing **potential shared-resource exposure** without any host switch. This is richer than “host ranges grow,” but **shared host != simultaneous local presence != resource competition**. An effect beyond a network-aware shared-host null has not been demonstrated.

**Receipt:** [shared-resource exposure and redundancy](../../provenance/reviewer_defenses/results/butterfly_resource_overlap_resilience_v0.1.json)

### E. A bigger resource network need not be locally redundant

- **8,574/14,553 (58.9%)** of newly opened butterfly × region cells depend on **one** contributing introduced host; the median number of contributing hosts is one.
- In a coarse regional removal stress test, suppressing introduced contributions from the 38 highest-credit plants removes **5,812/14,553 (39.9%)** added opportunity units, affects 123 butterflies, and erases all added opportunity for 16. The difference between 50.6% fractional contribution and 39.9% unique opportunity loss reflects non-identical credit and union-loss accounting.
- The comparison to random deletion of 38 hosts is **not** an independent significance test of an unexpectedly fragile network: the targeted 38 were chosen based on observed contributions.

**Ecological reading:** “More geographic opportunity” can coexist with **thin host redundancy**, yielding a testable vulnerability hypothesis. This does **not** justify retaining invasive hosts, predict local plant-removal outcomes, or quantify butterfly persistence.

**Receipt:** [shared-resource exposure and redundancy](../../provenance/reviewer_defenses/results/butterfly_resource_overlap_resilience_v0.1.json)

### F. There is a modest bridge to REAL butterfly distribution, but not to demographic success

- In the secondary climate-stratified, 32-species panel, **23** species had observed outside-native resource cells; **66/115 (57.4%)** outside-native butterfly × region occurrences were inside added introduced-host opportunity.
- Region-composition-matched structural null median: **49/115** recovered, against **66/115** observed. Species-level null and leave-one-species-out checks preserve the excess.
- Species-cluster bootstrap 95% interval of the pooled recovery fraction: **0.349–0.774**. Individual failures are real: **Erynnis tristis 0/11** outside-native units recovered.
- The panel was *not* selected to prospectively validate butterfly colonization and GBIF presence does not prove local larval host use.

**Ecological reading:** The reconstructed plant-resource changes align with observed butterfly occurrence beyond simple region composition, which is more informative than an area increase. However, presence is still not evidence of local breeding on introduced hosts, successful reproduction, or host-caused range expansion.

**Receipts:**
- [occurrence overlap](../../provenance/reviewer_defenses/results/butterfly_occurrence_resource_validation_v0.1.json)
- [region-matched overlap null](../../provenance/reviewer_defenses/results/butterfly_occurrence_overlap_null_v0.1.json)
- [species robustness](../../provenance/reviewer_defenses/results/butterfly_occurrence_species_robustness_v0.1.json)

## 2026-10-15 direct follow-up: species-level redundancy is not taxonomic breadth

The frozen figure source \`data/frozen/figure_sources/host_contribution_metrics.csv\` independently supports a **within-species / local-redundancy distinction**:

- Across all **206** butterflies with added host-resource regions, the median species-level fraction of new regions supported by exactly one introduced plant host is **68.8%**. This species-unweighted figure differs from the **58.9%** (8,574/14,553) opportunity-unit-weighted number; both denominators are explicitly distinguished.
- The median single-host fraction declines across resolved host-species portfolio classes: **100%** (1 host, n=14), **98.8%** (2–4 hosts, n=67), **67.6%** (5–9, n=45), **46.5%** (10+, n=80).
- Yet **37/80** butterflies with 10 or more resolved hosts have >50% of added regions supported by only one introduced host. Even requiring **20 or more newly added regions**, it remains **35/78**, and all 78 pass the host-taxonomy adequacy floor. Notable source-conditional examples: *Apodemia mormo* (17 hosts, 52 new regions, 100% singly supported) and *Erynnis horatius* (13, 25, 100%).
- Exploratory rank correlations of **multiple** host support across these 206 species: with host-family breadth **rho=0.291**, with resolved host-species richness **rho=0.669**. Partial rank correlation of family breadth, conditional on richness, is only **0.015**; also adjusting native geography and added geography yields **0.035**. These numbers are descriptive post-hoc, not corrected phylogenetic or equivalence tests; richness creates a built-in combinatorial opportunity for redundancy.
- The intuitive interpretation is not “host-rich butterflies are necessarily fragile.” Rather, **recorded host richness is an imperfect predictor of local substitute resources** because candidate hosts may occur in different regions.

**Occurrence cross-check—demoted, not promoted:** In the climate-stratified occurrence data, 23 species have outside-native records, but two have **zero** reconstructed introduced-host gain and hence an **undefined**, not zero, multiple-host fraction. Among the remaining **21**, the simple rank association between multiple-host support and outside-native occurrence recovery is **rho=0.649**; conditioning on resolved host richness, native resource geography and added opportunity drops it to **rho=0.262** (add outside-native count: **0.219**). Shared input geography, small selected sample, no local larval use and observational confounding block a novel establishment claim.

Reproducibility:
- [Self-contained species-level analysis](../../scripts/analyze_butterfly_species_redundancy.py)
- [Frozen-data and synthetic tests](../../tests/test_butterfly_species_redundancy.py)

## New structural null: do not confuse a planned contrast with a supported effect

A distinct [plant-range swap null](../../scripts/analyze_butterfly_shared_host_structural_null.py), with [synthetic invariance tests](../../tests/test_butterfly_shared_host_structural_null.py) and a [pinned data-rebuild workflow](../../.github/workflows/butterfly-plant-range-structural-null.yml), randomizes *introduced plant × WGSRPD3 region* incidences **within plant-family × WGSRPD Level1 blocks** while preserving plant introduced-region count within each Level1, plant-family-specific regional richness, native host exclusions and the fixed butterfly–host network.

It recomputes (a) novel exact shared-host butterfly-pair × region units, (b) newly added butterfly × region opportunity, (c) their ratio, and (d) fraction supported by a single introduced host. The model includes a predeclared **mixing/evaluability gate**, and records null medians and intervals only if movement is adequate.

**Essential null boundary:** This plant-level randomization does *not* preserve the number of added regions for each butterfly, nor demonstrate ecologically plausible host introduction routes. It is a post-hoc **structural decomposition**, complementary to—not a replacement for—the already reported butterfly- and region-fixed margin null. **Until the empirical workflow's receipt is checked, there is no positive or negative finding from this new null.**

## What NOT to promote as findings

1. **“More host-plant range = more butterflies.”** Only potential larval resource opportunity is reconstructed; butterfly density, recruitment and conservation benefit were not estimated.
2. **“Butterfly communities have homogenized by 66.8%.”** This number is potential resource-assemblage similarity, not observed composition; nearly all the raw gain is fixed-margin-predictable.
3. **“Introduced hosts intensified competition.”** Exact shared-host exposure is a necessary opportunity for interaction, not observed competition.
4. **“Plant network prominence is a proven biological cause.”** The HOSTS degree measure carries unknown recording bias, and identity-null conclusions depend on weighting.
5. **“Diet breadth has exactly no influence.”** The near-zero point estimate does not pass the declared post-hoc tight equivalence diagnostic.
6. **“The global temporal host–larva mismatch is proven.”** WCVP status has no introduction dates or seasonal edible-leaf trajectories. Kyoto's native-host/introduced-butterfly cage system is a distinct mechanistic test bed, not a direct global introduced-host validation.
7. **“Macro-data identify whether anthropogenic hosts produce self-sustaining butterflies.”** Presence/resource overlap falls short of that stronger demographic question.

## Manuscript storyline if prioritizing non-obvious ecology

**Title direction (not yet a title decision):** *Widespread opportunity, sparse redundancy: butterfly resource geography under plant globalization*.

**Question 1:** Does taxonomic host specialization capture the geography of larval resource access? **Answer:** weakly and incompletely (A).

**Question 2:** Do introduced host ranges just enlarge resource geography or reshape who can share which plants, where? **Answer:** both; but margin-controlled homogenization is **small**, whereas exact shared-host exposure broadens considerably (C, D).

**Question 3:** Does expanded opportunity mean diverse backup resources or ecological establishment? **Answer:** often single-host-backed (E); it exceeds geometric expectations for contemporary butterfly occurrence in a limited panel (F), but breeding and demographic effectiveness remain unresolved.

**Cross-cutting descriptor:** Plant network prominence and inventory composition organize resource changes more plausibly than a simple specialist/generalist dichotomy, but their ecological versus observational origins remain unresolved (B).

### Figure priorities

- **Fig A — two dimensions:** butterfly host-family count versus native geographic resource breadth, highlighting the 18 one-family/high-resource butterflies, with the 0.008 expansion slope and CI explicitly labelled as a different estimand.
- **Fig B — observed vs null:** regional resource Jaccard before/after alongside fixed-margin and continent-constrained medians; label the **+0.003616 residual**, rather than graphically inflating raw +66.8%.
- **Fig C — shared versus vulnerable opportunity:** exact shared-host pair × region gain, between-family fraction, and 58.9% single-host-added cells; show these as **structural**, not null-validated mechanisms.
- **Fig D — external reality check:** species-level occurrence recovery vs region-matched null, with 23-species variation, confidence/uncertainty and clear failure examples.

These are **figure storyboards**, not instructions to silently replace already frozen GEB Figures 1–3.

## High-value next tests and decision criteria

1. **Strongest feasible macroecological test:** strengthen butterfly occurrence corroboration across independent species/regions with predeclared GBIF effort and sampling filters, rather than continue to count new host-area cells. Predeclare species-held-out sensitivity and report both recovered and unrecovered occurrences. A disappearing region-adjusted excess would demote the “resource opportunity tracks realized geography” argument.
2. **Crucial structural counterfactual:** preserve host-plant redistribution intensity, geographic opportunity margins and butterfly-host network degree while perturbing which shared hosts are placed where. Then test whether the **+127.1% shared-resource exposure** or **58.9% single-host support** exceeds a corresponding null; until then, do not call either an “excess”.
3. **Biological mechanism beyond these maps:** independent observations of larval acceptance, young/mature leaf availability, phenology, survival and completed adult emergence can establish the filter between host presence and butterfly production. Treat this as a targeted case-study/field hypothesis, not a result from the global static panel.

**Stop rule:** Do not generate an unlimited list of post-hoc predictors to obtain larger effects. Prefer a few controls that can *invalidate* the most important claims.

## External position and interpretation limit

- [Yang et al. (2021), *The global loss of floristic uniqueness*](https://doi.org/10.1038/s41467-021-27603-y) already demonstrates global *plant* homogenization. The distinct proposed contribution here is the rearrangement of **butterfly larval-resource geography under fixed host links**.
- [Singer & Parmesan (2018), *Lethal trap created by adaptive evolutionary response to an exotic resource*](https://doi.org/10.1038/s41586-018-0074-6) establishes in one system that an introduced host can support larvae yet create management-sensitive dependence; it cannot validate worldwide ecological traps from mapped host occurrences.

**Recommended claim ceiling:** *Potential butterfly resource opportunity becomes more widely shared yet frequently relies on single introduced hosts. Partial alignment with observed butterfly occurrences supports its relevance, but population establishment and ecological outcomes remain outside the global estimand.*
