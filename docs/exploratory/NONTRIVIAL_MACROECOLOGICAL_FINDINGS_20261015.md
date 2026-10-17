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

## 2026-10-17 host-list source-completeness challenge: revised evidence priority

A **new source-pinned reconstruction** [run 38988757959](https://github.com/zuizui0223/chocho/actions/runs/38988757959) applied the *exact previously archived*, evidence-filtered European BCE/Clarke (2024) accepted-host candidate relationships to the same frozen WCVP and HOSTS host-geographic sidecars. It revalidated the original **239 butterflies, 26,530 native cells, 41,083 contemporary cells, 14,553 added resource cells and H=75,346 candidate introduced-host incidences** and independently reproduced the prior Supplementary Table S9 totals after augmenting all 1,027 candidate host links (81 butterflies, 679 WCVP accepted plant IDs). The [new audited receipt](../../provenance/reviewer_defenses/results/butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json) and [fixed analysis protocol](../../provenance/reviewer_defenses/protocols/BCE_CLARKE_RESOURCE_NOVELTY_SENSITIVITY_20261017.md) establish this as a **source-completeness sensitivity, not a novel biological null**.

| Fixed resource accounting, 239 species | HOSTS original | BCE-augmented partial list |
|---|---:|---:|
| Candidate introduced accepted-host × butterfly × region incidences, `H` | 75,346 | 103,349 |
| Already native-resource-supported incidences, `A` | **46,387 (61.6%)** | **67,740 (65.5%)** |
| Redundant extra supports within newly opened cells, `R` | 14,406 (19.1%) | 20,642 (20.0%) |
| Distinct newly opened butterfly×region cells, `U` | **14,553 (19.3%)** | **14,967 (14.5%)** |

The candidate denominator **H itself changes** (more known host species), so 65.5% is **not an independent global biological prevalence estimate** or a precise new "true" ecological fraction. It is an example of the uncertainty induced by adding European literature-derived accepted-WCVP candidate host links. The **+3.98 percentage point descriptive inventory shift** in A/H is on a much larger numerical scale than the **+~0.36 point old-source conditional spatial-null excess**, but these compare *different interventions and estimands*. In particular **the spatial null was not rerun after the host-list augmentation**. Do not call the ~0.36-point excess inventory-robust on that basis.

The net number of new cells changes by only **+414**, concealing substantial *cancellation*: **1,817** newly added butterfly×region cells versus **1,403** formerly introduced-added cells now reclassified to the *native* resource envelope. Among 81 augmented butterfly species, the introduced-added area count rises for **36**, falls for **35** and remains unchanged for **10**. This is a nontrivial *measurement/interpretation* lesson: geographic novelty of potential resource supply is strongly conditional on which known larval hosts define the baseline.

**Revised claim hierarchy:** (1) Resource presence is not equivalent to novel opportunity; the A+R+U decomposition is source-conditional but reproducible. (2) Catalogue incompleteness demonstrably changes the partition while sometimes cancelling in aggregate. (3) Small climate- and locality-constrained spatial overlap residuals remain exploratorily interesting **only for the frozen original HOSTS graph**, not a robust globally causal ecological mechanism. (4) Local feeding, larval development and demographic fitness remain unknown. **Do not promote the source sensitivity to a second stand-alone ecology paper or silently overwrite the GEB manuscript.**

## 2026-10-17 individual-butterfly host-inventory turnover: geographic identity matters more than aggregate net gain

After the source-completeness reconstruction [successful run 38988757959](https://github.com/zuizui0223/chocho/actions/runs/38988757959), an independent source-checked follow-up [run 38992305236](https://github.com/zuizui0223/chocho/actions/runs/38992305236) computes exact set-identity agreement in added butterfly×WGSRPD3 resource regions under the original HOSTS vs BCE/Clarke-augmented partial host catalogue. The original source includes per-butterfly counts for 81 European evidence-augmented butterflies. The complete derived species table and plot are generated from the original archived JSON without manually selecting regions; the compact [frozen receipt](../../provenance/reviewer_defenses/results/butterfly_bce_added_region_identity_turnover_v0.1.json) contains source, mathematical identities and auditable extreme cases.

Across those **81 modified butterflies**:
- `old newly accessible cells` = **4,628**, `newly accessible cells after supplementation` = **5,042**; net growth **+414** (+8.9%).
- **1,403** old "introduced-added" cells become *native-supported* once another recorded native host is recognized. **1,817** other cells become introduced-added under the expanded host catalogue.
- Shared classification `intersection = 4,628-1,403 = 3,225`; possible new classification `union = 4,628+1,817 = 6,445`; **Jaccard = 3,225/6,445 = 0.5004**, so **49.96% of the union is not shared between source classifications**.
- **54/81** modified butterflies simultaneously gain other newly supported cells and have some former introduced-added cells reclassified as native. Net added-region counts increase in **36**, decrease in **35**, remain unchanged in **10**.
- Across the **77** butterflies with at least one newly supported cell under at least one host inventory (excluding 4 with no added cells either way), the per-butterfly median resource-region Jaccard is **0.5705**. The *full-239* Jaccard is artificially higher at **0.8033**, because the other **158 species** have unchanged host inventories.

A particularly instructive **hidden-turnover example** is *Thymelicus lineola*: **65→65** potential introduced-added resource regions with **27 prior additions reclassified native and 27 newly added**, yielding within-species Jaccard **38/92 = 0.413** despite zero change in net regional count. Other source-extreme examples: *Lycaena dispar* **0→146** and *Minois dryas* **157→47**. These species are deterministic source-extreme/hidden examples, **not experimental demonstrations** of biological range change.

**Scientific meaning:** the apparent *location* of host-driven novel opportunity is strongly sensitive to the completeness of the larval-host relation inventory, even when global totals change relatively little. This result is a classification/measurement insight, not a temporal ecological map of larval feeding or colonization, and the added BCE candidate hosts are Europe-biased rather than independent global validation. It is also not a statistical refutation of the earlier original-HOSTS matched null, since that null was not rerun under augmented links.

**Updated hierarchy for one ecology paper:** First distinguish original resource geography, duplicated new-cell supports and unique newly accessible regions; second demonstrate that the geographic identities of these "new" resources are catalogue-conditional; third report the **small** original-HOSTS geography/climate-matched excess native overlap as a bounded supporting result. [The 250-word GEB structured abstract and proposed revised Results](BUTTERFLY_GLOBAL_RESOURCE_NOVELTY_PAPER_STORY_20261016.md) now state the 50% source-Jaccard outcome explicitly. The two-panel figure plots species-level gain count against augmented count and partitions 6,445 source-union cells into retained, reclassified and newly discovered cells. **Do not convert this into a second paper or claim demographic consequences.**

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

The frozen figure source `data/frozen/figure_sources/host_contribution_metrics.csv` independently supports a **within-species / local-redundancy distinction**:

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

## 2026-10-15 empirical plant-range null: result hierarchy changed

The initial [successful plant-family × Level1 range-shuffle run](https://github.com/zuizui0223/chocho/actions/runs/38681928507) is preserved **without modifying its outputs** as [v0.1 frozen result](../../provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_v0.1.json). It used **199** post-hoc swaps-chain draws, seed 20261015, holding recorded butterfly-host links, each plant's introduced range size per Level1, and the number of introduced plants of each family in each WGSRPD3 region. Of 30,373 introduced plant×region incidences, 28,670 belonged to viable swap groups; the burn-in changed 31.1% of total incidences relative to the initial arrangement. *This is a movement check, not proof of sampling convergence.*

| Effect (different geographic opportunity units) | Observed | Structural null median | Central implication |
|---|---:|---:|---|
| Newly supported butterfly × region cells | **14,553** | **15,383** (95% null interval 15,311–15,468) | Fewer than expected, 830-unit deficit; one-sided Monte Carlo p=0.005 |
| New exact shared-host butterfly-pair × region cells | **79,412** | **81,550** (95% interval 81,290–81,821) | Fewer than expected, 2,138-unit deficit; p=0.005 |
| Fraction of newly supported cells with one introduced supporting host | **0.5892** | **0.5924** (95% interval 0.5875–0.5968) | No resolved unusual single-host dependency; lower-tail p=0.12 |
| New exact shared-host pair×region / newly supported butterfly×region | **5.4567** | **5.3013** (95% interval 5.2796–5.3186) | Positive post-hoc cross-unit ratio residual; upper-tail p=0.005, **NOT mean co-users per newly added cell** |

The number of added resource cells is not fixed in this plant-level null, which makes the first row a meaningful but bounded comparison. The first two rows disprove the intended simple claim that host redistribution produces *more absolute shared resource geography* than expected by the plant range and link architecture. The third demotes any claim that 58.9% is unexpectedly fragile. The fourth is an **interpretation trap**: a region in which two butterflies already had *different native host plants* may acquire an introduced plant shared by both without either butterfly acquiring any *new resource region*. Therefore its two terms differ in incidence universe. It is a descriptive contrast rather than a clean competition or co-use density estimand.

**Replacement outcome and decomposition: completed two-chain check (v0.2).**

The [follow-up source workflow](https://github.com/zuizui0223/chocho/actions/runs/38685083853) successfully reconstructed the same HOSTS/WCVP resource data and ran *two independent 499-draw RNG chains*, seed 20261015 and 20261016, with doubled between-draw swap spacing. Unaltered JSON receipts and both chains' diagnostics are frozen at [plant-range structural null v0.2](../../provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_two_chains_v0.2.json). Burn-in changed **31.1%** and **30.7%** of introduced plant×region incidences relative to initial placement, respectively; this checks movement, *not* convergence.

**Key ecological partition and exact accounting identity:**

Let `H` be fixed introduced-host × known-focal-butterfly × region incidences, `A` those incidences already inside the butterfly's *native*-host resource envelope, `N` incidences in newly supported cells, `R` redundant supporting host incidences beyond one per new butterfly×region cell, and `U` unique newly supported butterfly×region cells. By construction, `H = A + N` and `U = N - R`. The plant-level null keeps `H` fixed.

| Quantity | Observed | Null median: seed 20261015 / 20261016 | Direction |
|---|---:|---:|---|
| `H`: total candidate incidences | **75,346** | 75,346 / 75,346 | Fixed exactly |
| `A`: absorbed within native-supported resource cells | **46,387** | 45,301 / 45,298 | **+1,086 / +1,089** more absorbed |
| `N`: candidate incidences supporting new cells | **28,959** | 30,045 / 30,048 | −1,086 / −1,089 |
| `R`: extra introduced host support within novel cells | **14,406** | 14,662 / 14,659 | −256 / −253 |
| `U`: unique added butterfly×region cells | **14,553** | 15,384 / 15,385 | **−831 / −832** (−5.4% vs median) |
| Mean introduced host supports *per newly added cell* | **1.9899** | 1.9532 / 1.9529 | Small positive difference |
| Mean other focal host-sharing butterflies *per newly added cell* | **5.7126** | 5.7483 / 5.7472 | Small negative difference |
| Single-host-backed fraction of new resource cells | **0.5892** | 0.5921 / 0.5920 | No clear unusual deficit or excess |

The first strong unexpected **conditional** result is the *native-envelope absorption surplus*: observed 46,387 vs null about 45,300; observed fraction **61.6%** of the fixed candidate incidence budget versus null **60.1%**, a difference of about **1.4 percentage points**. Both chains place the observed absorption and 5.4%-lower added-cell count outside all 499 respective random draws (nominal one-sided Monte Carlo `p = 0.002` per chain; not two independent 499-permutation hypothesis tests). Multiple exploratory endpoints were inspected, so do not turn this into a prospective global confirmatory p-value.

**Mechanistic decomposition:** Compared with matched plant-range permutations, more introduced host incidences fall within already-native-supported butterfly geography, reducing their ability to open *new* butterfly resource cells. The **lower** absolute redundant incidences within the newly opened cells partly offset, rather than cause, the observed 831–832-cell deficit. This isolates **native-resource saturation/overlap** as the main arithmetic contributor under this null, *not* extraordinarily high introduced-host redundancy.

**Equally important negative results:** 58.9% single-host-backed cells are *not unusually dependent* under this null; the mean number of other focal butterflies sharing exact hosts in newly added cells is slightly **lower**, not higher, than null (~5.713 observed vs ~5.747–5.748). Geographic shared-host pair×region counts increase from native to contemporary, but their absolute additions are **below** plant-range null expectation. Therefore neither “novel opportunity is unusually fragile” nor “more interspecific concentration than expected” is supported as a global ecological effect here.

**Hard inference ceiling:** This source-conditional excess absorption is plausible under common plant/insect climatic and geographic filtering, unequal habitat, historic introduction routes, or unevenly recorded host links. The null preserves plant family×Level1 margins and link identity but *not* plant climatic preferences, butterfly-specific added-row margins, or time since introduction. It cannot assign cause to ecological sorting, local host use, or the success/failure of colonization; confirmation would need a climate-conditioned constraint or independent observed larval-use data. The GEB submission paper remains frozen until the user chooses whether a new biological framing is warranted.

### Existing literature that bounds novelty

- Fricke & Svenning (2020), *Nature*, DOI [10.1038/s41586-020-2640-y](https://doi.org/10.1038/s41586-020-2640-y), had already shown global biogeographic homogenization of **observed plant–frugivore interaction networks**, so we cannot claim ecological interaction homogenization as an entirely new concept.
- A global 74-network analysis in *Science of the Total Environment* (2024), DOI [10.1016/j.scitotenv.2024.173949](https://doi.org/10.1016/j.scitotenv.2024.173949), examines exotic plants and plant–herbivore network topology; our distinct question is **potential larval resource geography at fixed butterfly–host identities**, with no measured local network topology.
- Galmán and colleagues (2025), *Ecology Letters*, DOI [10.1111/ele.70196](https://doi.org/10.1111/ele.70196), compare actual herbivory across native and non-native plants. This reinforces that resource presence cannot, by itself, establish realized local usage or successful butterfly development.

## Non-obvious pairwise result: expanding shared-host geography can dilute normalized niche overlap

An independent historical workflow [October 7, source run 37612135885](https://github.com/zuizui0223/chocho/actions/runs/37612135885) computed pairwise Jaccard on **exact host-plant ID × region resource cells** (not on geographic envelope alone). Its original numbers are preserved in [an auditable derived receipt](../../provenance/reviewer_defenses/results/butterfly_exact_host_pairwise_overlap_v0.1.json).

Across the **986 pairs with an exact shared host relationship** (from the same 239-butterfly panel), **616** pairs increased normalized exact-host resource-niche overlap, **366** decreased, and **4** did not change. Of those 986 pairs, **920** gained at least one new exact shared-host × region cell after plant redistribution. By set intersection alone, **at least 920+366−986 = 300 butterfly pairs** must therefore *gain shared-host resource geography while losing normalized pairwise niche similarity*. This is a rigorous lower bound; the exact cross-tab has not been independently extracted, so do not report 300 as an exact observed count.

A complementary descriptive comparison conditions the all-pairs means on the 986 sharing-host pairs, yielding mean Jaccard approximately **0.1163 → 0.1409**, versus all-pairs **0.00403 → 0.00489**. Denominators must be identified: **27,455** pairs have no exact shared host and necessarily contribute zero throughout.

**Ecological significance:** increasing the *number or geography* of resources shared by two consumers does not imply a uniformly increasing *proportion* of their resource niches that is shared. This supports explicit separation of a species-pair opportunity count from niche overlap, although changing normalized denominators explains part of the discordance and no tailored null has yet been used to establish whether the 366 declines are unexpectedly frequent.

## Direct correction of the cross-unit normalization

The same plant-range null, now with a **properly aligned butterfly×new-resource-region denominator**, ran successfully from [source run 38684924181](https://github.com/zuizui0223/chocho/actions/runs/38684924181); its unchanged numbers and source are frozen in [aligned-cell v0.1 result](../../provenance/reviewer_defenses/results/butterfly_added_cell_peer_null_v0.1.json).

- Observed mean number of **other focal butterflies sharing an exact host per newly added butterfly×region cell** = **5.71257**; null median **5.74938** (95% null interval 5.71275–5.78614), lower-tail exploratory Monte Carlo p=**0.030** (199 draws).
- Observed fraction with **no other focal butterfly** sharing an exact host = **0.09895**; null median **0.09973**, no clear residual (p=**0.205** lower tail).
- Therefore the earlier cross-unit ratio 5.4567 vs 5.3013 **must not be interpreted as a positive co-user concentration result**. When units align, the direction is slightly *negative*. The magnitude is small (about **−0.64%** relative to the null mean), and the completed two-chain 499-draw comparison reproduces the same slight negative direction (nominal p=0.032 and 0.028, before post-hoc/multiplicity cautions).

**Current claim decision:** Reject a global *excess co-user concentration* narrative and the inference that the 58.9% single-host-supported fraction is unusually high. Retain these only as descriptive structural properties. The absolute new opportunity deficit and its partition into native-resource overlap versus added-host redundancy remain the primary candidate for a testable, non-mechanical macroecological result.

## 2026-10-15 climate-conditioned null: residual native-resource absorption survives

**Empirical computation completed:** [GitHub Actions run 38696007445](https://github.com/zuizui0223/chocho/actions/runs/38696007445), with complete original JSON artifact and [compact frozen receipt](../../provenance/reviewer_defenses/results/butterfly_climate_matched_absorption_v0.1.json). The [post-hoc climate protocol](../../provenance/reviewer_defenses/protocols/CLIMATE_CONDITIONED_RESOURCE_ABSORPTION_20261015.md) was committed before this result was examined. This does **not** make the climate analysis prospectively confirmatory.

The source-independent TTF climate run 36270581743 provided CHELSA V2.1 1981–2010 BIO1/7/12/15 region climate vectors. **5,300** repeated species rows resolved to **365 unique regions**, covering **all 365** resource regions with no missing climate; repeated region vectors agreed (climate input SHA256 `0087117e0915a0847abaae8d11762b75f35b8eac94434d50913c2781bea727e2`).

Each climate-matched swap retains the exact recorded butterfly–host network, all native host ranges, each plant's introduced-region count **within its continent and climate stratum**, and each destination region's introduced-host count **within plant families**. Two separate seeds each generated 299 serial swap-chain null worlds. All three controls passed predefined climate coverage, viable-edge and swap-movement checks (81–92% of 30,373 introduced plant×region edges in viable groups; burn-in moved ~14–22% of all introduced plant-range edges).

### Key absorption result

`H = 75,346` focal introduced host×butterfly×region incidences, `A_observed = 46,387` incidences already in native-host-supported butterfly resource geography, `U_observed = 14,553` unique newly opened butterfly×region cells.

| Structural null | Median A in two chains | Observed excess of A | Null median U in two chains | Observed U deficit |
|---|---:|---:|---:|---:|
| Plant family × Level1 (prior baseline) | 45,301 / 45,298 | **+1,086 / +1,089** | 15,384 / 15,385 | −831 / −832 |
| Family × Level1 × CHELSA BIO1/BIO12 median bins (2×2) | 45,835 / 45,829 | **+552 / +558** | 14,979 / 14,984 | −426 / −431 |
| Family × Level1 × BIO1/BIO12 tertile bins (3×3) | 46,099 / 46,096 | **+288 / +291** | 14,817 / 14,813 | −264 / −260 |
| Family × Level1 × median bins BIO1/7/12/15 (2×2×2×2) | 46,006 / 46,007 | **+381 / +380** | 14,832 / 14,832 | −279 / −279 |

Across all six climate-matched chains, the observation exceeds null absorption and falls below null new-opportunity count. In each 299-draw chain no null sample reached the observed A (nominal one-sided Monte Carlo **p = 0.0033 per chain**, the minimum resolution, **not** a corrected confirmatory inference). Short-lag correlation in the A null draws ranges approximately 0.12–0.25, indicating serial dependence; **the 598 draw rows per mode are not treated as 598 independent random permutations**.

### Ecological interpretation and scope

The earlier apparent `+~1,090` native-envelope absorption is **partly attenuated** when coarse climate categories are preserved; the residual `+288–558` corresponds to about **0.38–0.74 percentage points of H**, compared with ~1.44 points in the original null. Therefore basic climatic sorting is a relevant alternative explanation but does not eliminate the geographic overlap residual in these particular constrained rearrangements.

The **arithmetically distinct** new-resource-area deficit is ~260–431 butterfly×region cells (about 1.8–2.9% of the matched null U median) after climate stratification. In all three climate comparisons, the absolute added shared-host pair×region count also remains *below* structural expectation, and single-host-supported share is close to null; neither exaggerated interspecific concentration nor unusually fragile new opportunity is established.

**Important:** Climate strata are not nested. Comparing their residuals is not a valid estimate of a causal percentage 'explained by climate', or proof that four variables are always more rigorous than two. Regional medians and quartile/tertile matching do not control continuous climate preferences, land-use, plant introduction dates, spatial autocorrelation, phylogeography, sampling bias, or site-level host coexistence.

**Recommended scientific framing:** *Anthropogenic plant redistribution reinforces already available geographic larval-host opportunity more than expected under constrained plant-range and broad climate allocation, while also adding new resource regions.* This is a **potential resource-geography structure**, not colonization, competition, host preference, or demographic benefits. Retain the exploratory status and do not silently replace the separate frozen GEB manuscript.

## Final numerical qualification: 1,500-km climate-and-geography-tethered null

The original spatial null [run 38699068001](https://github.com/zuizui0223/chocho/actions/runs/38699068001) was formally **demoted** because its 199 serial draws had substantial lag-1 autocorrelation (~0.78 for the 1,500-km primary endpoint), only ~9–10 approximate effective observations per seed, and evidence of nonstationarity. The initial spatial null's one-sided `p=0.005` is **not** defensible as an exact permutation p-value.

A **single finite qualification attempt**, with a frozen post-hoc [mixing protocol](../../provenance/reviewer_defenses/protocols/SPATIAL_NULL_EXTENDED_MIXING_GATE_20261015.md), was successful in [run 38700697916](https://github.com/zuizui0223/chocho/actions/runs/38700697916). We used the previously declared 1,500-km original-introduced-footprint-tethered geography plus plant family × Level1 × CHELSA annual-temperature/precipitation classes. Two *new*, independent-seed chains used 30×E burn proposals, 20×E spacing proposals between worlds and 99 recorded worlds each. The chain-specific diagnostic thresholds were satisfied (primary lag-1 = 0.105 and 0.064; approximate ESS = 84.6 and 75.0; early-versus-late quarter and between-chain median checks also pass). This is a **numerical reliability gate**, not prospective ecological hypothesis confirmation or proof of MCMC ergodicity.

The [qualified frozen result](../../provenance/reviewer_defenses/results/butterfly_spatial_null_longchain_qualified_v0.2.json) yields:

| Same fixed 239-butterfly resource reconstruction | Observed | Two null-chain medians | Conditional effect |
|---|---:|---:|---:|
| Introduced host×butterfly×region incidences, `H` | 75,346 | 75,346 / 75,346 | Exact conservation |
| Incidences falling within each butterfly's native resource envelope, `A` | **46,387** | 46,116 / 46,111 | **+271 / +276**, only ~**0.36 percentage points of H** |
| Distinct added butterfly×region resource cells, `U` | **14,553** | 14,817 / 14,820 | **−264 / −267**, about **1.8% below null** |

This is a **small residual coincidence** under the stated constraints, not evidence that plant transport targeted butterflies or that added hosts are locally used. Unmeasured land use, region topology, plant climatic niches within coarse bins, host-link reporting intensity and local environmental differences remain. Crucially, the null preserves each plant's introduced-area margins but **does not** fix every butterfly's total potential new regions. Independent host-use and field records are needed to infer any demographic ecological mechanism.

**Practical manuscript decision:** Keep the central nontrivial question as *novel vs overlapping resource geography*, but do not inflate the ~0.36-percentage-point spatial deviation into a large biological phenomenon. The 61.6% native overlap and the bookkeeping partition describe how the **candidate host incidence budget** is arranged; they are not a measured proportion of ecological benefit. Do not reuse the invalid short-chain p-values.

## What NOT to promote as findings

1. **“More host-plant range = more butterflies.”** Only potential larval resource opportunity is reconstructed; butterfly density, recruitment and conservation benefit were not estimated.
2. **“Butterfly communities have homogenized by 66.8%.”** This number is potential resource-assemblage similarity, not observed composition; nearly all the raw gain is fixed-margin-predictable.
3. **“Introduced hosts intensified competition.”** Exact shared-host exposure is a necessary opportunity for interaction, not observed competition.
4. **“Plant network prominence is a proven biological cause.”** The HOSTS degree measure carries unknown recording bias, and identity-null conclusions depend on weighting.
5. **“Diet breadth has exactly no influence.”** The near-zero point estimate does not pass the declared post-hoc tight equivalence diagnostic.
6. **“The global temporal host–larva mismatch is proven.”** WCVP status has no introduction dates or seasonal edible-leaf trajectories. Kyoto's native-host/introduced-butterfly cage system is a distinct mechanistic test bed, not a direct global introduced-host validation.
7. **“Macro-data identify whether anthropogenic hosts produce self-sustaining butterflies.”** Presence/resource overlap falls short of that stronger demographic question.

## Single-paper decision: elevate geographic novelty, not a second post-hoc paper

**Recommended question:** *When host plants are redistributed, where do butterflies gain genuinely new potential larval-resource geography, rather than another candidate host in a region they could already exploit?*

**Working title:** *Geographic novelty and redundancy of butterfly larval-host resources under plant globalization*.

The earlier global plant/consumer network homogenization literature already establishes that introductions alter interaction geography. Our distinction is a **fixed-interaction** butterfly resource reconstruction that partitions 75,346 introduced-host candidate incidences into (a) already-native-covered locations, (b) multiple introduced supports in newly available cells, and (c) 14,553 distinct butterfly×region gains. The additional small overlap beyond constrained nulls is a **secondary** result: it is robust to the named broad climate and spatial rearrangements, but source-conditional and post-hoc.

The companion [single-paper author draft](BUTTERFLY_GLOBAL_RESOURCE_NOVELTY_PAPER_STORY_20261016.md) contains one candidate English title, ~200-word abstract, biological questions, Results order, Discussion and appropriate inferential limits. **It is an alternative to the existing GEB presentation, not an extra manuscript to submit in parallel.** Author approval is needed before revising PR #38.

### Revised figure sequence

- **Fig 1 — resource accounting and analysis scale.** Show **239 butterflies, 365 WGSRPD3 regions and H=75,346 host×butterfly×region candidate incidences**; partition `H=A+R+U` (46,387 / 14,406 / 14,553), explaining that U is one first supporting incidence per distinct new butterfly×region cell.
- **Fig 2 — matched expectations with effect sizes.** Plot observed-minus-null native resource absorption and null-minus-observed new regions under plant-family/continent margins, + climate, and + numerically qualified 1,500-km geography. Use two visible chain-specific medians, and explicitly label post-hoc status. The [reproducible renderer](../../scripts/render_butterfly_resource_novelty_figure.py) and its [input-integrity test](../../tests/test_butterfly_resource_novelty_figure.py) create the combined exploratory evidence illustration.
- **Fig 3 — ecological heterogeneity and external check.** Taxonomic host-family count versus native resource geography (18 geographically wide one-family specialists in the more adequately documented subset), plus the secondary occurrence panel (66/115 outside-native butterfly presence cells aligned with resource gain vs region-matched structural null median 49). Mark 32 selected species, 23 informative for the occurrence summary, sampling and host-use ambiguity.
- **Supplement — other potential-resource structure.** Native/current resource-assemblage Jaccard 0.277→0.462 versus fixed-margin null 0.453, with small continental excess ~0.0036. Keep exact shared-host co-opportunity and 58.9% single-host support **descriptive**; neither excess competitive concentration nor excess fragility is supported by the corresponding plant-range null. Explain the >=300 pair discordance between new exact-host sharing geography and declining normalized host×region niche overlap.

These are exploratory storyboards, **not silent replacements** for the frozen GEB Figures 1–3.

### Stop rule and reviewer priorities

The partial **independent-host-inventory sensitivity has now been completed** with a Europe-derived BCE/Clarke checklist, but globally representative completeness and verified local larval use remain unmeasured; do not search for new significance across dozens of habitat covariates. The BCE host-list augmentation already showed that resource-envelope numbers change with additional literature-derived host associations (relative global expansion 54.9%→50.1%), but the new **native-resource absorption residual has not been rerun under a globally representative alternative host inventory**. State this as a source-dependence limitation.

The spatial null controls climate classes and polygon representative-point proximity, **not** elevation, fine-scale land use or introduction timing. Changing those one by one in the same data in pursuit of a lower p-value would not identify a mechanism. A direct field or host-rearing comparison is required to tell whether an extra mapped host region produces additional butterflies, and therefore lies beyond this paper's evidence.

## External position and interpretation limit

- [Fricke & Svenning (2020), *Nature*](https://doi.org/10.1038/s41586-020-2640-y): introduced species homogenize observed plant–frugivore network interactions; this topic is **not** an unoccupied novelty gap.
- [Araújo et al. (2024), *Science of the Total Environment*](https://doi.org/10.1016/j.scitotenv.2024.173949): exotic plants and land use alter some aspects of observed plant–herbivore network structure; our mapped larval-resource geography is a **different response variable**, not an observed interaction network.
- [Singer & Parmesan (2018), *Nature*](https://doi.org/10.1038/s41586-018-0074-6): a single butterfly/introduced-host case informs hypotheses about ecological traps and dependence but **does not** demonstrate worldwide trap effects from host maps.

**Recommended final claim ceiling:** *Introduced host distributions frequently overlap resource geography already supported by recorded native larval hosts, and the observed configuration adds modestly fewer new butterfly×region resource units than carefully specified climate- and spatially constrained alternative configurations. This is not a demonstration of realized host use, demographic success or human-caused evolutionary specialization.*
