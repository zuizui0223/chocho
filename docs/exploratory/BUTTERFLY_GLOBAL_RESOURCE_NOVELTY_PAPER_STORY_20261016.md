# Manuscript-ready ecological narrative — introduced larval-host geography

**Status (16 October 2026): ALTERNATIVE SINGLE-PAPER FRAMING / DISCUSSION ONLY.** This file proposes a revised scientific narrative for the existing 239-butterfly GEB paper; it is **not** a second manuscript, a change to the frozen submitted draft, a preregistered hypothesis or proof of population establishment. The earlier public manuscript and original sources remain unchanged.

## Biological question and novelty

**Question:** When humans move the known larval host plants of butterflies into new regions, how much of that geographic resource supply is genuinely *new to the butterfly*, and how much falls within the butterfly's already available native-host resource envelope?

The natural expectation “more host plant range produces more mapped opportunity” is essentially set-theoretic. The nontrivial test is whether the introduced host×consumer×region arrangement is **geographically more coincident with existing native resource access than a matched spatial rearrangement would produce**, and whether that remains true when plant taxonomic group, destination regional richness, broad climate and geographical proximity are constrained.

This explicitly concerns *resource geography*, not observed plant–larva edges in the introduced region or butterfly population demography. Before asserting novelty, distinguish from (i) plant–frugivore meta-network homogenization (Fricke & Svenning 2020, DOI 10.1038/s41586-020-2640-y), and (ii) local empirical plant–herbivore network topology under exotic plant prevalence (Araújo et al. 2024, DOI 10.1016/j.scitotenv.2024.173949). The distinguishing design here is keeping known butterfly–host link identities fixed and reconstructing **where potential resources overlap or expand**.

## Recommended English title and abstract

**Working title:** *Plant globalization often overlaps existing butterfly larval-resource geography rather than creating novel regional opportunity*

**Alternative (more neutral):** *Geographic novelty and redundancy of butterfly larval-host resources under plant globalization*

### Abstract — structured draft for GEB Research Article

**Aim:** Introduced host plants need not create geographically novel resources for consumers. We distinguished new butterfly larval-resource regions from overlapping candidate host supply and evaluated the sensitivity of this distinction to host-inventory completeness.

**Location:** Global; 365 WGSRPD3 geographic regions.

**Time period:** Native and introduced WCVP v13 plant distributions (without introduction dates); CHELSA 1981–2010 climate sensitivity.

**Major taxa studied:** 239 butterfly species and recorded larval host plants.

**Methods:** We joined accepted host–butterfly associations to native and contemporary plant regions and partitioned introduced host×butterfly×region incidences into native-supported, multiply introduced-supported and uniquely new butterfly×region resources. Plant-range nulls preserved host links, plant-family and regional margins, climate bins and local proximity. A separately curated European checklist added 1,027 candidate host links to 81 butterfly species.

**Results:** Of 75,346 introduced-host candidate incidences, 46,387 (61.6%) occurred inside butterfly native-host resource envelopes; 14,406 (19.1%) added support within new cells, leaving 14,553 distinct regional opportunities (19.3%). Under a numerically qualified 1,500-km climate-matched null, native-resource overlap exceeded expectation by just 271–276 incidences (0.36 percentage points). Host-list augmentation changed the native-overlap share to 65.5% and uniquely new share to 14.5%. Among 81 augmented butterflies, only 3,225 of 6,445 resource cells classified as new under either catalogue agreed between catalogues (Jaccard 0.500), despite a net increase of only 414 cells.

**Main conclusions:** Potential larval-host redistribution often overlaps pre-existing regional opportunity. More importantly, which places qualify as *new* depends strongly on recorded host relations; mapped supply does not establish realized local feeding, colonization or fitness.

**Word-choice caveat:** “Introduced hosts therefore mostly supplied geographically overlapping opportunities” describes the union-map reconstruction and incidence accounting. “More than expected” is *specific to the named post-hoc spatial randomization*. It must **not** imply that plants prefer butterfly-rich regions, or that individuals feed successfully on introduced plants.

## Main Results: proposed order

### 1. Widespread larval-resource increase masks considerable overlap

**Panel scope:** 239 butterflies; 365 geographically delineated regions; fixed HOSTS records and WCVP v13 native/introduced statuses; 670 plant species supply at least one distinct new resource unit.

Observed native butterfly×region resource cells **26,530**; contemporary **41,083**, an increase of **14,553 (+54.9%)**. Species gaining any extra region: **206/239 (86.2%)**. These are the scale and context, **not** the lead surprise.

To unpack the geographic novelty of candidate introduced host distribution:

\[
H = A + N,\qquad N = U + R
\]

where `H=75,346` introduced host×butterfly×region incidences, `A=46,387` incidences falling in regions already within the butterfly's native larval-resource envelope, `N=28,959` incidences falling into novel butterfly×region resource cells, `U=14,553` distinct novel butterfly×region cells and `R=14,406` extra supporting host incidences beyond one per new resource cell. **H=A+R+U is a bookkeeping identity with one incidence counted as the first support of each new cell, not a ratio of population growth or resource biomass.**

The last two terms are not independent: highly overlapping introduced plants can supply multiple supporting incidences within one new cell. None of these counts measures edible plant biomass, insect acceptance or reproductive success.

### 2. Introduced ranges are more likely to coincide with native-resource geography than the constrained structural null would predict — but only slightly

The initial null preserves planted-host degree and continental plant-family region margins. Of the fixed 75,346 candidate incidences, observed native-envelope overlap is **46,387**, versus **45,301/45,298** in two 499-draw null chains (**+1,086/+1,089**).

CHELSA temperature/precipitation-bin matching reduces the excess to **+552/+558**. A more stringent *locally tethered* comparison (1,500 km, original introduced host footprint, the same plant-family × continent × 2×2 climate bins, two independent longer chains) retains **+271/+276**, or **+0.36 percentage points** of the fixed candidate incidence budget. Independent-chain simulation diagnostic thresholds were satisfied: autocorrelation lag-1 ≈ 0.06–0.10; estimated ESS ≈ 75–85 for the primary A response; chain median difference 5 units. Distinct new resource additions remain **−264/−267** lower than the respective null medians (**14,817/14,820**).

These are **small conditional differences** in a coarse mapped resource system and **have not been re-estimated under the augmented BCE host inventory**. The inventory-sensitivity change in the native-overlap share (+3.98 percentage points) is numerically much larger than the original-inventory local-null excess (~0.36 percentage points), but they are different contrasts and must not be compared as if they were exchangeable effect estimates. Long-chain numerical qualification does not convert a post-hoc comparison into a confirmatory permutation test or into identification of climate, geography or introduction mechanisms.

The short-chain local null (199 serial draws per seed) had lag-1 around 0.78 at 1,500 km and ESS ≈ 9–10, and its nominal p-values were rejected as evidentially unreliable. The longer-chain sensitivity was run once with a finite frozen mixing gate; results were stable, not chosen by repeated threshold fishing.

### 3. Specialization does not capture geographic novelty well

Across focal butterflies, the association between host-family breadth and proportional introduced resource gain is **rho=0.008 (bootstrap 95% CI −0.111 to 0.128)**. This does **not** prove a mathematically exact zero effect; a post-hoc narrow ±0.10 equivalence margin was not met (approximate TOST p=0.078). Host-family breadth and native geographic resource breadth are only partially correlated (**rho=0.276** full sample). In the more host-taxonomy-adequate 215-species subset, **18 one-family host specialists** nevertheless occupy the top quartile of geographic resource opportunity.

For the introduced-added cells themselves, **8,574/14,553 (58.9%)** are supported by a single introduced host, but the proportion is not unusually high relative to a plant-range null with relevant margins. More host species typically allow more redundant added opportunities; **taxonomic host-family breadth by itself does not capture that spatial portfolio structure**.

### 4. Host-inventory completeness materially changes what counts as new opportunity

In a **source-completeness sensitivity**, we re-used the previously archived BCE-derived version of Clarke (2024) rather than newly selecting host taxa by result. The source had 92 exact matching butterfly names, 83 with source foodplant records, and 1,027 *candidate* additional accepted-WCVP-ID larval-host associations for 81 focal butterflies (679 plant species). Candidate relations were added **together** to the frozen HOSTS list and the original WCVP v13 native/contemporary plant maps were rebuilt. The original 239-species baseline reproduced exactly.

Across all 239 butterflies, reconstructed *native* resource cells increased **26,530→29,852**, contemporary cells **41,083→44,819**, and net new resource cells **14,553→14,967**. The fixed-host-list candidate introduced host×butterfly×region incidence pool `H` changed **75,346→103,349**. Its partition changed substantially:

| Accepted-host inventory | Native-overlap `A` | Extra supports in newly opened cells `R` | Unique newly opened cells `U` |
|---|---:|---:|---:|
| Frozen HOSTS | 46,387 (**61.6%**) | 14,406 (**19.1%**) | 14,553 (**19.3%**) |
| HOSTS + 1,027 BCE candidate links | 67,740 (**65.5%**) | 20,642 (**20.0%**) | 14,967 (**14.5%**) |

Percentages are proportions of `H` **within each source inventory**, not independent observations, ecological production or intrinsic probabilities of plant utilization. While the new source records add **1,817** previously unsupported butterfly×region cells, they also reclassify **1,403** previously introduced-added cells as natively supported, giving a **net +414** unique additions. Among 81 augmented butterflies, **36 gain**, **35 lose**, and **10 do not change** in the number of regions *classified* as introduced-additional. This shows that aggregate geographic gain can appear stable while the regional classification within butterfly species shifts substantially.

**Geographic identity turns over much more than the total count.** Among the **81** butterflies with additional BCE candidate host links, the original and augmented inventories classified **4,628 and 5,042** butterfly×WGSRPD3 region cells as introduced-added. Their *intersection* contains only **3,225** cells, with **1,403** old additions reclassified as native and **1,817** newly identified introduced additions; the union contains **6,445** cells. Thus the *set-level identity agreement* of new-resource regions is only **Jaccard = 3,225 / 6,445 = 0.5004**, with almost half of the union classified differently. The per-species median Jaccard is **0.570** among the 77 butterflies with at least one new-resource cell under either catalogue; four have no such cells, for which Jaccard is undefined. Of all 81 butterflies, **54** simultaneously lose old introduced-new classifications through native reclassification and gain other newly supported cells. Do not report the high **full-239 Jaccard = 0.803** without noting that the other 158 butterflies have identical inventories by design.

Species cases make the ecological interpretation concrete without supplying independent proof: *Lycaena dispar* changes from **0 to 146** candidate new-resource regions, *Minois dryas* from **157 to 47** (125 old-new cells reclassified as native, 15 previously unsupported newly discovered), and *Thymelicus lineola* remains **65 to 65** but has **27** old cells reclassified and **27** different new cells added (set Jaccard **0.413**). These are **deterministically selected extreme/hidden-turnover examples**, not independently sampled effects or observed range shifts. All values were rederived from the [source-validated 81-species archive](../../provenance/reviewer_defenses/results/butterfly_bce_added_region_identity_turnover_v0.1.json).

The source sensitivity is **Europe-biased and only one-directional**, and BCE/Clarke evidence comes from documented European feeding/oviposition records, not independent worldwide larvae sampled at every WCVP locality. It cannot establish a globally corrected 65.5% or demonstrate that the earlier small *conditional spatial-null* excess remains after the added host links. Its importance is instead to quantify a real constraint on the **interpretation and completeness** of host-resource geography.

### 5. Evidence that resource opportunity is biologically relevant is partial, and conservation claims are deliberately bounded

A **secondary, non-prospective** 32-species climate-stratified GBIF panel included **115** butterfly outside-native resource-region occurrence units from 23 informative species; **66** coincided with introduced host-resource gains versus a geography/region-matched null median **49**. Species sensitivity and uneven GBIF coverage warrant careful reporting; the result is **not** proof of introduced-host larval use or successful colonization.

Among the reconstructed regions, regional butterfly **potential resource assemblage** Jaccard similarity increased **0.277→0.462** (+66.8%), but a null fixing both butterfly-specific and region-specific additions predicted about **0.453**; the residual is only **+0.0086**, falling to **+0.0036** under continental margin controls. It is inaccurate to transfer the raw +66.8% to *observed butterfly community homogenization*, or claim the new shared-host exposure caused competition. Likewise, 920 of 986 host-sharing butterfly pairs acquire new shared-host regions, even though **at least 300 pairs** must have gained shared-host geography while their normalized host×region overlap declined; these are distinct quantities.

## Discussion — central implication

The primary ecological distinction is *access versus additional access*, not area versus population growth. Human redistribution creates many new potential plant occurrences, but the same introduced plant can (a) lie in a region with existing native host resources, (b) supply a new region alongside several other introduced plants, or (c) be the first introduced host to supply an otherwise unmatched butterfly×region resource cell. Conventional resource-envelope expansion reports only the union, making these different ecological situations look alike. The data show that most candidate supply maps to (a), with the observed arrangement producing a modest excess of (a) over plant-range alternatives that also preserve broad climate and source-geographic locality.

The source-level change of +3.98 percentage points in the overlap fraction is also a clear warning that recorded host-link completeness materially affects this geographic accounting. We should lead the paper with the transparent distinction between existing access, duplicate support and unique new access; the +0.36-point old-source spatial-null deviation remains conditional and secondary until an independent source-aware null is obtained.

The remaining overlap could reflect spatial autocorrelation, climatic niches at finer scales, land use, introduction history, botanical recording, host-database incompleteness or correlated native plant richness. These are not statistically eliminated. Importantly, no population success or reduction in fitness can be inferred without local plant abundance, larval performance, resource phenology, and verified host use.

## Final figure story (proposed, *not* currently adopted in submission)

- **Figure 1:** Global analysis overview, exact 239/365/75,346 scope and the **A+R+U** resource novelty accounting. Native, added and contemporary species×region totals shown as scale, not as headline inference.
- **Figure 2:** Reconstructed observed-minus-null native-resource absorption and new-cell deficit under plant-family/continent, plus climate, plus qualified 1,500-km locality controls. Show each chain separately; source receipts and the method's ESS/gates in caption. This figure should make the decreasing residual visible rather than exaggerate a small effect.
- **Figure 3:** Butterfly heterogeneity and biological relevance: specialist versus geographic resource-breadth discordance, plus outside-native observed butterfly occurrence recovery against a region-matched structural null. State separate sample sizes and biases.
- **Supplement:** Source-inventory comparison **H=A+R+U** for frozen HOSTS and augmented European BCE/Clarke host records, including the 1,817 newly identified cells and 1,403 reclassified native cells; plus regional assemblage Jaccard and matched null, with negative-result interpretation of fragility and competition claims.

### Rendered exploratory illustration and candidate figure legend

The [source-checked figure renderer](../../scripts/render_butterfly_resource_novelty_figure.py) generates a PNG, an editable vector SVG and a JSON values file in the isolated [exploratory figure workflow](../../.github/workflows/butterfly-resource-novelty-figure.yml). The source integrity test ties the plotted values to existing frozen receipts. The resulting figure is **not yet a GEB main figure**.

**Candidate combined figure legend (not automatically inserted in the manuscript):**

*Figure X. Most introduced larval-host resource incidences overlap geographically with pre-existing butterfly resource opportunity.* **(A)** Of 75,346 introduced accepted-host-plant × focal-butterfly × WGSRPD3-region incidences, 46,387 (61.6%) fall within the butterfly's already native-host-supported geographic resource envelope. A further 14,406 (19.1%) are additional introduced-host incidences in butterfly×region cells already newly supported by at least one introduced plant; one first-support incidence per novel cell yields 14,553 (19.3%) distinct added butterfly×region opportunities. **(B)** Excess native-resource overlap, calculated as observed count minus null median for a plant-family/continent fixed-margin rearrangement, additional coarse CHELSA temperature/precipitation matching, and a numerically qualified 1,500-km geographically tethered rearrangement. **(C)** Deficit in distinct new butterfly×region resource units against the same null constraints (null minus observed). Bars are means of two chain medians and points represent the individual medians (499, 299 and 99 simulated worlds per seed for respective constraints). The local-distance chains passed a separately declared finite simulation-mixing criterion after an earlier less-well-mixed run was rejected. These are descriptive, post-hoc structural opportunity comparisons, not independent confirmatory permutation p-values or evidence of realized plant–larva interactions.

### Supplementary figure: within-butterfly geographic identity

The [separate audited source-specific identity-turnover workflow](../../.github/workflows/butterfly-bce-source-turnover.yml) generates a 2-panel figure and the full 81-butterfly CSV from the previously frozen source-inventory sensitivity. Panel A compares original versus supplemented regional *count* per butterfly and marks the deterministic largest positive, largest negative and no-net-change-but-reclassified examples. Panel B compares the **geographic identities** of newly supplied cells, not just their totals: 3,225 retained as new under both inventories; 1,403 formerly new now supported by a newly recognized native host; and 1,817 other newly supplied cells. Sum of union = 6,445, intersection/union = **0.5004**.

*Candidate Supplementary Figure Sy.* **Nominally stable gains in butterfly larval-resource geography can conceal source-dependent shifts in which regions are classified as new.** All 81 butterfly species received one or more accepted-ID candidate host links from the Europe-derived BCE/Clarke evidence-filtered checklist. The two source scenarios are not chronological vegetation changes or independent observations of larval feeding. The classification changes include recoding already-introduced candidate resource cells as natively available when additional host records enter the baseline. The original and supplemented resource-region sets coincide exactly for 158 unchanged-inventory species outside the displayed subset. Accordingly, the union-normalized discordance illustrated here pertains **only to the source-augmented 81 butterfly species**.

### Supplementary figure: host-inventory source dependence

The [separate source-checked figure generator](../../scripts/render_bce_resource_novelty_sensitivity_figure.py) draws **two comparable 100%-partitioned candidate host-incidence compositions**, one for original HOSTS and one for original HOSTS plus the published evidence-filtered BCE/Clarke candidate host links. Its source is [the independently rebuilt, frozen numeric receipt](../../provenance/reviewer_defenses/results/butterfly_bce_resource_novelty_inventory_sensitivity_v0.1.json). Bar colours are consistent across the two source scenarios, and each bar prints its **own changed denominator** to prevent interpreting equal widths as equal total incidence counts.

*Candidate Supplementary Figure Sx. Recorded host relationships change the estimated novelty of globally mapped butterfly larval-host resources.* The 239-butterfly original HOSTS reconstruction yields `H=75,346` introduced accepted-host×butterfly×WGSRPD3 incidences, partitioned into native-supported opportunities (61.6%), additional supports within new butterfly×region resource cells (19.1%), and first supports for distinct new cells (19.3%). Adding 1,027 BCE/Clarke literature-derived candidate accepted host links for 81 butterflies (92 exact focal binomial matches, 83 with source foodplant records) changes the denominator to `H=103,349` and the corresponding partition to 65.5%, 20.0%, and 14.5%. A total of 1,403 previously introduced-added cells become native-supported, while 1,817 previously unsupported cells are now introduced-added (net +414 unique additions). This is a **Europe-biased partial host-inventory stress test**, not a globally corrected estimate, observed local feeding, demographic effect or rerun of the original-HOSTS spatial randomization null. The sensitivity should not be confused with uncertainty around the 0.36-percentage-point original-inventory spatial excess.

Submission-ready title/abstract and plot storyboard are a *candidate author decision*. Do not change the GEB submission bundle or silently merge PR #43 into PR #38.

## Source hierarchy and checks

- All reconstruction and null inputs remain pinned to the exact commits and source artifacts documented in `docs/REPRODUCIBILITY.md`.
- Base geographic null: `provenance/reviewer_defenses/results/butterfly_plant_range_structural_null_two_chains_v0.2.json`.
- Climate sensitivity: `provenance/reviewer_defenses/results/butterfly_climate_matched_absorption_v0.1.json`.
- Numerically qualified 1,500-km local null: `provenance/reviewer_defenses/results/butterfly_spatial_null_longchain_qualified_v0.2.json`.
- Figure generator: `scripts/render_butterfly_resource_novelty_figure.py`.
- Plant–frugivore global interaction homogenization: Fricke & Svenning (2020), *Nature*, DOI 10.1038/s41586-020-2640-y.
- Exotic plants and empirical plant–herbivore network topology: Araújo et al. (2024), *Science of the Total Environment*, DOI 10.1016/j.scitotenv.2024.173949.

**Final stop rule:** There is no need to keep searching new predictors until a larger post-hoc result appears. The supported conclusion is source-conditional, effect-size-qualified resource geography. If the paper requires a mechanism, it needs *independent plant/life-stage or temporal evidence*, not more re-shuffled versions of the same host map.
