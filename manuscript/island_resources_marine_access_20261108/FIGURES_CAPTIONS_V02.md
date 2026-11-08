# Figure architecture and source-safe captions — JBI integration v0.2

**Status:** *Figure specifications and submission-ready draft captions, not a claim that all five final figures exist.* The uniform-marine-current map in Fig. 4A exists, was visually inspected, and is reproduced by successful [PR #150 CI](https://github.com/zuizui0223/chocho/actions/runs/42439384640). Remaining composite panels require rendering and inspection.

The paper uses **five** figures, not five disconnected analyses. The logical progression is: resource gain → botanical assembly at ports → island pair composition and marine topology → can the apparent maritime effect be distinguished from geographic structure? The *resource* and *access* layers are intentionally distinct.

## Figure 1 — Resource versus access: conceptual framework and genuine sampling units

**Panels:** (A) two-layer conceptual diagram: candidate larval-host plant occurrence on islands vs potential over-water access, without arrows asserting actual plant or butterfly transport; (B) source-unit flow: 126 strictly insular and 237 mainland WGSRPD3 regions (six mixed excluded), 123 resource-active insular units, 98 botanical-effort-complete port-flora units; (C) separate original network sample: 51 WPI port-indexed landmass groups, 194 original ports, 182 qualified marine terminals, 1,275 non-independent pairs.

**Caption:** *Figure 1. Resource and access are separately measured layers of island butterfly biogeography.* The resource layer indexes the existence of documented larval-host plants for a fixed 239-butterfly panel within WGSRPD3 botanical regions. The access layer compares a separately assembled original 225-candidate butterfly panel across 51 port-indexed island/landmass groups. Islands and island archipelagos in the botanical data are not equated with marine network terminals, and none of the arrows represents observed shipments or colonisation.

**Audit gate:** Never combine the 239 and 225 source panels into one denominator. Do not portray a port as a known plant arrival route.

## Figure 2 — Disproportionate potential host opportunities in island botanical regions

**Panels:** (A) native-only and contemporary butterfly × botanical-region cells for island regions: 4,829→10,570 (+5,741, **+118.9%**) and mainland regions: 21,215→29,901 (+8,686, **+40.9%**); (B) contemporary added-resource shares: islands 54.3%, mainland 29.0%; (C) **matched** Level 1 botanical area-strata share difference **+10.7 pp** (0.107; descriptive permutation p=0.0002), distinct from broad Level 1-only +14.4 pp; mark the **ten** insular botanical regions with zero native-only opportunities in the fixed panel; (D, optional small source-coverage inset) **13/920 (1.413%)** candidate added island opportunity cells with an overlapping butterfly record in an independently archived **31-recorded-species GBIF sample** (32-species validation panel), DOI 10.15468/dl.pp5nc9.

**Caption:** *Figure 2. Introduced host-plant geography disproportionately expands reconstructed larval-host opportunity on islands.* Values count possible butterfly-by-WGSRPD3-region combinations, not butterflies, plant biomass or locally proven feeding. Pooled proportional increases use native opportunities as denominator, whereas added-resource shares use contemporary opportunities. The 10.7-point contrast is a weighted island–mainland comparison restricted to shared broad-biogeographic × area strata; the permutation value is descriptive under geographical dependence. Ten source units with no native-panel opportunity are botanical regions, not verified host-free physical islands. The optional archival validation inset records only whether a butterfly occurrence was indexed in a selected sample of 920 added cells (13 had records); the remaining cells are **not** confirmed ecological absences and 1.4% is **not** an estimated biological realisation rate.

**Audit gate:** Do not show a +14-point effect as the stricter adjusted result; do not call potential expansion observed abundance. Do not call 13/920 the fraction of all 239-species opportunities colonised.

## Figure 3 — Ports track entire introduced floras, without additional positive host-specific enrichment

**Panels:** three separate common-effort (n=98 insular WGSRPD3 regions) model coefficients with labels specifying *different responses*:
- Port indexed count → all introduced vascular plants: **β +0.2445**, descriptive p 0.0012.
- Same exposure → distinct introduced larval-host plants: **β +0.1639**, descriptive p 0.011.
- Port count → distinct introduced hosts **conditional on all introduced plants**: **β −0.0415**, descriptive p 0.093.

Plot intervals **only if computed from saved model/permutation outputs**; do not invent confidence bars. If intervals cannot be reproduced, show labelled point values and p annotations without uncertainty bars, and explain their descriptive nature.

**Caption:** *Figure 3. Contemporary indexed port density corresponds to an entire introduced flora that includes candidate butterfly host plants.* Each estimate controls for broad botanical geography, area, the corresponding native plant/resource baseline and all-plant GBIF indexing effort. The third estimate also conditions on total introduced flora; disappearance of the positive host association is not a formally identified causal mediation or selective transport test. The beta coefficients are on log(1 + species richness) response scales per log₂(1 + port count), not the number of butterflies per port.

**Audit gate:** No “port selectivity against hosts” causal interpretation from the negative conditional coefficient.

## Figure 4 — Island host-community similarity, marine-current geometry and complementary prediction

**Panels:** (A) the *existing*, verified uniform-ocean current map; (B) original island pair host versus butterfly Jaccard correspondence with later **ρ +0.2998**, earlier **ρ +0.4178** and descriptive partial later **ρ +0.2154**; (C) later-butterfly LOIO MSE improvement for hosts beyond the marine model **+1.7312%** and AIS/wind marine predictors beyond the host model **+6.4880%**. If space permits, annotate the separate **both-endpoints-unseen** check (**+0.8808% hosts**, **+4.4744% marine**) as a sensitivity in a subordinate panel or Supplement; do **not** draw them on the same reference-error axis without clear labels.

**Caption:** *Figure 4. Candidate host composition and modelled marine access describe overlapping but distinct parts of island butterfly biogeography.* (A) Uniform electrical current is calculated on a GSHHG-derived 0.5° sea graph (189,550 nodes; 370,996 undirected edges) using 182 marine-qualified WPI port terminals across 51 island groups. One unit of current is passed between each of 1,275 unordered island pairs. Mean absolute edge currents are accumulated to water nodes for display; bright corridors **are not observed shipping, insect movements or stepping-stone hubs**. The model has no through-island transport edges. (B) Original host and butterfly indexing similarities across 1,275 dependent island pairs; later records denote 2011–2022 first indexing rather than demonstrated arrival. (C) Exploratory 51-fold leave-one-island-out conditional prediction of later butterfly Jaccard with a shared joint model; improvement is percentage decrease in held-out mean squared error relative to *each model's own* reference. Hosts pool 2001–2024 and can postdate butterflies; all pair/fold results are dependent.

**Panel A rendered source:** [successful original pipeline](https://github.com/zuizui0223/chocho/actions/runs/42439384640), artifact `island51-figure4-uniform-marine-current-density`. Do not copy its map into Git as though this specification had assembled panels B–C. The previous conditional-model map/plots already exist in `results/island51_host_circuit_atlas_20261106/` and `results/island51_host_circuit_joint_20261106/`.

**Audit gate:** Panel A is current concentration averaged over unit-current endpoint pairs, with over-sampled source/sink geography. It does not identify actual intermediate islands. Never relabel high-current coastline pixels “relay hubs.”

## Figure 5 — Does apparent maritime accessibility identify ships?

**Panels:** (A) later butterfly double-endpoints-unseen latitude-adjusted effects: AIS only **+1.0020% MSE**, wind only **−4.6618%**, both **−4.0918%** under explicitly different model comparison labels; (B) true AIS compared with 199 latitude-band relabeling placebos (**5/199** at least as predictive) and the *more geographically constrained* 199 local ≤750 km placebos (**100/199** at least as predictive; median +1.0600% versus true +1.0020%); (C) historical 2006–2010 PLSCIQ on 17 island groups/136 pair predictions, full cohort **−2.9233% MSE**, prespecified 16-group sensitivity **+8.1461%** with training-set composition caveat.

**Caption:** *Figure 5. Maritime predictors contain geographically contingent retrospective information but do not establish ship-mediated transport.* Modern AIS vessel-position density sometimes improves later butterfly index-similarity prediction after latitude controls. That improvement is not exceptional under a close-distance graph-label placebo, and earlier port connectedness indices do not improve the full 17-island historical-sensitivity target. Placebo ranks are *not* causal p-values; the two null constructions preserve different aspects of geography. The 2006–2010 PLSCIQ source is an index of liner connectedness, not observed voyages or plant cargo. Comparisons involve different populations or baselines and must not be plotted as directly additive percentages.

**Audit gate:** Emphasize the positive absolute AIS model increment alongside its negative specificity test, not one to the exclusion of the other.

## Supplementary and exclusions

S1. Source taxonomy, WGSRPD classification, area and plotting geometries; S2. GBIF botanical-indexing failures versus true zeroes; S3. Additional native-versus-overview AIS and full circuit sensitivity; S4. All adverse earlier-butterfly and latitude-conditioned predictions; S5. Colonial/VOC network and plant trait screens with unresolved or negative outcomes; S6. Missing host chronology and historical-network source coverage; S7. Full figure source receipts, licences and overlap statement for GEB.

**Out of scope:** HydroWASTE, wastewater/WWTP, nitrogen-discharge exposure, and nutrient-dose mechanisms, now in `gesui`.

## Final editorial gates

Only Fig. 4A is already a newly rendered/verified output in the current-figure PR; other panels and captions are **design-level** until re-rendered against source files and visually checked. Journal submission must also wait for independent numeric audit, source licence checks, verified bibliography, overlapping GEB paper disclosure and coauthor approval.
