# Strict 78-port × 35-island nutrient-affinity comparison: source-locked observational audit

2026-11-02 · Draft PR #122 · **exploratory; GEB manuscript unchanged**

## Clear result

**Global source correlation retained:** among port-indexed candidate non-native butterfly larval hosts, 183/687 had an exact German Ellenberg N indicator and the cross-species association between N and the *number of ports with a GBIF record* was **Spearman rho +0.340819**. Genus-restricted permutation (24 genera / 77 species) was not significant (two-sided p=0.38). These facts are unchanged.

**Local matched-island mechanism not established:** same-island GBIF botanical speciesKey facets do **not** show a robust, consistent port-specific increase in the mean nutrient affinity of recorded butterfly host plants. This does not negate the global cross-species association, the GEB +54.9% potential-resource result, 155-butterfly host forecast MSE improvement (4.12% in temporal evaluation), or 51-island host–butterfly compositional correlation (rho +0.300). These represent different response variables, spatial grains, samples and observational questions.

## Full execution and data lineage

The original fixed first-paper candidate set was 36 GSHHG islands, 79 WPI port IDs from the source-locked 2026-10-31 audit. WPI locations are retrospective 2019 proxies; historic shipping cannot be inferred.

1. **Geometry:** [first WPI79 paired coastal construction](https://github.com/zuizui0223/chocho/actions/runs/41399982186) recovered 79 site pairs with nonport coastal candidate points, but original 2021 remote-coast point coordinates were not in this Git repository. The new locations were selected on GSHHG 2.3.7 level-1 polygons, 250 m inland, blind to plant and butterfly outcomes.
2. **1 km GBIF plant pilot:** [first 12-island pilot](https://github.com/zuizui0223/chocho/actions/runs/41400433588) found **0** comparable port–nonport pairs in 2001–2010 and **1** in 2011–2024. "Comparable" required general plant records >=30 and >=2 exactly trait-matched candidate host species **on both sides**. Exact GBIF accepted-species keys were resolved for 182/183 original N-scored hosts. Lack of botanical records is not plant absence.
3. **First 5 km pilot and critical rejection:** [79-site 5km pilot](https://github.com/zuizui0223/chocho/actions/runs/41401009854) produced data for 79 pairs, but an [independent overlap audit](https://github.com/zuizui0223/chocho/actions/runs/41401838417) established that **48 of 79 circular port/nonport 5km search footprints overlapped**. In the actually disjoint subset of that first design, only 1 island (early) and 3 islands (late) were usable. Those original 79-site 5km gap estimates must **not** be interpreted as independent port contrasts.
4. **Corrected strict source-blind geometry:** [strict GSHHG/WPI rerun](https://github.com/zuizui0223/chocho/actions/runs/41402076182) required nonport locations ≥12 km from **any indexed WPI port**, 2*5 km +0.25 km separation of actual port/nonport botanical site centers and within-island controls separated. **78 pairs / 35 islands** passed. One planned WPI port 1439 on island ID 1577 was excluded, rather than relaxed after inspecting responses. Minimum actual port/nonport separation = **11.443 km**. This deliberately differs from original 2021 nonport vegetation sample coordinates.
5. **Corrected same-radius 5 km independent GBIF queries:** [completed strict botanical run](https://github.com/zuizui0223/chocho/actions/runs/41402214941), 78 ports×2 coastal site types×2 periods = **312 site-period index queries**, source-locked location receipt, strict accepted-name speciesKey mapping, 2001–2010 vs 2011–2024, general plant effort summaries and candidate host N trait summaries. **No butterfly outcomes entered this site selection.** Uploaded source receipts as GH Actions artifacts.
6. **Offline independent verification (no new GBIF API calls):** [successful strict source receipt audit](https://github.com/zuizui0223/chocho/actions/runs/41402555134), reproducible using `scripts/exploratory_port_botanical_offline_audit.py`, exact source receipt URLs and code. Verifies 78 pairs, 35 islands, no 5km site-circle overlap, recomputes island-level means, botanical observation effort screens and identical-pair time comparison.

## Accepted observations in strict disjoint comparison

Paired metric: `mean_N(recorded exactly matched candidate non-native butterfly host species, unweighted across species) at the port − at the within-island non-port coast`. **Not** tissue N, local N deposition, amount of host plant biomass, butterfly reproduction or number of individual hosts.

| Summary | 2001–2010 | 2011–2024 |
|---|---:|---:|
| Complete surveyed port / nonport pairs | 78 | 78 |
| Pairs with both sides >=30 botanical GBIF records | 41 | 52 |
| Pairs with >=30 botanical records AND >=2 N-scored candidate hosts per side | **30** | **49** |
| Independent islands with at least one such usable port pair | **9** | **22** |
| Island-weighted mean port − nonport N | **+0.112912** | **−0.012278** |
| Descriptive island-bootstrap 95% interval | **[−0.0340, +0.3020]** | **[−0.3166, +0.2811]** |
| Median island mean difference | +0.0767 | +0.0230 |
| Positive island means | 7/9 | 12/22 |
| Pairs after additional general-record observation-effort ratio ≤5× | 29 | 45 |
| Independent islands after ≤5× general-record filter | 9 | 20 |
| Island-weighted mean after ≤5× filter | **+0.115005** | **−0.054838** |
| Island-bootstrap 95% after ≤5× filter | [−0.0303, +0.3037] | [−0.3771, +0.2493] |

All reported bootstrap intervals are descriptive over *islands*, not independently observed port replicates. None cleanly excludes zero. A different composition of eligible islands across years means year-specific means are not a within-island temporal trend by themselves.

### Truly identical port-pair temporal sensitivity

Only **30** of the strictly disjoint WPI port/nonport pairs across **9 independent islands** qualified for a matched N comparison *in both periods*. Keeping identical pairs fixed:

- Early port−nonport N: **+0.112912** (island-weighted)
- Late port−nonport N on **those same exact pairs**: **−0.049285**
- Late−early change: **−0.162196**; island-bootstrap 95% **[−0.59123, +0.17213]**
- 5/9 island means changed positively and 4/9 negatively, also consistent with no uniform shift.

This is variation in *GBIF-indexed plant composition* conditional on detection and eligibility, not directly observed spread of host plants away from ports. Site environments, observer effort, GBIF indexing, historical local introduced status and reporting completeness can differ.

### General botanical background vs candidate N-scored host richness

Among pairs with both sites >=30 general botanical records (41 early / 52 late):

- Island-weighted log(port/nonport) **all indexed plant species richness**: **−0.04627** early and **+0.41030** late.
- Island-weighted log(port/nonport) **N-scored candidate host species richness**: **−0.16893** early and **+0.51467** late.

These pairs have sufficient *general records* but **not necessarily adequate host trait richness on both sides**. Such ratios are observation proxies; do not read their raw positive signs as an independent host-specific mechanism, nor as host availability or biological population success. A true non-host N trait background is not yet represented.

## Interpretation and decision gate

- **Retained positive ecological association:** globally indexed candidate non-native butterfly larval hosts with high nutrient-affinity scores are recorded at more ports; different trait panels (L/F/S) produce substantially weaker correlations, conditional on their own exact-match subsets.
- **Unresolved local mechanism:** the larger, strict same-island port–nonport comparison offers no consistent port-specific *mean nutrient-affinity enrichment* pattern with reliable island uncertainty. Before the spatial repair the original 5km effect was artificially difficult to interpret because of overlapping botanical footprints.
- **Next discriminating requirement:** measure or externally annotate nutrient exposure, local botanical introduced history, urbanized vs undeveloped nonport coastline (same source/year/land cover), independent recording effort, plus **all plant** Ellenberg N background and independently observed butterfly occurrence/use. Distinguish broad globally successful ruderal taxa versus specifically port-mediated transport, and pair port/not-port sites by environment and historical risk set before interpreting a mechanism.
- **Publication:** no promotion to manuscript Results or mechanistic claim for port sewage, nutrient enrichment, larval performance, ship dispersal, butterfly colonization or stepping stones. This is one complete exploratory source audit; PR remains Draft, manuscript unchanged.

## Executable paths / frozen receipts

- Site selection and WPI / GSHHG verification: `scripts/exploratory_port_nonport_coast.py`, workflow `.github/workflows/chocho-strict-disjoint-coast.yml`.
- GBIF botanical facets: `scripts/exploratory_port_nonport_botany.py`, workflow `.github/workflows/chocho-strict-5km-botany.yml`.
- Pure offline source audit and same-site/time quality: `scripts/exploratory_port_botanical_offline_audit.py`, workflow `.github/workflows/chocho-strict-botany-offline-audit.yml`.
- Existing source-locked N plant trait pilot: `scripts/exploratory_port_nutrient_ellenberg.py`; source inputs and run notes in [README](README.md).
- Underlying original GEB result retained exactly in `main` (base SHA `e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660`) and untouched by this draft branch.
