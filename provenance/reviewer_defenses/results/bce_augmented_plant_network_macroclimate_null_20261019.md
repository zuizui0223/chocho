# Independent augmented-catalogue macroclimate test — source-verified, MCMC qualified

**Decision date: 19 October 2026. Scope: GLOBAL MACROECOLOGICAL POTENTIAL HOST RESOURCES; no new published-paper reanalysis, no butterfly recruitment.**

## Main question

Can the residual climatic sorting of newly available larval host resources reported from the original HOSTS-only 239-butterfly × WGSRPD3 analysis withstand (a) **1,027 additional documented WCVP-accepted butterfly–host links** for 81 species and (b) an independent plant-side randomization that preserves original accepted butterfly→host plant identities, introduced-host taxon families, native host geography, continent-scale plant introduction, coarse climate, and region×plant-family degrees?

**Qualified answer:** Under the Europe-biased *partial* BCE/Clarke host supplementation and the two prespecified native-host climate windows, a **small positive resource-climate alignment residual persists** even after the plant-network/plant-climate-constrained null. This is additional evidence of structure in **potential larval resources**, *not* causal butterfly population expansion, larval fitness or an evolutionary host preference mechanism.

## Exact source and denominator audit

- Original frozen source 239 butterflies; HOSTS-only geographic totals: native **26,530**, contemporary **41,083**, introduced-added **14,553** butterfly×WGSRPD3 cells.
- BCE/Clarke partial 1,027 accepted host links affecting 81 butterflies, 679 unique additional accepted plant IDs (of which 157 were already used by some other originally focal butterfly). Fully rebuilt augmented geography: native **29,852**, contemporary **44,819**, introduced-added **14,967**.
- Added accepted-host plant family assignments: **323** resolved by original HOSTS plant ID, **304** by a unique original HOSTS genus-family assignment, and exactly **52** by pinned **GBIF botanical family backbone** matches. One source spelling mismatch, WCVP `Galactites tomentosus` versus GBIF `Galactites tomentosa`, was corroborated at **genus Galactites = Asteraceae**, not silently auto-matched for host use.
- Across the 365 original climate-covered WGSRPD3 codes, **37,869 introduced host plant×region edges** are represented after BCE host supplementation, versus **30,373** in the original HOSTS focal graph. The 2 additional introduced regional codes `ANT`, `HMD` for *Coenonympha oedippus* lack the frozen CHELSA region climate vectors, and are explicitly **kept in geographic totals but omitted only from climate calculations**.
- Climate-evaluable native-host baseline species **234 → 235** because *Satyrium pruni* moves from 2 to 66 native-host resource regions when candidate hosts are included. Augmented climate-evaluable added opportunities are **14,961** (full added geography=14,967, 2 lack climate; 4 further cells belong to other species with insufficient native climate support).
- Compare BOTH the **augmented eligible 235 species** and a **fixed same original 234 eligible species** (exclude newly eligible *S. pruni* only), using the **SAME NULL DRAWS**. Matched observed added-cell denominator **14,903**. This avoids conflating an added host inventory with a new eligible butterfly in the climate panel.

Original pinned author-only HOSTS/WCVP/CHELSA archive: [prior source-reconstruction run 39006554609](https://github.com/zuizui0223/chocho/actions/runs/39006554609), recovered with exact fixed SHA; all original source manifest members verified. No original copyrighted plant geography is redistributed in the outputs.

## Results: difference from null median, percentage points

Each plant-null run includes **two independent seeds 20261019/20261020, 199 draws each, 10× plant-link-count burn attempts, and 4× between samples**. Invariant checks preserve plant native-range exclusions, plant×area(/climate stratum) added-range counts, and recipient region×plant-family introduced host degree. Both nulls pass every ESS≥30, chain agreement, burn movement and accepted swap gate, for **all four** predeclared statistics (both climate windows × all 235 / matched 234 eligible species).

| Macro estimand | Observed alignment | Plant family + continent null residual | Plant family + continent + coarse BIO1/BIO12 null residual |
|---|---:|---:|---:|
| 235 augmented eligible, 5–95% native-host climate | 9,949/14,961 = **66.500%** | +3.581 to +3.645 pp | **+0.972 to +0.977 pp** |
| 235 augmented eligible, 10–90% native-host climate | 7,370/14,961 = **49.261%** | +4.405 to +4.428 pp | **+1.571 to +1.599 pp** |
| **Matched original 234**, 5–95% | 9,903/14,903 = **66.450%** | +3.595 to +3.650 pp | **+0.982 to +0.993 pp** |
| **Matched original 234**, 10–90% | 7,331/14,903 = **49.191%** | +4.393 to +4.427 pp | **+1.566 to +1.592 pp** |

For the *augmented 235*, under the **most conservative climate-stratified plant network null**, 5–95% observed 66.4996% versus null medians **65.5278% / 65.5225%** (two-chain ESS ≈158.4 / 155.0); at 10–90%, observed 49.2614% versus null **47.6627% / 47.6901%** (ESS ≈168.8 / 163.7). The matched original 234 show the same sign and comparable size. **Null added-region denominators are variable**, median ~15,614/15,609 for 235 and ~15,550/15,544 for matched 234, while observed denominators are respectively 14,961 and 14,903. Thus we compare fractions `O/A` within each draw, *never raw aligned cell counts*. 

This pattern is directionally concordant with **original HOSTS-only** source-verified plant-network-climate null from [run 39212635930](https://github.com/zuizui0223/chocho/actions/runs/39212635930): 5–95% residual +0.861 to +0.903 pp and 10–90% +1.390 to +1.416 pp. No threshold-selection based on sign; both were fixed before augmented MCMC execution. The continental-only nulls differ between source catalogues in their sizes and should not be interpreted as additive causal decompositions.

### Qualification and provenance

**[GitHub Actions run 39221327128](https://github.com/zuizui0223/chocho/actions/runs/39221327128)**: successful **exact-source MCMC** including a nine-draw source/margin smoke check and full two-mode × two-seed × 199-draw execution. Derived machine-readable JSON uploaded as [artifact 12172690082](https://github.com/zuizui0223/chocho/actions/runs/39221327128/artifacts/12172690082). Every one of the eight independent outcome×mode qualification checks passed, including source family crosswalk 323/304/52 and all species/geography counts. **199 sampled MCMC states are serially correlated**; the minimum per-chain add-one nominal right tail was 0.005 in this run, **not** a valid independent permutation p-value and **not** a test of butterfly population increases.

An earlier 9-draw smoke attempt showed near-zero ESS and correctly remained NOT_EVALUABLE; the full 2×199 source/mixing gate passed. A previous first workflow attempt failed only on a `viable_edge_fraction` vs `viable_fraction_of_total` source diagnostic naming difference, which was explicitly harmonized before the verified full run. Do not substitute the failed smoke's numbers for qualified inferential results.

- Analysis: `scripts/analyze_bce_augmented_plant_network_climate_null_20261019.py`.
- Source-only botanical family backfill: `data/exploratory/bce_additional_plant_family_gbif_crosswalk_20261019.json`.
- Isolated source-verified CI: `.github/workflows/bce-augmented-plant-network-climate-null.yml`.
- Prior original-HOSTS result hierarchy: `docs/exploratory/BUTTERFLY_MACRO_CLIMATE_SORTING_RESULT_GATE_20261019.md`.

## What is ecologically nontrivial, and what is *not* identified?

1. The increasing *resource-envelope union* is a set-theoretic certainty; **climate sorting of newly supplied larval hosts and persistence of a small residual across two different host inventories** is not that certainty. A mostly plant-driven climate sorting pattern is quantitatively resolved.
2. Most of the large +8.28 pp continent-uniform climatic excess is explainable with much richer source, network, plant-family and climate-bin structure. After controlling broad plant climatic occupancy, the source-robust *potential resource* residual is approximately **one percentage point** under 5–95% native-host climate and **1.4–1.6 points** under 10–90% climate. Residuals are small in magnitude even when numerical MCMC tails are extreme.
3. The effect is conditional on coarse *plant* climatic bins and the particular HOSTS/BCE and WCVP fixed-source mappings. **A finer plant spatial-and-climatic matching null, geographically representative alternative host source, unequal-region-area correction, and true butterfly local occupancy** remain absent. These are possibilities for future macro empirical falsification, not claims supported by the current result.
4. The 81 BCE-augmented butterflies are European-evidence-biased, not 81 global independent experimental replications. **Do not claim worldwide butterfly host quality, host switching, reproductive success, abundance increase, butterfly community homogenisation, or conservation gain.**

**Scientific judgement:** A credible, independently source-recomputed *macroecological structural* result emerges: anthropogenic redistribution of plant resources preferentially duplicates known butterfly-host opportunity into native-host-like climates; the bulk is inherited from plant geography while a small, catalogue-robust conditional residual persists. The next novel biological leap would require independent butterfly occurrence/abundance, phenology and local host-use data; none was measured here. Original [GEB frozen PR #38](https://github.com/zuizui0223/chocho/pull/38) remains unchanged; alternative exploratory [PR #46](https://github.com/zuizui0223/chocho/pull/46) is still Draft.
