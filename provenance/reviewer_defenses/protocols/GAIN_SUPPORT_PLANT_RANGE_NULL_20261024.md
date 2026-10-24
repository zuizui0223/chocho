# Frozen exploratory plant-range null: is the gain–redundancy slope more than resource-set geometry?
**2026-10-24 | one primary structural falsification test | frozen before null results.**
Do not modify the GEB manuscript on `main` or PR #38. This follows the observed negative species-level gain–single-host slope from `analysis/gain-fragility-20261024`: **β = −0.2207187290** for 206 gained-resource butterflies (adequate-host β = −0.2711287143 for 191). It is a *posthoc null* designed after seeing this pattern, not a preregistered independent study.

## Test question and direction
Can the observed negative across-butterfly association between reconstructed geographic gain and the **fraction of added resource cells supported by one documented introduced host** be distinguished from a structural randomized geography in which all recorded host identities and major plant/regional opportunity margins are unchanged?

The original prior analysis only genus-bootstrapped a slope fitted to mathematically coupled responses from **the same plant-host union map**. A positive/negative within-source genus interval cannot answer whether plant-side geography alone produces the pattern. We now compare the signed adjusted OLS slope against its own constrained plant-range swap null.

## Fixed sources and exact reconstruction
- HEAD-snapshot source inputs: 239 butterflies from three git-blob SHA-locked files `anthropogenic_species_metrics.csv` (9f83d2d9ac8e46f1cfd4ded572c7a023273503dd), `host_contribution_metrics.csv` (4cf4851499ca2233ff2adadf1c4ee7f0e863f23d), `s1_resource_descriptors.csv` (09541b5894529cc8629c7515ea9022ded9edead2); do not refit or change hosts.
- Independently reconstruct plant native/introduced host records from pinned `globalbioticinteractions/HOSTS` commit `808e0b869f9ec1adf8efff87cf6a395adda103e0` and WCVP `matildabrown/rWCVPdata` commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`, using existing `scripts/build_wcvp_hosts_null_sidecars.R`.
- WGSRPD3→Level1 mapping from `data/frozen/wgsrpd3_level1_lookup_v0.1.json`, TDWG commit 52da7828aba9d461dd133c27b3bd7a4407161f54. Reuse deterministic, tested `build_model` / `attempt_range_swaps` from `scripts/analyze_butterfly_shared_host_structural_null.py`.
- Must reproduce **239 butterflies, 206 observed with additional resource regions, 26,530 native butterfly×region units, 14,553 introduced-added units and 8,574 single-host-supported added units**, and for each species match the frozen `host_contribution_metrics.csv` and `anthropogenic_species_metrics.csv` rows. Fail rather than silently mis-join or adjust any source.

## Margin-preserving graph null
Swaps two introduced `host plant × WGSRPD3` incidences **within the same plant family × WGSRPD Level1 group**, rejecting self-pairs, duplicate plant-area pairs and any new introduced position that is a native area for that same plant.

Preserved by construction: all butterfly–host identities and host breadth; all native plant range cells; each introduced plant's number of regions within each continent; each recipient region's number of introduced plants per plant family; no native overlap and global introduced plant-edge total. NOT preserved: exact native climatic analogues, plant species distribution continuity, introduced region identity for each butterfly, or actual feeding. **Changing butterfly newly added region totals is allowed.**

For each source state, reconstruct per-butterfly U (new cells) and S (cells supported by *exactly one* distinct introduced recorded host), use the same ecological estimand as original: `gain=log((native+U)/native)`, `single_support_share=S/U` when U>0. Equal-species OLS `single_support_share ~ 1+gain+log(native_count)+log1p(recorded_host_species)+host_family_count`.

Primary cohort rule: ALL 239 candidate butterflies are evaluated with the same `U>0` rule in observed and each null draw; record number of qualifying butterflies and host-catalogue adequately recorded count per state. This can alter the number of species across null states, so this is a **dynamic opportunity-panel null**, not a fixed-population counterfactual. Secondary fixed-original-panel statistic is estimable ONLY if all 206 originally gain-positive butterflies stay U>0 in each null draw; otherwise label it NOT EVALUABLE, not zero.

Secondary robustness: repeat the adjusted slope within adequately recorded host species (`host_taxonomy_lower_bound_adequate=True`), same data and null conditions. Do not search other taxonomic traits/thresholds to rescue an effect.

## Numerical procedure & predeclared gates
- Independent random swap chains seeds **20261024** and **20261025**.
- **5× E** proposed swaps to burn in, **2× E** between observed states, **199 sampled states per chain**; E must exactly be **30,373 introduced plant×WGSRPD3 edges** and at least 28,000 of them must be eligible for grouping.
- Both chains require burn accepted swaps >=100, original-edge movement >=10%, and per-draw accepted swaps min >=100.
- For adjusted β, approximate positive-sequence autocorrelation ESS >=30 per chain, and absolute difference between chain medians <=2×max(chain standard deviations). Test ESS also on adequate-host β, and log the null cohort size fluctuations. If gates fail, label NOT EVALUABLE and **do not interpret numerical tails**.
- Compare observed signed β with each chain's 2.5/50/97.5 percentiles, and summarize the observed-minus-null median (Δβ) and serial empirical left-tail fraction `(1+#null β<=observed)/(n+1)` **as an exploratory MCMC diagnostic, not an exact independent permutation p-value**.
- Observed β from reconstructed originals must match −0.2207187290 and adequate β −0.2711287143 within 1e−8.
- Reject unexpected source changes, total incidence drift, failed-native geography exclusions, zero-denominator handling, inadequate chain mixing and unverified 239-source join.

## Interpretive outcomes
1. If observed β lies comfortably within both null intervals, the nominal gain–support slope is **consistent with structural degree/geographic constraints** and must NOT be elevated to a new ecological law or conservation mechanism.
2. If outside both qualified null ranges on the same side, call it a **source-conditional additional geographic arrangement signal**, not evidence of competition, resource insurance, establishment or increased population abundance. A subsequent climate-and-locality matched null and independent natural-history evidence are then necessary.
3. If chains fail/contradict, return NOT EVALUABLE and stop, rather than reporting whichever chain supports the hypothesis.

Preserve machine JSON results, provenance, code and CI. The scientific object is only potential known-host resource geography. No GloBI source missingness or prior monarch/Kyoto experiment is used as biological replication.
