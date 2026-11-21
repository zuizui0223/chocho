# GEB Fig.3A — accepted-host × geographic region source-rebuild protocol (source audit only)

**2026-11-21; manuscript PR179 immutable at `99b876421af2b4f51e58569862994622a02717d7`. This document and the Actions workflow belong to independent Draft PR189.**

## Original source provenance, code and sample

Original source method is **not** a simple rank of the 670 plants in the figure; its original input consists of the accepted butterfly–plant relationships, native and contemporary WCVP species×WGSRPD3 ranges, and fixed original butterfly resource descriptors. Exact original code:

- `scripts/build_wcvp_hosts_null_sidecars.R` — Git blob `4adb65979333b9c2914f0f15f7b8b37169edf493`; accepts pinned R objects `data/wcvp_names.rda` and `data/wcvp_distributions.rda`, and HOSTS `resource.csv`, selecting nonextinct, nondoubtful species, accepting only exactly and unambiguously resolved scientific names.
- `scripts/analyze_butterfly_host_contribution_concentration.py` — Git blob `b49356fd51ce6313adef6ee2ea079fadbeb0ff35`; accepts original resolved plant IDs, both native/contemporary plant×region layers, and source `s1_resource_descriptors.csv`.
- Kew WCVP v13 original pinned source: `matildabrown/rWCVPdata@65bed76bae9d644ccb6ad200c05f9f5071d89e05` (cite `10.34885/0yex-xv26`).
- NHM HOSTS original pinned source: `globalbioticinteractions/HOSTS@808e0b869f9ec1adf8efff87cf6a395adda103e0` (cite `10.5519/havt50xw`).
- The original 339-butterfly S1 descriptor file is `data/frozen/figure_sources/s1_resource_descriptors.csv`, Git blob `09541b5894529cc8629c7515ea9022ded9edead2`, SHA256 **`894f48dbca1760fc4fa75bfee8f663540ab4b9380f8b9daf2bc09440e8bb0cdc`** (not the independent 239-species figure table).

The original archived source workflow `.github/workflows/geb-reviewer-defenses.yml` explicitly executes those **same two scripts**. Original successful historical GitHub Actions run #36511815412 (2026-09-29) logged `pairs=53449 native_rows=148254 contemporary_rows=222384`; the original `host_contribution_concentration_v0.1.json` source receipt recorded **206/239 butterfly source-eligible expanded species, 14,553 added butterfly×region units, 670 plants, 38 contributing half, Gini 0.777394331967**, and genus sensitivity **431 genera / 25 contributing half**.

## Exact source algorithm

For each eligible butterfly with at least one original documented native-resource region, form `N` = union of all native documented host-plant WGSRPD3 regions and `C` = union of all contemporary host-plant regions, ensuring `N ⊆ C`.

For each newly supplied butterfly–region resource unit in `C − N`, identify the subset of that butterfly's host plants whose contemporary presence in the region is outside its own native range. Split **exactly one unit** equally between these contributing plant species. For genus sensitivity, first deduplicate participating genera within the same butterfly–region unit and then split one unit equally among those genera. Both plant- and genus-level contributions must sum to **14,553**, not plant record counts or independent field sightings.

Rebuild acceptance gate is NOT just seeing 670 or 38 in a cached result. It independently reruns the R creation of the **three required sidecar CSVs** from fixed upstream repositories, then executes the Python fractional-credit allocation, and verifies:
1. Fixed source commits, R/Python script Git blobs and descriptor SHA256.
2. `206` expanded butterflies and `14,553` added units in the new Python output.
3. `670` accepted contributing plant species, `38` accounting for half, contribution Gini, top-k shares and HHI matching the **original frozen source receipt** within `10^−9`.
4. Independent genus model: `431` genera, `25` covering half, Gini and top-k shares matching.
5. All original Python output CSVs are generated **within the Action workspace**, but only a non-restricted **audit receipt containing source input/output SHA256 and aggregate statistics** is archived. Do not automatically distribute upstream HOSTS/WCVP raw tables or derived per-plant credits before source-specific licence/author review.

## Status / publication boundary

[Full original-source rerun workflow](../../../.github/workflows/rebuild-GEB-Fig3A-original-host-credit-20261121.yml), [Actions run #44847332630](https://github.com/zuizui0223/chocho/actions/runs/44847332630). **Until successful complete original-source re-execution is confirmed, treat Fig.3A as frozen-source numerical receipt only.**

Even if the numerical core passes, it demonstrates the original species×region fractional-credit algorithm can be reproduced from public pinned source repositories; it does not establish actual plant introductions, butterfly larval feeding, source legal redistribution, anonymous editor access, or a new citable data DOI. No original manuscript numerical result or PR179 source is changed.
