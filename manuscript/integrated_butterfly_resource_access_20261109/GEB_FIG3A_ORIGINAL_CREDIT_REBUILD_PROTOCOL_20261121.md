# GEB Fig.3A — accepted-host × geographic region source-rebuild protocol (source audit only)

**2026-11-21; manuscript PR179 immutable at `99b876421af2b4f51e58569862994622a02717d7`. This document and the Actions workflow belong to independent Draft PR189.**

## 実際の原データ再計算結果：SUCCESS（2026-11-21）

**[元のHOSTS・WCVP・339種記述子から実行したGitHub Actions #44847332630](https://github.com/zuizui0223/chocho/actions/runs/44847332630) が成功し、植物別・属別の再計算値が元の固定報告に一致した。** R中間データは原記録と同じ **解決済み蝶–植物ペア53,449、在来植物×地域148,254、現代植物×地域222,384**。Pythonで再計算した **206蝶種、14,553追加機会、670植物種、上位38植物で半分、Gini 0.7773943319671659、上位10植物25.10105136%、431属、上位25属で半分、属別Gini 0.7737917033318613** が元の報告と一致した（gini / concentration / top-k / hhi は `10^-9` 以内）。

[固定された元計算結果・入力/中間/出力SHA256（JSON）](../../../results/geb_v05_private_source_recovery_20261120/ORIGINAL_FIG3A_WCVP_HOSTS_CREDIT_REBUILT_SOURCE_QA_20261121.json)。結果ファイル自体のSHA256 `af5d7327bcb8883c85be94250a218348b40906e8e35b5f8636a87cea21960eff`。元の2つのRDataファイル、HOSTS `resource.csv`、生成した3つの中間CSV、植物・属のフルクレジット表のSHA256を収録し、**実ファイル本体は権利未確認につきGitHubへ追加していない**。厳密な元データの再配布・完全査読アーカイブ・間接匿名性の承認とは別である。

この結果をもって **Fig.3Aの植物別寄与計算は「数値サマリーとの照合のみ」から「元ソースから計算再実行済み」に昇格**。ただし以前の七ファイルFig.1/Fig.3B/C限定ZIP自体には信用配分の入力は含まれないため、そのZIPの再現範囲は遡及的に拡大しない。

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

[Full original-source rerun workflow](../../../.github/workflows/rebuild-GEB-Fig3A-original-host-credit-20261121.yml), [Actions run #44847332630](https://github.com/zuizui0223/chocho/actions/runs/44847332630). **Original-source re-execution succeeded in Actions #44847332630. Fig.3A's numerical core is reproduced from the pinned original HOSTS×WCVP interactions and regions; repository rights and reviewer deposit are still unapproved.**

Even if the numerical core passes, it demonstrates the original species×region fractional-credit algorithm can be reproduced from public pinned source repositories; it does not establish actual plant introductions, butterfly larval feeding, source legal redistribution, anonymous editor access, or a new citable data DOI. No original manuscript numerical result or PR179 source is changed.
