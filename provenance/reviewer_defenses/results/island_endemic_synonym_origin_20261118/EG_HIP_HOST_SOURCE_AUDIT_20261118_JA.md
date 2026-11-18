# 島固有蝶の幼虫寄主植物の直接証拠：Euchloe grancanariensisと近縁5種の未同定境界（2026-11-18）

## 新しい実測・文献確認（原典を追跡したsource-typed audit）

**独立したカナリア固有蝶 `Euchloe grancanariensis` の野外幼虫寄主植物は `Hirschfeldia incana` と種まで記録されている。** Acosta-Fernández（2010, SHILAP 38(149), pp. 121–128）は、2009年3月、グラン・カナリア島El Batánで`H. incana`の果実上から幼虫15個体を採集し、その後の飼育で **9個体の寄生・1個体の原因不明死亡・5個体の蛹化** を記述。蛹化したうち2個体は異常があり死亡し、残る3個体は休眠蛹と記載される。**成虫羽化成功はこの記録から確認できない。** 比率を野外個体群の寄生率や特定寄主の成績に一般化しない。

原著の名称は `Euchloe belemia grancanariensis`（亜種扱い）である。独立のカナリア政府BIOTAは `Euchloe grancanariensis` を有効な固有蝶種として分類している。系統・種/亜種概念を一致させず、原著の寄主を他島の `E. belemia` 全体へ流用しない。
- 昆虫野外採集・飼育・原著：https://www.redalyc.org/pdf/455/45514996013.pdf （原典 p. 123–125）
- 独立種・在来の分類：https://www.biodiversidadcanarias.es/biota/especie/E06769

別の2014年1月21日のグラン・カナリア島Cruz Verdeでは、雌が **`Descurainia preauxiana`** に産卵した野外写真と卵が記録されている。これは**産卵植物の直接証拠**であり、孵化後の幼虫がこの植物で生育・蛹化したことまでは示さない。
- 野外産卵・卵の画像と植物種名：https://lepiforum.org/wiki/page/Euchloe_grancanariensis
- 独立植物の同定・分布：https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:60460465-2

元GEBと同じ**WCVP v13** `matildabrown/rWCVPdata` 固定commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05` に戻り、植物受理名とCNY地域の在来/外来を再実行した。**`Hirschfeldia incana`・`Descurainia preauxiana` 両方で受理植物種が一意・CNY native** と確認。これは外来寄主利用ではなく、固有蝶と在来植物種の幼虫利用／産卵関係である。
- [原WCVP固定照合・15幼虫のアウトカム CI #44213255146 成功](https://github.com/zuizui0223/chocho/actions/runs/44213255146)

## 現生の島固有蝶の寄主属収録への影響

前ターンまでの独立42種（現生41種＋歴史的絶滅分類群1種）では、特定のデータベースと外部自然史資料の**文献寄主属記録を持つ現生種35/41種**を確認していた。今回の `E. grancanariensis` は以前の未記録6種の一つであり、受理植物種による野外幼虫観察が新たに見つかったため、

| 観測分母 | 以前のソース統合集計 | 今回の確認済み下限 |
|---|---:|---:|
| 現生の独立蝶41種 | 35/41 | **36/41（87.8%）** |
| カナリア固有蝶13種 | 7/13 | **8/13（61.5%）** |
| カナリア固有蝶の`Hipparchia`以外8種 | 7/8 | **8/8** |
| カナリア固有の`Hipparchia`5種 | 0/5 | **0/5：寄主属が未同定** |

これらは文献探索と分類名照合を行ったうえでの**出典上の寄主属記録の確認済み下限**であり、全種を無偏に調べた生態学的寄主既知率ではない。歴史的に絶滅した `Pieris wollastoni` の過去の寄主報告は別台帳に保持し、現生種の割合に入れない。

## Hipparchiaで何がまだ不明なのか

Spilani et al. (2025), *Insect Systematics and Diversity* 9(5), ixaf043 はカナリアの `Hipparchia` 5島固有種の独立分化を系統データで支持するが、野生で利用される植物種は不明とする。イネ科植物の利用は主に飼育観察から類推される。

独立の野外自然史記録（https://www.pyrgus.de/Hipparchia_wyssii_en.html ）では、**グラン・カナリア島の `H. tamadabae` 幼虫が2016年12月に株立ちの草上で複数観察された**が、植物は同定されていない。同所の`Hyparrhenia hirta`には幼虫を確認していない。したがって、寄主について記録可能な最小階級はPoaceae（イネ科）であって、特定属・特定種にはできない。この1例で他4島の近縁分類群すべての寄主植物種も決定しない。
- 種分化の独立源：https://doi.org/10.1093/isd/ixaf043
- 野生幼虫と未同定草本：https://www.pyrgus.de/Hipparchia_wyssii_en.html

**ここでの新しい判断は、島固有性に一律の欠落が生じているというより、現時点で寄主植物同定が進んでいない一つの近縁放散群に未解決が集中している**こと。元GEB239種の外来寄主資源拡大+54.9%（全球）/島嶼+118.9%（源集団）は変更しない。

## 別の正の比較：同じ寄主植物が近縁蝶に必ず使われるわけではない

同じAcosta-Fernández（2010）のフエルテベントゥラ島 `Euchloe belemia hesperidum`（現在の独立カナリア分類 `Euchloe hesperidum`）に関する**野外採集1幼虫の飼育観察**では、

- 提示した`Hirschfeldia incana`の果実を拒否。
- **`Sisymbrium erysimoides`と`S. irio`を飼育下で摂食し、その後に蛹化**。

一方グラン・カナリア島 `E. grancanariensis` の幼虫15個体は野外で `H. incana` 上から採集された。**同属・近縁分類群でも、在来寄主植物が存在することと実際の利用は同義でない**という非自明な源資料上の事例。ただし同じ設計の反復した相互受け入れ試験ではないため、種間の食性分化・選択・適応の統計的実証とはしない。

原WCVP v13に戻った照合では、`Hirschfeldia incana`、`Sisymbrium erysimoides`、`S. irio`はいずれもカナリアCNYで受理在来植物種（[再現CI #44213967074](https://github.com/zuizui0223/chocho/actions/runs/44213967074) 成功）。`S. irio`を外来植物として扱わない。3件の提示植物は**同じ幼虫1個体**で試されたので、3つの独立した実験標本として検定しない。

## 再現ファイル

- 原資料から源別に転記した4行（3 `E. grancanariensis`、1 `H. tamadabae`）：[E_grancanariensis_H_tamadabae_source_evidence_20261118.csv](E_grancanariensis_H_tamadabae_source_evidence_20261118.csv)。
- 原独立チェックリスト・原産区分・寄主・9/15寄生・5/15蛹化・未同定イネ科を検証するコード：`scripts/audit_Egrancanariensis_primary_host_20261118.py`。
- 元WCVP固定版の植物種在来性：`scripts/audit_Egrancanariensis_WCVP_v13_20261118.R`。
- 二種の野外15幼虫と1幼虫の飼育受け入れ対照：[Euchloe_two_island_host_use_field_vs_lab_2010.csv](Euchloe_two_island_host_use_field_vs_lab_2010.csv)。
- 擬似反復を拒否する照合：`scripts/audit_Euchloe_native_host_acceptance_field_lab_20261118.py`、`scripts/audit_Euchloe_field_lab_WCVP_v13_20261118.R`。

**推論境界**：固有蝶に外来寄主利用の事例がある（別の`Pieris cheiranthi`と`Tropaeolum majus`）一方、`Euchloe`のこの2種について確認した植物種は在来だった。存在する在来・外来植物相の属だけでは、実際の幼虫寄主植物種・成績・港を通じた定着を予測しきれない。全島固有蝶の平均寄主資源拡大率は引き続き**非同定**。
