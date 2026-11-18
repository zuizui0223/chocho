# 島固有蝶の外来幼虫寄主利用：植物種・個別島までの実証境界（2026-11-18）

## 追加の科学的結果

**カナリア固有蝶 `Pieris cheiranthi` は、野生化した外来植物 `Tropaeolum majus` を幼虫期に利用することが、少なくとも2島の別々の写真記録から確認できる。** 2015年2月17日テネリフェ島 Puerto de la Cruz のL4幼虫、2011年5月20日ラ・パルマ島 Las Manchas の庭園の幼虫（植物の一般名Kapuzinerkresse＝キンレンカ；種同定は原資料の文脈付き）をそれぞれ別の資料行として凍結。

- 昆虫源：[Lepiforum, Pieris cheiranthi](https://lepiforum.org/wiki/page/Pieris_Cheiranthi)。
- 外来植物由来：カナリア州[EXOS Tropaeolum majus](https://www.biodiversidadcanarias.es/exos/especie/F00988)（Introducido Seguro Invasor、主要7島で野生化）。庭園内植栽だけでなく**当該植物種が地域で帰化している**ことの独立した根拠。
- 元のGEBと同じWCVP v13（rWCVPdata commit `65bed76bae9d644ccb6ad200c05f9f5071d89e05`）では、受理植物種 `T. majus` の `CNY` 分布は **introduced**。植物種・由来・幼虫利用の3層が別々のソースで交差し、[照合GitHub Actions #44206655030](https://github.com/zuizui0223/chocho/actions/runs/44206655030) が成功。

## 対照：同じ固有蝶の在来寄主植物種

同じ固定CNY植物相で `Crambe strigosa`（テネリフェ島での寄主文献）と `Descurainia millefolia`（ラ・パルマ島での寄主文献）が **native** として受理植物種の分布に照合できた。外来植物のみを利用する種ではない。

- 植物・蝶の独立文献：[Acosta-Fernández ほか、2009年カナリア蝶類現地記録と保全状況](https://www.redalyc.org/pdf/455/45512170002.pdf)。
- `Crambe santosii`（`C. gigantea` に近い分類処理）および原著 `Brachypodium arbusculum`（地域名 `B. arbuscula` との照合）は、**受理植物種名の整合が未確認**。WCVP v13で機械的に同物異名へ変換しない。
- マデイラの `Gonepteryx maderensis` や歴史的絶滅 `Pieris wollastoni` の植物名には、カナリアの`CNY`の原産区分を適用しない。マデイラ地域コードの独立植物由来が必要。現在の絶滅種による利用/定着推定は行わない。

## 異なる二つの資源拡大を混同しない

*P. cheiranthi* の**同じ島**に在来寄主の文献記録がある一方で、外来 `T. majus` に幼虫が見つかっている。この資料からいえる最小限は、

| ソースの記録単位 | テネリフェ | ラ・パルマ |
|---|---:|---:|
| 在来の幼虫寄主植物種の文献記録 | `Crambe strigosa` | `Descurainia millefolia` |
| 外来植物上の幼虫の写真記録 | `Tropaeolum majus` | `Tropaeolum majus` |
| 文献で確認できる新規導入寄主植物種 | **+1** | **+1** |
| 以前から在来寄主が記録される島で、新たに増えた「寄主を持つ島」という記録上の単位 | **0** | **0** |

**生態学的な含意**：外来植物の利用が、もともと寄主がある島での**幼虫寄主の選択肢／資源ポートフォリオを追加**しうる。これは「新しい地理的島への進出」そのものではない。**ただし、植物が同じ場所・季節に共存したこと、幼虫成功率、利用頻度、個体群増加は測っていない**。在来・外来双方の利用は異なる年月・資料にまたがる。

この例は、従来のGEB239種「蝶×WGSRPD3の潜在寄主資源地理拡大」（全球 +54.9%、島サンプル相対 +118.9%）と**異なるestimand**である。既存の肯定結果は削除・変更しない。元の数値をこの島固有蝶へ代入したり、「島固有蝶の地理的資源が2倍」と主張したりしない。

## 制約：利用は成功ではない

2009年の調査論文は、テネリフェ島で1984年に採集された約12幼虫を `T. majus` により飼育したが、蛹化前に寄生蜂で全個体死亡した例を記録。これは**当該植物への産卵や幼虫利用と、その幼虫が成虫になる成功とを区別する実例**。植物毒性や「外来植物では繁殖できない」ことを証明したわけではない。野生個体群の持続や個体数増加・外来植物の因果効果の主張は禁止する。

## 実行と再現

- 原文・写真・植物由来・地理・確実度を独立に記録した13ソース行：[species-level evidence ledger](independent_butterfly_host_plant_species_events_v1.csv)。外来2件は画像・日付・島・幼虫段階あり。historically extinct 2件と植物分類名未解決3件は隔離。
- **原蝶母集団との島名一致**、報告された植物種/野外幼虫/歴史絶滅の区分を独立再監査：`scripts/audit_focal_island_larval_plant_species_20261118.py`。
- 固定WCVP v13の種受理名・カナリアCNY native/introducedで **植物種** を確認：`scripts/audit_focal_island_plant_WCVP_v13_20261118.R`。
- **幼虫寄主植物種の追加と「新しい寄主を持つ島」の追加の違い**を再演算：`scripts/audit_Pcheiranthi_host_portfolio_vs_geography_20261118.py`。
- 計算・データ保存：`results/focal_island_larval_plant_species_20261118/` をGitHub Actionsアーティファクトへ保存。

## 原資料ベースの全球拡張データ取得の制約

Clarke (2024) 欧州蝶類幼虫寄主総覧（Dryad DOI `10.5061/dryad.1vhhmgr12`）については、**18ファイルの公開メタデータ、版ID 274312と原ファイル名一覧まで確認済み**。しかしAPI上の実ファイル取得はHTTP **401 Unauthorized**（[失敗を記録したsource run #44204959410](https://github.com/zuizui0223/chocho/actions/runs/44204959410)）。[メタデータのみの成功run #44206778717](https://github.com/zuizui0223/chocho/actions/runs/44206778717) で**寄主行を1件も取得・解析していない**ことを明記。よってこの研究の「証拠レベル別ランク」を実測した結果として使ってはいない。今回の証拠ランクは、独立文献・写真の各sourceを研究者が区別したもの。

**判定**：島固有蝶における外来寄主利用は特定の1種・2島について確実性を上げられたが、全島嶼固有蝶の平均寄主資源地理拡大量は未推定。これが現段階で科学的に守れる境界。
