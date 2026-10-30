# 先行寄主記録の情報不足を分解する：固定した51島ペアの植物観測機会監査（2026-10-30）

**作業の位置付け：** [海上輸送→外来寄主→長距離記録→島の再記録→次島](maritime_host_rescue_steppingstone_20261030.md)を検証するためのデータ適格性調査。既存の2001–2010植物データを新規に取得せず、PR #104に公開済みの集計から、先行する植物記録が島の両側に存在する比較の範囲を固定する。これは登録済みの確証的仮説検定ではなく、既存の5/1,163一致を見た後の探索的QA。新しいp値を算出しない。

## 固定入力

- 194港の初期（2001–2010）と全年（2001–2024）植物インデックス件数：`provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_plant_availability.csv` Git blob `0a01831056c2cd59a4bd41bd1dbd4c84d77cf922`
- 51島×51島の全1,275ペアと2001–2024 butterfly data中の2011–2022新規共通記録種数：`provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_pairs.csv` Git blob `5663d2ffe42419f5c4496b92c9f3bbf76175907b`
- 早期植物の港取得適格島：`provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_port_early_inventory_islands.csv` Git blob `143f9956ce4f871c7a97b4777ab69d416963c2bb`
- 5件の年代付き寄主対応：`provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_links.csv` Git blob `88db773c74dea9cd2d58bde0b67a394d93226869`
- 以前の固定適格51島：`provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_islands.csv` Git blob `18cecde5bf3066f276380167f451f1fa5cc1caa1`

## 解析の定義（数値を見る前に固定）

(1) `landmass_id`と`port_id`は文字列。194港の各`early_plant_records`、`all_plant_records`は非負整数で`early<=all`。全51島に完全取得フラグがあり、全ペア1,275行の島ID集合と一致することを検証。

(2) 同じ島に複数港がある場合、**単一港内で最も多い早期植物記録数**`max_early`を島の*観測機会の代理量*にする（重なり合う5 km港範囲の件数を足さない）。補助として島内全港のearly件数合計とall件数合計を記録するが、独立した標本数・植物量とは呼ばない。`max_early=0`でも植物が生物学的に存在しないとはしない。

(3) 全島ペアを`both_early>=1`、`both_early>=10`、`both_early>=100`の固定階段に分類。1は実際の記録あり、10/100は記録の乏しさに関する記述的感度区分であり、観察完備の認定ではない。しきい値を結果の方向で変更しない。

(4) 各区分で対象島ペア数、分母となる**共通新規蝶種×島ペア数**（`shared_new_species`の合計）、対応する先行寄主一致の5行のうち何行が入るかを報告。各比較の`island_a/island_b`は順序なしの正規化キーとして照合。分母は同じ種が複数島ペアに現れ得るため独立nではない。

(5) 1,163件・5リンクの既存値と行単位で完全に照合。`early>0`は植物一般の観測機会であって、具体的な寄主植物の当時の不在／存在ではない。早期の寄主記録が両島に見つからない場合の上限推定や、捕捉補正した寄主効果係数を出さない。

(6) この監査は地域・年による観測努力と植物カタログの齟齬、前期・後期の植物種構成変化を直接説明しない。島間の標高・距離・港の大きさを追加調整しない。対象集団を事後に入れ替えない。最後に海域別・植物の年別/原座標データと研究対象蝶の先行寄主を結べるかの取得優先順位だけ提示。

## 合格・停止

公開派生CSVの全キー・合計が一致し、固定4条件（全群、1、10、100）の分母、5リンクの所属、植物記録零の別状態を再現できればPASS。全群ペアの合計が1,163でなければ即停止。部分集合で元の生態学的効果を再計算せず、既存の年代付き帰無モデルの抽出群を勝手に再定義しない。結果は #104 内の探索文書のみで、独立GEB原稿は変更しない。
