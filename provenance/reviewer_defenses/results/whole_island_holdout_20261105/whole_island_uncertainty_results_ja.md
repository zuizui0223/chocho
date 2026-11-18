# 島外予測の不確実性と港の追加価値

既存225候補・155種の主解析を変更せず、固定済みの25島・95種の追加評価を再標本化した。点推定を確認した後の追加診断であり、事前登録解析とは区別する。種と島を独立に復元抽出する交差再標本化2,000回。予測値は固定し、モデルの再学習はしていない。

|対象|重み|比較|MSE差|95%再標本化区間|
|---|---|---|---:|---|
|all_islands|species_equal|Natural geography|-0.001614|[-0.006728, +0.001833]|
|all_islands|island_equal|Natural geography|-0.002699|[-0.006394, +0.000109]|
|all_islands|species_equal|Geography + ports|-0.002210|[-0.007634, +0.001443]|
|all_islands|island_equal|Geography + ports|-0.002529|[-0.006296, +0.000574]|
|all_islands|species_equal|Ports added to geography|-0.000595|[-0.001608, +0.000193]|
|all_islands|island_equal|Ports added to geography|+0.000169|[-0.000888, +0.001208]|
|without_GB_Ireland|species_equal|Natural geography|-0.002680|[-0.007377, +0.000508]|
|without_GB_Ireland|island_equal|Natural geography|-0.002914|[-0.006654, +0.000257]|
|without_GB_Ireland|species_equal|Geography + ports|-0.003117|[-0.008333, +0.000218]|
|without_GB_Ireland|island_equal|Geography + ports|-0.002717|[-0.006613, +0.000756]|
|without_GB_Ireland|species_equal|Ports added to geography|-0.000436|[-0.001471, +0.000339]|
|without_GB_Ireland|island_equal|Ports added to geography|+0.000197|[-0.000998, +0.001334]|

MSE差は追加モデル−比較モデル。負なら改善。種等重みは各種内の全評価行を平均、島等重みは各島内の全評価行を平均する。再標本化では反対側のクラスタ重みも適用し、記録が残らない種／島はその反復の平均から除く。

区間は固定された予測誤差の標本構成への感度を表す。島間の地理依存、種間の系統依存、学習データやモデル推定の不確実性を全て解消する区間ではない。過去の記録は説明変数に保持される。港の輸送指標の被覆不足や後年の港位置情報も残る。自然地理群の改善を特定の海峡・山地の因果効果と同一視しない。

図：whole_island_uncertainty.png。島を一つずつ除外した全結果：whole_island_leave_one_out.csv。

今回の全12比較で95%再標本化区間はゼロを含んだ。自然地理の点推定での改善は保持するが、島や種を越えた確実な改善とまでは判断しない。港追加の点推定は種等重みと島等重みで符号が異なり、港の追加価値も未確定。既存の植物・港・チョウの組成相関を否定する検定ではない。
