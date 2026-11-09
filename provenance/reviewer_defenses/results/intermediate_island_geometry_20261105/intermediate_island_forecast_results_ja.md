# 中間島情報の追加予測価値

原225候補・適格155種を保持し、B5（食草・気候・土地利用・過去の記録接続性・海・セル内地形）に、幾何的な中間島数と判定不明フラグだけを追加したB6を比較した。RF設定と評価行を固定し、比較B5/M5のキー・応答・予測値の一致、特徴量ハッシュ、MSE/MAEの再計算を検証。

|評価|対象|比較対象|MSE改善率|MSE差の95%再標本化区間|
|---|---|---|---:|---|
|temporal|all155|M5_host_composition|-0.608%|[-0.001571, +0.001863]|
|temporal|all155|B5_sea_terrain|+0.539%|[-0.000382, +0.000079]|
|temporal|area_supported_islands|M5_host_composition|+7.753%|[-0.005532, +0.001825]|
|temporal|area_supported_islands|B5_sea_terrain|+0.422%|[-0.000924, +0.000809]|
|temporal|islands_without_GB_Ireland|M5_host_composition|+7.343%|[-0.006135, +0.001425]|
|temporal|islands_without_GB_Ireland|B5_sea_terrain|+0.548%|[-0.001063, +0.000774]|
|spatial|all155|M5_host_composition|-1.890%|[-0.001465, +0.002539]|
|spatial|all155|B5_sea_terrain|-0.049%|[-0.000305, +0.000303]|
|spatial|area_supported_islands|M5_host_composition|+11.135%|[-0.006778, +0.001114]|
|spatial|area_supported_islands|B5_sea_terrain|+1.152%|[-0.000680, +0.000270]|
|spatial|islands_without_GB_Ireland|M5_host_composition|+10.626%|[-0.007438, +0.000467]|
|spatial|islands_without_GB_Ireland|B5_sea_terrain|+0.873%|[-0.000714, +0.000268]|
|spatial|at_least100km|M5_host_composition|-0.594%|[-0.001822, +0.001999]|
|spatial|at_least100km|B5_sea_terrain|-0.082%|[-0.000306, +0.000328]|

MSE差はB6−比較対象で、負が改善。改善率は正が改善。区間は固定予測誤差の種×空間ブロック再標本化2,000回で、モデル再学習の不確実性や独立した島の反復ではない。

時間外評価の島部分集計は、島全体を学習から除いた評価ではない。中間島は直線経路の幾何情報で、中間島の占有・食草利用・実際の中継移動を示さない。海上距離などを既に含むB5に対する追加価値と、食草モデルM5との差を分ける。

評価状況：{"temporal": "COMPLETE", "spatial": "COMPLETE"}. 未完了の評価は結果として扱わない。

現時点の完成済み比較では全区間がゼロを含み、中間島情報の追加効果は未確定。小さな改善の点推定を保持するが、島による中継機構の確認には使わない。
