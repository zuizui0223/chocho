# 自然地理・輸送モデルの観測努力と提供者感度

225候補を維持し、各条件で適格な155種を全て評価した。同じ条件のM5とT1の予測値・評価キー・応答の一致を確認し、B5/J1を固定RF設定で再学習した。値は種等重みMSEの改善率。正なら比較モデルより誤差が小さい。

|条件|追加モデル|比較対象|評価行|改善率|
|---|---|---|---:|---:|
|temporal|B5_sea_terrain|M5_host_composition|208482|-1.152%|
|temporal|J1_natural_transport|M5_host_composition|208482|-1.134%|
|temporal|J1_natural_transport|B5_sea_terrain|208482|+0.018%|
|temporal|J1_natural_transport|T1_port_connectivity|208482|-1.362%|
|temporal_effort20|B5_sea_terrain|M5_host_composition|159084|+0.571%|
|temporal_effort20|J1_natural_transport|M5_host_composition|159084|+0.989%|
|temporal_effort20|J1_natural_transport|B5_sea_terrain|159084|+0.421%|
|temporal_effort20|J1_natural_transport|T1_port_connectivity|159084|+1.241%|
|temporal_exclude_top_provider|B5_sea_terrain|M5_host_composition|191004|-2.163%|
|temporal_exclude_top_provider|J1_natural_transport|M5_host_composition|191004|-2.808%|
|temporal_exclude_top_provider|J1_natural_transport|B5_sea_terrain|191004|-0.631%|
|temporal_exclude_top_provider|J1_natural_transport|T1_port_connectivity|191004|-1.836%|

条件ごとに応答となる標準化記録確率が異なるので、条件間の絶対MSE差を生物学的効果と解釈しない。最大提供者除外は観測努力・応答の変更であり、過去距離・候補領域・食草カタログ・地理の全入力を提供者なしで再構築した検証ではない。

本表は時間外の全155種の比較。島保持・空間外・種外の感度を追加実施したものではない。港距離と指数は実際の船舶運搬や食草利用を識別しない。

## 固定予測誤差の不確実性

種×空間ブロックを独立に再標本化した2,000回の95%百分位区間。通常解析の点推定を見た後の診断で、モデル再学習・系統や地理依存を完全に扱う区間ではない。負のMSE差が改善。

|条件|追加モデル|比較対象|MSE差|95%区間|
|---|---|---|---:|---|
|temporal|B5_sea_terrain|M5_host_composition|+0.000354|[-0.001335, +0.001983]|
|temporal|J1_natural_transport|M5_host_composition|+0.000349|[-0.001456, +0.002039]|
|temporal|J1_natural_transport|B5_sea_terrain|-0.000006|[-0.000272, +0.000321]|
|temporal|J1_natural_transport|T1_port_connectivity|+0.000418|[-0.001215, +0.002005]|
|temporal_effort20|B5_sea_terrain|M5_host_composition|-0.000233|[-0.003509, +0.002603]|
|temporal_effort20|J1_natural_transport|M5_host_composition|-0.000404|[-0.003344, +0.002217]|
|temporal_effort20|J1_natural_transport|B5_sea_terrain|-0.000171|[-0.000667, +0.000377]|
|temporal_effort20|J1_natural_transport|T1_port_connectivity|-0.000508|[-0.003575, +0.002134]|
|temporal_exclude_top_provider|B5_sea_terrain|M5_host_composition|+0.000694|[-0.001190, +0.002426]|
|temporal_exclude_top_provider|J1_natural_transport|M5_host_composition|+0.000900|[-0.000924, +0.002640]|
|temporal_exclude_top_provider|J1_natural_transport|B5_sea_terrain|+0.000207|[-0.000085, +0.000576]|
|temporal_exclude_top_provider|J1_natural_transport|T1_port_connectivity|+0.000594|[-0.001389, +0.002413]|

全12比較で区間はゼロを含んだ。全155種で一貫した追加改善は確認できず、点推定の符号も観測条件によって変わる。島部分集団の改善を広域全体に一般化しない。この結果は、独立した食草情報の予測改善や島間の組成相関を否定する検定ではない。
