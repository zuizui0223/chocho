# 島全体を学習応答から除外した予測検証

225候補・155種の主解析を保持した追加評価。対象は25島・95種・41,812行。2019–2022年の標準化された記録確率を予測する。各島の面積に交差する全セルを学習から除き、M5（食草・気候・土地利用など）、B5（M5＋過去の記録接続性・海・地形）、J1（B5＋港・輸送指標）を各foldで再学習した。

|対象|重み|モデル|MSE|M5比改善率|
|---|---|---|---:|---:|
|all_islands|species_equal|M5_host_composition|0.024289|+0.00%|
|all_islands|species_equal|B5_sea_terrain|0.022675|+6.65%|
|all_islands|species_equal|J1_natural_transport|0.022079|+9.10%|
|all_islands|island_equal|M5_host_composition|0.035184|+0.00%|
|all_islands|island_equal|B5_sea_terrain|0.032485|+7.67%|
|all_islands|island_equal|J1_natural_transport|0.032655|+7.19%|
|without_GB_Ireland|species_equal|M5_host_composition|0.029829|+0.00%|
|without_GB_Ireland|species_equal|B5_sea_terrain|0.027149|+8.99%|
|without_GB_Ireland|species_equal|J1_natural_transport|0.026713|+10.45%|
|without_GB_Ireland|island_equal|M5_host_composition|0.035922|+0.00%|
|without_GB_Ireland|island_equal|B5_sea_terrain|0.033008|+8.11%|
|without_GB_Ireland|island_equal|J1_natural_transport|0.033205|+7.56%|

島等重みは各島内の行平均誤差を平均する。種等重みとは異なる評価量であり、相互に代用しない。英国・アイルランドを除く結果も全対象を保持して集計した。

これは過去に未観測の島への予測ではない。予測時点までの島内記録は説明変数に利用する。港指標の被覆は英国・アイルランドに偏り、港位置・DEMの年代整合にも制約がある。自然地理群は記録接続性・地形・海を合わせたもので、個別の障壁効果や船による持ち込みを同定しない。改善率は点推定で、不確実性の検証は別途必要。

検証：固定テストキー、元パネルの応答、fold別の面積交差による学習除外行数、予測値からの種等重みMSE/MAEを独立に再計算した。
