# 地形と初記録時間差の現状

2026-10-30。海上距離のモデル比較、港・輸送・独立した導入年代の検証を含む目標は継続中。

## 地形追加：予測改善は確認できなかった
同じ155種・208,482評価組で、陸塊までを含むB2に標高中央値・起伏指標・地形不明フラグを足したB4を比較した。

時間外の種平均MSE差はB4−B2で+0.000297（固定予測の95%再標本化区間+0.000008〜+0.000515）。点推定・区間とも悪化方向だった。空間外は+0.000075（−0.000253〜+0.000383）で、追加的な改善は確認できなかった。モデル再学習の不確実性は含まない。

この結果は山が移動を妨げないという証拠ではない。今回の説明変数は候補セル内25地点から作った標高・P90−P10であり、過去の分布から候補セルまでの峠・山脈や、実際の飛翔経路を測ったものではない。

Copernicus90 DSMを648タイルから取得。海上サンプルを除くと、候補5,282セル中5,027セルで10地点以上の陸上標高を得た。残り255セルは不明フラグを付けて保持した。DSMには樹冠・建物が含まれ、当時に発行された予測ではなく、現代の地形代理変数を使った遡及比較である。

## 時間差：地域不整合を除外した現行集計

2026-11-05再照合。独立初記録資料の265組、チョウ9種・9地域・17種×地域が対応した。植物先行178組、チョウ先行68組、年代未確定16組、区間重複3組。Canada/Greveの地域不整合35組を隔離しており、旧300組を現行値として使わない。

種×地域単位では13/17群に植物先行記録、8/17群にチョウ先行記録があり、7群は両方を含む。食草が複数あるため、265組を独立した到来イベントとみなさない。最も早い対応食草との時間差図は14群が年代比較可能、3群が年代不足。線は原資料の年代幅であり、統計的信頼区間ではない。

時間差は初記録の差で、実際の導入・到来・定着の待ち時間ではない。未収録の在来食草や栽培植物、記録開始の差が残る。原典全行の確認は未完了。この9種を広域225候補／155解析種の代表値として扱わず、予測セルへの時間差説明変数としても使用していない。

出典：Alien Species First Records v3.1、DOI 10.5281/zenodo.10039630。FirstRecord_origの年代区間を使用。元資料の地域名は予測セル・評価島と同一視しない。メタデータとワークブックの利用条件の差は出典監査に保持する。

DEM出典：Microsoft Planetary Computerのcop-dem-glo-90コレクション。利用条件： https://dataspace.copernicus.eu/sites/default/files/media/files/2025-06/copernicus_contributing_mission_data_access_v2_cop_dem_licenses.pdf

produced using Copernicus WorldDEM™-90 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved.

The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM™-90.
