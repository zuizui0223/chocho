# GBIF Backbone照合：2006年スリランカ固有蝶18学名 × 元の225チョウ（2026-11-17）

**実行完了。** [GBIF公式 API v1：18学名の照会・応答SHA記録・原225種との比較・否定例テストがすべて成功](https://github.com/zuizui0223/chocho/actions/runs/44082761825)。この解析は**歴史的固有チョウの名前の解析選択枠**を確認するもので、**現代の種固有性、在来地、繁殖定着、個体数を調べた結果ではない**。

## 問いと結果

「2006年の歴史的E（スリランカ固有）学名が、単なる名前の違いで元のHOSTS選択225種に紛れ込んでいないか」を直接検査した。元の[18 E付きの学名と印刷頁リスト](SRILANKA_2006_HISTORICAL_ENDEMIC_E_18_NAME_MATCH_20261117.csv)と、固定の[225蝶×51陸塊GBIF source](../../data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv)は元Git blob SHAで検証。2026-11-17付のGBIF API /v1/species/match を Animalia/SPECIES/strict=true で18回独立照合し、EXACT、SPECIES、ACCEPTED/SYNONYM、confidence≥90の条件を満たす場合にのみ受理分類キーを225のspeciesKeyと比較した。

| 照合結果 | 2006年の18学名 |
|---|---:|
| GBIFの高信頼受理分類キーに対応し、そのキーは元225種にない | **12** |
| あいまい・不一致・確証不足：分類再精査が必要 | **6** |
| 受理分類キーを介して元225種に一致 | **0** |

**高信頼で元225にない12学名と受理分類キー：**
Troides darsius (1936992), Pachliopta jophon (1939171), Kallima philarchus (6230765), Lethe daretis (1918458), Mycalesis rama (8899821), Ypthima singala (1907507), Arhopala ormistoni (1927134), Tajuria arida (1922082), Nacaduba ollyetti (1932012), Nacaduba sinhala (1931902), Prosotas noreia (1926456), Halpe decorata (1942704).

**追加分類精査が必要な6学名：** Lethe dynaste, Elymnias singala, Spindasis nubilus, Jamides coruscans, Udara lanka, Halpe egena。歴史学名が現行分類で異名化・曖昧になっている可能性を残すため、無根拠の属/亜種変換・高次階級への代入をしない。

「18個の学名が単純に一致しない」という結果を、**そのうち少なくとも12は現在のGBIF受理分類キーで照合しても元225種の外側にいる**と強化できた。2006出版当時の固有種分類を2026年の固有種分類と同一視せず、18中6は未解決。元論文の2006年の総数20固有種も、現在の完全な固有蝶数や18種の網羅性を意味しない。

- [18学名とソースを記録したGBIF accepted key判定表](SRILANKA2006_E18_GBIF2026_TAXONKEY_RESOLUTION_20261117.csv)
- [実行したソースピン付きAPI照合コード](../../scripts/reconcile_SriLanka2006_historic_endemics_to_GBIF2026_20261117.py)
- [成功 GitHub Actions #44082761825](https://github.com/zuizui0223/chocho/actions/runs/44082761825)、アーティファクト14327897992、SHA256 4268c6ee8447b3c9b0df7fd9b2aa4afc336ac53671b6f765328c88be414a058f。GBIFの18生レスポンスと元応答SHAも含まれる（期限2027-02-15）。
- 原典：[Perera & Bambaradeniya 2006, pp53–64, in *The Fauna of Sri Lanka*](https://lk.chm-cbd.net/sites/lk/files/2022-06/2006-030.pdf)；[GBIF species API仕様](https://techdocs.gbif.org/en/openapi/v1/species)。

## 独立した生態学的な補強：2006年の寄主栽培・植生改変の記述

同じIUCN2006章の**印刷p.55**は、過去約20年間の土地改変と寄主植物の増加に伴うチョウの増加を、6種について**著者の知見として記述**している。これは標準化した個体数の長期時系列や交絡調整済み因果推定ではないため、**質的な歴史的地域説明・事例**としてのみ用いる。

| 2006年で増加に言及されたチョウ | 植物資源・人為環境（原資料で言及） | 元225種に掲載 |
|---|---|---|
| Papilio polymnestor | 家庭・農耕景観のミカン類やその他Rutaceaeの供給 | **なし** |
| Catopsilia pomona | 街路植栽の Cassia fistula | **あり** |
| Eurema hecabe | 庭園の Sesbania grandiflora、攪乱地の Cassia tora | **あり** |
| Eurema blanda | 同上の植物例 | **なし** |
| Euploea core | 生垣としての Nerium oleander | **あり** |
| Zizina otis | 攪乱地におけるイネ科草本の増加 | **あり** |

これらの寄主植物が全て外来かどうかはp.55に同じ地域原産フラグとして記されていない。**人の植栽・攪乱による寄主供給 ≠ 必ずしも外来植物の移入**。また6種が島に新たに到着したとも書いていない。人口増の定量や反実仮想はなく、既存GEB全球+54.9%、島+118.9%の種×植物地域の**潜在資源**という直接測定を外部事例が質的に支持する位置づけに限る。特に *Eurema hecabe* と *Zizina otis* は今回のSri Lanka陸塊でGBIF後期に初索引されたにもかかわらず、2006年著者が既に寄主植栽に伴う局所的増加を論じており、**GBIF初索引時期と生物学的な資源利用の時間軸のずれ**を明確にする。

本レポートは、[4島33文献照合／後期GBIF初索引317の監査](FOUR_ISLAND_33_HISTORY_VS_GBIF_INDEX_JA_20261117.md)に追加する補助証拠であり、元の原稿の潜在的資源拡大・港湾と外来植物相・海洋抵抗の主張は変更しない。
