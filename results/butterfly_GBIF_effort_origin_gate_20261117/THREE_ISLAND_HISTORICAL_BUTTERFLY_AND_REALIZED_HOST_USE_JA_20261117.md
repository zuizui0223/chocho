# 蝶のGBIF初索引と独立歴史資料・幼虫利用の照合：3島（2026-11-17）

**達成：** [一本化GEB原稿 PR #157](https://github.com/zuizui0223/chocho/pull/157) の既存225蝶種×51陸塊×3期間（11,475セル）を固定し、後期初索引317組のうち**3島・19組**が、2011年以前の独立島内資料に明記されていたことを照合した。**残り298組は未照合であり、新規移入・真の不在とは判定しない。** 栄養や港による定着機構の実証ではない。

## 1. 各島の裏付け

| WPI陸塊と代表港 | 元GBIF初索引の確認範囲 | 古い独立記録で同種を確認 | 主な原典 |
|---|---:|---:|---|
| フーコック島 471 / DUONG DONG (VN) | 2018–2022初索引12 | **12** | Wildlife At Risk 2006現地リスト9種＋2007図鑑3種 |
| モーリシャス島 218 / PORT LOUIS (MU) | 後期初索引5 | **3** | Davis & Barnes 1993島内蝶類総説、同名3種。2種は同定・分類未解決 |
| ケープブレトン島 79 / SYDNEY (CA) | 2011–2017初索引17 | **4** | Lepidopterists' Society 1948シーズン報告と1984年（1983年現地）ニュースレター |

ケープブレトン島で確認した4種は、*Colias philodice*（1948年の島内記録）、*Erynnis icelus*（1983年、Clyburn Brook/Freshwater Lake/North Mountain）、*Lycaena phlaeas*（1983年、McKenzie Mountain、原文では属名省略、2010年島内Atlas記録で別途学名を照合）、*Glaucopsyche lygdamus*（1983年、Grande Falaise/South Pond/Baddeck）。根拠表：[Cape Breton 4行](CAPE_BRETON_1948_1984_INDEPENDENT_FIELD_MATCHES_20261117.csv)。
- [1984 Lepidopterists' Society newsletter, 26(2), p.32](https://images.peabody.yale.edu/lepsoc/nls/1980s/1984/1984_v26_n2.pdf) には、1983年の具体的な島内地名がある。
- [1948 Lepidopterists' News field season summary, Nova Scotia](https://images.peabody.yale.edu/lepsoc/jls/1940s/1948/1948-2%28S%29-SEASON_SUMMARY.pdf) はCape Bretonで1948年6月20日に*Colias philodice*の存在を報告。刊行年の保守的上限を1950年として扱う。
- *Lycaena phlaeas*については1984年報告の属名が略されるため、[Maritimes Butterfly Atlasの2010年Cape Breton County同名・日付付き現地記録](https://accdc.com/mba/downloads/tabular-results/2010Lyca.htm)で学名を独立確認。少なくとも2010年より以前には島内にいた。厳格な「完全な原文二名法だけ」検定なら1984年記録を弱い名寄せとして別扱いする。

## 2. 新しい野外寄主利用の直接証拠

1983年、*Glaucopsyche lygdamus*の**幼虫**が同じケープブレトン島で2種の植物上に記録されていた。

| 幼虫が記録された植物と地点 | Kew POWO上のNova Scotiaでの区分 | 出典 |
|---|---|---|
| *Vicia cracca*（Grande Falaise） | **外来（introduced）** | [POWO: Vicia cracca](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A30263827-2) |
| *Lathyrus japonicus*（South Pond、Baddeck） | **在来（native）** | [POWO: Lathyrus japonicus](https://powo.science.kew.org/taxon/30381438-2) |

[1984年の独立現地観察原典](https://images.peabody.yale.edu/lepsoc/nls/1980s/1984/1984_v26_n2.pdf)の幼虫記録と、[地域植物由来の別ソース監査JSON](CAPE_BRETON_1983_LARVAL_HOST_ORIGIN_SOURCE_BRIDGE_20261117.json)を接続した。「外来の既知寄主が実際に幼虫に使われうる」という具体的な**実現利用例**。ただし、HOSTSが同じ公刊現地報告を引用している可能性はあるため、**HOSTSネットワークとは完全に独立した寄主相互作用データとは呼ばない**。幼虫生存、成虫羽化、個体数、港経由の伝播、長期の繁殖定着やチョウ自身の在来性は確認していない。植物の由来は**Nova Scotia地域単位**のPOWO分類であり、1983年の各野外地点で個別の来歴を測定したものではない。

## 3. 原225種×51陸塊・1,275島ペアを固定した感度分析

後期GBIF累積1,053種×島セルは不変。**前期に記録済みとする文献ソース集合だけ**段階的に補強した。

| 前期の既知種ソース | 前期セル | 平均前期Jaccard | 同じ後期0.177773との差 |
|---|---:|---:|---:|
| GBIF 2001–2010のみ | 736 | 0.120963 | **+0.056810** |
| フーコック島2006–07の12種を追加 | 748 | 0.124359 | +0.053414 |
| モーリシャス島1993の3種も追加 | 751 | 0.126495 | +0.051278 |
| **ケープブレトン島1948/1983の4種も追加** | **755** | **0.127075** | **+0.050698** |

元のGBIFのみの増分から **−0.006112 Jaccard**。異なる古い調査とGBIF累積索引を組み合わせた**記述的な先行存在のソース補強**であり、真の蝶群集均質化の不偏補正・潜在的移入割合・固有蝶の変化量ではない。サンプリングした3島は無作為標本ではないため19/317を全球の誤認率に変換しない。

## 4. 1987年の追加原資料と取得可能性

[Lafontaine et al. (1987) カナダ政府のCape Breton Highlands National Park昆虫報告](https://publications.gc.ca/site/eng/9.826269/publication.html)は、当該公園の蝶類を種別に記述することは書誌・目次から確認できる。しかし**外部のGitHub Actionsランナーでは政府PDF URLがHTMLを返し、昆虫学会ミラーは可用性が変動**しており、元資料本文の全ページ抽出は成功していない（[診断CI](https://github.com/zuizui0223/chocho/actions/runs/44079613150)）。未取得部分の種を、印刷本文の同定を確認したものとして追加しない。試験した抽出コードは[こちら](../../scripts/audit_capebreton_original_1987_butterfly_names.py)。継続的な無駄な自動再取得は行わない。

## 5. 研究上の主張と再現性

- 成果は**真の新規定着ではなくGBIF初索引に以前から存在した種が含まれる**ことを、3島と幼虫現地記録で示すこと。317件のうち独立歴史資料で検証済み19、未照合298。未照合≠不在。
- 報告の在来・移入・迷行、繁殖・長期定着、単一島・群島固有性はすべて別。文献に島内幼虫があってもそのチョウが本来そこに自然分布したかは分からない。
- [SHA固定の全317照合スクリプト](../../scripts/audit_three_island_independent_butterfly_history_20261117.py)で元11,475行、原51島、後期317初索引、3島19の厳密キー照合、幼虫の植物名・分類未昇格、Jaccardと島単位固定を検証。
- [CI #44079725593 SUCCESS](https://github.com/zuizui0223/chocho/actions/runs/44079725593)：数値検証とわざと偽の在来・固有性・寄主名を書いた負のテストも合格。
- 既存GEB本文（潜在資源の全球+54.9%、島+118.9%、港湾と外来植物相、海洋距離/抵抗の肯定結果）は保持。Figure 3Cは「**GBIF累積記録リスト**」、Supplementary S.Fに歴史的存在と現地の両寄主幼虫利用を追加し、野外寄主利用の一般性や真の蝶定着へ飛躍しない。
