# 2島の独立歴史記録によるGBIF初索引監査（2026-11-17）

**状態：再現計算CI成功。** 既存フーコック島12組にモーリシャス島1993年総説の同名種3組を独立追加。固定した225蝶種×51港関連陸塊（11,475組）の2011–2022年GBIF初索引317組のうち、**2島・15組の2001–2010年以前の島内記録**を確認した。残り302組は、過去の存在も真の不在も判定されていない。15/317は抽出2島のチェック率であり、全球の初索引誤認率の推定量ではない。

## 1. データ

- [元GBIF 225×51×3期間](../../data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv)。Git blob SHA b6a3c98eb76864052d45d520a55b7d830816126c。
- [元51陸塊のWPI代表港](../../results/island51_host_circuit_atlas_20261106/island_geography.csv)：フーコック=陸塊471, DUONG DONG (VN); モーリシャス=陸塊218, PORT LOUIS (MU)。国の記録を島内記録に代用しない。
- [フーコック島12種・2006/2007文献](HISTORICAL_FIELD_CHECKLIST_GBIF_FIRST_INDEX_AUDIT_JA_20261117.md)、[12行の根拠](PHUQUOC_PRE2011_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv)。
- 新規モーリシャス文献：Davis, P.M.H. & Barnes, M.J.C. (BHL書誌刊行日1993-12-01), *The Butterflies of Mauritius*. Journal of Research on the Lepidoptera 30(3–4):145–161. [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/part/332210)、[DOI 10.5962/p.332210](https://doi.org/10.5962/p.332210)。論文本体には1991の印刷表記もあるが、保守的な先行存在年の上限は書誌日1993とした。元研究者による1976–1980の現地観察と古い文献の引用も含まれる。
- [新規モーリシャスの3行根拠表](MAURITIUS_1993_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv)。Git blob SHA e06fb0032b8dd2a9bd6344e2e88db0f201f3a4ee。

## 2. 1993年モーリシャス総説とGBIF初索引5種

| GBIF元の初索引種名（すべて前期2001–2010は0） | GBIF初索引窓 | 独立歴史資料 |
|---|---|---|
| *Lampides boeticus* | 2011–2017 | **1993年以前に島内記録**。原典p.154、Trimen (1866) の島内記録にも言及、園芸寄主・作物との関連を記述 |
| *Zizula hylax* | 2011–2017 | **1993年以前に島内記録**。原典p.154、Manders (1907) の記録に言及。小型で見落とされやすいとの注記 |
| *Danaus chrysippus* | 2011–2017 | **1993年以前に島内記録**。原典p.147、島内で広く分布、1976–1980の観察を記述 |
| *Zizina otis* | 2011–2017 | **分類体系未照合**。原典の*Z. antanossa*は、Yago et al. (2008)で*Z. otis*の亜種として扱われるが、元GBIF speciesKeyに自動的に置き換えない |
| *Zerynthia polyxena* | 2018–2022 | **原記録未確認**。欧州分布で知られる種。元GBIFの同定・座標・原標本を点検するまで真の生物学的出現とも誤記録とも断言しない |

- 分類概念の出典：[Yago et al. 2008, Zootaxa 1746](https://doi.org/10.11646/zootaxa.1746.1.2)。
- *Zerynthia polyxena* の地域情報：[Naturalis, Butterflies of Europe](https://butterflies-europe.linnaeus.naturalis.nl/linnaeus_ng/app/views/species/taxon.php?id=114545)。

## 3. 元51陸塊×1,275ペアで再計算

GBIF後期1,053累積種×陸塊セル、後期317初索引セルと原51陸塊は変更しない。

| 前期2001–2010既知種リスト | 種×陸塊セル | 平均Jaccard | 2022累積GBIFとの差 |
|---|---:|---:|---:|
| 元GBIFのみ | 736 | 0.120963255 | +0.056809620 |
| フーコック2006/2007の12組を補強 | 748 | 0.124359023 | +0.053413852 |
| さらにモーリシャス1993の3組を補強 | 751 | 0.126494594 | **+0.051278280** |

15組の文献追加による増加量の差は **−0.005531340**。これは選択島の古い記録を異質な資料から追加した**記述的な前期ソース拡張**。真の蝶群集の変化量の不偏推定、全球のGBIF初索引誤判定率、定着や移入の割合ではない。2島の文献は統計的に独立な15個の島嶼調査ではない。

## 4. 在来性・固有性を測る母集団が違う

1993年のモーリシャス総説には*Euploea euphon*が同島の固有種と明記され、*Cyclyrius mandersi*も同島に限られると記述されている（p.148、p.154–155）。この2学名は元の225蝶種に含まれていない。前回のクリスマス島*Polyura andrewsii*の欠落と合わせて、**HOSTSに対応する選択225種は島固有蝶の完全な代表標本ではない**。これらの1993年の歴史的分類を2026年の受理種・現存分布の確定値に直接置き換えず、固有種の全体数や捕捉率は計算しない。

同じように*Zizina otis*と1993年*Z. antanossa*の分類上の包含関係は、種同定を精査せず一律にコード化しない。モーリシャスの元GBIF初索引5種のうち、**確証3、分類未確定1、原記録未点検1**を別々に残す。

独立に判定すべき軸は、(a) GBIFに索引されたか、(b) 独立した島内の過去の種記録、(c) チョウ自身の自然／移入／迷行の由来、(d) 野外繁殖・定着、(e) 単一島固有・群島固有・非固有の在来分布。前の2種から後ろ3種を機械的には推定しない。

## 5. 再現と論文への扱い

- [ソースSHA固定コード](../../scripts/audit_two_island_independent_butterfly_history_20261117.py)で原225×51、過去2島の独立文献表、各学名、元GBIF索引年を突合。
- [GitHub Actions #44076800265 成功](https://github.com/zuizui0223/chocho/actions/runs/44076800265)：317件の注釈、数値整合、分類概念を勝手に同一化する入力、根拠なしの在来・固有ラベルへの変換を拒否。
- 一本化GEB原稿のFigure 3CはGBIFの記録種リストを意味し、真の在来種群集均質化ではない。Supplementary Results S.Fの歴史資料監査として使用。既存HOSTS×WCVP潜在資源の全球+54.9%、島+118.9%と港・海洋抵抗の結果は維持する。


## 6. 公式GBIFで未判定2件の分類・データセット来歴を照合

[2026-11-17公式GBIF 7 endpoint小規模照会（成功）](https://github.com/zuizui0223/chocho/actions/runs/44077684191)、[監査JSON](MAURITIUS_OFFICIAL_GBIF_TAXON_OCCURRENCE_SOURCE_AUDIT_20261117.json)。元の225×51分析表・317初索引・確証15組は変更しない。

**Zizinaの分類体系**：GBIF現行バックボーンは *Z. antanossa* を**ACCEPTED speciesKey 8177775**、*Z. otis* を**別のACCEPTED speciesKey 7346444**と返した。Yago et al. (2008) が提案した亜種扱いを、公式GBIFバックボーンは一律には採用していない。したがって歴史的な*Z. antanossa*の島内記録を、元の*Z. otis*の前期記録へ付け替えることはまだできない。

**Zerynthiaの唯一の2022年モーリシャス記録**：GBIF occurrence ID **5023737393**、2022-10-11、緯度−20.038742、経度57.545460、座標不確実性25m、HUMAN_OBSERVATION。出版データセットは**NABU|naturgucker**（DOI [10.15468/uc1apo](https://doi.org/10.15468/uc1apo)）、原参照URLは [enjoynature.netの掲載元](https://enjoynature.net/?gebiet=-509720845&datum=11.10.22&art=1999567171)。GBIFのメディア配列は0、identifiedByは空。座標の丸め・WGS84仮定・座標由来の大陸判定という付随フラグがあるが、これだけでは誤同定を証明しない。**確認できたのは1件の市民科学の掲載記録と来歴であり、ヨーロッパ分布が基本の蝶がモーリシャスへ自然分散・人為移入・定着したことではない。** 画像・標本・原報告者の同定根拠を照合できるまで、要原同定検証のまま保持する。

監査ZIP artifact ID **14324344278**、SHA256 **4d1cb09306533df1e28d9212db9699d7cfec3964dd828f13944798658d7678ec**（保存期限2027-02-15）。[実行コード](../../scripts/audit_mauritius_GBIF_taxonomy_and_firstindex_20261117.py)は公開GBIF APIの小規模7照会の記録であり、真の繁殖・在来性・固有性を推定しない。
