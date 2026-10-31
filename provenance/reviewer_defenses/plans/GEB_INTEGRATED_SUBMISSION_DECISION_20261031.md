# GEB統合投稿の設計・誌面判断（2026-10-31）

**推奨：GEB (Global Ecology and Biogeography) Research Article に一報として投稿する場合、239種の資源機会+54.9%と155種の後年・空間外予測改善を主結果とし、51島・225種の島港相関は機構を主張しない地理的な探索結果として扱う。第一候補GEB、第二候補Ecography。旧v0.2は上書きせず保存。ここに示す新稿は未投稿。**

## 1. 投稿誌と主張の相性

| 優先順位 | 雑誌 | 主張の編集上の重心 | 現在の要注意点 |
|---|---|---|---|
| **第一** | **Global Ecology and Biogeography** | 相互作用相手の人為的再分配がチョウの潜在資源地理を変える、という全球マクロ生態学 | 港で蝶の定着が促進されたとは言えない |
| 第二 | Ecography | 寄主地理による時空間的な観測記録予測の改善 | 既知種の後ろ向きホールドアウト。種外移転は未支持 |
| 条件付き | Journal of Biogeography | 地理的障壁や島の生物地理的過程 | 港外対照、独立歴史的輸送・在来性・供給源の検証が不足 |
| 条件付き | Diversity and Distributions | 外来寄主植物の移入・侵入生態学や保全上の分布変化 | 現地の幼虫摂食や新規定着への独立した影響が未確認 |

公式範囲：
- https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html
- https://onlinelibrary.wiley.com/journal/16000587
- https://onlinelibrary.wiley.com/page/journal/13652699/homepage/forauthors.html
- https://onlinelibrary.wiley.com/page/journal/14724642/homepage/forauthors.html

GEB Research Article の代表的な目安は**本文約5,000語、合計6〜8点程度の図表、構造化要旨300語以下**。GEBは生態学的な一般性を重視し、double-anonymous査読を採用。Aim / Location / Time period / Major taxa studied / Methods / Results / Main conclusions の要旨見出しを使用する。

## 2. 統合論文の唯一の中心的な問い

**人為的に再配置された幼虫寄主植物はチョウの資源機会の地理をどれほど変え、その資源地理は実際の蝶の記録分布とどの程度対応するか。**

3つの段階を混ぜない。
1. 全球の潜在的な資源機会を再構築する：239種中206種で拡大。蝶×WGSRPD3地域の延べ機会が26,530→41,083（+14,553、**+54.9%**）。38植物種が増加の半分を担う。非作物の感度でも+46〜47%。これは蝶の実生息域面積ではない。
2. 別の蝶記録データによって予測上の意味を検証する：155種・208,482組、気候+土地被覆のM3に寄主資源と科組成を加えたM5で種等重みMSEは時間外**4.12%**減、未学習空間**3.00%**減、100km以上の事後感度で**2.74%**減。4感度でも同方向。ただし静的な現代植物地理を使う後ろ向き検証であり、植物導入の因果効果・真の占有ではない。
3. 島/港で記録分布と寄主候補の地理的な関係を評価する：51島・1,275非独立島ペアで港寄主候補相×新規蝶相 Spearman **ρ+0.300**、島間距離と以前の蝶相調整後**+0.215**。以前からの蝶相への相関はむしろ**ρ+0.418**。共通225種でGEBの種別追加資源機会と寄主候補の記録島港数**ρ+0.652**。植物の元の記録は2001–2024年プールなので時系列の先行効果は主張できない。

**強い主張として残す：** 人為的な植物再分配は既知寄主との潜在資源地理を拡大し、寄主資源の地理情報は後年・他地点の蝶記録の予測に価値がある。島・港は寄主候補の地理と蝶の記録分布が対応する空間的な接点である。

**まだ主張しない：** 港で輸送された植物が幼虫の餌になって定着に成功した、GEBの55%のうち何%が島/港で実現した、港間が実際に中継源になった。17島・113件の先行寄主一致は、植物索引頻度・地理を保った帰無に対し独立した追加効果が残らない。港の接続性と寄主候補数の単調相関は39島港でρ+0.023。

## 3. 編集順と主図の構成

**今のGEB v0.2は本文約4,943語、要旨約292語、主図3枚、SI図1枚。追加結果をそのまま後ろに足すと散漫になる。既存の感度解析・気候フィルタ・専門性の細部をSIに移し、結果の順番を組み替える。**

| 節 | 新しい主結果 | 主な図 |
|---|---|---|
| R1 潜在資源の全球再配置 | +54.9%、206/239、食性の広狭を問わず拡大 | Fig 1：在来・導入分布、gainと食性 |
| R2 資源再配置の植物側偏り | 半分の追加分を38種が担う、作物除外後も46–47% | Fig 2：寄主寄与集中・頑健性 |
| R3 既知チョウの記録予測 | 155種・208,482、4.12%・3.00%・2.74%。既存GEBの66/115種×地域の資源対象復元と概念的につなぐ | Fig 3：M3/M5、時間外・空間外 |
| R4 島港の観測的対応 | 51島ρ+0.300、既往蝶相ρ+0.418、調整後+0.215。港225種ρ+0.652は補助的 | Fig 4：島間の新旧蝶相と寄主相。港種相関と厳密な帰無はSI |

GEBの旧Fig 2、3や細かな予測・植物分類・港接続性の感度図はSIへ。原本の結果を削除するのでなく本文の階層を変える。サンプル数の単位は239蝶種×地域、155種×記録評価組、51島×1,275島ペア、225種×港集計、17島×リスクセルで異なるので合算禁止。

## 4. 統合論文のタイトルと300語以内の要旨案

**推奨タイトル案：**
Redistributed larval hosts expand butterfly resource opportunity and improve occurrence prediction

**元GEB向けのより保守的なタイトル案：**
Anthropogenic host redistribution expands butterfly resource geography across the specialization spectrum

**Aim:** To ask whether redistribution of known larval host plants reorganizes butterfly resource geography, and whether host geography helps predict recorded butterfly distributions.

**Location:** Global, with an exploratory island–port subset.

**Time period:** Contemporary host-plant ranges and source-specific butterfly records from 2001–2026, including temporal evaluations for 2019–2022.

**Major taxa studied:** Butterflies and their larval host plants.

**Methods:** We compared native-only and contemporary host distributions for 239 butterflies, holding host identities fixed across botanical regions. We apportioned additional resource opportunities among plant species and tested crop-host exclusion. For 155 butterflies, we assessed whether host-resource predictors improved occurrence-record models beyond climate and land cover under temporal and spatial holdouts. We also compared port-associated candidate-host similarity and butterfly recording assemblages across 51 sampled islands.

**Results:** Introduced host ranges added 14,553 butterfly–region resource-opportunity units (+54.9%), expanding opportunity for 206 of 239 species; 38 host plants accounted for half the added opportunity. Adding host predictors reduced species-equal mean squared error by 4.12% in temporal holdouts (208,482 evaluated pairs) and 3.00% in withheld spatial blocks. Across 1,275 nonindependent island pairs, recorded port-associated candidate-host similarity covaried with 2011–2022 first-indexed butterfly assemblage similarity (Spearman ρ = +0.300; partial rank correlation after distance and earlier butterfly assemblages = +0.215). More restrictive plant-geography and butterfly-recording nulls did not establish an additional host-indexing excess.

**Main conclusions:** Redistribution of known interaction partners can expand potential consumer resource geography even when recorded diet identities are held constant. Host geography contributes predictive information about butterfly recordings, and is positively associated with island assemblage structure. These patterns do not establish shipborne host transport, actual larval use or butterfly establishment.

## 5. 提出前のゲートと分割案

- **再現性が最優先：** 現行GEBの匿名データ・コードZIPはv0.2の資料を対象に作られている。155種予測と島港解析の再計算に必要な入力は、別途GEB査読者が再現できる形の匿名版を作る。非公開著者用HOSTS原辺・元Parquet・植物原座標を許可なしに公開しない。現状で集計JSONしか存在しないprivate入力に依存する図を、そのまま投稿本論文の主結果として増やさない。
- **誌面：** 5,000語程度に収め、300語以内の構造化要旨、British Englishを統一、図表6–8点の目安、引用の一般性を再検証。長くなる場合、GEBのCover Letterで超過の必要性を説明する。
- **匿名査読・著者：** 新規本文、claim map、SI、Cover Letter、匿名DOCX、添付用データコードZIPを同一コミットで凍結し再CI。共著者全員の承認と著者欄・所属・CRediT・資金・COI・連絡先が必要。現在のGEB旧v0.2は科学的準備完了だがこのメタデータは未完。
- **迅速投稿ルート：** 時間優先なら、既存の独立GEB v0.2（+54.9%、38寄主種、66/115の観測対応）を著者情報確定後にGEBへ投稿。予測+島港の続編は、独立した島外対照と外来寄主の時間的分類が整ってからEcographyまたはJournal of Biogeography向けに構成。
- **重複投稿回避：** +54.9%のGEB元結果を2報で別々の新規発見として使わない。二報に分ける場合は別の問い・別の予測対象・別の明確な新規証拠と原GEB引用を確保する。
- **投稿前提：** この文書は編集判断と未査読の要旨・図の設計であり、統合GEB本文や正式な投稿ではない。現行の投稿用v0.2本文・匿名パッケージ・図は変更せず保持。

**最終の雑誌判断：広域の資源地理+予測を1本に統合するならGEBを第一投稿先にする。島と港の観測相関は正面から示すが、港による繁殖・定着という未識別機構を論文のタイトルや主要結論に書かない。**
