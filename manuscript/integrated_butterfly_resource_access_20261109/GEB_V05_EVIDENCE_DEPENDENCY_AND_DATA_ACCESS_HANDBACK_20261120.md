# GEB v0.5 — 主張ごとのデータ依存性と査読用データ公開案（2026-11-20）

**判定：科学的原稿は共著者レビュー可能。完全な査読者用匿名データ・コード公開パックは未作成。**
対象の原稿は [PR179](https://github.com/zuizui0223/chocho/pull/179) の固定 HEAD \`99b876421af2b4f51e58569862994622a02717d7\`。ここでの出典監査は [PR189](https://github.com/zuizui0223/chocho/pull/189) の独立Draft上に追加し、本文・図の結果、未検証のDOI、著者情報を変更しない。

## 1. 現在の原稿の主結果と出典の依存関係

| 本文の部分 | 科学的応答と集計単位 | 原資料・出典の状態 | 1,596港の植物GBIF原キャッシュに依存するか |
|---|---|---|---|
| Fig.1、結果3.1 | 239種×126島・237本土植物地域の**潜在的な既知幼虫寄主資源**：島の増加率118.9%、全球14,553追加単位 | HOSTS \`10.5519/havt50xw\` + WCVP v13 \`10.34885/0yex-xv26\`；数値と図のソースレシート固定 | **いいえ** |
| Fig.2A、結果3.2 | 126植物地域のWPI港数と**在来/外来全植物相の種数**の関係（港×外来全植物相 β+0.2854） | WPI2019 + WCVP v13；面積・Level1調整、HC3、同一島の起源対比 | **いいえ**。GBIF港バッファの観測努力調査と別データ |
| Fig.2B–C、結果3.2 | 在来/外来相の島間類似性、104島海洋抵抗、22島EuDiS形質 | WCVP+原海洋幾何+EuDiS、source-only回帰と感度解析 | **いいえ** |
| Fig.2D、結果3.2 | 126島×239蝶種の本土距離と原/現在の潜在寄主機会、123島の外来のみ寄主割合 | WCVP+HOSTS+WPI、島ジオメトリ、cluster/Moran/MEM10；数値ソース再実行CI済 | **いいえ** |
| Fig.3A–C、結果3.3の主要部分 | 670既知寄主植物の寄与集中、239蝶種の寄主地理、97種の単一寄主科群 | WCVP+HOSTSの固定種・地域集計；original/source audit済 | **いいえ** |
| Fig.3の**補助的な港湾観測文** | 上位10寄主の8種を17島・38港で早期GBIF植物照会検出、より広い1,596港中1,150港で外来候補既知寄主を索引 | 2001–2010年GBIF植物検索と、元の1,596港\`speciesKey\` facet＋当時のWKTバッファ。**元の個別GBIF ID/\`datasetKey\`は未保持** | **はい**。但し直接の実食、定着、船舶輸送の証拠ではない |
| Fig.4、結果3.4 | 独立蝶記録の155種・208,482事例、および224種×169地域の群集記録 | **同じ元GBIFダウンロードDOI \`10.15468/dl.umqmsf\`**、14,637,152件；公式GBIF提供元**1,716件**を照合・保存済。出現原個票から派生表への完全復元は未公開 | **いいえ** |
| 225×51島の初索引・83,573行の歴史年表 | 島間記録年順の状況的根拠（初索引≠定着） | 異なる新規GBIF facetと非公開の古い集計Parquet。元のGBIF提供元は未復元 | **いいえ**。ただしこちらも別の出典問題あり |
| 独立の幼虫事例 | Canary固有蝶の幼虫による移入植物の直接利用 | 原写真・植物同定・権利確認が必要；個体群増加は未実証 | **いいえ** |

**結論：** 元の1,596港GBIF植物キャッシュの提供元帰属が欠落していても、Fig.1、Fig.2A–D、Fig.3の全球潜在資源主要結論、Fig.4の元ダウンロード出典はそれに**直接依存しない**。一方、Fig.3の**港湾現地照会補助文**と対応するSupplementの原観測については未完成のGBIF出典を正直に記載しなければならない。必要に応じ、共著者はこの補助文の表現・位置を検討できるが、ソースが無いことを理由に主要な肯定結果を消したり別のDOIを転用しない。

## 2. 2026-11-20の新しいGBIF出典監査

- **1,596/1,596原WPI港ID**が、新しいGBIF \`datasetKey\` facetにより照会完了。植物索引あり1,461港／なし135港、**1,482提供元データセット**。5km港バッファ索引ヒット総和**3,138,675**は港間重複を含む。原2019WPI中心座標は公開ミラーと**小数6桁まで一致**。
- 最新資料：[新1596港×datasetKey出典表](../../results/geb_v05_gbif_source_attribution_20261120/NEW_GBIF_ORIGINAL1596_PORTS_DATASETKEY_INDEX_20261120.json)と[新スナップショット監査](../../results/geb_v05_gbif_source_attribution_20261120/README_NEW_ALL1596_PORT_GBIF_PROVIDER_SOURCE_20261120.md)。GBIFの元\`speciesKey\`WKT照会とは**取得時点と検索ジオメトリが異なる**ため、旧出典の帰属には使えない。
- **新しい品質診断（PR189のみ）**：元の1,596港の「在来既知寄主」「現代に存在するが在来でない候補既知寄主」の絶対収録種数は、現在GBIF植物索引件数と強く連動（Spearmanρ **+0.842 / +0.821**）。国・領域コード142群の固定効果／国クラスタ誤差で、log1p記録数に対する絶対収録種数の係数はそれぞれ**+0.392**（95%CI +0.324〜+0.461）と**+0.414**（+0.360〜+0.467）。しかし2分類の**対数比**の係数は**+0.0275**（−0.0370〜+0.0919；p=.404）。感度4条件では比率係数は小さく、0を含む。これは観測努力の原理的な交絡を示す**探索的な異時点ソース品質診断**であり、国固定効果は観測機会を因果的に標準化したわけではない。
- 比較している元港湾候補は**在来の反証でない「contemporary-not-native-in-catalogue」**であり、専門家が独立に外来定着を認定した植物の集合ではない。新しい提供元数と記録数にもSpearmanρ **0.876** の強い共線性があり、提供元数を同時投入した弱い係数は部分群で不安定。
- **重要**：この品質診断はFig.2Aの126島×WCVP全植物相×WPI港数回帰の**補正解析ではない**。主張を後から取消したり、逆に主モデルの頑健性を実証したことにしてはいけない。

## 3. 次の論文データ公開上の決定

| 問題 | 現在決定できること | 残る外部入力・操作 |
|---|---|---|
| 155/224蝶の同一元GBIF DOI | **閉鎖：** \`10.15468/dl.umqmsf\`、1,716 original \`datasetKey\` の記録数が計14,637,152一致 | raw→155/224派生解析全体の再実行可能データ・処理条件 |
| 1,596港の元GBIF植物DOI | **未閉鎖：** 元の\`speciesKey\` facetの各発行者\`datasetKey\`は収録されていない。新索引の1,482件を元の帰属として申請しない | 元の港別GBIF原個票・元のWKT照会キャッシュ・元のダウンロード記録。無ければ新取得データを**新日付で独立に登録・引用**し、該当結果を再検証 |
| 225×51蝶初索引・古い83,573歴史年表 | 元の集計・構造・初年並びを監査済、GBIF\`datasetKey\`なし。新索引は旧データと明確に分離 | 元原個票／provider由来、又は分離した別サンプルとして限定使用 |
| 論文データ/コードの査読アクセス | **内部用**証拠ZIPあり、2026 v0.5源は検査済。古いv0.2の匿名ZIPとは別 | 実在する匿名化・著作権対応済の査読アクセスURL／永続アーカイブ先、著者の承認 |
| 著者・権利 | 元WPI2019版、WCVP v13、HOSTS引用は固定済 | 写真・元記録の再配布条件、著者のタイトルページ、COI/資金、共著者承認 |

**推奨停止線：** 当時の\`datasetKey\`が無い元GBIF speciesKey集合から元提供者一覧を一意に復元することは統計的に不可能（多数のdatasetKeyへの割り振りが同じ種×港カウントを生む）。今の元照会条件から新しいProvider集計を作っても過去の同一性は証明できない。今後は**旧ソース原記録の探索／別日付の代替収録を正直に分ける**。主論文の完成はこの補助解析だけを待たず、正確な出典・データアクセス・著者承認で進める。

## 4. 実際の査読用Data and Code Availability草案（著者レビュー用、公開済を偽らない）

The larval-host source used for potential resource reconstruction is the Natural History Museum HOSTS catalogue (doi:10.5519/havt50xw), with plant ranges from the World Checklist of Vascular Plants v13 (doi:10.34885/0yex-xv26), both using version-pinned transformations. The independent recorded-butterfly analyses of 155 and 224 modelled species draw on the **same original GBIF occurrence download** (doi:10.15468/dl.umqmsf); its 1,716 contributing source datasets reconcile with the archived 14,637,152 download records. Port locations are based on the World Port Index 2019, Publication 150, 27th edition, with geographic source identity and transformations documented. Additional original botanical-search and island-index sources were obtained via the GBIF occurrence-search API; their time-specific index-query receipts and derived aggregates are distinguishable from later re-queries. The original cached port and island facet sources do **not** retain individual occurrence identifiers or all contributing dataset keys and must not be treated as an independently DOI-registered GBIF download. Source-coded numerical results, query constraints, and independently rerunnable portions are retained in a version-pinned analysis record; original record-level data and third-party media are subject to source-specific access and reuse terms. **A reviewer-accessible anonymized repository link and complete final reproducibility inventory have not yet been approved or deposited.** This statement must be replaced with the actual approved anonymous reviewer access route and durable archive identification before journal submission.

**この英語は審査済の投稿文ではない。** 実際の査読者用URL、元個票へのアクセス、権利、著者承認が揃った時点で本文のData Availabilityに反映する。匿名化された公的ZIPや実在しないDOIを記した資料として提出しない。
