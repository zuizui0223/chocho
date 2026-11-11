# 統合原稿とリポジトリの成果の対応（2026-11-11）

GitHubの公開・終了済みを含む全157 PRのタイトル・説明と182リモートブランチを照合した。PRに対応しない現存ブランチは末尾に最新コミットを示す。これは成果・出典の棚卸しで、全解析の独立再実行ではない。

| 成果群 | 統合原稿での扱い |
|---|---|
| #1–38　原稿・公開準備／全球資源 | 主結果は§3.1。整理、引用、匿名化は投稿実務。同じ再集計を反復証拠に数えない。 |
| #39–45　実験案・独立個体データ | 京都の飼育設計、既存幼虫データは野外島嶼の繁殖検証と別。 |
| #46–46　撹乱・気候適性 | 条件付き資源推定。島の環境改変の因果推定から除外。 |
| #47–68　氷期の気候・古地理・化石 | 21 ka の解析および個別の地域解析。現代の港・蝶の時系列の分母に含めない。 |
| #69–70　資源の冗長性・均質化 | 補足。固定次数帰無で資源増加に特異な多重食草による緩衝は支持されない。 |
| #71–95　幼虫・寄主利用・侵入史・気候 | 補足／別研究。§3.1の資源指標と野外の幼虫利用、宿主転換、移入、定着を区別する。 |
| #96–103　予測・英国長期モニタリング | §3.4の155種解析と英国調査は別の目的変数。英国の気象バックテストから全球将来予測へ外挿しない。 |
| #104–116　港と食草の記録照合 | §3.2–3.5。記録時期の先行と地理条件付き帰無を併記。投稿先の旧提案はアーカイブ。 |
| #117–124　港の植物相・歴史・形質 | §3.1–3.3／補足。資源の収斂、VOC航海、植民地史、下水関連、生活形の独立性はそれぞれ別の検定。 |
| #125–125　島外予測・実航路・個別事例 | 本文Box 1と§4.4に統合。距離対応7組、同種内0組。DEMは補足感度。 |
| #126–140　海の地形・AIS・風と植物相 | §3.6・図5。51島と両端未知の990組を別々に評価。下水の継続研究は別リポ。 |
| #141–150　航路照合・過去の定期船 | 補足／§4.4。端点の照合は移動証拠と分ける。歴史的指数の符号反転を追記。 |
| #151–158　統合原稿・形質・時間順序・沿岸対照 | §3.2、3.5–3.6・図5。#151は元原稿。#158の二つの期間を本文に明記。 |

## 生物学的な主張への影響

- 資源機会の増加を、植物の追加によるチョウ種×地域の変化として要旨に定義した。206/239種で増加したが、複数食草による例外的な緩衝作用は次数固定帰無で裏付けられない（#69）。
- #75–86は食草利用報告、輸送経路、独立に記録された定着の区別に有用。地域内で導入食草と定着種が共存していても、港で食草と幼虫が一緒に運ばれたことにはならない。
- #117–124の植民地史、栄養塩・下水、植物生活形は代替説明。下水関連は在来・外来の両方の植物数と共変し、島での一年生草本の効果も植物ネットワーク帰無で p=0.495。本文の港数を、特定の歴史的経路や下水効果に読み替えない。
- #125の25島外検証と#149の過去の定期船指数には重みづけ・地域依存性がある。§4.4に加筆。#158の17島内対照は前期・後期とも弱く、§3.2に両期間を記載。
- #125の2024年事例は元51島の時系列を延長した値ではない。英国・アイルランドを除く別の港島抽出。個別事例を本文に入れ、分母と距離対応の対照を併記した。

## 全PR索引

各行は同じ成果の再検証や公開作業を含む。PR件数は独立した研究数を表さない。

- [#1 Prune chocho to the GEB paper dependency surface](https://github.com/zuizui0223/chocho/pull/1) — 原稿・公開準備／全球資源
- [#2 Prune chocho to the ecology-paper surface](https://github.com/zuizui0223/chocho/pull/2) — 原稿・公開準備／全球資源
- [#4 Harden GEB inference, anonymity, and release CI](https://github.com/zuizui0223/chocho/pull/4) — 原稿・公開準備／全球資源
- [#5 Restore frozen input and durable figure provenance](https://github.com/zuizui0223/chocho/pull/5) — 原稿・公開準備／全球資源
- [#6 Lead GEB abstract with ecological magnitude and complete reproducibility gate](https://github.com/zuizui0223/chocho/pull/6) — 原稿・公開準備／全球資源
- [#7 Build privacy-scrubbed blinded GEB review DOCX](https://github.com/zuizui0223/chocho/pull/7) — 原稿・公開準備／全球資源
- [#8 Fix v1 release identifiers and add administrative preflight](https://github.com/zuizui0223/chocho/pull/8) — 原稿・公開準備／全球資源
- [#9 Build de-identified anonymous review code bundle](https://github.com/zuizui0223/chocho/pull/9) — 原稿・公開準備／全球資源
- [#10 Build de-identified anonymous review code bundle](https://github.com/zuizui0223/chocho/pull/10) — 原稿・公開準備／全球資源
- [#11 Synchronize release readiness after review-package completion](https://github.com/zuizui0223/chocho/pull/11) — 原稿・公開準備／全球資源
- [#12 Document DOI-first final release sequence](https://github.com/zuizui0223/chocho/pull/12) — 原稿・公開準備／全球資源
- [#13 Lead butterfly manuscript with anthropogenic resource-portfolio novelty](https://github.com/zuizui0223/chocho/pull/13) — 原稿・公開準備／全球資源
- [#14 Expose dataset scale and ecological magnitude directly in figures](https://github.com/zuizui0223/chocho/pull/14) — 原稿・公開準備／全球資源
- [#15 Align Results and figure order with novelty-first story](https://github.com/zuizui0223/chocho/pull/15) — 原稿・公開準備／全球資源
- [#16 Summarize nested analysis scale in Methods overview](https://github.com/zuizui0223/chocho/pull/16) — 原稿・公開準備／全球資源
- [#17 Keep manuscript independent while preserving hypothesis provenance](https://github.com/zuizui0223/chocho/pull/17) — 原稿・公開準備／全球資源
- [#18 Clean repository layout and remove legacy package namespace](https://github.com/zuizui0223/chocho/pull/18) — 原稿・公開準備／全球資源
- [#19 Finish repository navigation and paper-facing layout](https://github.com/zuizui0223/chocho/pull/19) — 原稿・公開準備／全球資源
- [#20 Archive superseded climate contracts and complete repository indexes](https://github.com/zuizui0223/chocho/pull/20) — 原稿・公開準備／全球資源
- [#21 Archive superseded climate contracts outside active analysis directories](https://github.com/zuizui0223/chocho/pull/21) — 原稿・公開準備／全球資源
- [#22 Trim active evidence surfaces and archive audit history](https://github.com/zuizui0223/chocho/pull/22) — 原稿・公開準備／全球資源
- [#23 GEB v0.2: resource-side biogeography under host redistribution](https://github.com/zuizui0223/chocho/pull/23) — 原稿・公開準備／全球資源
- [#24 Separate climate execution provenance from active evidence](https://github.com/zuizui0223/chocho/pull/24) — 原稿・公開準備／全球資源
- [#25 Group scripts by paper, resource, climate, and diagnostics](https://github.com/zuizui0223/chocho/pull/25) — 原稿・公開準備／全球資源
- [#26 Reviewer-response reserve: realized-realm generality audit](https://github.com/zuizui0223/chocho/pull/26) — 原稿・公開準備／全球資源
- [#27 Add crop-host sensitivity before GEB submission](https://github.com/zuizui0223/chocho/pull/27) — 原稿・公開準備／全球資源
- [#28 Align final GEB title with specialization-spectrum framing](https://github.com/zuizui0223/chocho/pull/28) — 原稿・公開準備／全球資源
- [#29 Integrate v0.2 manuscript with crop-host sensitivity](https://github.com/zuizui0223/chocho/pull/29) — 原稿・公開準備／全球資源
- [#30 Fix final v0.2 title in release preflight](https://github.com/zuizui0223/chocho/pull/30) — 原稿・公開準備／全球資源
- [#31 Deepen crop-independent ruderal interpretation](https://github.com/zuizui0223/chocho/pull/31) — 原稿・公開準備／全球資源
- [#32 Final GEB submission polish](https://github.com/zuizui0223/chocho/pull/32) — 原稿・公開準備／全球資源
- [#33 Separate initial GEB submission from public archive release](https://github.com/zuizui0223/chocho/pull/33) — 原稿・公開準備／全球資源
- [#34 Complete final manuscript citation audit](https://github.com/zuizui0223/chocho/pull/34) — 原稿・公開準備／全球資源
- [#35 Add manual GBIF DOI mint workflow](https://github.com/zuizui0223/chocho/pull/35) — 原稿・公開準備／全球資源
- [#36 Sharpen invasion ecology and conservation framing](https://github.com/zuizui0223/chocho/pull/36) — 原稿・公開準備／全球資源
- [#37 Finalize abstract logic and conservation scale](https://github.com/zuizui0223/chocho/pull/37) — 原稿・公開準備／全球資源
- [#38 Reframe paper around resource homogenization and shared-resource exposure](https://github.com/zuizui0223/chocho/pull/38) — 原稿・公開準備／全球資源
- [#39 Kyoto temporal host-window hypothesis: theory plus source-grounded no-go diagnostics](https://github.com/zuizui0223/chocho/pull/39) — 実験案・独立個体データ
- [#40 Experiment design: test invader food rescue before species-identity mechanisms](https://github.com/zuizui0223/chocho/pull/40) — 実験案・独立個体データ
- [#41 Kyoto four-arm invader-food experiment: 3-cage greenhouse incomplete-block solution](https://github.com/zuizui0223/chocho/pull/41) — 実験案・独立個体データ
- [#42 Use all 30 Kyoto greenhouse cages: contrast-optimal near-balanced food-rescue design](https://github.com/zuizui0223/chocho/pull/42) — 実験案・独立個体データ
- [#43 Butterfly larval-resource geography: full alternative GEB paper and source-robustness audit](https://github.com/zuizui0223/chocho/pull/43) — 実験案・独立個体データ
- [#44 Original monarch adult output: plot-level temperature × host species audit (prior-art calibration)](https://github.com/zuizui0223/chocho/pull/44) — 実験案・独立個体データ
- [#45 Exploratory: four-butterfly original-cohort host-quality audit (3,458 larvae; no recruitment claim)](https://github.com/zuizui0223/chocho/pull/45) — 実験案・独立個体データ
- [#46 Exploratory macroecology: disturbance affinity and climate-filtered butterfly resource opportunity](https://github.com/zuizui0223/chocho/pull/46) — 撹乱・気候適性
- [#47 Establish native-WCVP-screened host occurrences for 392 insect LGM comparisons](https://github.com/zuizui0223/chocho/pull/47) — 氷期の気候・古地理・化石
- [#48 Hindcast 416 WCVP-native host niches and 513 insects' potential climate-resource turnover at 21 ka](https://github.com/zuizui0223/chocho/pull/48) — 氷期の気候・古地理・化石
- [#49 Use verified PMIP4 ICE-6G_C 21 ka land and ice barriers in native host climate analysis](https://github.com/zuizui0223/chocho/pull/49) — 氷期の気候・古地理・化石
- [#50 Qualify PMIP4–CHELSA 0.5° paleoland climate coverage and audited native-host potential topology](https://github.com/zuizui0223/chocho/pull/50) — 氷期の気候・古地理・化石
- [#51 Complete 229-taxon country macro analysis and document unresolved association](https://github.com/zuizui0223/chocho/pull/51) — 氷期の気候・古地理・化石
- [#52 Measure 0 BP–21 ka potential native-host habitat fragmentation across 513 insects](https://github.com/zuizui0223/chocho/pull/52) — 氷期の気候・古地理・化石
- [#53 Audit 513 insect LGM potential habitat graphs against eight-neighbor and 349 distinct host combinations](https://github.com/zuizui0223/chocho/pull/53) — 氷期の気候・古地理・化石
- [#54 Audit LGM host climate-graph fragmentation against novel-climate exclusion](https://github.com/zuizui0223/chocho/pull/54) — 氷期の気候・古地理・化石
- [#55 Distinguish LGM host-climate patch splitting from loss of contact at fixed modern insect sampling pairs](https://github.com/zuizui0223/chocho/pull/55) — 氷期の気候・古地理・化石
- [#56 Map original modern insect sampling anchors to potential host-climate patches at 0 BP and 21 ka](https://github.com/zuizui0223/chocho/pull/56) — 氷期の気候・古地理・化石
- [#57 Separate native-host LGM climate turnover from ice-free terrain gain and loss across 513 insects](https://github.com/zuizui0223/chocho/pull/57) — 氷期の気候・古地理・化石
- [#58 Audit 292-species real-other-plant LGM patch null for missing shared edges and statistical identifiability](https://github.com/zuizui0223/chocho/pull/58) — 氷期の気候・古地理・化石
- [#59 Evaluate known-host 21ka patch splitting against preselected climate-matched real plants](https://github.com/zuizui0223/chocho/pull/59) — 氷期の気候・古地理・化石
- [#60 Qualify independent Neotoma fossil evidence sources for LGM native host occurrence without inferring ancient occupancy](https://github.com/zuizui0223/chocho/pull/60) — 氷期の気候・古地理・化石
- [#61 Verify Neotoma plant fossil taxon filters with impossible-name control before 21 ka occupancy inference](https://github.com/zuizui0223/chocho/pull/61) — 氷期の気候・古地理・化石
- [#62 Inspect source-authenticated Neotoma taxon-ID fossil sample ages before host LGM occupancy claims](https://github.com/zuizui0223/chocho/pull/62) — 氷期の気候・古地理・化石
- [#63 Identify 18–24ka Alnus pollen source candidates with exact Neotoma taxon ID and sample age checks](https://github.com/zuizui0223/chocho/pull/63) — 氷期の気候・古地理・化石
- [#64 Audit six Neotoma Alnus pollen datasets for mismatched LGM ages and radiocarbon calibration scales](https://github.com/zuizui0223/chocho/pull/64) — 氷期の気候・古地理・化石
- [#65 Audit real-other-plant LGM host specificity null across 292 insects and 173 host combinations](https://github.com/zuizui0223/chocho/pull/65) — 氷期の気候・古地理・化石
- [#66 Audit full Neotoma downloads for exact six Alnus pollen sample chronology provenance](https://github.com/zuizui0223/chocho/pull/66) — 氷期の気候・古地理・化石
- [#67 Audit raw six Neotoma pollen datasets for Alnus-positive samples with explicit 18–24ka age models](https://github.com/zuizui0223/chocho/pull/67) — 氷期の気候・古地理・化石
- [#68 Audit 21 LGM-labelled Alnus pollen samples against original Neotoma chronology models and age-depth controls](https://github.com/zuizui0223/chocho/pull/68) — 氷期の気候・古地理・化石
- [#69 Falsify proposed butterfly gain–buffering mechanism with fixed-degree plant-range null](https://github.com/zuizui0223/chocho/pull/69) — 資源の冗長性・均質化
- [#70 Calibrate GEB ecological homogenization claims and verify submission-only blockers](https://github.com/zuizui0223/chocho/pull/70) — 資源の冗長性・均質化
- [#71 Audit monarch larval growth-to-pupation stage conversion and unresolved source label swap](https://github.com/zuizui0223/chocho/pull/71) — 幼虫・寄主利用・侵入史・気候
- [#72 Test whether indexed butterflies enter climate-novel versus analogous introduced-host opportunities](https://github.com/zuizui0223/chocho/pull/72) — 幼虫・寄主利用・侵入史・気候
- [#73 Verify butterfly and introduced larval-host point proximity beyond country-level niche opportunity](https://github.com/zuizui0223/chocho/pull/73) — 幼虫・寄主利用・侵入史・気候
- [#74 Reconstruct and validate all 73 climate-novel vs analogue butterfly positive countries for point-level niche testing](https://github.com/zuizui0223/chocho/pull/74) — 幼虫・寄主利用・侵入史・気候
- [#75 Curate invasion histories and recover source HOSTS Location for host adoption vs resource redistribution](https://github.com/zuizui0223/chocho/pull/75) — 幼虫・寄主利用・侵入史・気候
- [#76 Map HOSTS Location host-use reports into WCVP native versus introduced plant regions for 239 butterflies](https://github.com/zuizui0223/chocho/pull/76) — 幼虫・寄主利用・侵入史・気候
- [#77 Test whether butterfly diet breadth predicts introduced host-use records beyond generic HOSTS reporting bias](https://github.com/zuizui0223/chocho/pull/77) — 幼虫・寄主利用・侵入史・気候
- [#78 Test native vs introduced HOSTS host-pair documentation overlap against fixed-source nulls](https://github.com/zuizui0223/chocho/pull/78) — 幼虫・寄主利用・侵入史・気候
- [#79 Audit independent invasive butterfly larval host evidence versus HOSTS record geography](https://github.com/zuizui0223/chocho/pull/79) — 幼虫・寄主利用・侵入史・気候
- [#80 Test whether introduced larval-host plants become shared butterfly resource hubs beyond HOSTS reporting bias](https://github.com/zuizui0223/chocho/pull/80) — 幼虫・寄主利用・侵入史・気候
- [#81 Bridge 239 documented butterfly host portfolios with independent field caterpillar detection](https://github.com/zuizui0223/chocho/pull/81) — 幼虫・寄主利用・侵入史・気候
- [#82 Separate seed/live-plant butterfly arrival vectors from host-resource redistribution; correct all-observed 115 cohort](https://github.com/zuizui0223/chocho/pull/82) — 幼虫・寄主利用・侵入史・気候
- [#83 Connect chocho butterfly host-resource opportunities to independent Lepidoptera establishment and insect interception records](https://github.com/zuizui0223/chocho/pull/83) — 幼虫・寄主利用・侵入史・気候
- [#84 Cross independently documented alien-insect establishment records to chocho 239 host-resource opportunities](https://github.com/zuizui0223/chocho/pull/84) — 幼虫・寄主利用・侵入史・気候
- [#85 Map 35 independently established alien butterflies to recipient-specific WCVP hosts beyond Hawaii](https://github.com/zuizui0223/chocho/pull/85) — 幼虫・寄主利用・侵入史・気候
- [#86 Measure local wild-host bottlenecks among 24 alien butterflies and verify plant-contaminant vs native-host pathway cases](https://github.com/zuizui0223/chocho/pull/86) — 幼虫・寄主利用・侵入史・気候
- [#87 Audit species-paired butterfly climate-recording contrasts and replication attrition across 73 cases](https://github.com/zuizui0223/chocho/pull/87) — 幼虫・寄主利用・侵入史・気候
- [#88 Test geographic climate opportunity against the 73 butterfly climate-recording contrasts](https://github.com/zuizui0223/chocho/pull/88) — 幼虫・寄主利用・侵入史・気候
- [#89 Test whether cultivated native milkweeds alter seasonal monarch larval host-record composition](https://github.com/zuizui0223/chocho/pull/89) — 幼虫・寄主利用・侵入史・気候
- [#90 Project CMIP6 warming effects on monarch larval thermal opportunity with explicit milkweed leaf-persistence scenarios](https://github.com/zuizui0223/chocho/pull/90) — 幼虫・寄主利用・侵入史・気候
- [#91 Measure real milkweed live-leaf phenology and winter sampling gates before coupling monarch CMIP6 projections](https://github.com/zuizui0223/chocho/pull/91) — 幼虫・寄主利用・侵入史・気候
- [#92 Test CMIP6 winter thermal opportunity against directly monitored native milkweed leaf months at the same sites](https://github.com/zuizui0223/chocho/pull/92) — 幼虫・寄主利用・侵入史・気候
- [#93 Test whether annual weather predicts directly observed milkweed leaves at held-out USA-NPN sites](https://github.com/zuizui0223/chocho/pull/93) — 幼虫・寄主利用・侵入史・気候
- [#94 Falsify milkweed-leaf climate transport using temperature-only, rainfall-only and future-year holdouts](https://github.com/zuizui0223/chocho/pull/94) — 幼虫・寄主利用・侵入史・気候
- [#95 Project conditional milkweed leaf-status response to real CMIP6 monthly warming anomalies with out-of-domain gates](https://github.com/zuizui0223/chocho/pull/95) — 幼虫・寄主利用・侵入史・気候
- [#96 Pilot historical European butterfly range-tracking hindcast with source-frozen host and effort identification gates](https://github.com/zuizui0223/chocho/pull/96) — 予測・英国長期モニタリング
- [#97 Recover original UK butterfly survey visits, 1976–2021 matched transects and source-verified historical detection outcomes](https://github.com/zuizui0223/chocho/pull/97) — 予測・英国長期モニタリング
- [#98 Backtest UK butterfly detections using pre-1970 BSBI larval host geography versus actual historical climate](https://github.com/zuizui0223/chocho/pull/98) — 予測・英国長期モニタリング
- [#99 Test actual pre-1987 BSBI larval-host geography against 2005–2021 UK butterfly monitoring outcomes](https://github.com/zuizui0223/chocho/pull/99) — 予測・英国長期モニタリング
- [#100 Build global H3 r4 × WGSRPD3 fractional crosswalk for butterfly–host–landscape year panels](https://github.com/zuizui0223/chocho/pull/100) — 予測・英国長期モニタリング
- [#101 Validate host-resource gains in temporal and spatial butterfly forecasts](https://github.com/zuizui0223/chocho/pull/101) — 予測・英国長期モニタリング
- [#102 Island barrier tests: exact source-distance validation and original230 host recovery](https://github.com/zuizui0223/chocho/pull/102) — 予測・英国長期モニタリング
- [#103 Use published pre-2000 British plant atlas to test historical host-mediated butterfly climate tracking](https://github.com/zuizui0223/chocho/pull/103) — 予測・英国長期モニタリング
- [#104 Document landmass validation and exploratory island port associations](https://github.com/zuizui0223/chocho/pull/104) — 港と食草の記録照合
- [#105 Prepare historically dated host/nonhost comparison targets for island-port relay hypothesis](https://github.com/zuizui0223/chocho/pull/105) — 港と食草の記録照合
- [#106 Resolve pre-butterfly larval-host plant identity and dated island-port observations](https://github.com/zuizui0223/chocho/pull/106) — 港と食草の記録照合
- [#107 Reconstruct pre-2011 plant recording composition across all 17 near-island nodes and test butterfly similarity under reporting-quality gates](https://github.com/zuizui0223/chocho/pull/107) — 港と食草の記録照合
- [#108 Test pre-2011 source-verified larval-host records against 2011–2017 island butterfly re-recording](https://github.com/zuizui0223/chocho/pull/108) — 港と食草の記録照合
- [#109 Test original230 host-indexing excess under geographic island fixed margins](https://github.com/zuizui0223/chocho/pull/109) — 港と食草の記録照合
- [#110 Geography-matched original230 HOSTS null removes apparent historical host co-indexing excess](https://github.com/zuizui0223/chocho/pull/110) — 港と食草の記録照合
- [#111 Test whether GEB's 54.9% potential host-resource expansion is expressed through island ports](https://github.com/zuizui0223/chocho/pull/111) — 港と食草の記録照合
- [#112 Integrate 155-species forecast improvements with positive island and port correlations in GEB synthesis](https://github.com/zuizui0223/chocho/pull/112) — 港と食草の記録照合
- [#113 Lead with observed GEB–island–port correlations while keeping mechanisms unclaimed](https://github.com/zuizui0223/chocho/pull/113) — 港と食草の記録照合
- [#114 Set GEB-first integrated submission route for global host gains, forecast validation and island associations](https://github.com/zuizui0223/chocho/pull/114) — 港と食草の記録照合
- [#115 Island- and port-first JBI submission gate: lock 40-island coastal control and 36-island botanical joins](https://github.com/zuizui0223/chocho/pull/115) — 港と食草の記録照合
- [#116 Butterfly island expansion: test independent migration and alien-host resources using source-verified record histories](https://github.com/zuizui0223/chocho/pull/116) — 港と食草の記録照合
- [#117 Island port networks and butterfly host geography: genuine 17th-century VOC voyages and botanical cargo](https://github.com/zuizui0223/chocho/pull/117) — 港の植物相・歴史・形質
- [#118 Global plant-to-butterfly resource convergence and historical trade legacy: source-locked pilot](https://github.com/zuizui0223/chocho/pull/118) — 港の植物相・歴史・形質
- [#119 Global colonial history × plant/butterfly potential-resource beta: spatial crosswalk and negative test](https://github.com/zuizui0223/chocho/pull/119) — 港の植物相・歴史・形質
- [#120 Global host-resource convergence versus real butterfly communities: 14.6M GBIF records and colonial-history controls](https://github.com/zuizui0223/chocho/pull/120) — 港の植物相・歴史・形質
- [#121 Island butterfly resources and three gates: host, climate, historical access (Q1–Q4)](https://github.com/zuizui0223/chocho/pull/121) — 港の植物相・歴史・形質
- [#122 Draft: nutrient affinity, port-indexed butterfly hosts, and matched all-flora controls](https://github.com/zuizui0223/chocho/pull/122) — 港の植物相・歴史・形質
- [#123 SOURCE ARCHIVE (one-paper PR #157): Draft JBI: port-linked introduced-flora assembly and island butterfly resource expansion](https://github.com/zuizui0223/chocho/pull/123) — 港の植物相・歴史・形質
- [#124 Draft test: no positive introduced annual-herb filtering with island port richness](https://github.com/zuizui0223/chocho/pull/124) — 港の植物相・歴史・形質
- [#125 Validate whole-island geography and port forecasts with uncertainty](https://github.com/zuizui0223/chocho/pull/125) — 島外予測・実航路・個別事例
- [#126 Draft Fig.4 gate: 51-island circuit resistance fails distance LOIO; verify commercial AIS access](https://github.com/zuizui0223/chocho/pull/126) — 海の地形・AIS・風と植物相
- [#127 Run observed marine circuits for 51 islands with official and native-AIS validation](https://github.com/zuizui0223/chocho/pull/127) — 海の地形・AIS・風と植物相
- [#128 [island51] Joint host and marine circuit prediction with independent audit](https://github.com/zuizui0223/chocho/pull/128) — 海の地形・AIS・風と植物相
- [#129 Map and account for fixed island host–circuit prediction gains](https://github.com/zuizui0223/chocho/pull/129) — 海の地形・AIS・風と植物相
- [#130 Test whether candidate-host identity predicts butterflies beyond host-list richness](https://github.com/zuizui0223/chocho/pull/130) — 海の地形・AIS・風と植物相
- [#131 Test 51-island host/circuit prediction against pre-2011 botanical recording effort](https://github.com/zuizui0223/chocho/pull/131) — 海の地形・AIS・風と植物相
- [#132 Evaluate island host and marine-circuit predictors with both island-pair endpoints unseen](https://github.com/zuizui0223/chocho/pull/132) — 海の地形・AIS・風と植物相
- [#133 Separate AIS vessel-density and wind circuit increments on two-unseen-island predictions](https://github.com/zuizui0223/chocho/pull/133) — 海の地形・AIS・風と植物相
- [#134 Do AIS and wind prediction gains persist after latitude and hemisphere controls?](https://github.com/zuizui0223/chocho/pull/134) — 海の地形・AIS・風と植物相
- [#135 Does latitude-controlled AIS prediction depend on two port-heavy island groups?](https://github.com/zuizui0223/chocho/pull/135) — 海の地形・AIS・風と植物相
- [#136 Falsify AIS-island alignment using 199 same-hemisphere latitude-band graph relabelings](https://github.com/zuizui0223/chocho/pull/136) — 海の地形・AIS・風と植物相
- [#137 Test AIS network specificity against only geographically local island-label swaps (≤750 km)](https://github.com/zuizui0223/chocho/pull/137) — 海の地形・AIS・風と植物相
- [#138 Reacquire actual 2001–2010 GBIF plant species identities at all original 194 ports; test on 51 unseen islands](https://github.com/zuizui0223/chocho/pull/138) — 海の地形・AIS・風と植物相
- [#139 Document island/port/AIS/wind marine research scope and move wastewater ownership to gesui](https://github.com/zuizui0223/chocho/pull/139) — 海の地形・AIS・風と植物相
- [#140 Does maritime connectivity predict both 2001–2010 port flora and early/later butterfly communities?](https://github.com/zuizui0223/chocho/pull/140) — 海の地形・AIS・風と植物相
- [#141 Audit public 2020 AIS vessel–port triplet archive for direct routes among original 194 ports](https://github.com/zuizui0223/chocho/pull/141) — 航路照合・過去の定期船
- [#142 Can even all 27 geographically nearby but unverified AIS port aliases rescue island-route coverage?](https://github.com/zuizui0223/chocho/pull/142) — 航路照合・過去の定期船
- [#143 Audit independently published global ~1,300-port maritime network for original WPI194 compatibility](https://github.com/zuizui0223/chocho/pull/143) — 航路照合・過去の定期船
- [#144 Locate PortWatch's public AIS-derived global Maritime Connections table before WPI194 matching](https://github.com/zuizui0223/chocho/pull/144) — 航路照合・過去の定期船
- [#145 Pilot 12 directed PortWatch model links among four original WPI-verified ports](https://github.com/zuizui0223/chocho/pull/145) — 航路照合・過去の定期船
- [#146 Recover original GLSN2018 weighted network edges and historical UNCTAD PLSCIQ source bytes](https://github.com/zuizui0223/chocho/pull/146) — 航路照合・過去の定期船
- [#147 Independently geoverify 17 original WPI ports shared by GLSN2018, historic UNCTAD and AIS2020](https://github.com/zuizui0223/chocho/pull/147) — 航路照合・過去の定期船
- [#148 Geoverify remaining 13 original WPI ports shared by GLSN2018 and UNCTAD but not indexed in AIS2020](https://github.com/zuizui0223/chocho/pull/148) — 航路照合・過去の定期船
- [#149 Does historical 2006–2010 port liner connectedness predict later butterfly communities on 17 original islands?](https://github.com/zuizui0223/chocho/pull/149) — 航路照合・過去の定期船
- [#150 Draft Fig.4: map source-locked marine current corridors without claiming island relays](https://github.com/zuizui0223/chocho/pull/150) — 航路照合・過去の定期船
- [#151 SOURCE ARCHIVE (one-paper PR #157): Draft JBI v0.2: island host-resource reservoirs, whole-flora ports and marine access](https://github.com/zuizui0223/chocho/pull/151) — 統合原稿・形質・時間順序・沿岸対照
- [#152 Draft NC first-shot: freeze temporal host × marine access test and verify historical-source bottleneck](https://github.com/zuizui0223/chocho/pull/152) — 統合原稿・形質・時間順序・沿岸対照
- [#153 Draft NC: reconstruct 225×51 source visibility without false non-event labels](https://github.com/zuizui0223/chocho/pull/153) — 統合原稿・形質・時間順序・沿岸対照
- [#154 Draft NC: locate original Library Parquet and recover 4 verified later first-indexed island cases](https://github.com/zuizui0223/chocho/pull/154) — 統合原稿・形質・時間順序・沿岸対照
- [#155 Draft NC: independent public 225×17 historical butterfly risk versus ocean circuit access](https://github.com/zuizui0223/chocho/pull/155) — 統合原稿・形質・時間順序・沿岸対照
- [#156 Butterflies first: wing-size moderation of marine resistance across island GBIF records](https://github.com/zuizui0223/chocho/pull/156) — 統合原稿・形質・時間順序・沿岸対照
- [#157 One GEB butterfly paper: island resource amplification, port floras and marine geography](https://github.com/zuizui0223/chocho/pull/157) — 統合原稿・形質・時間順序・沿岸対照
- [#158 Test matched island port vs nonport coast: historical candidate hosts and later butterfly indices](https://github.com/zuizui0223/chocho/pull/158) — 統合原稿・形質・時間順序・沿岸対照

## PRのhead名に対応しない現存ブランチ

- [analysis/butterfly-macro-hostgain-abundance-20261026](https://github.com/zuizui0223/chocho/commit/fc702c0a65714c87e2a262869331d0cf22f0dd01) — Fix missing-as-zero versus recorded-positive adult trend sensitivity before full macro population model outcomes
- [analysis/disturbance-filter-macro-20261024](https://github.com/zuizui0223/chocho/commit/b6513108eee2e194d886be4f8446cac5fdc47f89) — Archive posthoc adequacy-group imbalance explaining trait association sensitivity
- [analysis/gain-fragility-20261024](https://github.com/zuizui0223/chocho/commit/8d24a511a610655a2a2e111e448c6b4e569f64dc) — Link visually inspected SVG/PNG and bounded manuscript addendum from source-locked report
- [analysis/geb-integrated-submission-track-20261031](https://github.com/zuizui0223/chocho/commit/169c1133ef4fa661f68ad4c911f56f7357e4af48) — Tighten integrated GEB structured abstract to leave word-count margin and preserve cohorts
- [analysis/globi-larval-host-source-pilot-20261024](https://github.com/zuizui0223/chocho/commit/adfb8beb28325fb7b7af99aec80e786d84e502a9) — Close realized-use macro evidence gate with 19/20 pagination audit and proven upper bound
- [analysis/island51-observed-AIS-conductance-20261106](https://github.com/zuizui0223/chocho/commit/0f739cf6f502112699302aae46635b039361ba19) — Verify source GSHHG 1-degree ocean topology and all 51 authentic WPI port-centre attachments before fitting AIS/wind resistance
- [analysis/larval-evidence-20261011](https://github.com/zuizui0223/chocho/commit/5d70550cf315d9fd15083653e188bc0df1090150) — Separate actual larval evidence from host omissions and wrong-host geographic agreement
- [analysis/larval-host-seasonal-realization-20261026](https://github.com/zuizui0223/chocho/commit/7f75f51b96e46f755e78f5a83d5d4d5b4f68d86a) — Close butterfly-specific source identification gate, maintain moth/butterfly distinction, archive final QA and figure
- [analysis/monarch-host-stage-realization-20261026](https://github.com/zuizui0223/chocho/commit/185d9ceb6c7855d4fe2927cb1949cc028173a467) — Archive verified public Dryad v10 original monarch juvenile and host plant dataset source quality only
- [analysis/remote-sensing-20261029](https://github.com/zuizui0223/chocho/commit/30c43c884b521260572e12e2f117a9ea2be5e811) — Request multi-year coordinate source for landscape response analysis
- [analysis/resource-homogenization-v01](https://github.com/zuizui0223/chocho/commit/835d1e5d0f652d201a2ddb14f92d4144354151e8) — Freeze exact matching deterministic Level1 null receipts from two CI runs and four passing tests
- [analysis/resource-redundancy-20261011](https://github.com/zuizui0223/chocho/commit/fc89336a7d6f0a3bdecbcbc4777c5aeb68295162) — Bound source-g phenology label rarity with all-year flowering and clone-out controls
- [analysis/reviewer-defense-v01](https://github.com/zuizui0223/chocho/commit/0cca23e133cdc6260c8363cb046103ae618794af) — Index reviewer-defense evidence for manuscript v0.2
- [analysis/sbus-environment-shock-20261026](https://github.com/zuizui0223/chocho/commit/b917ec5249b102265ce09f0e92b5887a87cab7b0) — Finish Illinois weather-vs-niche falsification with fully source-verified figure and five green workflows
- [analysis/sbus-host-sharing-abundance-20261025](https://github.com/zuizui0223/chocho/commit/4dbd711288862f4d9a04add4ea59d64a05bfed25) — Prominently demote unqualified adult competition coefficient after 199 source-lock host-graph nulls
- [analysis/ukbms-niche-competition-20261025](https://github.com/zuizui0223/chocho/commit/58794816192ff57895ce0733f0f89eb79478ce17) — Close documentation-only UKBMS acquisition gate and specify independent US standardized-count alternative before outcomes
- [ecology/lgm-native-plant-niche-hindcast-v01](https://github.com/zuizui0223/chocho/commit/8a000a53a379c0bdf4254f7aae5503ffbc25ff6f) — Freeze WCVP native-screened host SDM eligibility and model source receipt
- [ecology/lgm-neotoma-Alnus-model-id-chronology-qc-v19](https://github.com/zuizui0223/chocho/commit/059962344f7ae3178dc5bd807161fc4ef7e08a4f) — Audit complete Neotoma original chronology IDs for 21 fixed 18-24ka Alnus source pollen sample metadata
- [ecology/lgm-neotoma-occurrence-sample-age-crosswalk-v17](https://github.com/zuizui0223/chocho/commit/175b2300245a730c31ca5a49356b4d6fd0bc286e) — Audit exact sample/chronology ID joins of fixed first Alnus fossil occurrence per six Neotoma sites
- [ecology/lgm-neotoma-six-site-chronology-v14](https://github.com/zuizui0223/chocho/commit/ded48dc2b597e27bee4de107e1bb2d95cb6a2393) — Archive source chronologies and references for six 21ka Alnus pollen candidate datasets
- [exploration/geb-generality-audit-v1](https://github.com/zuizui0223/chocho/commit/fc8f6b603dd2c76e6e65a63e3441b4234584c72f) — Freeze leave-one-group-out GEB robustness protocol
- [exploration/realized-host-window-v01](https://github.com/zuizui0223/chocho/commit/7924a02c5c1824102b65acbc594a7d2f8f01041d) — Surface vital 2020 native density-food life history prior art and total-abundance confounder in Kyoto integrated handoff
- [main](https://github.com/zuizui0223/chocho/commit/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660) — Clarify single integrated JBI island-port focus versus GEB-first split, preventing duplicate +54.9% result claims
- [ops/gbif-run-v01](https://github.com/zuizui0223/chocho/commit/7960fdf82ddfeb4aaa8b1e8339621fd1cea3936b) — Add manual GBIF DOI mint workflow
- [repo/archive-supporting-receipts](https://github.com/zuizui0223/chocho/commit/e2a108b2ac4b75d606cd8e0de6d0f363064b987c) — Add resource-development archive index
- [repo/archive-supporting-receipts-v2](https://github.com/zuizui0223/chocho/commit/4d930f0acd8641f44212b1699cd9318912816306) — Archive superseded climate contracts outside active analysis directories (#21)
- [repo/split-active-results-from-execution-history](https://github.com/zuizui0223/chocho/commit/e1611c8e50e52736ce4ce9147915ef8c02c68a76) — Move butterfly_climate_release_preclimate_result_v0.2.1.json into provenance
