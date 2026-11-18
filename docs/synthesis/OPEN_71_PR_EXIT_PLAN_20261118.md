# chocho 残り71 PRの出口計画（2026-11-18）

**現在地:** 当日時点の `open` 71件を完全に排他的に分類。分類は PR # と head SHA に結びつく静的監査であり、PRをマージした証拠ではない。**目標: 常時進行は単一のGEB原稿とその統合テストだけ。最終的に open を0〜2件へ削減。** 旧ブランチは削除しない。

## 実行の順番
1. **重複した旧研究／旧原稿を終了（C）**：すでに main の evidence ledger にある結果・限界へ索引でつなぎ、11件を閉じる。終了とは研究の否定ではない。
2. **mainの結果原本を科学的に固定（B→A）**：source-pinned scripts + CI assertions + 155/239/225/224/41 と各分母＋新しい幼虫証拠を、現在の main から作る **一つのクリーンな統合変更**へ移す。現在mainに保存された811資産は検証済み主張と限界を保持するが、元の未統合スクリプトを全て含んだわけではない。
3. **巨大な積み重ねPRを研究史へ戻す（D）**：22 claim register / 最新週の研究結果 / source hash へ各PRを対応づけ、独立していない重複・null・非代表性を確認した後、一括Git mergeせず閉じる。
4. **一つの論文を校閲・CI・査読可能な版にする（A）**：GEB v0.2 frozen main を維持しつつ、v0.3 integrated manuscript の採用は、単一代表データの科学的妥当性・引用・倫理・共著者確認・図と本文数値の一致に合格した場合だけ。PR #177 の固有蝶独立コホートは主結果へ未校正外挿しない。
5. **最終PR数を再監査**：GitHub の state=all / open 状態、ファイル差分、CI完走と本投稿の承認ステータスを記録して終了する。

## 具体的な科学的閉鎖基準
- A: 最新の統合原稿が対象コホートに忠実で、再現性がある、というより強い完成条件を要求。PR #157だけで論文の唯一の投稿候補とする。PR #177は外的妥当性・幼虫利用の検証であり、新規の別論文として提出しない。
- B: 各sourceの `input identity → exact-code → result receipt → CI → limited scientific claim`をmain新たな統合ブランチで突合する。サンプルの欠如を生物学的なゼロとして埋めない。
- C: 統合原稿と枝PRの結果の出典は残し、旧モデル・旧図・主張をmainへの生Gitマージで復活させない。**ここだけすぐクローズ可能**。
- D: 時系列上の最終マージから一括閉鎖する前に、各研究系列の最後のソース監査結果と反証が main にあるか点検し、必要な欠落は main 追加してから閉じる。現在mainにないコードが必要ならBへ昇格。
- **絶対条件:** 239種の+54.9%、島+118.9%、155種の予測改善、海洋・港の空間関連、固有蝶の移入寄主幼虫利用を保持。ただし船による植物輸送、蝶の新規到着・個体数増大、global endemic平均、港/海流の独立因果は断定しない。

## 71件の分類（完全一覧）

### A — 活かす2つの中心PR（2件）
#157 is the single GEB paper draft; #177 is the independent living-island-endemic external-validity audit. These are **not** two different manuscripts. Consolidate original source-typed results and manuscript into one tested main candidate before either can be merged/closed.

[#157](https://github.com/zuizui0223/chocho/pull/157)、[#177](https://github.com/zuizui0223/chocho/pull/177)

### B — 原データ・解析コードの統合まで一時保留する12件（12件）
Source-bearing PRs for global 239 resource reconstruction, 155-butterfly recording forecast, port/whole flora, island comparisons, marine geometry and metadata, historical source-index controls, primary larval report comparisons. On each: identify minimal reproducible scripts/tests + pinned original dataset IDs + expected result receipts; import selectively into one clean candidate branch based on current main. Close only after reproducible CI and source coverage, not simply after collecting prose.

[#117](https://github.com/zuizui0223/chocho/pull/117)、[#120](https://github.com/zuizui0223/chocho/pull/120)、[#121](https://github.com/zuizui0223/chocho/pull/121)、[#122](https://github.com/zuizui0223/chocho/pull/122)、[#125](https://github.com/zuizui0223/chocho/pull/125)、[#159](https://github.com/zuizui0223/chocho/pull/159)、[#160](https://github.com/zuizui0223/chocho/pull/160)、[#161](https://github.com/zuizui0223/chocho/pull/161)、[#162](https://github.com/zuizui0223/chocho/pull/162)、[#163](https://github.com/zuizui0223/chocho/pull/163)、[#175](https://github.com/zuizui0223/chocho/pull/175)、[#176](https://github.com/zuizui0223/chocho/pull/176)

### C — 主張をmainの証拠台帳へ引き継いだ旧PR／不要な並行設計：今閉じる11件（11件）
Archive source links and previous results, but do not merge code and alternative manuscripts. Caution: negative/null controls must remain in the main scientific evidence register. #118/#119 subsumed by #120; #116 by later island/port/historical audits; #124/#156/#158 are rebuttal/sensitivity layers not standalone effects; #69 is a fixed-degree null; #70 is older submission gate; #100/#102 map/provenance pilots; #150 outdated figure proposal.

[#69](https://github.com/zuizui0223/chocho/pull/69)、[#70](https://github.com/zuizui0223/chocho/pull/70)、[#100](https://github.com/zuizui0223/chocho/pull/100)、[#102](https://github.com/zuizui0223/chocho/pull/102)、[#116](https://github.com/zuizui0223/chocho/pull/116)、[#118](https://github.com/zuizui0223/chocho/pull/118)、[#119](https://github.com/zuizui0223/chocho/pull/119)、[#124](https://github.com/zuizui0223/chocho/pull/124)、[#150](https://github.com/zuizui0223/chocho/pull/150)、[#156](https://github.com/zuizui0223/chocho/pull/156)、[#158](https://github.com/zuizui0223/chocho/pull/158)

### D — 一連の旧ソース検証PR46件：検証台帳の出典だけ残して順に閉じる（46件）
These are not 46 independent discoveries: a large stacked maritime circuit/port provenance line (#126–#174) and an earlier larval/use/interception line (#75–#86). Save every original positive, changed interpretation, null, input provenance and relevant script/test with correct species/island/time denominator; check main file coverage and workflows. Once direct results are imported/linked and 'no causal port→arrival' caveats retained, close in dependency order WITHOUT deleting branches.

[#174](https://github.com/zuizui0223/chocho/pull/174)、[#173](https://github.com/zuizui0223/chocho/pull/173)、[#172](https://github.com/zuizui0223/chocho/pull/172)、[#171](https://github.com/zuizui0223/chocho/pull/171)、[#170](https://github.com/zuizui0223/chocho/pull/170)、[#169](https://github.com/zuizui0223/chocho/pull/169)、[#168](https://github.com/zuizui0223/chocho/pull/168)、[#167](https://github.com/zuizui0223/chocho/pull/167)、[#166](https://github.com/zuizui0223/chocho/pull/166)、[#165](https://github.com/zuizui0223/chocho/pull/165)、[#164](https://github.com/zuizui0223/chocho/pull/164)、[#149](https://github.com/zuizui0223/chocho/pull/149)、[#148](https://github.com/zuizui0223/chocho/pull/148)、[#147](https://github.com/zuizui0223/chocho/pull/147)、[#146](https://github.com/zuizui0223/chocho/pull/146)、[#145](https://github.com/zuizui0223/chocho/pull/145)、[#144](https://github.com/zuizui0223/chocho/pull/144)、[#143](https://github.com/zuizui0223/chocho/pull/143)、[#142](https://github.com/zuizui0223/chocho/pull/142)、[#141](https://github.com/zuizui0223/chocho/pull/141)、[#140](https://github.com/zuizui0223/chocho/pull/140)、[#138](https://github.com/zuizui0223/chocho/pull/138)、[#137](https://github.com/zuizui0223/chocho/pull/137)、[#136](https://github.com/zuizui0223/chocho/pull/136)、[#135](https://github.com/zuizui0223/chocho/pull/135)、[#134](https://github.com/zuizui0223/chocho/pull/134)、[#133](https://github.com/zuizui0223/chocho/pull/133)、[#132](https://github.com/zuizui0223/chocho/pull/132)、[#131](https://github.com/zuizui0223/chocho/pull/131)、[#130](https://github.com/zuizui0223/chocho/pull/130)、[#129](https://github.com/zuizui0223/chocho/pull/129)、[#128](https://github.com/zuizui0223/chocho/pull/128)、[#127](https://github.com/zuizui0223/chocho/pull/127)、[#126](https://github.com/zuizui0223/chocho/pull/126)、[#86](https://github.com/zuizui0223/chocho/pull/86)、[#85](https://github.com/zuizui0223/chocho/pull/85)、[#84](https://github.com/zuizui0223/chocho/pull/84)、[#83](https://github.com/zuizui0223/chocho/pull/83)、[#82](https://github.com/zuizui0223/chocho/pull/82)、[#81](https://github.com/zuizui0223/chocho/pull/81)、[#80](https://github.com/zuizui0223/chocho/pull/80)、[#79](https://github.com/zuizui0223/chocho/pull/79)、[#78](https://github.com/zuizui0223/chocho/pull/78)、[#77](https://github.com/zuizui0223/chocho/pull/77)、[#76](https://github.com/zuizui0223/chocho/pull/76)、[#75](https://github.com/zuizui0223/chocho/pull/75)

## 完全再照合用テーブル
| PR | 判定 | base | head | name |
|---|---|---|---|---|
| [#69](https://github.com/zuizui0223/chocho/pull/69) | C | `analysis/gain-fragility-20261024` | `a72fdc1f7ea2` | Falsify proposed butterfly gain–buffering mechanism with fixed-degree plant-range null |
| [#70](https://github.com/zuizui0223/chocho/pull/70) | C | `paper/ecological-resource-homogenization-v01` | `ea05c4ad97aa` | Calibrate GEB ecological homogenization claims and verify submission-only blockers |
| [#75](https://github.com/zuizui0223/chocho/pull/75) | D | `main` | `7ca3710c8799` | Curate invasion histories and recover source HOSTS Location for host adoption vs resource redistribution |
| [#76](https://github.com/zuizui0223/chocho/pull/76) | D | `analysis/hosts-location-pathways-20261027` | `2b315701066a` | Map HOSTS Location host-use reports into WCVP native versus introduced plant regions for 239 butterflies |
| [#77](https://github.com/zuizui0223/chocho/pull/77) | D | `analysis/hosts-location-wcvp-origin-20261027` | `19aae0f4bff9` | Test whether butterfly diet breadth predicts introduced host-use records beyond generic HOSTS reporting bias |
| [#78](https://github.com/zuizui0223/chocho/pull/78) | D | `analysis/hosts-location-bias-model-20261027` | `06de0d0b6dc6` | Test native vs introduced HOSTS host-pair documentation overlap against fixed-source nulls |
| [#79](https://github.com/zuizui0223/chocho/pull/79) | D | `analysis/hosts-native-introduced-overlap-null-20261027` | `ba54c205e71f` | Audit independent invasive butterfly larval host evidence versus HOSTS record geography |
| [#80](https://github.com/zuizui0223/chocho/pull/80) | D | `analysis/invasive-larval-host-evidence-20261027` | `80d5623322a9` | Test whether introduced larval-host plants become shared butterfly resource hubs beyond HOSTS reporting bias |
| [#81](https://github.com/zuizui0223/chocho/pull/81) | D | `analysis/introduced-host-local-sharing-20261027` | `0f9872c4e014` | Bridge 239 documented butterfly host portfolios with independent field caterpillar detection |
| [#82](https://github.com/zuizui0223/chocho/pull/82) | D | `analysis/introduced-host-local-sharing-20261027` | `3ab929f2d75b` | Separate seed/live-plant butterfly arrival vectors from host-resource redistribution; correct all-observed 115 cohort |
| [#83](https://github.com/zuizui0223/chocho/pull/83) | D | `analysis/arrival-resource-channel-decoupling-20261027` | `58e995003436` | Connect chocho butterfly host-resource opportunities to independent Lepidoptera establishment and insect interception records |
| [#84](https://github.com/zuizui0223/chocho/pull/84) | D | `analysis/arrival-interception-establishment-dryad-20261027` | `190aa7da88c5` | Cross independently documented alien-insect establishment records to chocho 239 host-resource opportunities |
| [#85](https://github.com/zuizui0223/chocho/pull/85) | D | `analysis/established-insects-zenodo-species-20261027` | `f4c1fb3ad66f` | Map 35 independently established alien butterflies to recipient-specific WCVP hosts beyond Hawaii |
| [#86](https://github.com/zuizui0223/chocho/pull/86) | D | `analysis/invasion-recipients-wcvp-host-opportunity-20261027` | `a46ad39b32a8` | Measure local wild-host bottlenecks among 24 alien butterflies and verify plant-contaminant vs native-host pathway cases |
| [#100](https://github.com/zuizui0223/chocho/pull/100) | C | `analysis/uk-public-transect-event-hindcast-20261028` | `8e016a97c228` | Build global H3 r4 × WGSRPD3 fractional crosswalk for butterfly–host–landscape year panels |
| [#102](https://github.com/zuizui0223/chocho/pull/102) | C | `analysis/h3-wgsrpd3-grid-bridge-20261029` | `0ed18840caf3` | Island barrier tests: exact source-distance validation and original230 host recovery |
| [#116](https://github.com/zuizui0223/chocho/pull/116) | C | `main` | `216dff8230da` | Butterfly island expansion: test independent migration and alien-host resources using source-verified record histories |
| [#117](https://github.com/zuizui0223/chocho/pull/117) | B | `main` | `bb9f7d08e094` | Island port networks and butterfly host geography: genuine 17th-century VOC voyages and botanical cargo |
| [#118](https://github.com/zuizui0223/chocho/pull/118) | C | `main` | `bb4d8f00fbdb` | Global plant-to-butterfly resource convergence and historical trade legacy: source-locked pilot |
| [#119](https://github.com/zuizui0223/chocho/pull/119) | C | `analysis/global-colonial-host-resource-convergence-20261031` | `0712619d0b5b` | Global colonial history × plant/butterfly potential-resource beta: spatial crosswalk and negative test |
| [#120](https://github.com/zuizui0223/chocho/pull/120) | B | `analysis/lenzner-colonial-region-source-20261031` | `dfd49cfcf9fb` | Global host-resource convergence versus real butterfly communities: 14.6M GBIF records and colonial-history controls |
| [#121](https://github.com/zuizui0223/chocho/pull/121) | B | `main` | `b70c43213b83` | Island butterfly resources and three gates: host, climate, historical access (Q1–Q4) |
| [#122](https://github.com/zuizui0223/chocho/pull/122) | B | `main` | `b4f331c4eacc` | Draft: nutrient affinity, port-indexed butterfly hosts, and matched all-flora controls |
| [#124](https://github.com/zuizui0223/chocho/pull/124) | C | `analysis/island-mainland-host-opportunity-20261102` | `facede8ec886` | Draft test: no positive introduced annual-herb filtering with island port richness |
| [#125](https://github.com/zuizui0223/chocho/pull/125) | B | `main` | `ff642c878e67` | Validate whole-island geography and port forecasts with uncertainty |
| [#126](https://github.com/zuizui0223/chocho/pull/126) | D | `manuscript/island-ports-flora-assembly-jbi-v01-20261105` | `eb94b2cbe81e` | Draft Fig.4 gate: 51-island circuit resistance fails distance LOIO; verify commercial AIS access |
| [#127](https://github.com/zuizui0223/chocho/pull/127) | D | `analysis/island51-circuit-steppingstone-20261106` | `a32975a45d37` | Run observed marine circuits for 51 islands with official and native-AIS validation |
| [#128](https://github.com/zuizui0223/chocho/pull/128) | D | `analysis/island51-ocean-circuit-run-20261106` | `17a4768cbb6c` | [island51] Joint host and marine circuit prediction with independent audit |
| [#129](https://github.com/zuizui0223/chocho/pull/129) | D | `analysis/island51-host-circuit-joint-20261106` | `a9c608a47bd7` | Map and account for fixed island host–circuit prediction gains |
| [#130](https://github.com/zuizui0223/chocho/pull/130) | D | `analysis/island51-host-circuit-atlas-20261106` | `6d7b923313c3` | Test whether candidate-host identity predicts butterflies beyond host-list richness |
| [#131](https://github.com/zuizui0223/chocho/pull/131) | D | `analysis/island51-richness-control-20261107` | `84ae978772e3` | Test 51-island host/circuit prediction against pre-2011 botanical recording effort |
| [#132](https://github.com/zuizui0223/chocho/pull/132) | D | `analysis/island51-pre2011-effort-control-20261107` | `c09b2897c805` | Evaluate island host and marine-circuit predictors with both island-pair endpoints unseen |
| [#133](https://github.com/zuizui0223/chocho/pull/133) | D | `analysis/island51-double-unseen-longitude-block-20261107` | `c4026cd2c101` | Separate AIS vessel-density and wind circuit increments on two-unseen-island predictions |
| [#134](https://github.com/zuizui0223/chocho/pull/134) | D | `analysis/island51-double-unseen-AIS-vs-wind-20261107` | `b06d1defd155` | Do AIS and wind prediction gains persist after latitude and hemisphere controls? |
| [#135](https://github.com/zuizui0223/chocho/pull/135) | D | `analysis/island51-double-unseen-latitude-controls-20261107` | `aee2499d5dd2` | Does latitude-controlled AIS prediction depend on two port-heavy island groups? |
| [#136](https://github.com/zuizui0223/chocho/pull/136) | D | `analysis/island51-top-port-training-influence-20261107` | `c364148606d1` | Falsify AIS-island alignment using 199 same-hemisphere latitude-band graph relabelings |
| [#137](https://github.com/zuizui0223/chocho/pull/137) | D | `analysis/island51-AIS-latitude-matched-node-null-20261107` | `eabfb4bdccc9` | Test AIS network specificity against only geographically local island-label swaps (≤750 km) |
| [#138](https://github.com/zuizui0223/chocho/pull/138) | D | `analysis/island51-AIS-local-750km-node-null-20261107` | `b2a50cb598fe` | Reacquire actual 2001–2010 GBIF plant species identities at all original 194 ports; test on 51 unseen islands |
| [#140](https://github.com/zuizui0223/chocho/pull/140) | D | `analysis/island51-pre2011-public-flora-direct-20261107` | `b0fa17f0944d` | Does maritime connectivity predict both 2001–2010 port flora and early/later butterfly communities? |
| [#141](https://github.com/zuizui0223/chocho/pull/141) | D | `analysis/island51-three-community-marine-comparison-20261107` | `9b4018c110b7` | Audit public 2020 AIS vessel–port triplet archive for direct routes among original 194 ports |
| [#142](https://github.com/zuizui0223/chocho/pull/142) | D | `analysis/island51-zenodo2020-port-voyage-audit-20261107` | `1a4b5434a556` | Can even all 27 geographically nearby but unverified AIS port aliases rescue island-route coverage? |
| [#143](https://github.com/zuizui0223/chocho/pull/143) | D | `analysis/island51-zenodo2020-geoonly-27-coverage-envelope-20261107` | `e27c6e2335b4` | Audit independently published global ~1,300-port maritime network for original WPI194 compatibility |
| [#144](https://github.com/zuizui0223/chocho/pull/144) | D | `analysis/island51-mendeley-global-port-network-audit-20261107` | `59e50fd5c9ef` | Locate PortWatch's public AIS-derived global Maritime Connections table before WPI194 matching |
| [#145](https://github.com/zuizui0223/chocho/pull/145) | D | `analysis/island51-portwatch-global-connections-source-20261107` | `217724d3da62` | Pilot 12 directed PortWatch model links among four original WPI-verified ports |
| [#146](https://github.com/zuizui0223/chocho/pull/146) | D | `analysis/island51-portwatch-four-port-directed-connections-20261107` | `7e02f3aeefdc` | Recover original GLSN2018 weighted network edges and historical UNCTAD PLSCIQ source bytes |
| [#147](https://github.com/zuizui0223/chocho/pull/147) | D | `analysis/island51-GLSN2018-UNCTAD-raw-source-recovery-20261108` | `360865fed8f2` | Independently geoverify 17 original WPI ports shared by GLSN2018, historic UNCTAD and AIS2020 |
| [#148](https://github.com/zuizui0223/chocho/pull/148) | D | `analysis/island51-17-portwatch-independent-geography-20261108` | `689a1334360b` | Geoverify remaining 13 original WPI ports shared by GLSN2018 and UNCTAD but not indexed in AIS2020 |
| [#149](https://github.com/zuizui0223/chocho/pull/149) | D | `analysis/island51-PortWatch13-GLSN-UNCTAD-source-geoverification-20261108` | `6c2164b25f05` | Does historical 2006–2010 port liner connectedness predict later butterfly communities on 17 original islands? |
| [#150](https://github.com/zuizui0223/chocho/pull/150) | C | `analysis/island17-historic2006PLSCIQ-LO2-20261108` | `8124a086d9eb` | Draft Fig.4: map source-locked marine current corridors without claiming island relays |
| [#156](https://github.com/zuizui0223/chocho/pull/156) | C | `analysis/ncomms-public17-access-holdout-20261108` | `6367a7fe5a09` | Butterflies first: wing-size moderation of marine resistance across island GBIF records |
| [#157](https://github.com/zuizui0223/chocho/pull/157) | A | `analysis/chocho-butterfly-trait-marine-20261108` | `365526cc70e5` | One GEB butterfly paper: island resource amplification, port floras and marine geography |
| [#158](https://github.com/zuizui0223/chocho/pull/158) | C | `analysis/butterfly-trait-public51-20261109` | `5119bdaf2464` | Test matched island port vs nonport coast: historical candidate hosts and later butterfly indices |
| [#159](https://github.com/zuizui0223/chocho/pull/159) | B | `analysis/butterfly-trait-public51-20261109` | `d465e09b30d3` | Historical GBIF validity and port-host × marine-access tests on 225 butterflies × 51 islands |
| [#160](https://github.com/zuizui0223/chocho/pull/160) | B | `analysis/port-host-ocean-interaction-20261114` | `bd0b17abe43a` | Draft: test host–marine predictive overlap and native/alien host record links on islands |
| [#161](https://github.com/zuizui0223/chocho/pull/161) | B | `analysis/host-marine-bridge-audit-20261116` | `e21cd107b312` | Draft: connect native/introduced port-host indexing with history-corrected butterfly records |
| [#162](https://github.com/zuizui0223/chocho/pull/162) | B | `analysis/native-introduced-temporal-bridge-20261116` | `b5b245aa4b8e` | Draft: test geographic HOSTS use records against introduced plants, FAO crops, island origin and ports |
| [#163](https://github.com/zuizui0223/chocho/pull/163) | B | `analysis/hosts-real-use-island-origin-ports-20261116` | `de8fd4f264f3` | Draft: test island alien-host literature signal against source-locked whole- and family-matched WCVP flora |
| [#164](https://github.com/zuizui0223/chocho/pull/164) | D | `analysis/hosts-use-vs-background-flora-20261116` | `1797d66ce98c` | Draft: source-gate EuDiS sea-dispersal traits against native/alien WCVP island floras |
| [#165](https://github.com/zuizui0223/chocho/pull/165) | D | `analysis/eudis-thalasso-island-source-gate-20261116` | `f01daddc4939` | Draft: verify historical OSCAR v2 ocean-current data access before directional seed-drift tests |
| [#166](https://github.com/zuizui0223/chocho/pull/166) | D | `analysis/oscar-v2-1993-2010-source-access-gate-20261117` | `21b56c379fda` | Draft: measure actual 2005 directed offshore currents between source-matched Atlantic island regions |
| [#167](https://github.com/zuizui0223/chocho/pull/167) | D | `analysis/atlantic-legacy-oscar-directed-2005-pilot-20261117` | `6bd7fbfe3b94` | Draft: test decadal persistence of directed Atlantic OSCAR flow alignment, 2001–2010 |
| [#168](https://github.com/zuizui0223/chocho/pull/168) | D | `analysis/atlantic-oscar-decadal-robustness-20261117` | `2df111bc7855` | Draft: crosswalk native-source thalassochorous plants to regional recipient flora along Atlantic corridors |
| [#169](https://github.com/zuizui0223/chocho/pull/169) | D | `analysis/atlantic-eudis-wcvp-source-dest-20261117` | `36ee18bc2c12` | Draft: test whether ocean-dispersal trait predicts Atlantic recipient flora beyond plant-family composition |
| [#170](https://github.com/zuizui0223/chocho/pull/170) | D | `analysis/atlantic-legacy-oscar-directed-2005-pilot-20261117` | `772d11f8c0d8` | Draft: advect real Atlantic 2005 NOAA OSCAR currents with source-verified GSHHG coast mask |
| [#171](https://github.com/zuizui0223/chocho/pull/171) | D | `analysis/atlantic-2005-oscar-lagrangian-pilot-20261117` | `1d272990e99d` | Draft: resolve NOAA particle coast contacts to original GSHHG landmasses, not false island arrivals |
| [#172](https://github.com/zuizui0223/chocho/pull/172) | D | `analysis/atlantic-gshhg-exact-landmass-contact-20261117` | `c092327d4fdd` | Draft: audit same-port AIS voyage × host × butterfly observation-effort identifiability |
| [#173](https://github.com/zuizui0223/chocho/pull/173) | D | `analysis/port-nonport-voyage-effort-identifiability-20261117` | `5008e8021711` | Draft: test why indexed ports correlate more strongly with alien flora on small/native-poor islands |
| [#174](https://github.com/zuizui0223/chocho/pull/174) | D | `analysis/port-alien-flora-island-moderators-20261117` | `830b5e246872` | Draft: test nonlinear and historical botanical-index confounding of island port–alien-flora moderation |
| [#175](https://github.com/zuizui0223/chocho/pull/175) | B | `analysis/port-alien-flora-nonlinear-historical-source-20261117` | `93730965d81b` | Draft: cross-test complete WCVP native vs alien island floras against sea geography and verified shipping coverage |
| [#176](https://github.com/zuizui0223/chocho/pull/176) | B | `analysis/126-island-native-alien-flora-crossed-marine-20261117` | `22e9f6a04bed` | Draft: measure island-endemic butterfly host reporting from an independent checklist denominator |
| [#177](https://github.com/zuizui0223/chocho/pull/177) | A | `analysis/independent-island-butterfly-host-coverage-20261117` | `417bbc77f273` | Draft: island-endemic butterfly external validity, host taxonomy and nonidentifiable population-level expansion |

## 対応するmainの結果台帳
- [科学的結果と限界](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md)
- [直近1週間の成果](LAST_WEEK_POSITIVE_EVIDENCE_20261118.md)
- [mainへ保存した811資産・Git blob manifest](SELECTIVE_RESULT_MERGE_20261118.md)
- [56件の関連しないPRクローズと別研究の保存](CLOSED_OUT_OF_SCOPE_PRS_20261118.md)

この一覧は決定・作業順の記録。PRの状態は実行後のGitHubライブ結果で別途更新する。

## 実行記録：第一波のクローズ完了（2026-11-18）

2026-11-18時点の最新GitHub PR状態を再検証。上記 **Cグループ11件すべて** を未マージのままクローズした（#69, #70, #100, #102, #116, #118, #119, #124, #150, #156, #158）。成功11件、エラー0件。コミット・元ブランチは削除していない。旧モデルの結果・境界は[全PR台帳](ALL_PR_AUDIT_20261118.md)・[主張台帳](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md)・対応PR本文に保持。

**最新のオープンは60件 = A 2件 + B 12件 + D 46件。** 全176PRの現状態はマージ46、クローズ未マージ70、オープン60。上の一覧「71」は今回の分類作業開始時点の固定母集団であり、現在オープン数と混同しない。

**次の主作業はBのコード再現性の統合。** Bグループ差分の概要（親枝との差分でありmainへの新規ファイル数とは異なる）：
- #120: 解析スクリプト13・テスト9、#121: スクリプト57、#122: スクリプト22・テスト19、#159: スクリプト22。
- #160〜#163: それぞれスクリプト2/2/6/6、#175: 15。これらを一括mainへ入れると旧CIや未認定モデルが混入するため、**本文5つの主張を再現する最小実行入口＋入力ハッシュ＋結果検証だけを選択**。
- #117/#125/#176は追加原ソース/テスト/限定データ中心。#157/177は原稿・最新出典・反証の中心として継続。

**Dの46件**は今回閉じていない。#126〜#174の海洋／港ネットワーク系と#75〜#86の記録／到着系から、最終的なソースと限界をすべて含むかをmain台帳で照合した後、コードを残す必要のあるもの以外は閉じる。未チェックのコードを完了と称して閉じない。
