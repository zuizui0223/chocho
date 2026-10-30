# Task 2 修復取得と寄主joinの独立レビュー

## 現在の判定

**植物取得修復・joinコード・実join出力: PASS。修正必須blockerなし。** 全2385候補、93蝶coverage、28target集計を別実装で独立に再構成し、実出力と全項目が一致した。

対象は`repair_island_plant_facets_20261030.py`、`analyze_island_dated_hosts_20261030.py`、各tests、修復planおよび実修復sourceである。共有repo・cache・source・remoteには変更を加えていない。既存テストの再実行や追加のネットワーク取得はせず、実source archiveを読み取り、独立した集計コードで検証した。本ファイルのみを作成・更新する。

## 1. 修復方針と実装

修復planは`provenance/archive/host_hindcast/ISLAND_PLANT_FACET_REPAIR_20261030.md`。元100-bin paginationで重複した3島だけを対象とする。元25 completeを再取得・置換しない。元結果、原応答archive、失敗statusを保存した上で、別manifestとして修復と統合を記録する。

単一pageはfacetLimit=10000。返却bin数がこの上限未満であることに加え、各facetについて返却実bin数nをoffsetとするlimit1のterminal probeを行う。4つのprobeが同じ全体countを維持し、空である場合にのみ先へ進む。その後、full requestを再取得し、全4facetのkey→count写像とcountを比較する。tieの順序には依存しない。

このため、APIが要求10000より小さいpageへ黙示的に制限しても、まだ残るbinがterminal probeに出ればmissingになる。元の重複binを単に削除してcompleteへ変える経路はない。単一page内のduplicate、上限到達、count変化、再取得時のmap変化、YEAR/DATASET margin不整合、期間外年を拒否する。

## 2. 修復sourceの実検証

元manifest SHA256は前回レビュー時から変わっていない。

- original manifest: `8de586ffb615fb70d1e5c2debefb5eba691db5d49402e2f14847fa58dda2fd48`
- original archive: `d93166643d05864f223f1488f389f46120bbf945a439e1cb9490c0326e69803e`
- repair script: `f9ab251c1815b243c0653b41993cf18d135fd438b5e0dc43fab56dc04a34659d`
- repair plan: `f6792c453cfea9bfb74f9820d799fdecc0ab54a93b30035ac7e7c582323acdd2`
- repair manifest: `b506d33d56b67d23b816579b9804c0539a68a79d7014eeb6ddccb0752b70097e`
- resolved manifest: `a38711cb91afbe7aae86d24477ed046db1101849b4e60bb9d1fdea6d899cf915`

修復lockは`2026-10-30T02:13:33.132158+00:00`に作成され、元manifest・元script・repair script・planと失敗島集合を固定している。現物の各hashと一致した。

| 島 | 全植物records | species bins | speciesに割り当てられたrecords | terminal offsets: species / basis / dataset / year |
| --- | ---: | ---: | ---: | --- |
| `gshhg237_0667315a23db9839` | 300,967 | 2,669 | 295,045 | 2669 / 5 / 78 / 8 |
| `gshhg237_4dbf385289bcf888` | 175,920 | 2,017 | 173,552 | 2017 / 5 / 81 / 8 |
| `gshhg237_e8ca9c42f6d3f947` | 12,994 | 712 | 7,431 | 712 / 3 / 18 / 8 |

3島とも6 requests（full + terminal4 + full recheck）、合計18 requests、retry0。上記12 terminalはすべて空で、全体countはfull・terminal・recheckを通じて同じだった。countは初回失敗取得の最初の応答とも一致していた。

### 実施した独立検査

- 全18 request/response/receiptのhashとresponse byte数を検証。
- 各実requestを元失敗島の実requestと比較し、変更はfacetLimit、facetOffset、terminal時のfacet選択だけと確認。predicate、checklist、limit0、期間、PRESENT、座標条件、geometryは同じ。
- 全responseがresults=[]であり、個別occurrence recordを含まないことを確認。
- 3島のfullとrecheckから独立にmapを構成し、duplicateなし・同map・各manifest facetと一致を確認。
- terminalの実offsetは要求limitではなく各実bin数であり、limit1、対象facetだけ、空であることを確認。
- species/basis residualは負でなく、YEARとDATASET合計は全植物countと一致。全YEARは2010–2017。
- 元25 completeのresult objectはresolved manifest内で完全同一。元3 missingだけが修復resultと一致する形で追加された。
- resolved対象は元固定28島と完全一致。全28島のgeometry SHA256は固定selectionと一致。
- source-route欄は25 original / 3 repairを正しく区別。
- resolved5 CSVのhashと、4 facet CSVの全bin・countを再計算して一致。

修復archiveは`early_plants_repair_raw_facets.tar.gz`、119,132 bytes、54 JSON member、SHA256 `583cba2c4e9092affb555671a83bfea4c0e7e3a096398bd599bc2dcbbf55b492`。全memberがmanifest参照18 triplesと正確に対応し、余分・重複memberはない。sort順、gzip/tar timestamp0、tar uid/gid0、mode0644も確認した。

統合結果の**28 complete / 0 missing / 0 pending、728,040 records、botanical-support16島**を独立に再計算した。元25島の238,159件に修復3島の489,881件を加えた値である。補われた3島のうちspecies-assignmentが低い島も残し、support screenで島を除外していない。

## 3. joinコードの判定

`analyze_island_dated_hosts_20261030.py`は元sourceのhash/count検証、taxonomy出力hash、1062 source-ID集合、固定28島とgeometry、植物query scope・期間・marginsを検証する。旧exact-null入力loaderを通じ、元risk mask・observed table・margins・距離scoreの固定hashと整合性も確認する。primary2385候補・66陽性以外では停止する。

`build_catalogues()`は蝶ごとに元source plant IDの集合を維持し、eligible行だけからaccepted keyのsetを作る。複数source IDsが同じaccepted conceptへ写っても、植物recordを足すときはそのconceptを一度だけ数える。source-ID数、accepted-concept数、未解決ID数、taxonomy request missing数を別に保持する。

CSVの`True`/`False`を明示的にparseし、文字列`False`を真とみなさない。primary/accepted-verification FAILEDか、keyのないeligible行は拒否する。taxonomy科学的拒否とsource-name未照合は未解決として保持し、確認済み別hostの陽性まで捨てない。

host record状態は次を区別する。

| 条件 | 出力状態 |
| --- | --- |
| 植物query失敗 | `PLANT_QUERY_MISSING`、host数/record数/confirmedはNone |
| 解決済みのcatalogue-compatible植物記録あり | `RECORDED_KNOWN_HOST`。別hostが未解決でも、その未解決数を併記した上で確認済み記録を保持 |
| matchなし、taxonomy未解決を含む | `NO_MATCH_UNRESOLVED_HOST_CATALOGUE` |
| matchなし、catalogue全IDを解決済み | `NO_MATCH_COMPLETE_RESOLVED_CATALOGUE` |

確認割合の分母は植物query completeの候補である。taxonomy未解決が残る場合も、これは「確認できたcatalogue-compatible植物記録の割合」であり、host存在率やabsence推定ではない。コードはこの記述的用途と期間facetの限界を明示し、host効果モデル、導入/在来の区別、摂食、定着推定を追加していない。

初読時にCLI defaultが元25/3 manifestだった点を親へ通知した。親はdefaultを`early_plants_resolved_manifest.json`へ変更し、実行時も明示指定した。この変更は取得結果・科学的条件を変えず、通常のreplayが完成した28島を使うための修正である。変更後のdefaultと実summaryのinput記録を確認した。

## 4. 実join出力の独立検証

実script出力を見る前に、元primary risk cells、元1440 catalogue edges、最終taxonomy bridge、verified28島species facetsから独立実装で期待値を作成した。joinコードの関数は使用していない。期待値の全行をkey-sortしたJSONに正規化し、3表のhashを保持した。その後、実CSVを同じ型（整数、浮動小数、明示boolean、None、文字列）へ正規化して照合した。

| 表 | 独立検証行数 | 期待値と実出力のcanonical SHA256 |
| --- | ---: | --- |
| primary risk cells | 2,385 | `099087e0c415281f155c117a96a31fad4b9dedc4eb790ac61e5ebaf4753473ec` |
| butterfly catalogue coverage | 93 | `105f29769fc477f50c6dc699dea6f0d64c84508890f2529b7d1b3001909726c9` |
| target support | 28 | `063fb6344b092f08819e1af4dabbd895fe714f8e7f0baf0f6fcec6eb2beaba3f` |

これは件数だけの一致ではない。全候補の蝶key、島、longitude block、後期陽性label、後期record数、元source距離、catalogue各数、全植物記録、support screen、host記録key集合・数・record合計、未解決/failure/confirmed各状態が一致した。

### 実集計の確認

| 指標 | 独立再計算と実summary |
| --- | ---: |
| 元primary候補 | 2,385 |
| 元後期陽性 | 66 |
| plant query complete候補 | 2,385 |
| 前期のcatalogue-compatible植物記録を確認 | 336 / 2,385 = 14.0881% |
| 後期陽性内で前期のcatalogue-compatible植物記録を確認 | 20 / 66 = 30.3030% |
| `RECORDED_KNOWN_HOST` | 336 |
| `NO_MATCH_UNRESOLVED_HOST_CATALOGUE` | 1,961 |
| `NO_MATCH_COMPLETE_RESOLVED_CATALOGUE` | 88 |
| `PLANT_QUERY_MISSING` | 0 |
| 元catalogueリンク | 1,440 |
| 解決済みsourceリンク | 865 |
| 未解決sourceリンク | 575 |
| 少なくとも1つの寄主conceptが解決した蝶 | 84 / 93 |
| source名がすべて解決した蝶 | 4 / 93 |
| taxonomy source requestがmissingの蝶 | 0 |

taxonomyの科学的な未解決575リンクは残る一方、primary/accepted verificationの取得failureはない。この二つを混同していない。陽性66候補のうち確認がない46候補も、これだけでhost不在とは読めない。実scriptとsummaryは記録ベースの記述的結果として保持している。

summary内のprimary集計とcatalogue集計は、独立に計算した辞書と完全一致した。3 CSVの実byte hash、summaryが記録する全入力hash、source pins、plan hash、analysis script hashも一致した。実inputはresolved manifestであり、元25/3 manifestを誤用していない。

### 最終identity

- join script SHA256: `1008675a8948f514c7bfcbbcd280bd16cb3c1d4ebf3b77e7afa17b9fe92f7c05`
- `dated_host_recording_summary.json`: `9146d9eb58bcdc90beef9335284196ce8e18734dd45667d96d0d243b79f8f583`
- `plant_taxonomy_manifest.json`: `680ad64ecaff0ad4f9cb8dd30e428371f1a08c5cc05c6391934a446254ca3af9`
- `plant_taxonomy_bridge.csv`: `23ae32b4b09b4f38b5ce5357cd56ea3e2c880dffb52f2584becaed6f6da03b2b`
- `dated_host_primary_risk_cells.csv`: `cf34d03f3d01242bd182557114dca711f3efd47e0d1d066e6d7d0548c22ecbe8`
- `dated_host_butterfly_catalogue_coverage.csv`: `958592af6ee654fdf5567974338e0736b13e01ac7c5d96c6ac1739ac5f4528cb`
- `dated_host_target_support.csv`: `25e462ec3c328f5469e4b229fda0cbde54d7eb49a5f0bc5dafaed2e401aae9e8`

## 最終結論

取得修復は元の対象・期間・geometry・科学的predicateを維持し、providerのpagination問題を完全性確認付きの別取得で処理している。joinは固定2385/66と1440 sourceリンクを保持し、current conceptの重複、未解決taxonomy、query failure、確認済み記録を適切に区別する。**提出・統合を妨げる問題はない。** 親が報告した最終suite631 passed・5 subtests passedは再実行せず、ここでは実sourceと実derived tablesの独立検証を行った。
