# Task 2 取得実装の独立レビュー

## 判定

**取得実装に修正必須の blocker は見つからない。最初の取得結果は25島 COMPLETE、3島 MISSINGであり、全28島の完全取得とは記載できない。** 3島の失敗は保存されたGBIF応答のfacetページ重複で説明でき、部分的な種集計を採用せず `count=None` とした処理は仕様どおりである。取得コード・first-run結果・原応答archiveを維持し、別の明示した取得方式で3島だけを修復する方針は妥当。

レビュー対象は取得部分のみ。植物分類照合と寄主joinは別レビュー対象である。稼働中のスクリプトは編集せず、既存10テストを重複実行せず、外部APIへの追加読み取りも行っていない。親から提供された公式API schemaのローカルコピーを読んだ。

## 対象と固定identity

- 実装: `/workspace/scratch/10d2ad986b7b/chocho/scripts/acquire_island_early_plants_20261030.py`
- 実装SHA256: `b42b8cb7c2e39ec016c40cb03546e5286464346ab29b480e4cd01625aac9724a`
- 仕様: `/tmp/chocho-dated-workflow/task-2-brief.md`
- diff: `/tmp/chocho-dated-workflow/task2-acquisition-review.diff`
- first-run manifest: `data/exploratory/island_dated_host_20261030/early_plants_manifest.json`
- manifest更新日時: `2026-10-30T02:09:47.202264+00:00`
- manifest SHA256: `8de586ffb615fb70d1e5c2debefb5eba691db5d49402e2f14847fa58dda2fd48`
- run lock作成: `2026-10-30T02:05:47.822834+00:00`
- 最初の保存済みresponse receipt: `2026-10-30T02:06:00.379134+00:00`

現在のスクリプト、plan、target-support CSV、selection JSONの実bytesはrun lockのSHA256と一致した。

## 仕様との照合

| 検査項目 | 判定・根拠 |
| --- | --- |
| 固定28島 | source hashを検証し、元96島のsupport列から固定28島だけを選択。最終manifestの島集合はその28島と完全一致。botanical screenで島を除外していない。 |
| 日付 | 全保存requestでYEAR >=2010、YEAR <=2017。complete結果の全YEAR facetもこの期間内。 |
| 分類namespace | 全requestでPlantae `P` とCOL XR `7ddf754f-d193-4cc9-b351-99906754a03b`。POST predicateとparametersの双方を実requestで確認。 |
| 点を取得しない | 全152 requestはlimit=0、offset=0。全152 responseのresults=[]を確認。保存前にも同条件を検証する。 |
| 位置条件 | 全requestで座標あり、indexed geospatial issueなし、PRESENT。withinの実WKT SHA256はそれぞれの固定geometry SHA256と一致。simplification/rectangleへの代替はない。 |
| source geometry復元 | SHA-pinned GSHHG archiveを検証する既存helperを使用。全28 geometryを検証後に最初の植物queryへ進む。現在の実requestでも各geometry hashを独立照合した。 |
| facet完全性 | 100 bins/page、最大100 pages/island。count変化、重複bin、異常page、dataset/year margin不整合、first-page recheck不一致でfailする。25 complete島の全pageを独立に再集計しmanifestと一致。 |
| 欠測とzero | 失敗3島はcount=None、facetsなし。部分的な種facetは出力CSVに入らない。complete zeroは別扱い。今回のcomplete zero実績は0。 |
| botanical screen | >=20 records、species assignment >=90%、>=3 positive yearsの記述的screen。queryの選択、対象島の変更、後期取得に使用されていない。 |
| rate limit/retry | clientは各request startの前に最低1.5秒を確保し、1 query最大3 attempts。429はRetry-Afterを尊重。receiptにはretrieved時刻のみのため実際のstart間隔はartifactから直接再測定できないが、コードの実装は明確。 |
| archive | 全152組=456 JSONを原bytesのまま保存。member集合はmanifestが参照するrequest/response/receiptの集合と完全一致。sort順、mtime=0、uid/gid=0、mode=0644、gzipのmtime=0とfilename省略を確認。 |

## 実データによる独立確認

単体テストの再実行に代えて、完了したfirst-runの原応答を独立した読み取りコードで検証した。

- 152 request/response/receipt tripleの全hashとresponse byte数が一致。
- 25 COMPLETE島について4種類のfacetを原pageから再構成し、manifestの全bin・count・unassigned-record数と一致。
- 25島のfirst-page再取得が最初のpageと一致。
- 全complete島でDATASET_KEYとYEARのcount合計が全植物record countに一致。
- 最終4 facet CSVをmanifestから独立再構成した期待値と照合し一致。
- 最終5 CSVの全output SHA256が一致。
- 3 MISSING島はcount=Noneかつ部分facetsがないことを確認。
- archiveの456 memberは全てcacheの元bytesと一致し、余分なmember・duplicate memberはない。

最終集計は28 requested、25 complete、3 missing、0 pending、238,159 complete indexed plant records、14 botanical-support島。これは完全取得できた25島の集計であり、欠測3島を0として加えた結果ではない。

archive: `data/exploratory/island_dated_host_20261030/early_plants_raw_facets.tar.gz`

- bytes: `253559`
- files: `456`
- SHA256: `d93166643d05864f223f1488f389f46120bbf945a439e1cb9490c0326e69803e`

## 取得元の3件の不完全性

| 島 | 保存された重複 | 判定 |
| --- | --- | --- |
| `gshhg237_0667315a23db9839` | offset100と200に`4DMBB`、`58PNJ`。各count365。全体count300,967は一定。 | 同件数の並びがpage境界で安定しない挙動と整合的。原因の断定には追加source確認が必要。 |
| `gshhg237_4dbf385289bcf888` | offset100と200に`324SZ`、count205。 | 同様。 |
| `gshhg237_e8ca9c42f6d3f947` | offset100と200に`6HYB3`、count12。 | 同様。 |

これらを重複削除だけでcompleteへ昇格させると、同じ境界で脱落した別binを見落とせる。現実装がfailした判断を維持する必要がある。

## 公式API schemaと修復方式の条件

親が取得した公式source:

- URL: `https://techdocs.gbif.org/openapi/occurrence.json`
- local: `/tmp/chocho-dated-workflow/official_occurrence_api.json`
- bytes: `225700`
- SHA256: `ffde21586ee99495708d1fcb37cb88f792c50605fbb7a263217b3936f07fce39`

GET `/occurrence/search`の公式パラメータには以下がある。

| パラメータ | schemaにある型・制約 |
| --- | --- |
| `facet` | string。複数指定可。 |
| `facetMincount` | integer/int32。 |
| `facetMultiselect` | boolean。 |
| `facetLimit` | integer/int32。maximumの記載なし。 |
| `facetOffset` | integer/int32、minimum=0。 |

`facetSort`はschema全体に存在しない。したがって未確認のソート引数を加えて解決したとは扱えない。

POST `/occurrence/search/predicate`の説明はfacet/limitをPOST bodyに置くことを示す。一方、参照される`OccurrencePredicateSearchRequest` schemaは`predicate`だけを列挙しており、facet propertiesと上限の記載は不完全である。このschemaからfacetLimit=10000の実際の受理上限を保証できない。

**有界の単一pageへの修復は可能だが、返却bin数が要求limit未満であることだけではsilent server capを排除できない。** 次の条件を新しい修復仕様へ固定することを推奨する。

1. 元3失敗島だけ、同じgeometry・期間・分類・count-only条件で取得し、元失敗結果とarchiveを保持する。
2. 単一pageで各facet内のkey重複を拒否し、同じrequestの再取得でkey→count集合と全体countが一致することを確認する。件数が同じbinの順序そのものには依存しない。
3. 各facetについて実返却bin数nをoffsetにしたterminal probeが空であることを確認する、または別に検証したspecies-assigned record countとspecies facet合計を完全一致させる。要求limit未満という理由だけでcompleteにしない。
4. speciesのunknown assignmentを許す元queryではspecies合計が全植物count未満でも正当な場合がある。その差分だけでtruncationを否定しない。
5. 上限到達、terminal probeのnonempty、count変化、再取得不一致は引き続きmissingにする。

この修復は観測対象・時点の選び直しではなく取得の完全性を確認する変更として記録できる。修復結果自体はまだ本レビューの承認範囲に含めない。

## 軽微な記述上の注意

仕様の「three bounded retries」に対しコードは最大3 attempts、つまり初回+最大2 retriesである。公開説明は実装どおり「最大3回の試行」とするのが正確。今回transport retry実績に依存しない、データ完全性を損なわない相違であり、first-runの再取得やコード変更を要求するblockerではない。

`FINISHED`は処理が28島に到達して終了した状態を表す。scientific readoutでは必ず25 complete/3 missingを併記し、全28島の完全な植物情報が得られたと表現しない。
