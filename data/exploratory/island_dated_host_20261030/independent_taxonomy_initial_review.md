# Task 1 植物taxonomyの独立レビュー

## 判定とレビュー範囲

**仕様適合性または科学的判定を損なう修正必須blockerは見つからない。コードレビューはPASS。** 固定したsource rosterの保持、exact synonymの受容条件、競合概念の拒否、family/kingdom確認、重複したaccepted conceptの保持、request failureと科学的な不一致の区別は仕様に沿っている。

この判定は実装、テスト内容、および進行中cacheの少数の具体例に対するものである。最終1,040 primary照合、accepted-name検証、1,062-row bridge、最終collision inventory、最終archiveの完了・hash検証を代行したものではない。レビュー最終観測時は2026-10-30T02:16:32.938180Z、primary 600/1,040 SUCCESS、accepted-name verification receiptはまだ0件だった。

指示どおり実行中のscript/test、取得cache、source、remoteは変更していない。subagentは起動せず、追加API requestと既存テストの再実行もしていない。18 focused tests/613 full-suite testsの実行実績は担当者reportを証拠とし、内容を静的に確認した。

## 対象identity

- root: `/workspace/scratch/10d2ad986b7b/chocho`
- brief: `/tmp/chocho-dated-workflow/task-1-brief.md`
- 実装者report: `/tmp/chocho-dated-workflow/task-1-report.md`
- diff: `/tmp/chocho-dated-workflow/task1-review.diff`
- script: `scripts/bridge_island_plant_taxonomy_20261030.py`
- script SHA256: `5f26c60c1a2b743b4463e4e1744cee20f2ca8f4e923c1a7cd8f2bcb774862f75`
- tests: `tests/test_island_plant_taxonomy_20261030.py`
- tests SHA256: `717eccb734bf59b405e24e92c42d68d12cd66ec4b55b3ec430a72411e8ed006b`
- read-only cache: `/tmp/chocho-plant-taxonomy-20261030`

## 固定sourceとquery対象

source CSVを別の読み取りコードで照合し、3つのhashを実bytesから再計算した。

| source | 実SHA256 |
| --- | --- |
| `plant_taxon_query_candidates.csv` | `a1cc87b20fbd19415fa81ca4269b8414b0c32167bb6033ca0425354c89c1906f` |
| `original230_host_edges.csv` | `e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22` |
| `baseline93_host_coverage.csv` | `f702ce482ca6b6411bb4530aa8b829400300cb5ef2464e367eef4e5a365ce9aa` |

early93のoriginal butterfly nameを基準に元230 rosterから辺を抽出する操作は、93蝶・1,440辺・1,062 unique plant IDs・1,040 eligible query names・22 unqueried source namesを再現した。現取得cacheにも22 unqueried IDのprimary receiptは存在しなかった。source IDは文字列として保持・sortされ、数値変換・名称切り捨て・後期記録による選別はない。

`load_sources()`はこの3hashと固定件数の双方を検証してからrequestを生成する。元1,440辺を改変していない。全1,062 IDを最終bridgeへ出す設計になっている。

## 判定ロジック

| 仕様 | 実装での確認 |
| --- | --- |
| explicit current namespace | `/v2/species/match`へCOL XR checklist、scientificName、kingdom=Plantae、source family、taxonRank=SPECIES、strict=true、verbose=trueを明示。確認した実receiptも一致。 |
| exact canonical | `diagnostics.matchType=EXACT`に加え、usageのcanonicalNameと入力名の一致を要求。空白正規化だけを許し、fuzzy/higherrankを許容しない。 |
| usageとaccepted concept | usageはACCEPTED/SYNONYM、両者のrankはSPECIES、accepted keyは空でないstring、accepted statusはACCEPTEDが必要。synonym flag/statusの不整合、acceptedUsageのないsynonymも拒否。 |
| kingdom/family | classification内のPlantae/Pとsource familyの完全一致を要求。別familyへ移されたsource名もこのcontinuationでは未解決のままになる。 |
| ambiguity/homonym | ambiguity診断を拒否し、同じcanonical名・EXACT・species rankで別accepted keyを持つalternativeを拒否。低rankingやAMBIGUOUS_SYNONYMを理由に競合を無視しない。 |
| 同じconceptの別usage | alternativeが同じaccepted keyへ戻る場合は別conceptと扱わない。usageの重複だけで安全なaccepted conceptを失わない。 |
| 欠落したaccepted status | 欠落値を既定ACCEPTEDに置換しない。accepted canonical nameを同じnamespace/familyで別queryし、そのresponse自体が全gateを通り、usage.status=ACCEPTED、usage/accepted双方のkeyが元keyと一致する場合だけ受容。 |
| verificationの連鎖 | 別queryも同じ`assess_match()`で厳格に判定し、未確認statusをさらに暗黙に補う経路はない。primary側の競合・family等のrejection reasonはverification成功でも消えない。 |
| 重複concept | `annotate_collisions()`はeligible source行を削除せず、accepted keyごとのsource IDs/countを注記し、collision CSVを作る。 |

## 実responseの小標本

最初の8 source IDsのrequest/receipt/raw hash・主要taxonomy項目を確認した。追加で2つの実synonym例と数件の既知名を確認した。全件の結果集計や再取得は行っていない。

| source ID / source name | 原応答に見える条件 | 実装の判定 |
| --- | --- | --- |
| `107865` / Lantana achyranthifolia | EXACT、ACCEPTED、SPECIES、Verbenaceae/Plantae、競合なし、key `3S7KH` | eligible |
| `107386` / Lamium amplexicaule | primaryはconfidence100/EXACT/ACCEPTEDだが別exact concept `VL8C2`あり | `competing_exact_current_species_concepts`で拒否 |
| `114913` / Livistona chinensis | alternativeは別usageだが同じaccepted key `3VP7Y`へ戻る | eligible。usage数を概念数と混同しない |
| `117504` / Magnolia champaca | exact synonym alternativesはいずれもaccepted key `3XH7J` | eligible |
| `2369031` / Medicago sativa | distinct exact concepts `3YNKJ`, `VL7SM`, `VL7TC`あり | 競合として拒否 |
| `154769` / Phyllanthus reticulatus | SYNONYM → Kirganelia reticulata / `3R8FX`。acceptedUsageにstatusなし | secondaryなしではeligible=False、`accepted_status_not_verified` |
| `212400` / Verbena laciniata | SYNONYM → Glandularia laciniata / `6KH8R`。acceptedUsageにstatusなし | secondaryなしではeligible=False、`accepted_status_not_verified` |

最後の2例はraw応答とreceipt hashも一致した。親が承認したaccepted-name verificationが必要になる実API形状を実際に確認できた。verificationの実結果はまだ取得前だったので、将来の成功を先取りしていない。

## failure、unresolved、zeroの境界

`query_status=SUCCESS`はHTTP/JSON取得成功であり、科学的に利用可能なtaxonomyの意味ではない。以下の分離が実装されている。

| 状態 | primary query_status | accepted_verification_query_status | conservative_join_eligible |
| --- | --- | --- | --- |
| 原名が今回のspecies query対象外 | UNQUERIED_SOURCE_NAME | 空 | False |
| primary request失敗 | FAILED | 空 | False |
| primary取得成功だが競合/不一致 | SUCCESS | 空または該当verificationの状態 | False |
| primary synonymの必須verification失敗 | SUCCESS | FAILED | False |
| 全scientific gateを通過 | SUCCESS | 不要なら空、必要ならSUCCESS | True |

失敗行を削除せず、primary failureは`taxonomy_request_failed`、source-name未解決は専用理由、verification不成立は`accepted_status_not_verified`を保持する。候補accepted keyが残る科学的拒否行もあるため、join側はkeyが空でないことだけで採用せず、明示的なeligible=Trueを条件にする必要がある。CSVの文字列`False`をtruthyとして扱ってはならない。

本scriptはoccurrence countを生成しない。したがってtaxonomy failureから植物recordの0を生成する経路はない。join側には、全source辺を保持したままeligible accepted keyの集合を蝶ごとに作り、conceptを一度だけ数えるための情報が揃っている。その実際のdedupとmissing-state伝播はjoin実装のレビューで確認する。

## 品質・provenance

単一Clientがprimary・accepted verification・retryの全taxonomy requestを逐次処理し、開始間隔>=1.05秒を保つ。最大3 attempts、45秒timeout、8MB response boundを実装。既に保存されたFAILED receiptを再利用したときに追加の4回目へ進まない。cacheはURLとraw body hashを確認し、HTTP/JSON成功とscientific gateを混ぜない。

GETの完全なquery URL、method、各attempt時刻、response hashがreceiptに残る。accepted検証も独立したラベル・receipt・response hashを持ち、primary1,040件と別のrequest/attempt数で集計される。archiveはsort順、gzip filename省略/mtime=0、tar mtime=0/mode0644で固定している。scriptは最終bridge/collision/archiveのhashをmanifestへ記録する。

最終取得が未完了なので、archive member数・最終hash・全件のsource-ID集合・accepted verification失敗件数・eligible concept数はまだ本レビューで承認していない。これらは取得担当者の完了時検証と最終統合検証で確定させること。

## 結論

**現実装はそのまま取得を継続できる。修正要求なし。** 最終readoutではprimary SUCCESS数とconservative eligible数を分け、accepted検証が必要な同物異名を未検証のまま利用しない。多対一conceptは元sourceリンクを保持しつつ、occurrence joinで蝶ごとにdeduplicateする。
