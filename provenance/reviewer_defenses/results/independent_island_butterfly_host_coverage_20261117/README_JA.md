# 独立島嶼蝶類チェックリスト × HOSTS 収録監査 — 2026-11-17

## 問い：HOSTSに収録されていない種を分母から落とさない

島嶼蝶類の母集団を**独立した在島チェックリスト**で定義し、既存の固定HOSTSデータに幼虫寄主植物の報告がどれだけあるか計測する。これは**データベース収録率＋学名照合率**の検証であり、未照合は野外の寄主不存在・蝶の不在・幼虫の生存不能を意味しない。

HOSTSから事前選抜した既存239種の全世界潜在的寄主機会（+54.9%、島での相対+118.9%）とは異なる分母である。元の肯定的結果を削除せず**一般化可能性だけを監査**する。

## 原データと固定ID

- INBO [National checklists and red lists for European butterflies](https://github.com/inbo/eurobutt-checklist)、[GBIF DOI 10.15468/ye7whj](https://doi.org/10.15468/ye7whj)：
  - data/processed/taxon.csv（496分類群、Git blob 67f3126314c7facd5b133ca06f18f397b8538549）。
  - data/processed/distribution.csv（8,020分布行、Git blob 1e57f4cbf2ec78e84a7ba2d21b112eda49396130）。
- [GloBI/NHM HOSTS mirror](https://github.com/globalbioticinteractions/HOSTS)：固定commit 808e0b869f9ec1adf8efff87cf6a395adda103e0、resource.csv Git blob 06749c04bd5b7a754d6f7717a31c01c51357cae8（140,485報告行）。
- **独立した在来性・固有性**：Borges et al. (eds., 2008), *A list of the terrestrial fungi, flora and fauna of Madeira and Selvagens archipelagos*, [公開目録](https://ifcn.madeira.gov.pt/images/Doc_Artigos/Divulgacao/publicacoes/livros/ListfunffterrMadSel.pdf), 342–343頁。原資料の n（在来）、END（固有）、MAC（マカロネシア）、i（導入）、m（渡り）、†（絶滅）の区分をマデイラ諸島の該当分類群のみに適用。

### 分母・照合規則

1. 分母：INBO locality に個別の島名があるアゾレス、カナリア、マデイラの18島。諸島の一括集計行を除き、occurrenceStatus が present の蝶×島を採用し同じ島同じ種は重複除去。irregular 単独は不採用。*Issoria lathonia*（Madeira）にはpresent/irregular同時登録があり、競合を記録する。
2. 照合：蝶の二名法を正規化し、HOSTSの Insect Genus＋Insect Species に完全一致した行のみカウント。植物の属名および種小名が実際に記録された割合を別々に集計。亜種・同物異名・属替えは推測で結合しない。
3. 在来性：Borges et al. 2008でマデイラについて独立確認された状態を使う。他17島に在来・固有のラベルを推測で補わない。原資料は2008年の歴史的目録であり、現代の繁殖確認ではない。

## 得られたデータ原表集計

| 分母 | 総数 | HOSTS二名法一致＋寄主属 | HOSTS寄主種まで |
|---|---:|---:|---:|
| 蝶種×個別島 | **227** | **168 (74.0%)** | **162 (71.4%)** |
| 独立チェックリストのユニーク種 | **42** | **19 (45.2%)** | **18 (42.9%)** |

**23/42種**に、厳密な二名法一致でHOSTS寄主属が見つからなかった。種×島組の74.0%よりユニーク種の45.2%が低いのは、複数島に出る広域種を種×島で重複カウントするため。独立母集団の収録率は後者を主要指標とする。

| 諸島 | 個別島 | 種×島分母 | 寄主属完全一致 | 寄主種完全一致 |
|---|---:|---:|---:|---:|
| アゾレス | 9 | 61 | 61 | 61 |
| カナリア | 7 | 141 | 88 | 82 |
| マデイラ | 2 | 25 | 19 | 19 |

### マデイラ本島：独立在来・固有分類との照合

| 2008年目録による分類 | INBOのpresent に含まれる分類群 | HOSTS寄主属完全一致 |
|---|---:|---:|
| **在来・非固有（n）** | **7** | **7/7** |
| **現生固有種（END）** | **3** | **0/3** |
| マカロネシア固有域（MAC） | 1 | 0/1 |
| 固有亜種がある広域種（END at subspecies） | 1 | 1/1（ただし亜種固有の寄主証拠ではない） |
| 歴史的に絶滅した固有分類群（END†） | 1 | 0/1 |
| 導入（i） | 1 | 1/1 |
| 渡り（m） | 1 | 1/1 |

固有3種は *Hipparchia maderensis*、*Pararge xiphia*、*Gonepteryx maderensis*。同書本文には *H. maderensis* 幼虫の複数イネ科植物の利用と *G. maderensis* のMacaronesian Buckthorn利用が明記されている。**少なくとも2種は実際に食草に関する文献知識が存在するにもかかわらず、当該学名ではHOSTSに照合できていない**。これは寄主の生物学的不存在ではなく、データベースへの収録・名前照合過程の欠落を示唆する。

- *Pieris wollastoni* は旧目録では **Pieris brassicae wollastoni END†**。2008年の絶滅分類群を現在の固有現生種に加えない。母種 *P. brassicae* のHOSTS寄主を自動代入しない。
- *Lycaena phlaeas* の目録ENDは**固有亜種 phlaeoides**に付く。種全体へのHOSTS寄主報告は存在するが亜種固有の寄主利用は別途要調査。
- *Vanessa vulcania* は旧目録では **MAC**。島内限定固有種と扱わない。
- *Issoria lathonia* は旧目録 **m**。INBOではMadeiraのpresent/irregular両方に収載されているため、在来/固有の対比へ混ぜない。

### 文字列一致の保守性

HOSTSには *Leptotes webbianus* と *Pieris daplidice* の文字列があるのに対し、INBOは *Cyclyrius webbianus* と *Pontia daplidice* を採用している。これは**分類名の同物異名候補**として記録し、受理名・原記載と分類群同一性が確認されるまでHOSTSカバーに自動加算しない。これ以外にも類似の問題があり、45.2%は**受理名を統合した最終収録率ではない**。

## 科学的帰結・保留条件

- HOSTS掲載種のみで復元した島の+118.9%という大きな相対資源増は、**観測された寄主の分母に条件づけられた計算**である。この肯定結果を壊すわけではないが、固有蝶類を含む全島蝶相が同じだけ資源拡大した証拠ではない。
- 在来寄主が見かけ上少ない地域で、未収録の在来蝶・固有蝶の寄主情報を0扱いすると、外来寄主の「穴埋め」という生態的主張が強く見える可能性がある。**HOSTS掲載選択と寄主不足を別々の測定過程としてモデル化**する必要がある。
- マデイラ7/7対0/3は一島の小さな歴史的例であり、全球固有蝶のHOSTS収録率と呼べない。島間生物地理・調査史・分類名を統制した多島チェックリストが次の検証になる。
- 従来の港・植物相、155種のGBIF既知種予測、海の回路の肯定的結果は**元の別々のendpointを保って維持**する。今回の収録監査は港からの蝶の移入・野外定着の実証でもない。

## 再現ファイル

- audit_independent_checklist_hosts.py：Python標準ライブラリだけを使用、**HOSTSに存在しない蝶も分母に残す**。入力は上述の固定Git blobと一致するCSV。個別島の種×島表と集計JSONを再生成する。
- species_island_coverage.csv：227の全種×島ペア。学名、チェックリストの状態、HOSTS完全一致、植物の属/種、未検証別名、Madeiraの独立原産地・固有区分。
- summary.json：出典Git blob、明示的な分母、諸島別・原産地別の結果、および注意事項。

実行例（原資料CSVを上記固定版に揃えた場合）：

    python audit_independent_checklist_hosts.py --taxon taxon.csv --distribution distribution.csv --hosts resource.csv --out-dir replay

本監査では原稿の主結果を変更しない。先に分類名監査と独立チェックリストの多地域展開を進める。
