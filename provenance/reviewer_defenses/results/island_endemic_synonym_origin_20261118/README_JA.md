# 島嶼蝶類：学名統合後のHOSTS収録率と独立原産地分類の再現性監査（2026-11-18）

## WCVP v13で生態学的な一般化の最終ゲートを検査（ソース実行完了）

[オリジナルGEBと同じWCVP v13](https://github.com/matildabrown/rWCVPdata/tree/65bed76bae9d644ccb6ad200c05f9f5071d89e05) の植物 accepted species を、TDWG Level3 `CNY` (Canary Islands全域、個別島ではない) について native / introduced / extinct / doubtful 区分付きで取得。上の一律24蝶種の明記された**幼虫寄主52属（チョウ×属86件）**と交差した。

| 受理植物種を含む幼虫寄主属の分類 | 属数 |
|---|---:|
| 在来植物種が掲載され、導入植物種なし | 25 |
| **在来植物種と導入植物種の両方が掲載** | **18** |
| 導入植物種が掲載され、在来植物種なし | 6 |
| 該当するCNY受理植物種の記録なし | 3 |
| 合計 | **52** |

同52属内でWCVP CNYに収載される確認済み受理植物種は348件。52属中43属に在来植物種、24属に導入植物種の少なくとも1件あり、**18属には両方**がある。24蝶種のうち19種には少なくとも1つ導入植物種を含む幼虫寄主属がある。これを原産区分別に分けると、**源資料で幼虫寄主属が確認されたカナリア固有蝶6/6種すべて**に導入植物を含む属、かつ同じ属内に在来と導入の両方がある。在来非固有で幼虫寄主属確認18種中13種には導入植物を含む属、10種には在来・導入両方がある属。

**絶対に避ける誤解**：「その属に外来植物が存在する」≠「幼虫が実際にその外来植物種を食べる」。例えば `Rhamnus` 全体と実際に幼虫が利用する特定 `Rhamnus` 種を混同しない。WCVP `CNY` は全カナリア諸島の区域で、実際の島ごとの同所性を保証しない。庭で植栽されただけの植物はWCVPの「野生で導入」分布に必ずしも含まれない。これらを含め、結果は **寄主候補属に存在する植物相の記述**にとどまる。

再現コード：`scripts/audit_canary_host_genus_WCVP_v13_20261118.R`。出力は原種名・科・原産区分・属×CNY受理植物種の在来/外来件数・蝶ごとの属情報と、原産区分別記録を含む。 [WCVP原版照合 CI #44201479559](https://github.com/zuizui0223/chocho/actions/runs/44201479559) および[原産区分別出力CI #44201926869](https://github.com/zuizui0223/chocho/actions/runs/44201926869) 成功。後者のアーティファクト `sha256:a7b3b4086fb5a4be91ed6b2117cb86a17e81ce5ceacff3b23bcc44ced9eae622`。

全結果の不変な数値・出典台帳：`VERIFIED_GENERALIZATION_GATE_20261118.json`。

**投稿用の境界**：既存GEBの+54.9%/+118.9%は、既知寄主植物**種**の地理分布を固定資料から組み直した推定である。今回の独立島固有蝶の寄主は多くが**属**までしか特定できない。したがって、どちらかの大きさを島固有蝶の「資源拡大率」に代入することはできない。一方、外来植物候補の存在そのものが除外されないと測定できたため、次に必要なのは**植物種単位の幼虫利用関係と現在の個別島の分布・由来**を照合することである。


## 再現実行：一律文献ソースと35/42種の寄主属情報（2026-11-18）

**今までの二つの測定（DBへの収録と実際の寄主属の文献記録）は同じ意味ではない。** 第三の独立した記録ソースとして [Polinizadores de Canarias](https://polinizadoresdecanarias.org/repositorio/) の種別ページをINBOカナリア **全33種**に同じURI構築ルールで照会した。成虫の `PLANTAS VISITADAS` を使わず、幼虫についての `CICLO VITAL` 部分だけを採集し、文中の植物属名を86件（24種×属）、ユニーク52属、すべて原文に実在することをプログラムで検証した。種ページがない場合は **SOURCE_NOT_AVAILABLE** であり、寄主欠如とコードしない。

| INBOカナリア独立群 | 一律資料で幼虫寄主属を確認 | 元HOSTS厳密一致（別の資料） |
|---|---:|---:|
| 種レベル固有13種 | **6/13** | 0/13 |
| 在来非固有19種 | **18/19** | 13/19 |
| 導入推定1種 | **0/1**（資料ページなし） | 1/1 |
| 合計33種 | **24/33** | 14/33 |

**資料収録の比較だけ**：固有6/13対在来18/19、Fisher two-sided p=**0.0032221**、科内の固有/在来とページ数を固定する厳密分布に基づく片側p=**0.00410084**。だがカナリア島固有 `Hipparchia` 5種を除くと固有6/8 vs 在来18/19、Fisher two-sided p=**0.20103**（差は支持されない）。在来判定に独立の異論がある `Cacyreus marshalli` と `Pieris rapae` を除外しても、固有6/13 vs 在来16/17、p=**0.00936**。単一の近縁放散群の記録不足が大きく寄与し、科内pだけで普遍的な「固有性が原因」とは断定できない。

一律資料で新たにHOSTSの既知名や先行の固有10種の外部文献と重複せず記録できたのは **`Aricia cramera`、`Euchloe charlonia`、`Polyommatus celina`、`Vanessa vulcania`** の在来非固有4種。固定42種を分母として、
**元HOSTS19＋確認済み同物異名2＋以前の固有蝶外部文献10＋一律資料の追加4＝35/42種（83.3%）** に幼虫寄主属の文献記録がある（重複除去済み）。歴史的な絶滅分類群 `Pieris wollastoni` 1種を除く分母では **35/41種**。いずれも文献探索が全ソースで均等ではないため、**真の生物学的収録率の無偏推定ではない**。未確認の現生6種は、カナリアの `Euchloe grancanariensis` と近縁な `Hipparchia` 5種である。

**GEB種採用の別の学名誤差も補正**。独立INBOの `Colias crocea` と元GEB239の `Colias croceus` は、[NCBIのhomotypic synonym](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=72248)として受理同一性を独立検証。従来の厳密二名法交差 **12/42 → ソース検証済み最少13/42**。カナリア在来非固有19種のGEB採用は従来厳密9種→**10/19**。マデイラ在来非固有7種も **6/7→7/7**。 `Pararge xiphioides` と `P. aegeria` の分類学的衝突は除外したまま、島固有13種のGEB239採用は0/13で変わらない。まだ全239と全42の異名を手作業で完全照合した値ではない。

**ソース／コード／実行成功**
- 33種のHTTP200/404、原HTML SHA256、原文幼虫欄：`scripts/audit_uniform_canary_source_20261118.py`
- 原文に現れる86件の植物属名のみ受理・源資料反復抽出：`uniform_source_larval_genera_24.json`、`scripts/summarize_uniform_canary_hosts_20261118.py`
- 異なる科を固定した厳密条件付き検定と `Hipparchia`/原産性の感度検定：同スクリプト
- 文献重複除去と全42種源別フラグ：`scripts/audit_independent42_multisource_host_union_20261118.py`
- 239種の独立同物異名補正：`scripts/audit_independent42_original239_taxonomic_selection_20261118.py`
- [実行成功：一律24種/86属対応・分類名13/42・合計35/42](https://github.com/zuizui0223/chocho/actions/runs/44200903066)；[感度検定追加の再実行成功](https://github.com/zuizui0223/chocho/actions/runs/44201042984)

**残る同定限界**：サイトの24種はほとんど幼虫寄主を**属レベル**で記述している。属内に在来植物と外来植物が共存する場合、実際に利用する植物種がどちらかは同定できず、元GEBの+54.9%/+118.9%を固有蝶に移植する根拠にはならない。WCVP v13の原始CNY属内の在来/外来植物種の分布と照合する追加測定を別ワークフローで実行する（属内共存＝チョウの外来寄主利用とはしない）。


## 追加完了：カナリアの元INBO全33蝶種を独立原産地区分に照合

実際のINBOに出てくるカナリア蝶 **33種すべて**について、カナリア州政府BIOTAの種別掲載ページから、種レベルの固有性と現地原産区分を独立照合した。最初の20種以外の13種（`canary_BIOTA_remaining_13.csv`）を取得し、`13固有+19在来非固有+1導入推定=33`を検証。元の20種に固有13種・在来非固有7種が含まれており、その21/42統合結果は保持される。

| 原産区分（島別、BIOTA掲載） | INBO独立分母 | NHM HOSTS厳密二名法 | 分類名統合後HOSTS寄主属 | 元のGEB239種に採用 |
|---|---:|---:|---:|---:|
| **カナリア種レベル固有** | **13** | 0 | **0/13** | **0/13** |
| **カナリア在来・非固有** | **19** | 13 | **15/19 (78.9%)** | **9/19** |
| ├ 在来確実NS | 5 | 1 | 3/5 | 0/5 |
| └ 在来推定NP | 14 | 12 | 12/14 | 9/14 |
| カナリア導入推定IP | 1 | 1 | 1/1 | 0/1 |
| 合計 | **33** | **14** | **16/33** | **9/33** |

この分類は **現代BIOTAの地域別分類**に厳密準拠。例えばBIOTAは `Cacyreus marshalli` をNPに分類するが、同種の移入・侵入履歴に関する他の文献と矛盾しうる。これは別途感度評価し、単に生物学的な在来の真実として強制しない。BIOTAでは `Pieris rapae` がカナリアNPだが、マデイラBorges 2008では導入 `i` であり、原産性を地理間で勝手に統一しない。BIOTAの `Argynnis pandora` と `Maniola jurtina` の照合には亜種階級の掲載行を用いたため、種階級の区別を維持。GBIFの `Pararge xiphioides` 統合衝突は隔離。

**検証済み公開CI**： https://github.com/zuizui0223/chocho/actions/runs/44197938200 （全33個別名の分類・分母等の内部整合確認、成功）。属や科の構成による見かけの差を調べるため、各科×原産区分ごとの原表も `canary_full33_family_origin_host_coverage.csv` としてActions artifactへ追加している（ソースの種名・科・区分の組合せのみ、独立進化の反復数ではない）。

## 追加完了：HOSTS未収録の島固有蝶の寄主知識を独立文献から回収

**HOSTS正規化後21/42（50.0%） + 外部自然史資料で幼虫寄主属を確認した新規10種 = 重複なし31/42（73.8%）**。追加10種の内訳はマデイラ固有3/3、カナリア固有7/13。各種について出典URL・属名・一次観察か二次整理かの証拠種類を保存：`independent_larval_host_rescue_10.csv`。元42種との厳密照合、独立固有区分の確認、元HOSTSとの非重複を標準ライブラリの `scripts/audit_literature_host_rescue_20261118.py` で実行。**検証CI成功** https://github.com/zuizui0223/chocho/actions/runs/44197227988 。

重要：外部文献はHOSTSで欠落する固有種に**偏って意図的に**探索した。この31/42は資料上確認できた幼虫寄主**属名**の下限であり、左右の群を同じ努力で調査した真の生物学的収録率ではない。特に文献の記録の種類は現地幼虫観察、文献レビュー、自然史ウェブ資料など強さが一様でない。幼虫寄主属の知識が回復しても、対象島の在来植物の分布、植物が現在存在すること、幼虫成績、島嶼へのチョウ到達・定着は示せない。元GEBの+54.9%/+118.9%は保持し、独立固有蝶への一般化だけ未確定とする。


## 確認済み結果（2026-11-18、数値確定済みのデータベース監査）
**独立42種の厳密二名法一致19/42（45.2%）→受理名・別名の保守的統合後21/42（50.0%）**。追加でHOSTS寄主属報告を検出したのは `Argynnis pandora` → `Pandoriana pandora` と `Pontia daplidice` → `Pieris daplidice`。HOSTS原データの科一致を再点検して成功。5種には未解決の同物異名・分類階級問題があるため、21/42を最終の「真の野外寄主既知率」とは呼ばない。

| 独立分類の種群 | HOSTS元の二名法一致 | 保守的異名統合後 | 元の固定GEB239種採用 |
|---|---:|---:|---:|
| マデイラ在来・非固有 | 7/7 | 7/7 | 6/7 |
| マデイラ現生固有 | 0/3 | 0/3 | 0/3 |
| カナリアBIOTA固有13種（7属・4科） | **0/13** | **0/13**（未解決4） | **0/13** |
| カナリアBIOTA在来非固有7種 | **4/7** | **4/7**（未解決1） | **3/7** |
| アゾレス独立固有 `Hipparchia` | 2/2 | 2/2（分類階級未解決1） | 0/2 |

カナリアの元の「近縁 `Hipparchia` 5種 0/5」から、公式 BIOTA で種単位の固有性が確認された **13種・7属・4科の0/13** に検証範囲を拡張した。種別BIOTA ID、正式名、原産地区分、掲載URLは [canary_BIOTA_source_status_20.csv](canary_BIOTA_source_status_20.csv) に凍結。独立した在来非固有7種は `Nativo Seguro` 2種（1/2）と `Nativo Probable` 5種（3/5）に分けている。すべてのカナリア在来非固有蝶を列挙した分母ではない。

重要な分類学的衝突：GBIFでは `Pararge xiphioides` と `Pararge aegeria` が同じ受理種ID `8049830` に対応するが、BIOTAは前者を独立固有種とする。これを機械的に統合して、広域種の寄主報告を島固有種に流用する操作は**行っていない**。またGBIFで `H. bacchus`、`H. gomera`、`Thymelicus christi` が種同物異名と処理されていても、独立地域分類・分子系統との衝突があるため自動結合しない。

**再現証拠**：[名字・原産地の統合解析 CI](https://github.com/zuizui0223/chocho/actions/runs/44195494325) と [凍結原表の科一致も再監査した CI](https://github.com/zuizui0223/chocho/actions/runs/44195890870)。前者の原API応答ZIP digest: `sha256:0cf6a967cacf7b7de4e6577201e909a3815b1c1b1c7bd97a341458dedd807dc7`。後者のdigest: `sha256:3940eb4562e6d87b6353c0aebf9db74b9cbd667e6ed82ae1616894a5cff8c59e`。手元に永久保存した集計は [VERIFIED_SOURCE_OUTCOME_20261118.json](VERIFIED_SOURCE_OUTCOME_20261118.json) を参照。

**生態学的判定**：固有蝶への寄主資源拡大の外挿を阻害する収録偏りは「マデイラだけ・`Hipparchia`だけ」ではない。ただしアゾレスでは逆、カナリア在来群も未網羅で、固有性そのものの効果を断定できない。現時点で検証できたのは蝶種の寄主「報告の収録率」であり、実際の寄主植物有無・利用・蝶の定着確率ではない。


## 目的と明確な分母
2026-11-17時点の元の独立チェックリスト（INBO 18実島、42ユニーク蝶種、227種×島）から開始する。**HOSTSで見つかった蝶だけを分母にしない**。元のGEB解析で採用された239種とは別の選抜過程を測る。元の +54.9% 全体、+118.9% 島嶼の既知寄主機会の肯定結果と原稿は変更しない。

対象：上流の [独立分母監査 PR #176](https://github.com/zuizui0223/chocho/pull/176) と同じ凍結42種データ（Git blob `3c56054fa21e667eed422a67010059c8e198b1bc`）。NHM HOSTS は GloBI mirror `808e0b869f9ec1adf8efff87cf6a395adda103e0` の `resource.csv` 全140,485行を使用する。

- 第1指標：INBOで列挙されたユニーク42蝶種の **HOSTS寄主属報告**。元の二名法の完全一致は19/42 = 45.2%。
- 第2指標：同じHOSTS原表で**受理名・種階級の同物異名を、GBIF分類階級・科・受理speciesKeyの整合を確認して統合**した後に新たに見つかる寄主属報告。分類が不確定な名前は未照合として残す。
- 第3指標：独立42種から元のGEB固定239種へ入った割合12/42 = 28.6%。HOSTS寄主属があっても追加の7種が除外されており、理由は未同定。

生物学的な「真の寄主既知率」はHOSTSだけでは測定不能。和名・同種異名統合後に得るのは **分類名統合後のHOSTS報告収録率**。文献に寄主が存在してもHOSTSにはない場合、未収録のままである。

## 独立した在来・固有種の事前指定コホート

| 原産地区分・独立根拠 | 生物学的独立単位 | INBOとHOSTSの元二名法完全一致 | 元239種に採用 |
|---|---|---:|---:|
| マデイラ在来非固有7種、Borges et al. (2008) pp.342–343 | 本島の7種 | 7/7 | 6/7 |
| マデイラの現生固有3種、同書 | `Gonepteryx maderensis`, `Hipparchia maderensis`, `Pararge xiphia` | 0/3 | 0/3 |
| カナリアの島固有 `Hipparchia` 5種、Spilani et al. (2025) | 1つの近縁放散群、5島固有種 | 0/5 | 0/5 |
| アゾレスの固有 `Hipparchia` 2種、アゾレス州SIARAM固有分類 | `H. azorina`, `H. miguelensis`（階級・受理名監査対象） | 2/2 | 0/2 |

カナリア原著：Spilani et al. 2025, *Insect Systematics and Diversity* 9(5), ixaf043, https://doi.org/10.1093/isd/ixaf043 。5島の `H. wyssii`, `H. tamadabae`, `H. gomera`, `H. bacchus`, `H. tilosi` は異なる単一島の固有種として支持されるが、同じ放散群であり「独立の5反復」ではない。原著も野外幼虫寄主が十分に報告されていないことを記す。

アゾレス州の独立固有種リスト：https://siaram.azores.gov.pt/fauna/lista/artropodes.html （`H. azorina` と `H. miguelensis` に E 区分）。アゾレス大学のより包括的な分類・島別データ：https://www.gbif.org/dataset/2d91cfd8-0a48-4d80-8128-080e52a1e650 （GBIF DOI 10.15468/vzpxhx）。`H. azorina` 種群の種・亜種階級の不一致はさらに監査が必要。

マデイラ史資料：Borges et al. (eds.) (2008), *A list of the terrestrial fungi, flora and fauna of Madeira and Selvagens archipelagos*, pp.342–343, https://ifcn.madeira.gov.pt/images/Doc_Artigos/Divulgacao/publicacoes/livros/ListfunffterrMadSel.pdf 。2008年当時の地位であり現在の繁殖地位ではない。絶滅した `Pieris wollastoni` と固有亜種のみの `Lycaena phlaeas phlaeoides` は「現生固有種」コホートに混ぜない。

## 自動化された名前統合の境界
`scripts/audit_island_endemic_HOSTS_synonyms_20261118.py` は42種それぞれのGBIF *species/match*、受理種実体、受理種の *synonyms* を取得し、返答をタイムスタンプ付きで保存する。GBIFの分類名 API は https://techdocs.gbif.org/en/openapi/v1/species を参照。

同じ受理speciesKeyを異なるチェックリスト種に割り当てた場合は**衝突として隔離**し、属違いなどの受理名を勝手に固有種へ合算しない。科違い、曖昧一致、上位階級への一致、亜種のみの異名、不完全API応答は別に出力し、成功した別名だけ加算する。元の完全一致19種を削除しない。結果は `results/independent_island_endemic_alias_20261118/` に出し Actions artifact に原API応答を含める。

**現時点での代表的な注意例**：INBO `Pontia daplidice`、HOSTS `Pieris daplidice` はGBIFで同物異名の対応が確認できる（https://www.gbif.org/taxon/6W83L）。しかし対応名のHOSTS行に**有効な寄主属報告があるかは別途原表で検証する**。同様に `Cyclyrius webbianus` と `Leptotes webbianus` は候補であって、受理名同一性と記録行を調べてから加算する。

## 科学的判定
1. 元の単一島マデイラ7/7対0/3から、カナリアの別系統 `Hipparchia` 0/5まで「固有種の資料収録に偏りがある」観測を拡張。ただしカナリアは1属1放散群で、系統と研究史が完全に交絡。
2. アゾレスの固有 `Hipparchia` 2/2は重要な反例。3諸島を一つの普遍的固有種効果とせず、地域×固有性（加えて属・分類群）を区別する。
3. 同じ偏りが学名・異名を統合しても残るかを、このブランチの主要な反証的ゲートにする。**未照合を「真の在来寄主がゼロ」と代入しない**。
4. 広範囲の在来・固有種を推論する最終ゲートは、アゾレス大学の種・亜種階級区分を全チョウへ照合し、カナリアBIOTAの全蝶分類・在来性を取得し、同じ分母で比較すること。現在の5種選定を全カナリア蝶相と呼ばない。
5. 収録バイアスが残るなら、元239種の+54.9%/+118.9%は**既知HOSTS寄主を持つ選抜標本での潜在的機会**として保持し、島固有種全般の資源拡大率への一般化を保留する。

## 実行
```bash
python3 scripts/audit_island_endemic_HOSTS_synonyms_20261118.py \
  --hosts _hosts/resource.csv \
  --out-dir results/independent_island_endemic_alias_20261118
```
同じ結果の固定再生は出力した `gbif_source_receipts/` を元の出力ディレクトリに復元して `--offline` を指定する。HOSTS commit、元42 CSV blob、全GBIF応答を元の再現ソースとして記録する。原稿のマージ・投稿は行わない。
