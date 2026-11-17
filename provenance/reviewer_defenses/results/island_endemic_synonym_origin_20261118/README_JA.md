# 島嶼蝶類：学名統合後のHOSTS収録率と独立原産地分類の再現性監査（2026-11-18）

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
