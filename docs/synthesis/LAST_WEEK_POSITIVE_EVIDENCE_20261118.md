# この1週間の肯定結果・反証的再検証と一論文への統合（2026-11-18）
**対象期間:** 2026-11-11〜18（JST）の追加・更新研究。PR全体の検索では「2026-11-11 00:00 UTC以降に更新」の **21 PR** を直接照合。従来の +54.9% 全球、+118.9% 島嶼、+4.12% 予測は**既存の土台**であり、今週の新発見として数えない。  
**Source:** [全PR索引](ALL_PR_AUDIT_20261118.md)、[長期証拠台帳 R01–R22](POSITIVE_EVIDENCE_AND_ONE_STORY_20261118.md)。  
**Freeze type:** 解析の結果と後続反証を `main` に**非破壊的に保存**。各Draft PRのコードや未提出原稿の強制マージ・再解析・査読承認は行わない。

## まず結論（今週、新たに言えるようになったこと）

1. **「外来植物が潜在的な寄主資源を増やす」は野外の生態学的な使用例に到達した。** カナリア固有蝶 `Pieris cheiranthi` の幼虫がカナリア外来植物 `Tropaeolum majus` 上で報告された（テネリフェとラパルマ、後者庭園）。植物種・地域原産区分・昆虫幼虫を別々に照合。**同一島での新規定着や個体数増加までは証明していない**。[#177](https://github.com/zuizui0223/chocho/pull/177) [source replay](https://github.com/zuizui0223/chocho/actions/runs/44207292437)
2. **島の固有蝶にも拡大機会はあるが、均一ではない。** 寄主植物を種まで確定した固有蝶8種（独立現生固有種18の恣意的な部分集団）を元WCVPで再計算：4種で全球**潜在植物地域**増加（Pieris 2→80, +78; `Euchloe grancanariensis` 30→77, +47; `E. hesperidum` 60→111, +51—飼育1幼虫の摂食のみ; `Pararge xiphioides` 92→101, +9—文献のみ）。残り4種で追加0（**現在確認した不完全寄主リストに限定**）。8種すべてに出発諸島地域の在来寄主種があり、この8種には新たに寄主が存在するようになった出発諸島地域は0。直接幼虫が名付けた植物上で記録された個別島の組は22組中**4組**。[#177](https://github.com/zuizui0223/chocho/pull/177) [source replay](https://github.com/zuizui0223/chocho/actions/runs/44223652218)
3. **島固有蝶18種への一般化は元239種からはできない、という正確な理由が確認された。** 独立チェックリスト42種から歴史的絶滅1種を除く現生41種＝現生固有18＋非固有23。元339候補および239採用で固有種 **0/18**、239の非固有は**13/23**。外部資料で少なくとも寄主属が分かる現生種 **36/41**。これが証明するのは「HOSTSが知らない＝植物がない」と誤認してはならないこと。36/41は複数資料の非均等な追跡から得た**文献確認下限**であり、無偏な寄主知識の母集団収録率ではない。未知の寄主を追加すると導入由来増分は数学的に増減可能で、4/8を18種へそのまま外挿できない。[#176](https://github.com/zuizui0223/chocho/pull/176) [#177](https://github.com/zuizui0223/chocho/pull/177) [verified CI](https://github.com/zuizui0223/chocho/actions/runs/44234108322)
4. **島で外来植物を寄主とする文献報告が相対的に多いことを、全蝶群と同種比較で具体的に測った。** 原HOSTS地域資料では島13地域の外来寄主報告割合中央値42.4%、本土15地域9.1%。作物除外でも**20.9% vs 8.1%**。ただし厳密に同じ蝶67種を比較すると島−本土の生データ差 **+11.32 pp**（CI [-.89,+23.53]）、植物科別の利用可能な外来植物割合からの期待 **+15.19 pp**、期待を超える残差 **−3.87 pp**（CI [-15.85,+8.13]）。**島では外来寄主利用の報告が多いが、植物相を超える選好・幼虫適応の証明ではない**。同一蝶×同一植物の場所違い確認はFAO作物種除去後11組（ただし本土比較10/11がタイに集中）。[#162](https://github.com/zuizui0223/chocho/pull/162) [#163](https://github.com/zuizui0223/chocho/pull/163)
5. **寄主と海洋は分布記録に異なる追加情報を持つが、相乗作用を強弁できない。** 同じ51島の後期群集予測で、海洋AIS/風情報を寄主に追加 **+6.488% MSE削減**（[+3.338,+9.609]）、海洋情報に寄主を追加 **+1.731%**（[-1.305,+4.440]）、host×marine交互作用 **+0.343%**（[-2.798,+3.219]）。肯定的なのは**条件付き海洋予測情報**であり、蝶定着に対する港・航路と寄主の相乗的因果効果ではない。[#160](https://github.com/zuizui0223/chocho/pull/160)

## 一週間の追加正・追加制限の詳細

| 追加日 | 新しく確認・精査したソース結果 | 分類 | 正しい論文内での位置 |
|---|---|---|---|
| 11/14–15 | **2020・2022年 Santa Catalina の `Pieris rapae` 幼虫**を別資料・座標で再確認（2記録、近接31.56m）。**1927年から島にいた歴史資料を発見**、後期GBIF初記録は新規移入でない。 | 堅固な生物記録＋必須反証 | 幼虫発見を新規定着に言い換えない [#159](https://github.com/zuizui0223/chocho/pull/159) |
| 11/14–16 | 全**225種×51島＝11,475**の歴史GBIFを同じ定義で再照会。2001–2010不記録後の「初記録」317例の**145（45.7%）は2001年以前に記録済み**。 | 強い史料監査結果 | 最重要の因果解釈ガード、港→新定着の誇張を防ぐ [#159](https://github.com/zuizui0223/chocho/pull/159) |
| 11/16 | 2001–2010の港付近に**在来のみ**の適合寄主が記録された場合、2011–2022にチョウを再記録する確率の調整差 **+11.24 pp**（概算CI [+3.00,+19.48]）。過去1700–2000でも在来寄主と記録の正の関連 **+17.11 pp**がすでに存在。外来のみでは正確な独立利益を示さず。 | 肯定的な在来寄主-記録地理（因果なし） | 既存資源の役割も明示し、移入固有の定着効果と混同しない [#161](https://github.com/zuizui0223/chocho/pull/161) |
| 11/16 | 文献HOSTSで島の導入寄主報告高比率と、FAO除外後の正方向。後続で**植物科×全植物相**に合わせると強い残差なし。 | 報告上の正＋選好仮説反証 | 「供給された植物プールから実際に使う」生態学へ [#162](https://github.com/zuizui0223/chocho/pull/162) [#163](https://github.com/zuizui0223/chocho/pull/163) |
| 11/16 | 別の欧州・大西洋22島サブセットでは**在来植物に海洋散布形態が+2.95 pp多い**（17/22地域で高い、4低い、1同点）。 | 測定された形態構成の正差 | 海での天然種子移動を証明することではない [#164](https://github.com/zuizui0223/chocho/pull/164) |
| 11/16–17 | NOAA旧OSCARの**実際の海流**データ：AZO→CNY沿経路平均 +0.01292 m/s、18/18年の原流速サンプルが充足（13/18年の年平均が正）。CNY→CVIは長期14/18年で事前ゲート未達、AZO→CVIは不安定。全球モデル280本、別地域接近12方向だが**良い対照が2/12のみ、2/2背景期待どおり**。 | 実測海流の部分肯定、種子伝播は未同定 | マクロ物理学的背景として補足 [#166](https://github.com/zuizui0223/chocho/pull/166) [#169](https://github.com/zuizui0223/chocho/pull/169) [#175](https://github.com/zuizui0223/chocho/pull/175) |
| 11/17 | WCVP全植物相の**126島で在来／外来の組成を同じ島ペアで比較**。在来植物相Jaccard中央値 .004487、外来 .046030（7,875ペア）。洋上最短経路を加えると在来 **+2.234%**、外来 **+1.775%** の誤差減少。港の現代数proxy単独はほぼ0。 | 新しい同じ母集団の対照 | 「自然は海／外来は船」の単純二分を退ける [#175](https://github.com/zuizui0223/chocho/pull/175) |
| 11/17 | 港×小さな島（β−.1394）／在来貧弱な島（β−.1815）は単独なら正の港-外来相関が強い方向。ただし面積と在来植物数 r=.817、同時投入・非線形化・未学習生物区で**安定しない**。 | 探索的な正の空間異質性＋強い外的検証制約 | 補足に保存し「小島の港は独立して強い」と見出しにしない [#173](https://github.com/zuizui0223/chocho/pull/173) [#174](https://github.com/zuizui0223/chocho/pull/174) |
| 11/17 | 港近傍のCIRADモデル農業N勾配は外来／在来寄主候補の**記録組成比**と正（β+.08617、CI [.02375,.14859]、1,013港/40か国）。ただし厳しいラスタ被覆では弱まり、絶対外来種数や土壌への実測N効果を示さない。 | 探索的な環境相関 | 港周辺の人為土地利用の文脈に限る [#122](https://github.com/zuizui0223/chocho/pull/122) |
| 11/17–18 | 独立蝶類原産区分と正式受理名を照合：現生固有18種が元239に**0/18**、非固有23種で**13/23**。外部資料で少なくとも現生**36/41**種の寄主属を確認。 | 重要な正の情報回収＋母集団欠陥同定 | 一般化の限界を強い論点にして記載 [#176](https://github.com/zuizui0223/chocho/pull/176) [#177](https://github.com/zuizui0223/chocho/pull/177) |
| 11/18 | 固有 `Pieris cheiranthi` の移入 `Tropaeolum majus` 幼虫利用2島を出典照合。他方 `Euchloe grancanariensis` の在来 `Hirschfeldia incana` 上で**野外15幼虫、うち5蛹化**、羽化未確認。 | 最も生態学的に直接的な今週の肯定例 | 「利用」は実証、「局所個体群増大」は未実証 [#177](https://github.com/zuizui0223/chocho/pull/177) |
| 11/18 | 種階級まで確認した固有蝶8種の全球WCVP再構築は**4種でモデル上の追加地域**、4種は既知寄主リスト内で0。元18種には全く代表性がない。未知食草を加えると潜在的「追加」地域は減る場合もある。 | 不均一な潜在寄主地理＋非同定の解明 | 強い外的妥当性チェックとしてSI [#177](https://github.com/zuizui0223/chocho/pull/177) |

## 直近21 PRを漏らさず最終判断

| PR | 今週の役割 | 論文採用 |
|---|---|---|
| [#122](https://github.com/zuizui0223/chocho/pull/122) | 港周辺のモデルN負荷と外来・在来候補植物の相対記録比 | 探索的補足 |
| [#157](https://github.com/zuizui0223/chocho/pull/157) | 統合単一GEB論文、最新の航路と記録年の対照、Bornholm例 | 本文の候補。ただし年次港×初記録係数の解釈は本文・元データの整合確認が必要 |
| [#159](https://github.com/zuizui0223/chocho/pull/159) | 歴史的初記録の再分類、実幼虫写真、寄主×海洋交互作用 | 初記録反証を必須採用、相互作用は主要結果にしない |
| [#160](https://github.com/zuizui0223/chocho/pull/160) | 海洋と寄主の条件付き予測利得 | 海洋の残差予測情報を主結果、host×seaは未確定 |
| [#161](https://github.com/zuizui0223/chocho/pull/161) | 在来寄主と再記録・古い歴史の正の関連 | 生物学的文脈として採用、導入による新定着の推定でない |
| [#162](https://github.com/zuizui0223/chocho/pull/162) | 地域つきHOSTS寄主報告の島vs本土の正の差 | 本文小節。ただし独立の幼虫野外記録ではない |
| [#163](https://github.com/zuizui0223/chocho/pull/163) | 植物相と科で標準化して上記の過剰選好を限定 | #162と**対で**必須採用 |
| [#164](https://github.com/zuizui0223/chocho/pull/164) | 22島の海洋散布形質差 | 植物相供給の補足 |
| [#165](https://github.com/zuizui0223/chocho/pull/165) | NOAA原データ確認と島位置同定 | データ成立性、機構の証明ではない |
| [#166](https://github.com/zuizui0223/chocho/pull/166) | 歴史海流の方向の部分正 | 物理層の探索的補足 |
| [#167](https://github.com/zuizui0223/chocho/pull/167) | 2001–2010季節海流の継続性検査 | 符号と欠測含む補足 |
| [#168](https://github.com/zuizui0223/chocho/pull/168) | EuDiS植物源地域と受け手植物相の照合 | モデル出典との照合 |
| [#169](https://github.com/zuizui0223/chocho/pull/169) | 植物科別の海洋散布形質帰無 | 選択的伝播を主張しない |
| [#170](https://github.com/zuizui0223/chocho/pull/170) | NOAAで実軌跡計算・原港再接近 | 物理計算のみ |
| [#171](https://github.com/zuizui0223/chocho/pull/171) | GSHHG到達点の真正島同定、他諸島到着0 | 偽陽性を防ぐ重要なQA |
| [#172](https://github.com/zuizui0223/chocho/pull/172) | 同一24港で植物・2020船舶・蝶調査努力を統合、十分な同時比較は僅か3組/3島 | 輸送因果機構をまだ実行しない根拠 |
| [#173](https://github.com/zuizui0223/chocho/pull/173) | 小島／貧在来植物島の港-外来種相関異質性 | 補足のみ |
| [#174](https://github.com/zuizui0223/chocho/pull/174) | 非線形・歴史・地理外予測で異質性を反証 | #173と必ずセットで報告 |
| [#175](https://github.com/zuizui0223/chocho/pull/175) | 同じ島ペアの在来/外来植物相×海地理・港比較 | 植物相地理の独立サポート |
| [#176](https://github.com/zuizui0223/chocho/pull/176) | 独立蝶類分母・固有分類の収録差 | 一般化限界の必須SI |
| [#177](https://github.com/zuizui0223/chocho/pull/177) | 名前統合、実幼虫寄主、固有蝶8/18再構築、非同定 | **今週最重要の生物学的発見**、本文実例＋SI |

**監査要注意:** #157のPR説明には、年次調整係数 +0.354 pp の信頼区間 [+0.031,+0.676] と「すべての年次区間がゼロを含む」という文章が併存する。これらは同じ比較を指すなら矛盾するため、**原結果ファイルで整合性を確認するまで肯定的な独立した新規定着効果として主張しない**。

## どう一本のストーリーに接続するか — 単位をずらさずに「必要条件の連鎖」として扱う

**Ecological claim:** “Human-assembled plant communities expand **potential** larval resources unevenly, and islands magnify that opportunity because of a sparse documented native-host baseline; local butterfly use and recorded geography are then conditional on host identity, marine access, and independent observation/establishment filters.”

**論理橋1：供給の非一様性。** 港の多い島ほど広い外来植物相を記録するという主結果。これは寄主を船が運んだ因果証拠ではない。港に寄主植物だけが特異的に集まる説よりも、**全植物相の集積から既知の寄主植物が供給される**というモデルが現データに近い（#163, #173–175）。

**論理橋2：初期条件による増幅。** 在来寄主カバーが島16.0%・本土37.5%、残り空白の充足率が島22.7%・本土24.5%。**類似の空白充足でも相対利益が島で大きい**ため、港で植物が単純に多い、だけではない効果を得る（従来R02）。

**論理橋3：潜在的なホストと実際の摂食。** 固有 `P. cheiranthi` 幼虫の移入 `T. majus` 利用が**存在**する一方、`E. grancanariensis` は在来 `H. incana` で幼虫が見つかる。原239の効果を未知の18固有種へは外挿不可。**植物相→実際の摂食は一部で成立しうるが自動的ではない**。

**論理橋4：分布にはほかの条件も要る。** 全球155種の寄主地理が蝶記録予測を改善する（従来R11）。同じ51島の後期群集では海洋情報が寄主に独立して予測利得を足す（今週R10）。しかし、実際の湾港輸送、個体群の成立、全球実現群集の大きな転換は未同定。**「資源を利用できる」から「その島に新たに定着した」への飛躍をしない**。

**総合的な一文（査読耐性）:** *Human-redistributed floras create large but heterogeneous opportunities for familiar butterfly larval hosts, amplified on islands by low native-resource coverage. Port-linked floristic assembly and marine geography organize these opportunities, and documented endemic larval use shows that at least some are biologically usable, whereas independent butterfly records demonstrate strong additional spatial and ascertainment filters.*

## 推奨する5図
1. **島で増幅する算術とその生態学的意味**（16% vs 37.5%、gap-fill 22.7% vs 24.5%、+118.9% vs +40.9%）。
2. **島の外来植物相・港・供給**（全外来植物と候補食草の関連、島の地理交絡を記載）。
3. **既知食性のもとで誰がどれだけ潜在機会を得るか**（全球239、寄主植物寄与集中、null制御済み潜在共通化）。
4. **実際の幼虫寄主と外的妥当性**（`P. cheiranthi` vs `E. grancanariensis`、独立現生41/固有18への収録欠落・8再構築の母集団明示）。
5. **チョウ記録と海の制約**（155種予測、51島海洋予測、globals群集比較の小さな差、歴史GBIF first-index反証）。

## 原稿に昇格しない主張
港が実際にその植物種を輸送した／海流がその植物を他島へ漂着させた／小さい島の港が独立因果要因である／チョウが新しい島に定着した／外来食草が幼虫の適応度や個体数を増加させた／全球のチョウ群集が大規模な生物学的均質化を起こした、という主張はこのソース集合では未証明。

**Source notes:** Data and citations all originate from the listed immutable PR branch datasets and succeeded CI records; no full re-execution of all 176 PRs was performed in this documentation-only freeze. All 21 recent PRs remain separate Draft research branches. The original `main` GEB v0.2 manuscript is unchanged.
