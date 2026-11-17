# 独立した島の歴史資料でGBIF初索引を検証（2026-11-17）

**状態：実行完了・ソース固定CI成功。** これは統合GEB論文[PR #157](https://github.com/zuizui0223/chocho/pull/157)の一部。元の225種×51陸塊（11,475組）・2011–2022初索引317組は変更せず、独立現地資料が確実に先行する行だけを補助的に照合した。**真の移入・定着、チョウ在来性、島固有性は判定しない。**

## 1. フーコック島12組がすべてGBIFより前に報告済み

固定した[元GBIF種×陸塊×年窓表](../../data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv)の陸塊**471**はWorld Port Indexの**DUONG DONG（VN）**を代表港とし、同名港はフーコック島にある（[フーコック市観光案内局の港一覧](https://dulich.phuquoc.gov.vn/en/detailutilities/?child=1137&id=1130)）。地域一致は国名一致ではなく島・港の一致。

2001–2010と2011–2017には、以下12種×この島のGBIF索引がそれぞれ0件。**2018–2022年に12種すべてが初索引**された。一方で次の独立した調査資料では**12/12すべてが2007年までにフーコック島に報告済み**で、同一の学名に完全一致した。

| 最も古い独立報告の上限年 | 種 | 根拠 |
|---|---|---|
| 2006 | Melanitis leda; Ypthima baldus; Delias hyparete | 2006年、Wildlife At Risk、フーコック島73種現地調査リスト、印刷頁4 |
| 2006 | Papilio demoleus; Papilio polytes; Junonia almana | 同上 |
| 2006 | Junonia lemonias; Eurema hecabe; Leptosia nina | 同上 |
| 2007 | Euchrysops cnejus | 2007年図鑑印刷頁131、27 Dec 2007、An Thoi の資料 |
| 2007 | Danaus genutia | 2007年図鑑印刷頁48、20 Mar 2007、Eo Xoai/Nui Tuong 現地写真 |
| 2007 | Zizina otis | 2007年図鑑印刷頁129、23 Apr 2007、Bai Thom 現地写真 |

**独立文献（ユーザーの原データやGBIFの後期索引とは独立）：**
- Nguyen Vu Khoi & Bui Huu Manh / Wildlife At Risk (2006), *Butterflies of Phu Quoc National Park, Kien Giang, Vietnam*. Copyright/print license 2006、8日間の現地調査、73種のチェックリスト。原本：https://wildlifeatrisk.org/wp-content/uploads/2019/01/books.P6.-Initial-photographic-checklist-of-Butterfly-of-PQ.pdf
- Bui Huu Manh / Wildlife At Risk (2007), *Các loài Bướm ngày ở Phú Quốc*. 2007年12月印刷、約30日間の2007年調査の図鑑と種別現地写真。原本：https://wildlifeatrisk.org/wp-content/uploads/2019/01/books.P9.-PQ-Butterfly-book.pdf

これら2冊は**同じ保全団体の継続現地調査**であり、統計上独立な2調査や12独立島ではない。図鑑の写真、本文や地図は著作権保護対象なので、本リポでは**必要な学名・記録年・頁・出典URLなどの照合事実だけ**を保存し、画像やPDF本体は複製しない。

厳密な結論は**この12種の2018–2022 GBIF初索引は、種ごとの「フーコック島への初到達日」ではありえない**、ということ。3種が2006年チェックリストに載らないのは2006年の不在を意味しない。12種が2007年にすでにいたことも、その時点以前の自然分布、いつから在来か、野外繁殖・存続、または港からの侵入を証明しない。2018年以前のどの時点で初到来したかは不明。

## 2. 原51陸塊の歴史資料追加感度分析（観測指数の再定義に注意）

全51陸塊・1,275島間ペアを固定し、2018–2022初索引12セルのうち**2006–2007の独立記録で確認済みの12だけ**を前期の種リストへ加えた比較。これは異なる観測過程の資料を合わせた**下限的な文献補強ベースライン**であり、GBIF観測努力を補正した純粋な生態学的群集推定量ではない。

| 種組成の定義 | 2001–2010種×島索引 | 平均ペアJaccard | 2022年末までとの差 |
|---|---:|---:|---:|
| 元GBIF索引のみ | 736 | **0.120963** | **+0.056810** |
| 2006/2007独立文献を前期に足す | 748 | **0.124359** | **+0.053414** |
| 2022年末の同じGBIF索引 | 1,053 | **0.177773** | — |

唯一の事後注釈による全島ペアΔの変化は −0.003396。これは「真の変化量が減った」ではなく、**元のGBIF時系列が古い野外記録をすくえていなかったことが、数量的に示された**という意味である。残る305組は未精査であり、全球の初索引誤認率を12/317から推定しない。

## 3. 在来・固有性が元の225種に表現される保証はない

[Australian Government (Director of National Parks), Christmas Island Biodiversity Monitoring Program, report (2008), p.100–101](https://www.aph.gov.au/binaries/senate/committee/eca_ctte/estimates/add_0809/ewha/christmas_island_report_sept08.pdf) は次の2つの**異なる**情報を持つ。

- *Melanitis leda* ssp. *bankia* はクリスマス島の歴史的な迷行種として記載される。しかし元の225×51 GBIF索引ではこの種×陸塊**1210**の2001–2010、2011–2017、2018–2022記録数は**いずれも0**。記録ゼロを「島で知られていない・種が不在」とすることはできない。2008年の迷行状態を2026年の定着状態と取り違えない。
- クリスマス島固有のChristmas Island Emperor（2008年報告の表記 *Polyura andrewsii*）が記載されているが、元の**225蝶の候補リストにはこの学名/対応するandrewsi綴りの蝶がない**。これはHOSTSで寄主が知られる225種の選択セットが、少なくとも一つの確認済み島固有蝶を含んでいない例である。島固有蝶全体の捕捉率は算出できない。分類群の現行異名は未照合のため、taxon keyを独断で割り当てない。

これらを主解析の現行の在来・移入ラベルへ直接代入しない。定着／迷行／在来／島固有は観測時点、分類階級（種と亜種）、地理範囲、チェックリストの定義が違う。

## 4. ソース固定・再現性・論文への反映

- [12種の根拠表](PHUQUOC_PRE2011_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv)は、受理学名、GBIF speciesKey、陸塊471、資料年、印刷頁、独立資料リンクと由来/固有性のunknownラベルを残す。
- [source SHA固定の全317組の再集計コード](../../scripts/audit_independent_island_historical_butterfly_presence_20261117.py)は、元11,475行・51島代表港・2006 9種＋2007 3種・2018–2022の12新規索引・クリスマス島の対照を機械検証する。
- [GitHub Actions CI #44074103818（成功）](https://github.com/zuizui0223/chocho/actions/runs/44074103818) に独立数値検証・故意の誤ラベル検出を含み、全317組の候補注釈と結果JSONを保存する。元のHOSTS寄主データ・GEB本体の+54.9%資源拡大、島+118.9%資源拡大、港と外来植物相、海洋抵抗の肯定結果は変更しない。

**結論の採用基準：** Figure 3Cでは「GBIFに累積記録された蝶の種リストが似る」という記録上の結果を保持するが、GBIF初索引を蝶の真の新規定着、在来蝶の増減、島固有蝶の入れ替わり、または生態学的均質化と断言しない。今後は、**島別独立チェックリスト×年代×在来・迷行・繁殖・固有ステータス**が判明する種に限って比較。レビュー済みでない種×島を真の不在や非固有へ補完しない。
