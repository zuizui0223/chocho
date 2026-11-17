# スリランカ島固有チョウの分類と寄主制約：歴史18学名のうち未解決6学名を再照合（2026-11-17）

**状況：独立資料・公式GBIFの再検証完了。** 2006年IUCNのE（固有）印がある18学名は元のHOSTS由来225種と二名法完全一致0。以前、旧GBIF Backboneで12学名の受理分類キーが元225にないことを検証した。今回は地域蝶類保全団体・2024年査読研究・2021年寄主植物研究で残る6学名を追い、厳格な旧GBIF Backbone NUBキー照合を再実行した。

## 名前を確認した6件と同一性の確証レベル

| 歴史2006学名 | 現在の地域専門資料が使う学名 | 旧GBIF Backbone種キーの元225照合 |
|---|---|---|
| Lethe dynaste | **Lethe dynsate** | **1918426**、元225に不在（厳格照合） |
| Elymnias singala | **Elymnias singhala** | **1895537**、元225に不在（厳格照合） |
| Spindasis nubilus | **Cigaritis nubilus** | チェックリスト名は確認、旧Backboneへの種キー連結は未確定 |
| Jamides coruscans | **Jamides coruscans** | 同上 |
| Udara lanka | **Udara lanka** | 同上 |
| Halpe egena | **Halpe egena** | 同上 |

**分類キーの重要な訂正：** GBIFのspecies/searchは個別チェックリスト固有の使用IDも返す。検索された数字IDを、元の225種の旧GBIF Backbone種IDと直接比べてはならない。再実行では、旧GBIF v1のEXACT SPECIES ACCEPTED/SYNONYMか、確認済みBackbone NUB参照のみを比較可能とした。結果として**今回2学名の不在が新たに検証され、以前の12と合わせ14学名**が旧Backbone種キー比較で元225の外側。残り**4学名は現行地域名・固有種扱いが外部資料で確認できても、旧Backboneキー未解決**。それらの自然界での不在を意味しない。

実行根拠：[6学名の地域名・寄主・原著リンク](SRI_LANKA_2006_SIX_UNRESOLVED_NAME_CURRENT_SOURCES_AND_HOSTS_20261117.csv) ／ [厳格に比較可能な旧GBIF Backboneキーの判定](SRILANKA_SIX_MODERN_NAMES_STRICT_NUB_GATES_20261117.csv) ／ [公式GBIF API・元225種SHAの再現成功CI #44088382011](https://github.com/zuizui0223/chocho/actions/runs/44088382011)。ソース応答JSONを含むアーティファクト14331025522、SHA256 6520459bab9eb803400525ff492a31bad51664805c303a2a50ddd5397a06c161。対応コード：[source-locked six-name audit](../../scripts/reconcile_srilanka_six_regional_modern_to_GBIF_20261117.py)。

出典：Butterfly Conservation Society of Sri Lanka の現行 [Elymnias singhala](https://www.slbutterflies.lk/explore/butterfly/Elymnias%20singhala)、[Jamides coruscans](https://www.slbutterflies.lk/explore/butterfly/Jamides%20coruscans)、[Udara lanka](https://www.slbutterflies.lk/explore/butterfly/Udara%20lanka)、[Halpe egena](https://www.slbutterflies.lk/explore/butterfly/Halpe%20egena) 種別資料（2025–2026）。[2024年査読論文 Cigaritis nubilus はSri Lanka固有](https://www.taprobanica.org/Archives/Volume-9-14-2020-26/Volume-13-Number-2-November-2024/v13i2343.html)、[Spindasis nubilusとの名称対応を掲載する分類カタログ](https://bio.nic.funet.fi/index/Science/bio/life/insecta/lepidoptera/ditrysia/papilionoidea/lycaenidae/theclinae/cigaritis/)、[Lethe dynsateを扱う2021年査読論文](https://threatenedtaxa.org/JoTT/article/download/6875/7356?inline=1)。歴史学名から地域現行名への対応は外部文献の概念対応であり、全ての組について厳密な法的/分類学的異名関係を独立検証したわけではない。GBIF公式が2026年に導入したCOL XRは別体系なので、元225の旧Backbone数字キーと混用しない：[GBIF移行解説](https://data-blog.gbif.org/post/catalogue-of-life-taxonomic-backbone/)。

## 固有蝶の資源依存は同じではない

[Jayasinghe, Rajapakshe & Ranasinghe 2021, Journal of Threatened Taxa, DOI:10.11609/jott.6875.13.2.17731-17740](https://threatenedtaxa.org/JoTT/article/download/6875/7356?inline=1) はスリランカの幼虫食草の既報を統合し、**223蝶種、植物582種（地域固有70、在来351、外来161）、蝶–植物相互作用1,091件**と報告する。植物の分布・密度、実際の供給量や蝶の個体群成績を直接測った値ではない。

| 島固有蝶 | 2021年時点の既知幼虫寄主 | 解釈上の注意 |
|---|---|---|
| **Halpe egena** | 島固有植物 Davidsea attenuata のみを当時の確認済み寄主として記録 | 追加寄主が未発見の可能性 |
| **Lethe dynsate** | 島固有植物 Ochlandra stridula のみを当時の確認済み寄主として記録 | 他の竹類も利用し得ると著者が指摘 |
| **Elymnias singhala** | 島固有 Calamus ovoideus、Loxococcus rupicola、その他のヤシ類 | 固有植物に排他的に依存しない |

さらに同論文は、Lethe daretis（島固有）で当初の既知植物以外のイネ科植物への産卵と羽化までの飼育を記載する。したがって**島固有蝶＝必ず自島固有寄主だけを利用する**という単純化は支持されず、出発時の寄主相と記録の完全性が重要。現行地域専門資料では、Jamides coruscans に Humboldtia laurifolia、Udara lanka に Persicaria chinensis・Smithia blanda の幼虫寄主も掲載しているが、後2種の植物の地域原産区分はこの監査では未照合。

**本体との接続：** chochoの全球+54.9%、島+118.9%は、既存HOSTSで名前が分かっている蝶とその潜在的寄主植物の地理的供給量の結果であり、島固有蝶全体の資源拡大率ではない。今回の研究は、元225種へ採用されなかった島固有蝶にも幼虫寄主文献があり、寄主の分類群・地域固有性・未発見性が一様でないという**選択機構と生態学的対照**を追加する。今回の3種に外来植物利用が必ずないとは言えず、全球の比較効果を推定したものでもない。

**次の可検証な主目標：** 対象蝶の母集団をHOSTSから作るのでなく、独立の島嶼蝶類チェックリストで定義してからHOSTS被覆率を算出し、固有性の違いによる記録網の欠測と資源利用の記述的偏りを検証する。その際、GBIF未索引を生態学的不在や移入初到達へ変換しない。
