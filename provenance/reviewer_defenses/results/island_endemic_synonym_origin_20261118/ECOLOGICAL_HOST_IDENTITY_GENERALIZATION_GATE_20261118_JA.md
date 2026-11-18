# 1本のGEB蝶論文に入れる寄主利用の新証拠：分類名・寄主種・個別島・幼虫成績の境界（2026-11-18）

## 現在地と論文上の意味

**独立在島蝶の寄主属の文献記録は、現生41種中36種（87.8%）まで確認できた。** 元のHOSTSに蝶の二名法が完全一致するのは42種中19種、保守的な異名統合でも21種のみ。HOSTSに現れないことは寄主が存在しないことではない。

*Euchloe grancanariensis* の野外記録を追加した結果、独立に分類したカナリア固有蝶13種のうち**8種で幼虫寄主属の名前が出典付きで確認できる**。この8種はすべて非`Hipparchia`で、残る5種はカナリアの`Hipparchia`近縁放散群に一致。**全非`Hipparchia`固有8/8種、`Hipparchia`固有0/5種**の源資料収録率である（文献努力を揃えた生物学的な食草既知確率ではない）。

原41種ごとの元HOSTS・異名・先行寄主文献・統一カナリア文献・今回2010年野外論文の出典区分を、`scripts/rebuild_independent41_named_larval_host_20261118.py`で**1種ずつ再合成して成功**： [CI #44214691925](https://github.com/zuizui0223/chocho/actions/runs/44214691925)。新規記録1種（`E. grancanariensis`）と従来35種の重複を排除。歴史的絶滅`Pieris wollastoni`には過去の寄主記録が存在しても、現生41種の分母には入れない。

## どの蝶が何の植物を使ったのか

| 独立した蝶種と地理 | 寄主植物の正確な同定と在来性 | 直接証拠の階級 | この資料で言える結果 |
|---|---|---|---|
| *Pieris cheiranthi*（テネリフェ・ラパルマ） | 外来・野生化 *Tropaeolum majus* | 各島の外来植物上で幼虫を観察・撮影 | 外来植物が**幼虫寄主として実際に利用される**（2蝶×島単位） |
| *Euchloe grancanariensis*（グラン・カナリア） | 在来 *Hirschfeldia incana* | 2009年に幼虫15個体をその植物の果実から採集 | **野生の幼虫利用**が確認でき、飼育で5個体蛹化 |
| 同種・同島 | 島固有の在来 *Descurainia preauxiana* | 2014年の産卵行動と卵の写真 | **産卵植物**を確認。幼虫の成長・羽化成績は別 |
| *Gonepteryx cleobule*（ラ・ゴメラ） | 在来 *Rhamnus crenulata* | 2013年、卵・幼虫の野外観察 | 在来植物への直接的な卵・幼虫利用 |
| *Hipparchia tamadabae*（グラン・カナリア） | **未同定の株立ち草（Poaceae）** | 2016年に草上の幼虫を野外観察 | 野外幼虫の資源**科**は示唆されるが、寄主属・種は特定できない |

後ろの5種の`Hipparchia`で、2016年の`H. tamadabae`幼虫は植物名未同定である。`Hyparrhenia hirta`は近くに生育したが著者が幼虫未確認と記述したため、寄主として**登録しない**。Spilani et al. 2025は五つの島固有種が独立分化したことを支持するが、野外幼虫寄主植物種を同定したとは述べない。系統・調査史が同じなので「固有種5種に独立にゼロ効果が生じた」とも言わない。

## 野外利用と飼育成績を接続する際の限界

Acosta-Fernández (2010, SHILAP 38(149):121–128) は`Euchloe belemia grancanariensis`（現在の地域独立分類では`E. grancanariensis`）の幼虫を2009年3月に`Hirschfeldia incana`の果実上から**15個体**採集した。**9個体は寄生蜂に寄生、1個体は原因不明死亡、5個体は蛹化**した。蛹化した2個体は異常蛹で死亡、3個体は休眠蛹。成虫羽化成功は原文から証明できない。よって5/15蛹化を外来/在来植物の質の優劣として扱わず、9/15を島全体の寄生率と解釈しない。

同論文フエルテベントゥラ島の姉妹群`Euchloe hesperidum`では、野外採集した**幼虫1個体**が飼育下で`Hirschfeldia incana`の果実を拒否し、`Sisymbrium erysimoides`と`S. irio`を食べ、蛹化した。受理植物種名は原GEBの固定WCVP v13で**3種ともカナリア在来**。これは同じ`Hirschfeldia`が**野外で1蝶種の寄主**であり、**別種の一幼虫によって飼育下で拒否**されたという対比。反復・移植・相互比較実験ではなく、種間の適応的選択差・植物の質・一般則を導出しない。

- 原著：https://www.redalyc.org/pdf/455/45514996013.pdf （原著は独立種名でなく亜種名を使用。分類群概念を保守的に引き継ぐ）
- `E. grancanariensis` 野外産卵：https://lepiforum.org/wiki/page/Euchloe_grancanariensis
- WCVP同じ原固定版で2植物がnativeと照合：[CI #44213255146](https://github.com/zuizui0223/chocho/actions/runs/44213255146)
- 近縁2蝶種・3植物の源資料対比、3植物すべてnativeと照合：[CI #44213967074](https://github.com/zuizui0223/chocho/actions/runs/44213967074)
- `Hipparchia`の種認識・野外寄主の未同定：Spilani et al. 2025 doi:10.1093/isd/ixaf043、https://www.pyrgus.de/Hipparchia_wyssii_en.html

## 原研究の主張をどう保持するか

この進展は**元GEB239種での潜在的資源拡大全球+54.9%／島嶼サンプル相対+118.9%を否定しない。** 今回の独立在島蝶の生物学的な実際の利用は、元の「寄主植物種の固定された文献ネットワーク×WCVP地域分布」の別指標である。生物学的に起きる過程を段階で区別する：

1. **資源供給の可能性：** 在来植物・外来植物が島にあるか（WCVP v13・地域EXOS）。
2. **寄主としての利用：** その植物種に卵・幼虫があるか（原論文や野外写真）。
3. **繁殖・個体群応答：** 成長、寄生、蛹化、羽化、個体数・定着が起きるか（2010文献の限られた蛹化データ、人口増減の直接証拠なし）。

現時点で支持できる肯定結果は`Pieris cheiranthi`で**外来植物の幼虫利用が2島で存在する**こと、`Euchloe grancanariensis`で**在来植物への幼虫利用と蛹化が存在する**こと、そして資料収録の偏りが**一つの近縁放散群へ集中**していること。これを**島固有蝶全般の平均資源拡大率・外来植物による定着促進・海洋障壁の因果媒介**に拡張するのはまだ不適切。

## 原資料・再現ファイル

- `E_grancanariensis_H_tamadabae_source_evidence_20261118.csv`：新たな幼虫・産卵・未同定草のsource-typed 4行。
- `Euchloe_two_island_host_use_field_vs_lab_2010.csv`：2010年の2蝶種×3提示植物の野外/飼育証拠を分類した4行（飼育群の標本は幼虫1個体）。
- `scripts/audit_Egrancanariensis_primary_host_20261118.py`、`scripts/audit_Egrancanariensis_WCVP_v13_20261118.R`：蝶名/原産性/寄主種と採集成績。
- `scripts/audit_Euchloe_native_host_acceptance_field_lab_20261118.py`、`scripts/audit_Euchloe_field_lab_WCVP_v13_20261118.R`：擬似反復を拒否した姉妹蝶の在来寄主比較。
- `scripts/rebuild_independent41_named_larval_host_20261118.py`：独立41現生蝶を行単位で再集計し、源別フラグを保持。
- 結果は各成功CIのZIP artifact（元の論文PDF自体のバイト同一性を今回のCIで検証したわけではなく、公開原文を手動精査した数値の転記とその追跡可能性を保存）。

投稿上は主論文を1本の`Global Ecology and Biogeography`として維持し、現在の源資料追加は**サンプリング・寄主同定境界の補足**へ位置づける。マージ・原稿書換・投稿は未実施。
