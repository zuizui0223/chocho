# 島への寄主植物の時間的先行を種キー・原年表で再検証（2026-10-31）

**実測結果：*Medicago sativa* と *Alnus incana* のGBIFの受理種キー・2001–2010年の植物観測を再取得し、既存のチョウ×島年表と接続した。** 3蝶リンク・2植物種・4植物種×島端点、重複を含む6蝶種×島端点で、植物が14〜19年早く索引された。

| チョウ | 既知の適合寄主 | 島 | 寄主植物の初索引年 | 蝶の初索引年 | 年差 |
| --- | --- | ---: | ---: | ---: | ---: |
| *Colias croceus* | *Medicago sativa* | 169 | 2002 | 2019 | 17 |
| *Colias croceus* | *Medicago sativa* | 276 | 2001 | 2015 | 14 |
| *Colias hyale* | *Medicago sativa* | 169 | 2002 | 2019 | 17 |
| *Colias hyale* | *Medicago sativa* | 276 | 2001 | 2015 | 14 |
| *Nymphalis antiopa* | *Alnus incana* | 276 | 2001 | 2015 | 14 |
| *Nymphalis antiopa* | *Alnus incana* | 715 | 2002 | 2021 | 19 |

### ソースの識別
GBIF厳密照合は *M. sativa* = `speciesKey 9151957`、*A. incana* = `speciesKey 9148577`。いずれもrank SPECIES、status ACCEPTED、match EXACT、confidence=100。植物元集計：2001–2010年の9島・全23港、9,788植物種キー×港×年のGBIFファセット（gzip SHA256 `b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e`）。蝶の索引年：元Parquet SHA256 `8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f` から既存のstrict year resolverで導出。

植物種キー×島端点は4つだけ。*Colias* 2種は同じ2島と寄主を共有しており、**6独立事例とは数えない**。個々の植物記録が幼虫の摂食場所・植物の定着年や島の外来性を示すわけではない。2001/2002年の初索引は実際の到着日ではない。船舶・コンテナ輸送、蝶の移動経路、同島での繁殖、H2/H3の寄主効果は未実証。

**肯定的に昇格するのは「この2つの既知寄主植物が対象の蝶記録より以前に、該当する島の港周辺で実際に索引されていた」という独立再取得の分類群・年の確認だけ**。同じGBIF索引基盤からの再取得なので独立した現地実験・繁殖証拠ではない。

次の中心検定には、全島×種×年でのラーバル寄主対応、当時の地域外来性、港と非港の植物記録・調査努力を繋げた**未選択の比較群**が必要。「見つからない」を植物不在に変換せず、現段階で寄主による定着率を主張しない。

[GBIF種名照合結果](focal_historical_host_taxon_year_audit.json) · [年順の出力JSON](focal_host_prior_to_butterfly_20261031.json) · [解析コード](match_dated_host_names.py)。
