# HOSTSの所在地付き寄主利用報告と島の外来植物：正式解析プロトコル

**2026-11-16 | 原本取得・ソース監査済み。利用率の推定値はActions出力を別途確認すること。**

## 何を実測しているか

- 原データ: Natural History Museum HOSTS 原本（GloBIで固定された `globalbioticinteractions/HOSTS` コミット `808e0b869f9ec1adf8efff87cf6a395adda103e0`、NHM取得DOI `10.5519/qd.bsucrxdz`）。140,485件、Location非空140,484件、所在地表記285種類。原本の原記録 ID, Insect Genus/Species, Hostplant Genus/Species, Location, Lab Rearingを読む。
- 植物学名・分布: WCVP v13（2024-05-21）固定 `matildabrown/rWCVPdata` SHA `65bed76bae9d644ccb6ad200c05f9f5071d89e05`。HOSTS学名をWCVP受容名に厳密に対応し、地域の導入区分を適用。未一致・同一地域で両義的・複数植物地理区で区分不一致の記録は **unknown**（割合の分母に入れない）。
- 作物: FAO ICC Appendix 4（https://www.fao.org/4/a0135e/A0135E10.htm）。既存の `scripts/extract_fao_crop_taxa.py` を再使用し、学名完全一致による作物除外と、FAOが「spp.」と明示した属の保守的除外を分ける。作物を広く属全体とみなすことは、過剰除外の感度分析。
- 港: 既存 `results/three-gate-wpi-port-v02/wgsrpd3_port_counts_corrected_v02.csv` の369地域の2023 WPI港数（対数化）。船による運搬や港からの植物導入年代ではない。
- 面積: `results/island_resource_gap_filling/areas.csv` のWGSRPD3面積。観測個体数や実際の幼虫生息面積ではない。
- チョウ: 固定原稿の239種（`data/frozen/figure_sources/anthropogenic_species_metrics.csv`）。HOSTSの*全鱗翅目*を対象にした解析は**補足の比較**とし、チョウの実測と混ぜない。

## 位置の対応と島の分類

元の `Location` は自由記述。事前固定した [地域照合表](../../data/frozen/HOSTS_locality_WGSRPD3_strict_crosswalk_20261116.csv) に完全一致した記録だけを分類し、世界全域（Nearctic、Europe等）や地域未確定の記録を除外する。

**主解析**: 単独の植物地理区への確実な照合／複数区のすべてでWCVP在来・導入区分が一致する複合表記。Japan、Indiaなど行政範囲とWGSRPD3範囲が厳密には一致しないものは**拡張感度**に分ける。区分が混在する島々（Japan、Philippines、Cubaなど）は海洋島／大陸島の二値対比から除き、島全体の集計では区分未確定として保持。

*島の地質学的な操作定義*: oceanicには孤立火山島・島弧を含み、continentalには大陸棚島・大陸地殻片を含む。太平洋のタイプ整理には Nunn et al. (2016, Island Studies Journal, https://doi.org/10.1186/s40562-016-0041-8) の地域分類を参照。プエルトリコの島弧構造はUSGS等の独立地質資料に合致。個々の地域・広域名での帰属には不確実性があり、混合を二分しない。これは公開された全球の島由来属性データを完全に置き換えるものではないため、現段階では**手動監査付き操作分類**と呼ぶ。

## 統計単位と判定

1. **各HOSTS所在地の報告割合:** WCVPで地域外来と確認できた記録／（地域外来＋在来と確認できた記録）。異なる昆虫・植物ペア・出版記録が完全に独立とはみなさない。HOSTS `Lab Rearing=Y` を主解析から除外するが、空欄が**野外利用を保証するわけではない**。
2. **島 vs 本土**、**海洋島 vs 大陸島**: 地点中央値の差と個別地点の記録件数、Mann–Whitney（探索的）、所在地単位ブートストラップを併記。数値を報告するときにn地点を必ず示す。記録数で重み付けした割合は別列。
3. **港数を説明変数:** 所在地を等重み単位とした経験logitに、島ダミー、`log2(1 + WPI港数)`、植物地理区面積、HOSTS記録数を加える。海洋島・大陸島分類を追加したモデルも実行。小標本／共通支持不足では係数を採用しない。係数は観測上の関連で、港による因果媒介ではない。
4. **作物除外:** FAO学名完全一致、FAO種と属単位ワイルドカードの保守的除外。サンプル数も変わるので必ず同時出力。
5. **観測偏り対照:** 固定239チョウ種と全HOSTS鱗翅目を分離。元記録行ベースと同じ昆虫×寄主×所在地の重複除去ベースを分離。主張を左右する定義依存性を隠さない。

## 推論と投稿の境界

**「実際に使われている」には段階がある。** この解析で確かめられるのは「当該地域に紐づいた寄主利用の文献報告」の割合で、確認された野外の産卵・幼虫摂食・その植物での発育完了を全行で同定したことではない。HOSTSは刊行物・飼育記録の集成であり、文献収集の重点が農作物や害虫、研究の盛んな島に偏る。

この割合を既存の「島の潜在資源の54.3%が外来宿主により追加」（蝶×地域という**別の分母**）と同じ指標として扱わない。二つが近い場合も一致は独立再現ではない。港数とこの地域報告率の関連は供給・運搬の経路や新規蝶定着の証拠ではない。

**既存の肯定結果を維持する:** 全球の潜在資源+54.9%、島の相対増加+118.9%、寄主地理による155蝶種の記録予測+4.12%、港での広い移入植物相の記録と海洋地理の条件付き予測情報。新しい利用率は成功した正式集計が原本に再現された後でだけ、第三の結果候補として原稿に追加する。

## 再現コマンド

- `scripts/audit_HOSTS_location_source_20261116.py`: 原本の地名・飼育記録・FAO/港ソース監査。
- `scripts/build_wcvp_hosts_sidecar.R`: HOSTS→WCVP受容植物学名の厳密対応。
- `scripts/extract_HOSTS_location_origin_20261116.R`: WCVP地域固有の在来/導入フラグ抽出。
- `scripts/analyze_HOSTS_reported_use_island_fao_ports_20261116.py`: チョウ固定パネル／全鱗翅目、島の起源、作物除外、港の条件付きモデル。
- `.github/workflows/chocho-HOSTS-reported-use-oceanic-ports.yml`: source-pinned入力の取得から結果artifactの確認まで。
