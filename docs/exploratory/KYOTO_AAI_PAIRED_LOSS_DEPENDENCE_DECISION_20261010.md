# 京都アゲハ2種・AAI試験：欠測が葉ペア内で集中することを原データで確認

> **Follow-up source-format correction (2026-10-10):** All seven original paired joint losses are in source groups with simple `l.id` records; 11 other source groups have comma-composite IDs. Full-cohort baseline initial leaf area is smaller in joint-loss groups (6.09 versus 9.17; source-strain conditional p=0.0057), but within the 19 *single-ID source groups* it is 6.09 versus 6.82 and the same source-strain conditional test is **p=0.2207**. This is an exploratory, posthoc necessary confounding check; **do not claim leaf size caused joint loss**, nor use the whole-cohort p-value alone. See [final stop decision](KYOTO_AAI_PRETREATMENT_ATTRITION_FINAL_STOP_20261010.md). The overlap signal is a paired observation pattern, **not** an identified plant-chemical mechanism.

**2026-10-10 | 実際に実行した新しい原データ監査・探索的統計診断 | 新しい植物化学生態学の因果効果は未発見**

## 何がデータで識別できたか

[Hashimoto & Ohgushi (2023), Ecol Evol, DOI 10.1002/ece3.10164](https://doi.org/10.1002/ece3.10164) の原著データ（[Figshare 23170898](https://doi.org/10.6084/m9.figshare.23170898), original `bioassay_larvae.csv`, MD5 `982d4931da306a7ff8e8d550cffe5246`）には、3齢幼虫を同じ食草由来の半葉2枚にそれぞれ割り当てた、AAI添加とエタノール対照の24時間試験がある。合計120個体、*A. alcinous* と *S. montela* 各60個体・30対応葉ペア。原著では `loss=1` を haphazard な幼虫死亡として、成長・摂食量のモデルから除外している。

今回、欠測を結果から削除せず、同一の `l.id` ごとに対応を復元した。

| *A. alcinous*（元の30葉ペア） | 対照は非脱落 | 対照は脱落 |
| --- | ---: | ---: |
| AAIは非脱落 | 16 | 1 |
| AAIは脱落 | 6 | 7 |

AAI側の脱落率 **13/30 =43.3%**、対照 **8/30 =26.7%**。対応葉を単位にしたリスク差（AAI−対照）は **+16.7百分率ポイント**、**9,999回ペア単位ブートストラップ95%区間 [0.0,+33.3]百分率ポイント**。**対応あり正確McNemar両側検定 p=0.125**。したがって、この試験の `loss` コードから **AAIが死亡を増加させたという結論は出ない**（同時に無害・同等とも結論できない）。

*Sericinus montela* 60個体の原著 `loss` は両処理とも0である。種差は発達状態、採集系統、飼育・取り扱いの差などとも交絡するため、種間のAAI耐性差をこの脱落頻度だけから主張できない。

## 新しい探索的診断：処理差より葉ペア共通の要因に注目

*Atrophaneura* の30葉ペアでは、AAI側と対照側 **両方が脱落した葉が7組**。両側の脱落が葉間で独立である場合、元のAAI13件・対照8件の周辺合計を固定した期待値は **13×8/30=3.47組**。

- 葉ペア内の脱落関連のオッズ比：**(7×16)/(6×1)=18.67**。
- 周辺合計を固定した **探索的Fisher正確両側 p=0.00936**。
- しかし各葉ペアの幼虫は著者の `strain` 区分を共有している。全30ペアの半葉で `strain` は一致し、区分は4種類（a12, a29, a36, a5）。`strain` は家系/由来・飼育群の代理となりうる。
- **追加の事後的感度分析**：各 `strain` 区分内の葉数とAAI/対照の脱落件数をすべて固定し、両側脱落件数の厳密な超幾何分布を4区分で畳み込んだ。**期待両側脱落3.729組に対し実際7組、同側集中の片側尾確率 p=0.00711**。
- この事後的感度分析は **仮説生成上の依存性診断** であり、独立の確証ではない。Fisher両側 p と strata固定片側 p は異なる帰無仮説・方向規約であり、都合よく比較してどちらかを「主解析」とはしない。

著者の `strain` だけでは両側脱落の同時性を説明しきれないが、そこから直ちに「葉の化学物質が2幼虫を殺した」とは推論できない。共同の葉片の状態、採集・取り扱い・容器、同じ葉の2半片で相関する未測定の特性などが残る。**実際の死亡原因はデータにない**。この観察は原著の24時間AAI無効果を一般化できないことを支持する、主にデータ品質/依存構造の新しい監査結果である。

## 生存者だけに対する成長率比較の制約

同じ既報の原著表を対応ペアのまま再解析した結果、*A. alcinous* の成長率比（24時間後−初期体重）/初期体重のAAI−対照は、**生存・観測された16対応葉だけで +0.0744、95%区間 [−0.0778,+0.2364]**。*S. montela* は30葉の比較で **−0.0018、[−0.1948,+0.1927]**。両者とも広く、**有意差がない ≠ 生物学的な同等性**である。16葉のみの解析は治療後の脱落で対象を選択する。

この試験が測定するのは **AAI単一化合物を外から塗布した24時間の3齢幼虫摂食/成長**。先行 [Nishida & Fukami (1989), DOI 10.1007/BF01014731](https://doi.org/10.1007/BF01014731) が報告する複数のアリストロキア酸類と水溶性葉成分の相互作用や、実際に他の幼虫が食害した直後の植物誘導応答は測っていない。

## 次の研究の実験設計に反映する具体事項

1. 化学添加・葉の質試験は **同一元葉から分けた処理・対照**を必ず追跡し、原葉/採集系統/日付/容器/観察者のIDを事前記録する（すべて同じものと混同しない）。
2. **幼虫の脱落・死亡・事故・不明**を実験前に定義し、すべての初期割当を保持する。生き残りだけの成長率を機構の一次エンドポイントにしない。
3. 可能なら葉の物理的損傷、食害履歴、補充された若葉の質を分けて操作し、独立した葉/飼育ケージ間で真の反復を得る。
4. **本研究の未解決の主要因果問題**は、外来 *S. montela* 存在下での在来 *A. alcinous* 生存低下が、幼虫の高需要齢期の一時的資源制約か、摂食後の植物質か、直接的干渉かである。今回のAAI再解析だけではどれも識別できない。
5. 追加の「有意差が出るサブセット」探索を禁止し、独立の操作実験に進む。

## 再現可能性と manuscript の境界

- 原本照合、葉ペア内脱落率、McNemar、bootstrap、source-strain条件付き超幾何分布：`scripts/analyze_kyoto_aai_paired_loss.py`。
- 先に主効果の方法を固定：`KYOTO_AAI_PAIRED_LOSS_ANALYSIS_PROTOCOL_V01.json`。事後の `strain` 感度分析は別の `KYOTO_AAI_PAIRED_LOSS_SOURCE_BATCH_SENSITIVITY_V01.json` として明示。
- Source-integrity + test workflow [38040771219](https://github.com/zuizui0223/chocho/actions/runs/38040771219)：**9ソフトウェアテスト成功**、原データハッシュ一致、JSON receipt保存。
- 独立した Global Ecology and Biogeography 論文 PR #38 は **一切変更せず**。蝶の外来化と植物寄主の外来化は別の生態学的問い。
