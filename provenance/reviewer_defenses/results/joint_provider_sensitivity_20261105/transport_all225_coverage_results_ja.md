# 全225候補と自然地理・輸送検証の対象範囲

元の225候補を保持し、完成済みの予測Parquetから種名・taxon key・評価行数を再照合した。元の採否列の MODELED_ALL_PLANNED_VALIDATIONS は元予測解析に限ったラベルであり、今回の拡張目標の完了ではない。

|検証|評価種数|評価行数|
|---|---:|---:|
|temporal|155|208482|
|spatial|155|208482|
|species|155|208482|
|temporal_effort20|155|159084|
|temporal_exclude_top_provider|155|191004|
|joint_temporal|155|208482|
|joint_spatial|155|208482|
|joint_temporal_effort20|155|159084|
|joint_temporal_exclude_top_provider|155|191004|
|joint_island|95|41812|

## 島検証の採否

- EVALUATED_SECONDARY_ISLAND_SUBSET: 95種
- NOT_PRIMARY_ELIGIBLE: 70種
- PRIMARY_MODELED_NO_ELIGIBLE_ISLAND_TEST_ROWS: 60種

島検証に評価行がない主解析種は、チョウが島に不在であることを意味しない。原座標・セル面積の支持、期間をまたぐ観測努力、島対象範囲を満たす評価行がないという意味。全225種の行数と島数は transport_all225_evaluation_coverage.csv に保存した。

島検証は全25島・95種の追加評価。元の155種の主解析を島95種に置換しない。自然地理・輸送の種外比較は港単独追加で実施済みだが、統合モデルの未知種への一般化はこの表では検証していない。全目標は継続中。
