from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
path=W/'host_butterfly_first_record_pairs.parquet'
p=pd.read_parquet(path)
expected=json.loads((O/'first_record_lag_audit.json').read_text())
np.testing.assert_allclose(p.lag_lower,p.lower_butterfly-p.upper_plant,equal_nan=True)
np.testing.assert_allclose(p.lag_upper,p.upper_butterfly-p.lower_plant,equal_nan=True)
ordering=np.select([p.lag_lower.isna()|p.lag_upper.isna(),p.lag_lower.gt(0),p.lag_upper.lt(0)],['unresolved_date','plant_recorded_first','butterfly_recorded_first'],default='overlapping_intervals')
np.testing.assert_array_equal(ordering,p.ordering)
assert not (p.Region.eq('Canada') & (p.Source_butterfly.eq('Greve et al. (2017)')|p.Source_plant.eq('Greve et al. (2017)'))).any()
assert p.ordering.value_counts().to_dict()==expected['date_order_counts']
d=pd.read_csv(O/'transport_all225_evaluation_coverage.csv')
assert len(d)==225 and d.species.is_unique and set(p.species)<=set(d.species)
for name,q in [('first_record_pairs',p),('resolved_order_pairs',p[p.ordering.isin(['plant_recorded_first','butterfly_recorded_first'])])]:
    counts=q.groupby('species').size()
    d[name]=d.species.map(counts).fillna(0).astype(int)
groups=p.groupby('species').Region.nunique()
d['first_record_regions']=d.species.map(groups).fillna(0).astype(int)
d.to_csv(O/'date_forecast_all225_coverage.csv',index=False)
summary=[]
for label,mask in [('all225',pd.Series(True,index=d.index)),('primary_forecast155',d.temporal_evaluated_rows.gt(0)),('island_holdout95',d.joint_island_evaluated_rows.gt(0))]:
    subset=d[mask];q=p[p.species.isin(subset.species)]
    summary.append(dict(cohort=label,candidates=len(subset),species_with_date_pairs=q.species.nunique(),pairs=len(q),regions=q.Region.nunique(),species_region_groups=len(q[['species','Region']].drop_duplicates())))
audit=dict(status='PASS_INTERVAL_ARITHMETIC_AND_CANDIDATE_CROSSWALK',pair_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),summary=summary,ordering=p.ordering.value_counts().to_dict(),limits='Species overlap only; source regions are not equivalent to forecast cells or island landmasses. First-record interval is not introduction/establishment date. Multiple hosts and rows are not independent butterfly arrivals. No source geography repair inferred.')
(O/'date_forecast_coverage_audit.json').write_text(json.dumps(audit,indent=2))
lines=['# 独立初記録資料と予測対象種の被覆','', '既存の初記録組について、元の年代上下限から時間差区間と先行順序を再計算した。既知のCanada/Greve地域不整合が除外されていること、225候補との種名対応を確認。原典の全行確認や実際の導入年の確定ではない。','', '|対象|候補種数|年代組のある種数|植物×蝶×地域の組|地域数|蝶×地域群|','|---|---:|---:|---:|---:|---:|']
for x in summary:lines.append(f"|{x['cohort']}|{x['candidates']}|{x['species_with_date_pairs']}|{x['pairs']}|{x['regions']}|{x['species_region_groups']}|")
lines+=['','この表は種の重なりだけを示す。国・島・地方を含む原資料の地域名を、予測セルや25評価島へ同一視して結合していない。年代組がある種でも、評価島に独立した導入年代があることは保証しない。','', '全265組の内訳は植物の記録先行178、チョウ先行68、年代未確定16、区間重複3。複数食草を独立の到来イベントとして検定しない。導入／定着を直接裏付ける年代は別途原典確認が必要で、予測の時間差変数として使用済みとは数えない。']
(O/'date_forecast_coverage_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(audit))
