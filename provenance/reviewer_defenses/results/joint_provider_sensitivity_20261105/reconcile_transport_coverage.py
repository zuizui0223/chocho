from pathlib import Path
import json
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast'
d=pd.read_csv(O/'forecast_all_225_final_coverage.csv')[['species','taxon_key','stage_disposition','final_disposition']]
d=d.rename(columns={'final_disposition':'original_forecast_disposition'})
assert len(d)==225 and d.species.is_unique and d.taxon_key.is_unique
states={}
details={}
for prefix,run in ([('transport',m) for m in ['temporal','spatial','species','temporal_effort20','temporal_exclude_top_provider']]+[('joint',m) for m in ['temporal','spatial','temporal_effort20','temporal_exclude_top_provider','island']]):
    mode=run if prefix=='transport' else 'joint_'+run
    receipt=O/f'forecast_{prefix}_{run}_results.json'
    status=json.loads(receipt.read_text())['status'] if receipt.exists() else 'NO_RECEIPT'
    states[mode]=status
    # Never certify the first model in an incomplete two-model file as a completed comparison.
    if status!='COMPLETE':
        d[mode+'_evaluated_rows']=pd.NA
        continue
    p=pd.read_parquet(F/f'predictions_{prefix}_{run}.parquet')
    cols=['M5_host_composition','T0_port_location','T1_port_connectivity'] if prefix=='transport' else ['M5_host_composition','B5_sea_terrain','T1_port_connectivity','J1_natural_transport']
    if run=='island':cols=['M5_host_composition','B5_sea_terrain','J1_natural_transport']
    assert np.isfinite(p[cols+['target']].to_numpy()).all()
    assert p[cols+['target']].ge(0).all().all() and p[cols+['target']].le(1).all().all()
    assert not p.duplicated(['speciesKey','cell']).any()
    assert set(p.speciesKey).issubset(set(d.taxon_key))
    expected_names=d.set_index('taxon_key').species
    assert (p.species==p.speciesKey.map(expected_names)).all()
    counts=p.groupby('speciesKey').size()
    d[mode+'_evaluated_rows']=d.taxon_key.map(counts).fillna(0).astype(int)
    assert d[mode+'_evaluated_rows'].sum()==len(p)
    details[mode]=dict(rows=len(p),species=int(p.speciesKey.nunique()),models=cols)
    if run=='island':
        details[mode]['islands']=int(p.landmass_id.nunique())
        counts=p.groupby('speciesKey').landmass_id.nunique()
        d['joint_island_evaluated_islands']=d.taxon_key.map(counts).fillna(0).astype(int)
        d['joint_island_disposition']=np.where(d[mode+'_evaluated_rows'].gt(0),'EVALUATED_SECONDARY_ISLAND_SUBSET',np.where(d.temporal_evaluated_rows.gt(0),'PRIMARY_MODELED_NO_ELIGIBLE_ISLAND_TEST_ROWS','NOT_PRIMARY_ELIGIBLE'))
d.to_csv(O/'transport_all225_evaluation_coverage.csv',index=False)
(O/'transport_all225_coverage_audit.json').write_text(json.dumps(dict(status='PASS_COVERAGE_RECONCILIATION_NOT_GOAL_COMPLETION',
    candidates=225,mode_states=states,mode_coverage=details,original_dispositions=d.original_forecast_disposition.value_counts().to_dict(),island_dispositions=d.joint_island_disposition.value_counts().to_dict(),missing='Blank means no completed comparison yet; zero means absent from completed comparison, not biological absence. Original disposition labels refer only to original forecast validations, not completion of this extended goal.',
    checks=['225 unique candidate names and keys','prediction keys and species names agree with original candidates','all required model columns and targets finite and bounded','candidate row sums match completed prediction tables']),indent=2))
lines=['# 全225候補と自然地理・輸送検証の対象範囲','', '元の225候補を保持し、完成済みの予測Parquetから種名・taxon key・評価行数を再照合した。元の採否列の MODELED_ALL_PLANNED_VALIDATIONS は元予測解析に限ったラベルであり、今回の拡張目標の完了ではない。','', '|検証|評価種数|評価行数|','|---|---:|---:|']
for name,v in details.items():lines.append(f"|{name}|{v['species']}|{v['rows']}|")
lines+=['','## 島検証の採否','']
for name,n in d.joint_island_disposition.value_counts().items():lines.append(f'- {name}: {n}種')
lines+=['','島検証に評価行がない主解析種は、チョウが島に不在であることを意味しない。原座標・セル面積の支持、期間をまたぐ観測努力、島対象範囲を満たす評価行がないという意味。全225種の行数と島数は transport_all225_evaluation_coverage.csv に保存した。','', '島検証は全25島・95種の追加評価。元の155種の主解析を島95種に置換しない。自然地理・輸送の種外比較は港単独追加で実施済みだが、統合モデルの未知種への一般化はこの表では検証していない。全目標は継続中。']
(O/'transport_all225_coverage_results_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(states))
