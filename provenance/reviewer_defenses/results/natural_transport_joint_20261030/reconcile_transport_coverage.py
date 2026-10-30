from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast'
d=pd.read_csv(O/'forecast_all_225_final_coverage.csv')[['species','taxon_key','stage_disposition','final_disposition']]
d=d.rename(columns={'final_disposition':'original_forecast_disposition'})
assert len(d)==225 and d.species.is_unique and d.taxon_key.is_unique
states={}
for prefix,run in ([('transport',m) for m in ['temporal','spatial','species','temporal_effort20','temporal_exclude_top_provider']]+[('joint',m) for m in ['temporal','spatial','temporal_effort20']]):
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
    assert p[cols].notna().all().all()
    assert not p.duplicated(['speciesKey','cell']).any()
    assert set(p.speciesKey).issubset(set(d.taxon_key))
    counts=p.groupby('speciesKey').size()
    d[mode+'_evaluated_rows']=d.taxon_key.map(counts).fillna(0).astype(int)
    assert d[mode+'_evaluated_rows'].sum()==len(p)
d.to_csv(O/'transport_all225_evaluation_coverage.csv',index=False)
(O/'transport_all225_coverage_audit.json').write_text(json.dumps(dict(status='PASS_COVERAGE_RECONCILIATION_NOT_GOAL_COMPLETION',
    candidates=225,mode_states=states,missing='Blank means no completed comparison yet; zero means absent from completed comparison.',
    checks=['225 unique candidate names and keys','completed prediction keys within original candidates','all three model columns populated','candidate row sums match completed prediction tables']),indent=2))
print(json.dumps(states))
