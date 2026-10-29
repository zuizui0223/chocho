"""Reconcile model scoring support to the entire 225-species candidate census."""
from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parents[2];S=R/'work/temporal/forecast';O=R/'outputs'
d=pd.read_csv(O/'forecast_candidate_disposition.csv')
expected=pd.read_csv(O/'forecast_sensitivity_all_225_support.csv').set_index('species')
assert len(d)==225 and d.species.is_unique
modes=['temporal','species','spatial','temporal_effort20','temporal_radius250','temporal_exclude_top_provider','temporal_precision100']
checks=[];pending=[]
for mode in modes:
    path=O/f'forecast_{mode}_results.json'
    if not path.exists() or json.loads(path.read_text())['status']!='COMPLETE':pending.append(mode);continue
    r=json.loads(path.read_text())
    p=pd.read_parquet(S/f'predictions_{mode}.parquet',columns=['species','cell','M5_host_composition'])
    assert p.M5_host_composition.notna().all() and not p.duplicated(['species','cell']).any()
    counts=p.groupby('species').size()
    d[mode+'_test_pairs']=d.species.map(counts).fillna(0).astype(int)
    assert len(p)==r['test_rows'] and p.species.nunique()==r['test_species']
    assert d[mode+'_test_pairs'].sum()==len(p)
    if mode.startswith('temporal_'):
        variant=mode.removeprefix('temporal_')
        want=d.species.map(expected[variant+'_test_pairs']).to_numpy()
        assert (d[mode+'_test_pairs'].to_numpy()==want).all()
    checks.append({'mode':mode,'scored_pairs':len(p),'species':p.species.nunique(),'candidate_rows':len(d)})
d['final_disposition']=d.stage_disposition
eligible=d.stage_disposition.str.startswith('BASELINE_ELIGIBLE')
d.loc[eligible,'final_disposition']='MODELED_VALIDATIONS_PENDING' if pending else 'MODELED_ALL_PLANNED_VALIDATIONS'
d.to_csv(O/'forecast_all_225_final_coverage.csv',index=False)
(O/'forecast_coverage_reconciliation.json').write_text(json.dumps({'status':'COMPLETE' if not pending else 'PARTIAL','checks':checks,'pending_modes':pending},indent=2))
print(json.dumps({'verified_modes':len(checks),'pending':pending}))
