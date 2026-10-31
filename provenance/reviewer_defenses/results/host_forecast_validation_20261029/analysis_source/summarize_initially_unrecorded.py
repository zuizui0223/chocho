"""Predeclared subgroup: no focal recording in predictor period, not true absence."""
from pathlib import Path
import json
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs';S=R/'work/temporal/forecast'
rng=np.random.default_rng(20261029);rows=[]
for mode in ['temporal','species','spatial']:
    path=O/f'forecast_{mode}_results.json'
    if not path.exists() or json.loads(path.read_text())['status']!='COMPLETE':continue
    p=pd.read_parquet(S/f'predictions_{mode}.parquet');p=p[p.baseline_probability.eq(0)].copy()
    count=pd.crosstab(p.species,p.spatial_block)
    sw=rng.multinomial(len(count),np.repeat(1/len(count),len(count)),size=2000)
    bw=rng.multinomial(len(count.columns),np.repeat(1/len(count.columns),len(count.columns)),size=2000)
    den=bw@count.to_numpy().T;weights=sw*(den>0)
    for metric in ['mse','mae']:
        a=p.M5_host_composition-p.target;b=p.M3_both-p.target
        a=a.pow(2) if metric=='mse' else a.abs();b=b.pow(2) if metric=='mse' else b.abs()
        delta=a-b
        agg=p[['species','spatial_block']].assign(delta=delta).groupby(['species','spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=count.index,columns=count.columns,fill_value=0)
        means=np.divide(bw@agg.to_numpy().T,den,out=np.zeros(den.shape),where=den>0)
        bootstrap=(means*weights).sum(axis=1)/weights.sum(axis=1)
        low,high=np.quantile(bootstrap,[.025,.975])
        reference=b.groupby(p.species).mean().mean();full=a.groupby(p.species).mean().mean()
        rows.append({'mode':mode,'metric':metric,'species':p.species.nunique(),'rows':len(p),'reference_M3':reference,'full_M5':full,'difference':full-reference,'relative_improvement_percent':100*(reference-full)/reference,'ci_low':low,'ci_high':high})
pd.DataFrame(rows).to_csv(O/'forecast_initially_unrecorded_comparison.csv',index=False)
(O/'forecast_initially_unrecorded_scope.json').write_text(json.dumps({'status':'COMPLETE_FOR_AVAILABLE_MODES','definition':'Focal recording probability=0 in2015-2018, evaluated in2019-2022; predeclared subgroup in model lock; NOT never-recorded historically, true absence, new colonisation or breeding','interval':'2000 crossed species x spatial-block resamples, fixed fitted predictions'},indent=2))
print(pd.DataFrame(rows).to_string(index=False))
