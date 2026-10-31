"""Recompute saved prediction scores and geographic support without model fitting."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/temporal/forecast'
geo=pd.read_parquet(R/'work/barriers/cell_geography.parquet')
rows=[];checks=[];support=[]
for mode in ['temporal','spatial']:
 path=W/f'predictions_barrier_{mode}.parquet'
 if not path.exists():continue
 p=pd.read_parquet(path)
 receipt=json.loads((O/f'forecast_barrier_{mode}_results.json').read_text())
 assert not p.duplicated(['speciesKey','cell']).any()
 p=p.merge(geo,on='cell',validate='many_to_one')
 masks={'all':np.ones(len(p),bool),'classified':p.landmass_id.notna(),'noncontinental':p.landmass_area_km2.lt(7e6),'continental':p.landmass_area_km2.ge(7e6),'unclassified':p.landmass_id.isna(),'noncontinental_excluding_two_largest_sampled':p.landmass_area_km2.lt(7e6)&~p.landmass_id.isin(['15','26']),'initially_unrecorded':p.baseline_probability.eq(0)}
 for group,mask in masks.items():
  q=p[mask]
  support.append(dict(mode=mode,group=group,rows=len(q),species=q.speciesKey.nunique(),cells=q.cell.nunique(),landmasses=q.landmass_id.nunique()))
  for model in ['persistence','M5_host_composition','B1_history','B2_landmass']:
   if model not in q:continue
   assert np.isfinite(q[model]).all() and q[model].between(0,1).all()
   err=q[model]-q.target
   for metric,values in [('mse',err**2),('mae',abs(err))]:
    score=values.groupby(q.speciesKey).mean().mean()
    rows.append(dict(mode=mode,group=group,model=model,metric=metric,score=score))
    if group=='all':
     expected=receipt['results'][model]['species_equal_'+metric]
     assert abs(score-expected)<1e-12
     checks.append(dict(mode=mode,model=model,metric=metric,max_error=abs(score-expected)))
 old=pd.read_parquet(W/f'predictions_{mode}.parquet')
 joined=p.merge(old[['speciesKey','cell','target','M5_host_composition']],on=['speciesKey','cell'],suffixes=('','_old'),validate='one_to_one')
 assert len(joined)==len(old)==len(p)
 assert np.array_equal(joined.target,joined.target_old)
 assert np.array_equal(joined.M5_host_composition,joined.M5_host_composition_old)
 checks.append(dict(mode=mode,original_targets_and_M5_identical=True,prediction_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
pd.DataFrame(rows).to_csv(O/'barrier_recomputed_scores.csv',index=False)
pd.DataFrame(support).to_csv(O/'barrier_geographic_support.csv',index=False)
(O/'barrier_prediction_validation.json').write_text(json.dumps({'checks':checks,'scope':'Saved prediction arithmetic and row identity only, not validation of ecological mechanism or model fitting','status':'PASS'},indent=2))
print(pd.DataFrame(support).to_string(index=False))
print('PASS',len(checks),'checks')
