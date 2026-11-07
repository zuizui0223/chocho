"""Locked B5 plus intermediate-island geometry; unchanged paired evaluation."""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,hashlib,json
import numpy as np,pandas as pd,joblib
from sklearn.ensemble import RandomForestRegressor
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';F=R/'work/temporal/forecast';O=R/'outputs'
parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['temporal','spatial']);args=parser.parse_args();mode=args.mode
lock=O/'intermediate_island_geometry_lock.json'
assert json.loads((O/'intermediate_island_geometry_audit.json').read_text())['status']=='GEOMETRY_EXTRACTION_COMPLETE_NOT_MODEL_OR_RELAY_RESULT'
model_lock=json.loads((O/'forecast_model_lock.json').read_text())
base_features=joblib.load(F/'model_B5_sea_terrain.joblib')['features']
addon=['log_intermediate_islands','intermediate_identity_unknown'];features=base_features+addon
assert len(features)==len(set(features))
p=pd.read_parquet(F/'species_cell_period_panel.parquet')
source_specs=[('history_geography_features.parquet',['speciesKey','cell','period']),('land_dem_features.parquet',['cell']),('sea_route_features.parquet',['speciesKey','cell','period']),('intermediate_island_features.parquet',['speciesKey','cell','period'])]
for source,keys in source_specs:
    q=pd.read_parquet(W/source)
    cols=[c for c in features if c not in p.columns and c in q.columns]
    p=p.merge(q[keys+cols],on=keys,how='left',validate='many_to_one' if keys==['cell'] else 'one_to_one')
p['baseline_probability']=p.recording_probability10
p['log_days']=np.log1p(p.days);p['log_distance']=np.log1p(p.distance_km)
def transition(k):
    a=p[p.period.eq(k)].copy()
    b=p[p.period.eq(k+1)][['speciesKey','cell','recording_probability10']].rename(columns={'recording_probability10':'target'})
    q=a.merge(b,on=['speciesKey','cell'],validate='one_to_one')
    q=q[q.baseline_probability.notna()&q.target.notna()].reset_index(drop=True)
    assert np.isfinite(q[features+['target']].to_numpy()).all()
    q['spatial_fold']=q.spatial_block.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
    return q
train=transition(0);test=transition(1);del p
assert train.species.nunique()==test.species.nunique()==155
reference=pd.read_parquet(F/f'predictions_sea_{mode}.parquet')
np.testing.assert_array_equal(test[['speciesKey','cell','target']],reference[['speciesKey','cell','target']])
out=reference.copy()
pred=np.full(len(test),np.nan);completed=[]
for fold in ([None] if mode=='temporal' else range(5)):
    tr=train if fold is None else train[train.spatial_fold.ne(fold)]
    mask=np.ones(len(test),dtype=bool) if fold is None else test.spatial_fold.eq(fold).to_numpy()
    te=test[mask]
    if fold is not None:assert set(tr.spatial_block).isdisjoint(te.spatial_block)
    weights=tr.species.map(1/tr.species.value_counts()).to_numpy();weights*=len(weights)/weights.sum()
    reg=RandomForestRegressor(**model_lock['random_forest'])
    reg.fit(tr[features].to_numpy(dtype='float32'),tr.target.to_numpy(),sample_weight=weights)
    pred[mask]=reg.predict(te[features].to_numpy(dtype='float32'))
    completed.append(dict(fold=fold,train_rows=len(tr),test_rows=len(te)))
    (O/f'forecast_intermediate_{mode}_progress.json').write_text(json.dumps(dict(status='FITTING_NOT_COMPLETE',completed_folds=completed)))
    print(mode,'completed fold',fold,flush=True)
assert np.isfinite(pred).all();out['B6_intermediate_islands']=pred
out.to_parquet(F/f'predictions_intermediate_{mode}.parquet',index=False)
results={}
for name in ['M5_host_composition','B5_sea_terrain','B6_intermediate_islands']:
    sq=(out[name]-out.target)**2;ab=abs(out[name]-out.target)
    results[name]=dict(species_equal_mse=float(sq.groupby(out.species).mean().mean()),species_equal_mae=float(ab.groupby(out.species).mean().mean()))
receipt=dict(status='COMPLETE',mode=mode,candidates=225,train_species=155,test_species=155,train_rows=len(train),test_rows=len(test),results=results,settings=model_lock['random_forest'],features=features,completed_folds=completed,source_hashes={s:hashlib.sha256((W/s).read_bytes()).hexdigest() for s,_ in source_specs},lock_sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),reference_prediction_sha256=hashlib.sha256((F/f'predictions_sea_{mode}.parquet').read_bytes()).hexdigest(),limits='New post-primary-result comparison; geometric intermediate-island candidates only, not verified relay or occupied stepping stones; historical catalogue and coastal geometry remain retrospective.')
(O/f'forecast_intermediate_{mode}_results.json').write_text(json.dumps(receipt,indent=2))
(O/f'forecast_intermediate_{mode}_progress.json').write_text(json.dumps(dict(status='COMPLETE',completed_folds=completed)))
print(json.dumps(receipt['results']),flush=True)
