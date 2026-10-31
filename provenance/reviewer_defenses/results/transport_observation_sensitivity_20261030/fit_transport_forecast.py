from pathlib import Path
import os,json,time,hashlib,sys,argparse
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np,pandas as pd,joblib
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import BallTree
R=Path(__file__).resolve().parents[2];W=R/'work/temporal/forecast';O=R/'outputs'
parser=argparse.ArgumentParser()
parser.add_argument('mode',choices=['temporal','species','spatial'],default='temporal',nargs='?')
parser.add_argument('--sensitivity',choices=['effort20','radius250','exclude_top_provider','precision100'])
args=parser.parse_args();mode=args.mode
if args.sensitivity is not None:
 if mode!='temporal' or args.sensitivity not in {'effort20','exclude_top_provider'}:
  parser.error('Only locked temporal effort20 and exclude_top_provider transport sensitivities are supported')
 sensitivity_lock=json.loads((O/'transport_observation_sensitivity_lock.json').read_text())
 assert args.sensitivity in sensitivity_lock['variants']
runname=mode+('_'+args.sensitivity if args.sensitivity else '')
baseline_runname=runname
panel_path=W/(f'panel_{args.sensitivity}.parquet' if args.sensitivity else 'species_cell_period_panel.parquet')
lock=json.loads((O/'forecast_model_lock.json').read_text())
p=pd.read_parquet(panel_path)
if args.sensitivity:
 from scipy.stats import hypergeom
 n=20 if args.sensitivity=='effort20' else 10
 ok=p.days.ge(n)
 assert np.array_equal(p.recording_probability10.notna(),ok)
 np.testing.assert_allclose(p.loc[ok,'recording_probability10'],hypergeom.sf(0,p.loc[ok,'days'].astype(int),p.loc[ok,'success'].astype(int),n),atol=1e-12)
geo=pd.read_parquet(R/'work/barriers/transport_forecast_features.parquet')
extra=['log_nearest_indexed_port_km','indexed_port_support_unknown','log_nearest_port_connectivity']
p=p.merge(geo[['cell','period']+extra],on=['cell','period'],how='left',validate='many_to_one')
p['baseline_probability']=p.recording_probability10;p['log_days']=np.log1p(p.days);p['log_distance']=np.log1p(p.distance_km)
base=['baseline_probability','log_days','log_distance','latitude','longitude','provider_dominance']+['season_q'+str(i) for i in range(4)]
cl=['temperature','precipitation','temperature_range','warmest_month','coldest_month'];cl+= [x+'_trend_per_year' for x in cl]
lc=['cropland','crop_mosaic','natural_crop_mosaic','tree','tree_mosaic','herbaceous_mosaic','shrub','grass','wetland','built','bare','sparse','lichen','water','snow'];lc += [x+'_trend_per_year' for x in lc]
host=['available_hosts','available_families','global_hosts','available_share','native_host_count','introduced_host_count','introduced_share']
family=[c for c in p.columns if c.startswith('host_family_')]
full=base+cl+lc+host+family
models={'T0_port_location':full+extra[:2],'T1_port_connectivity':full+extra}
runname='transport_'+runname
def transition(period):
 a=p[p.period==period].copy();b=p[p.period==period+1][['speciesKey','cell','recording_probability10','days']].rename(columns={'recording_probability10':'target','days':'future_days'})
 a=a.merge(b,on=['speciesKey','cell'],validate='one_to_one');a=a[a.baseline_probability.notna()&a.target.notna()].copy()
 assert a[models['T1_port_connectivity']+['target']].notna().all().all()
 return a.reset_index(drop=True)
train=transition(0);test=transition(1);del p
for d in [train,test]:
 d['species_fold']=d.species.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
 d['spatial_fold']=d.spatial_block.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
assert not set(models['T1_port_connectivity'])&{'target','future_days','speciesKey','species','period'}
def evaluate(pred):
 q=test[['species','cell','region','spatial_block','baseline_probability','target']].copy();q['prediction']=pred;q=q[np.isfinite(q.prediction)]
 q['sqerror']=(q.prediction-q.target)**2;q['abserror']=abs(q.prediction-q.target)
 sp=q.groupby('species').agg(mse=('sqerror','mean'),mae=('abserror','mean'),rows=('cell','size'))
 new=q[q.baseline_probability==0]
 return dict(rows=len(q),species=len(sp),species_equal_mse=float(sp.mse.mean()),species_equal_mae=float(sp.mae.mean()),pooled_mse=float(q.sqerror.mean()),initially_unrecorded_rows=len(new),initially_unrecorded_mse=float(new.groupby('species').sqerror.mean().mean()) if len(new) else None),sp
out=test[['speciesKey','species','cell','region','spatial_block','climate_id','baseline_probability','target','days','future_days','species_fold','spatial_fold']].copy()
if mode=='spatial':
 out['distance_to_training_cell_km']=np.nan
 for fold in range(5):
  trcells=train[train.spatial_fold!=fold][['cell','latitude','longitude']].drop_duplicates('cell')
  tecells=test[test.spatial_fold==fold][['cell','latitude','longitude']].drop_duplicates('cell')
  tree=BallTree(np.radians(trcells[['latitude','longitude']].to_numpy()),metric='haversine')
  distances=tree.query(np.radians(tecells[['latitude','longitude']].to_numpy()),k=1)[0].ravel()*6371.0088
  lookup=dict(zip(tecells.cell,distances));mask=out.spatial_fold==fold
  out.loc[mask,'distance_to_training_cell_km']=out.loc[mask,'cell'].map(lookup)
 assert out.distance_to_training_cell_km.notna().all()
out['persistence']=test.baseline_probability.to_numpy();results={};results['persistence'],_=evaluate(out.persistence)
old=pd.read_parquet(W/f'predictions_{baseline_runname}.parquet')
assert np.array_equal(old[['speciesKey','cell']].to_numpy(),out[['speciesKey','cell']].to_numpy())
assert np.array_equal(old.target.to_numpy(),out.target.to_numpy())
out['M5_host_composition']=old.M5_host_composition.to_numpy();results['M5_host_composition'],_=evaluate(out.M5_host_composition)
settings=lock['random_forest'];folds=[None] if mode=='temporal' else list(range(5));foldcol=mode+'_fold'
for name,features in models.items():
 started=time.monotonic();prediction=np.full(len(test),np.nan)
 for fold in folds:
  tr=train if fold is None else train[train[foldcol]!=fold]
  mask=np.ones(len(test),dtype=bool) if fold is None else test[foldcol].eq(fold).to_numpy();te=test[mask]
  if not len(te):continue
  if fold is not None:
   key='species' if mode=='species' else 'spatial_block';assert set(tr[key]).isdisjoint(set(te[key]))
  weights=tr.species.map(1/tr.species.value_counts()).to_numpy();weights*=len(weights)/weights.sum()
  reg=RandomForestRegressor(**settings)
  reg.fit(tr[features].to_numpy(dtype='float32'),tr.target.to_numpy(),sample_weight=weights)
  prediction[mask]=reg.predict(te[features].to_numpy(dtype='float32'))
  if mode=='temporal':joblib.dump({'model':reg,'features':features,'settings':settings},W/f'model_{name}{"_"+args.sensitivity if args.sensitivity else ""}.joblib',compress=3)
  print(mode,name,'fold',fold,'train',len(tr),'test',len(te),flush=True)
  del reg
 assert np.isfinite(prediction).all()
 out[name]=prediction;results[name],sp=evaluate(prediction)
 sp.to_csv(O/f'forecast_{runname}_{name}_species_metrics.csv')
 out.to_parquet(W/f'predictions_{runname}.parquet',index=False)
 receipt=dict(status='IN_PROGRESS',mode=mode,train_rows=len(train),train_species=int(train.species.nunique()),test_rows=len(test),test_species=int(test.species.nunique()),results=results,settings=settings,model_lock_sha256=hashlib.sha256((O/'forecast_model_lock.json').read_bytes()).hexdigest(),interpretation='Predicts standardized recording probability; temporal holdout2019-2022; not verified true occupancy or causal drivers')
 receipt.update(extension_lock_sha256=hashlib.sha256((O/'transport_forecast_extension_lock.json').read_bytes()).hexdigest(),transport_feature_sha256=hashlib.sha256((R/'work/barriers/transport_forecast_features.parquet').read_bytes()).hexdigest(),sensitivity=args.sensitivity,input_panel=str(panel_path),rarefaction_n=20 if args.sensitivity=='effort20' else 10)
 if args.sensitivity:
  receipt.update(sensitivity_lock_sha256=hashlib.sha256((O/'transport_observation_sensitivity_lock.json').read_bytes()).hexdigest(),input_panel_sha256=hashlib.sha256(panel_path.read_bytes()).hexdigest(),baseline_predictions_sha256=hashlib.sha256((W/f'predictions_{baseline_runname}.parquet').read_bytes()).hexdigest())
 (O/f'forecast_{runname}_results.json').write_text(json.dumps(receipt,indent=2))
 print(name,results[name],'seconds',round(time.monotonic()-started,1),flush=True)
receipt['status']='COMPLETE';(O/f'forecast_{runname}_results.json').write_text(json.dumps(receipt,indent=2));print('COMPLETE',runname,flush=True)
