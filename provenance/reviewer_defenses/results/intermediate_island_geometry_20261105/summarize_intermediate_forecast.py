from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,hashlib,sys
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';F=R/'work/temporal/forecast'
mode=sys.argv[1]
r=json.loads((O/f'forecast_intermediate_{mode}_results.json').read_text());assert r['status']=='COMPLETE'
path=F/f'predictions_intermediate_{mode}.parquet';p=pd.read_parquet(path)
refpath=F/f'predictions_sea_{mode}.parquet';ref=pd.read_parquet(refpath)
assert hashlib.sha256(refpath.read_bytes()).hexdigest()==r['reference_prediction_sha256']
np.testing.assert_array_equal(p[ref.columns],ref)
assert hashlib.sha256((O/'intermediate_island_geometry_lock.json').read_bytes()).hexdigest()==r['lock_sha256']
for name,digest in r['source_hashes'].items():assert hashlib.sha256((W/name).read_bytes()).hexdigest()==digest
models=['M5_host_composition','B5_sea_terrain','B6_intermediate_islands']
assert np.isfinite(p[models+['target']].to_numpy()).all()
for m in models:
    p[m+'_mse']=(p[m]-p.target)**2;p[m+'_mae']=abs(p[m]-p.target)
    for metric in ['mse','mae']:np.testing.assert_allclose(p.groupby('species')[m+'_'+metric].mean().mean(),r['results'][m]['species_equal_'+metric],rtol=1e-12)
geo=pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']]
strict=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
area=pd.read_csv(O/'cell_land_overlap.csv')[['cell','area_supported']]
p=p.merge(geo,on='cell',validate='many_to_one').merge(strict,on='cell',how='left',validate='many_to_one').merge(area,on='cell',validate='many_to_one')
islands=p.strict_single.fillna(False)&p.area_supported&p.landmass_area_km2.lt(7e6)
cohorts={'all155':p,'area_supported_islands':p[islands],'islands_without_GB_Ireland':p[islands&~p.landmass_id.isin(['15','26'])]}
if mode=='spatial':cohorts['at_least100km']=p[p.distance_to_training_cell_km.ge(100)]
rows=[]
for cohort,q in cohorts.items():
    sp=sorted(q.species.unique());blocks=sorted(q.spatial_block.unique());rng=np.random.default_rng(20261105)
    counts=pd.crosstab(q.species,q.spatial_block).reindex(index=sp,columns=blocks,fill_value=0).to_numpy()
    sw=rng.multinomial(len(sp),np.full(len(sp),1/len(sp)),size=2000);bw=rng.multinomial(len(blocks),np.full(len(blocks),1/len(blocks)),size=2000)
    den=bw@counts.T;weights=sw*(den>0);assert (weights.sum(axis=1)>0).all()
    for reference in models[:2]:
        delta=q.B6_intermediate_islands_mse-q[reference+'_mse']
        sums=q[['species','spatial_block']].assign(delta=delta).groupby(['species','spatial_block']).delta.sum().unstack(fill_value=0).reindex(index=sp,columns=blocks,fill_value=0).to_numpy()
        means=np.divide(bw@sums.T,den,out=np.zeros(den.shape),where=den>0)
        samples=(means*weights).sum(axis=1)/weights.sum(axis=1)
        low,high=np.quantile(samples,[.025,.975]);point=delta.groupby(q.species).mean().mean()
        old=q.groupby('species')[reference+'_mse'].mean().mean();new=q.groupby('species').B6_intermediate_islands_mse.mean().mean()
        rows.append(dict(mode=mode,cohort=cohort,reference=reference,rows=len(q),species=len(sp),landmasses=q.landmass_id.nunique(),mse_difference=point,improvement_pct=100*(1-new/old),interval_low=low,interval_high=high))
pd.DataFrame(rows).to_csv(O/f'intermediate_{mode}_comparison.csv',index=False)
(O/f'intermediate_{mode}_validation.json').write_text(json.dumps(dict(status='PASS_PAIRED_PREDICTION_REPLAY',mode=mode,prediction_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),checks=['all reference columns exactly match','feature and lock hashes','MSE MAE receipt replay'],limits='2000 fixed-prediction crossed species/spatial-block percentile intervals; not model-refit or causal uncertainty; retrospective geometric proxy.'),indent=2))
print(pd.DataFrame(rows).to_string(index=False))
