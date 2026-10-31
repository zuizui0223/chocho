from pathlib import Path
import json,sys
import numpy as np,pandas as pd,joblib
R=Path(__file__).resolve().parents[2];O=R/'outputs';F=R/'work/temporal/forecast';W=R/'work/barriers'
mode=sys.argv[1]
receipt=json.loads((O/f'forecast_joint_{mode}_results.json').read_text())
assert receipt['status']=='COMPLETE'
p=pd.read_parquet(F/f'predictions_joint_{mode}.parquet')
models=['M5_host_composition','B5_sea_terrain','T1_port_connectivity','J1_natural_transport']
for prefix,column in [('',models[0]),('sea_',models[1]),('transport_',models[2])]:
    if receipt.get('sensitivity') and prefix=='sea_':continue
    old=pd.read_parquet(F/f'predictions_{prefix}{mode}.parquet')
    np.testing.assert_array_equal(p[['speciesKey','cell','target']],old[['speciesKey','cell','target']])
    np.testing.assert_array_equal(p[column],old[column])
natural=joblib.load(F/'model_B5_sea_terrain.joblib')['features']
transport=joblib.load(F/'model_T1_port_connectivity.joblib')['features']
assert set(receipt['fitted_features'])==set(natural)|set(transport)
assert len(receipt['fitted_features'])==len(set(receipt['fitted_features']))
if receipt.get('sensitivity'):
    assert set(receipt['fitted_model_features']['B5_sea_terrain'])==set(natural)
    assert set(receipt['fitted_model_features']['J1_natural_transport'])==set(natural)|set(transport)
geo=pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']]
strict=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
p=p.merge(geo,on='cell',validate='many_to_one').merge(strict,on='cell',how='left',validate='many_to_one')
island=p.strict_single.fillna(False)&p.landmass_area_km2.lt(7e6)
cohorts={'all155':np.ones(len(p),bool),'strict_islands':island,'strict_islands_without_GB_Ireland':island&~p.landmass_id.isin(['15','26'])}
if mode=='spatial':cohorts['at_least_100km']=p.distance_to_training_cell_km.ge(100)
for m in models:
    assert np.isfinite(p[m]).all();p[m+'_sq']=(p[m]-p.target)**2
rows=[]
for name,mask in cohorts.items():
    q=p[mask];scores=q.groupby('species')[[m+'_sq' for m in models]].mean().mean()
    for ref in models[:-1]:
        rows.append(dict(cohort=name,reference=ref,rows=len(q),species=q.species.nunique(),landmasses=q.landmass_id.nunique(),
                         reference_mse=scores[ref+'_sq'],joint_mse=scores['J1_natural_transport_sq'],
                         improvement_pct=100*(1-scores['J1_natural_transport_sq']/scores[ref+'_sq'])))
    if name=='all155':
        for m in models:np.testing.assert_allclose(scores[m+'_sq'],receipt['results'][m]['species_equal_mse'])
summary=pd.DataFrame(rows);summary.to_csv(O/f'joint_{mode}_comparison.csv',index=False)
(O/f'joint_{mode}_validation.json').write_text(json.dumps(dict(status='PASS_JOINT_COMPARISON',mode=mode,
    checks=['all reference predictions keys and targets exactly match original runs','joint feature set is exact union of saved B5 and T1 features','independent species-equal MSE matches receipt'],
    limits='Fixed settings, retrospective heldout prediction. No causal identification, no uncertainty estimated here.'),indent=2))
print(summary.to_string(index=False))
