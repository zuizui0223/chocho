from pathlib import Path
import json,sys
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers';F=R/'work/temporal/forecast'
mode=sys.argv[1] if len(sys.argv)>1 else 'temporal'
receipt=json.loads((O/f'forecast_transport_{mode}_results.json').read_text());assert receipt['status']=='COMPLETE'
p=pd.read_parquet(F/f'predictions_transport_{mode}.parquet');old=pd.read_parquet(F/f'predictions_{mode}.parquet')
np.testing.assert_array_equal(p[['speciesKey','cell','target']],old[['speciesKey','cell','target']])
np.testing.assert_allclose(p.M5_host_composition,old.M5_host_composition)
geo=pd.read_parquet(W/'cell_geography.parquet')[['cell','landmass_id','landmass_area_km2']]
strict=pd.read_parquet(W/'forecast_cell_original_island_overlap.parquet')[['cell','strict_single']]
p=p.merge(geo,on='cell',validate='many_to_one').merge(strict,on='cell',how='left',validate='many_to_one')
p.strict_single=p.strict_single.fillna(False)
models=['M5_host_composition','T0_port_location','T1_port_connectivity']
for m in models:
    assert np.isfinite(p[m]).all();p[m+'_sq']=(p[m]-p.target)**2
island=p.strict_single&p.landmass_area_km2.lt(7e6)
cohorts={'all155':np.ones(len(p),bool),'strict_islands':island,'strict_islands_without_GB_Ireland':island&~p.landmass_id.isin(['15','26']),
         'strict_continental':p.strict_single&p.landmass_area_km2.ge(7e6),'initially_unrecorded':p.baseline_probability.eq(0)}
rows=[]
for name,mask in cohorts.items():
    q=p[mask];by=q.groupby('species')[[m+'_sq' for m in models]].mean();scores=by.mean()
    for m in models:
        rows.append(dict(cohort=name,model=m,rows=len(q),species=len(by),cells=q.cell.nunique(),landmasses=q.landmass_id.nunique(),
                         species_equal_mse=scores[m+'_sq'],improvement_vs_M5_pct=100*(1-scores[m+'_sq']/scores['M5_host_composition_sq'])))
    by.to_csv(O/f'transport_{mode}_{name}_species_errors.csv')
s=pd.DataFrame(rows);s.to_csv(O/f'transport_{mode}_comparison.csv',index=False)
overall=s[s.cohort.eq('all155')].set_index('model')
for m in models:np.testing.assert_allclose(overall.loc[m,'species_equal_mse'],receipt['results'][m]['species_equal_mse'])
v=dict(status='PASS_EXISTING_FORECAST_EXTENSION',mode=mode,rows=len(p),species=p.species.nunique(),
       t1_vs_t0_improvement_pct=100*(1-overall.loc['T1_port_connectivity','species_equal_mse']/overall.loc['T0_port_location','species_equal_mse']),
       checks=['identical heldout keys targets and M5 predictions','all finite predictions','independent error aggregation matches fit receipt'],
       limits='Retrospective2019 port geocoding; species/space external comparisons required separately. Subgroups are not new model fits.')
(O/f'transport_{mode}_validation.json').write_text(json.dumps(v,indent=2),encoding='utf-8')
print(s.to_string(index=False));print(json.dumps(v))
