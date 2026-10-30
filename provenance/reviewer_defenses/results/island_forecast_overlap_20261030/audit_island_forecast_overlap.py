from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from pyproj import Transformer
R=Path(__file__).resolve().parents[2];O=R/'outputs';W=R/'work/barriers'
files=[O/'island_forecast_overlap_lock.json',W/'original_coordinate_geography.parquet',
       W/'global_cell_landmass_2001_2024.parquet',R/'work/temporal/forecast/predictions_temporal.parquet',
       O/'island_host_network_islands.csv',O/'forecast_candidate_disposition.csv']
lock=json.loads(files[0].read_text())
g=pd.read_parquet(files[1]);t=Transformer.from_crs(4326,6933,always_xy=True)
x,y=t.transform(g.decimalLongitude.values,g.decimalLatitude.values)
g['cell']=np.floor(x/10000).astype(int).astype(str)+np.char.add(':',np.floor(y/10000).astype(int).astype(str))
g['known']=g.landmass_id.notna()
group=g.groupby('cell').agg(coordinate_points=('known','size'),all_known=('known','all'),
                            landmass_count=('landmass_id','nunique'),observed_landmass=('landmass_id','first'))
centres=pd.read_parquet(files[2]);centres['cell']=centres.cell_x_10km.astype(str)+':'+centres.cell_y_10km.astype(str)
assert not centres.cell.duplicated().any()
group=group.merge(centres[['cell','landmass_id']],left_index=True,right_on='cell',how='outer',validate='one_to_one')
group['strict_single']=group.all_known.fillna(False)&group.landmass_count.eq(1)&group.observed_landmass.eq(group.landmass_id)
islands=pd.read_csv(files[4],dtype={'landmass_id':str});target=set(islands.loc[islands.eligible,'landmass_id'])
group['target_centre']=group.landmass_id.isin(target)
group['target_strict']=group.strict_single&group.target_centre
group.to_parquet(W/'forecast_cell_original_island_overlap.parquet',index=False)
p=pd.read_parquet(files[3]);p=p.merge(group[['cell','strict_single','target_centre','target_strict','observed_landmass','landmass_count','all_known']],on='cell',how='left',validate='many_to_one')
for col in ['strict_single','target_centre','target_strict']:p[col]=p[col].fillna(False)
p['M3_squared_error']=(p.target-p.M3_both)**2;p['M5_squared_error']=(p.target-p.M5_host_composition)**2
p['M3_residual']=p.target-p.M3_both;p['M5_residual']=p.target-p.M5_host_composition
rows=[]
for label,q in [('all_existing_test_rows',p),('target51_cell_centre',p[p.target_centre]),('target51_strict_coordinate_agreement',p[p.target_strict])]:
    means=q.groupby('species')[['M3_squared_error','M5_squared_error','M3_residual','M5_residual']].mean().mean()
    row=dict(cohort=label,rows=len(q),cells=q.cell.nunique(),species=q.species.nunique(),islands=q.observed_landmass.nunique())
    row.update(means.to_dict());row['relative_MSE_improvement_pct']=100*(1-row.get('M5_squared_error',np.nan)/row.get('M3_squared_error',np.nan))
    rows.append(row)
s=pd.DataFrame(rows);s.to_csv(O/'island_forecast_overlap_summary.csv',index=False)
eligible=p[p.target_strict]
byisland=eligible.groupby('observed_landmass').agg(rows=('species','size'),cells=('cell','nunique'),species=('species','nunique')).reset_index().rename(columns={'observed_landmass':'landmass_id'})
islands[islands.eligible][['landmass_id']].merge(byisland,on='landmass_id',how='left').to_csv(O/'island_forecast_overlap_all51.csv',index=False)
c=pd.read_csv(files[5])[['species','taxon_key']];assert len(c)==225
c['forecast_rows']=c.species.map(p.groupby('species').size()).fillna(0).astype(int)
c['strict_island_rows']=c.species.map(eligible.groupby('species').size()).fillna(0).astype(int)
c['disposition']=np.where(c.forecast_rows.eq(0),'NOT_IN_EXISTING_TEST',np.where(c.strict_island_rows.eq(0),'NO_STRICT_TARGET_ISLAND_ROWS','HAS_STRICT_TARGET_ISLAND_ROWS'))
c.to_csv(O/'island_forecast_overlap_all225.csv',index=False)
eligible.to_parquet(W/'forecast_predictions_target51_strict.parquet',index=False)
audit=dict(status='GEOGRAPHIC_OVERLAP_AND_EXISTING_MODEL_SUBSET_NOT_NEW_VALIDATION',
           source_coordinate_points=len(g),forecast_rows=len(p),strict_target_rows=len(eligible),
           strict_target_islands=eligible.observed_landmass.nunique(),target51_centre_rows=int(p.target_centre.sum()),
           all225_accounted=True,limits=lock['limits'],
           input_sha256={str(f.relative_to(R)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
(O/'island_forecast_overlap_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(s.to_string(index=False));print(json.dumps({k:v for k,v in audit.items() if k!='input_sha256'}))
