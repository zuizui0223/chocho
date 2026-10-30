from pathlib import Path
import numpy as np,pandas as pd,json
from sklearn.neighbors import BallTree
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';F=R/'work/temporal/forecast'
p=pd.read_parquet(F/'species_cell_period_panel.parquet',columns=['speciesKey','cell','period']).query('period<2')
c=pd.read_parquet(W/'cell_geography.parquet')
h=pd.read_parquet(W/'species_history_landmass.parquet')
h=h[np.isfinite(h.latitude)&np.isfinite(h.longitude)]
result=[]
for (sp,period),q in p.groupby(['speciesKey','period']):
 q=q.merge(c,on='cell',validate='many_to_one')
 a=h[(h.speciesKey==sp)&(h.year<=2014+4*period)].drop_duplicates(['cell_x_10km','cell_y_10km'])
 assert len(a)
 t=BallTree(np.radians(a[['latitude','longitude']]),metric='haversine');x=np.radians(q[['latitude','longitude']])
 q['log_history_distance']=np.log1p(t.query(x,k=1)[0].ravel()*6371.0088)
 for radius in [50,100]:q['log_history_count'+str(radius)]=np.log1p(t.query_radius(x,r=radius/6371.0088,count_only=True))
 q['landmass_unknown']=q.landmass_id.isna().astype(int)
 q['log_landmass_area']=np.log1p(q.landmass_area_km2.fillna(0))
 q['same_landmass_recorded']=0;q['log_same_landmass_distance']=np.log1p(20000.);q['log_same_landmass_count100']=0.
 for lm,idx in q.groupby('landmass_id').groups.items():
  b=a[a.landmass_id==lm]
  if b.empty:continue
  t=BallTree(np.radians(b[['latitude','longitude']]),metric='haversine');x=np.radians(q.loc[idx,['latitude','longitude']])
  q.loc[idx,'same_landmass_recorded']=1
  q.loc[idx,'log_same_landmass_distance']=np.log1p(np.minimum(t.query(x,k=1)[0].ravel()*6371.0088,20000))
  q.loc[idx,'log_same_landmass_count100']=np.log1p(t.query_radius(x,r=100/6371.0088,count_only=True))
 result.append(q)
out=pd.concat(result,ignore_index=True)
assert len(out)==len(p) and not out.duplicated(['speciesKey','cell','period']).any()
out.to_parquet(W/'history_geography_features.parquet',index=False)
print('COMPLETE',out.shape,'species',out.speciesKey.nunique(),flush=True)
