from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,geopandas as gpd
from sklearn.neighbors import BallTree
from pyproj import Transformer
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs';F=R/'work/temporal/forecast'
land=gpd.read_file(next(W.rglob('GSHHS_h_L1.shp')))
land['landmass_id']=land.sibling_id.astype(str)
def join(q):
 g=gpd.GeoDataFrame(q.copy(),geometry=gpd.points_from_xy(q.longitude,q.latitude),crs=4326)
 j=gpd.sjoin(g,land[['landmass_id','area','geometry']],how='left',predicate='intersects')
 assert not j.index.duplicated().any(), 'Overlapping landmass assignments'
 return pd.DataFrame(j.drop(columns=['geometry','index_right'])).rename(columns={'area':'landmass_area_km2'})
c=join(pd.read_parquet(F/'candidate_cells.parquet',columns=['cell','latitude','longitude']))
ports=json.loads((W/'ports_2019.geojson').read_text())['features']
xy=np.array([p['geometry']['coordinates'][:2] for p in ports],dtype=float)
assert np.isfinite(xy).all() and (abs(xy[:,1])<=90).all()
tree=BallTree(np.radians(xy[:,::-1]),metric='haversine')
c['distance_port_2019_km']=tree.query(np.radians(c[['latitude','longitude']]),k=1)[0].ravel()*6371.0088
c.to_parquet(W/'cell_geography.parquet',index=False)
c.to_csv(O/'island_port_candidate_cell_audit.csv',index=False)
p=pd.read_parquet(R/'work/temporal/species_year_cells_10km.parquet',columns=['speciesKey','year','cell_x_10km','cell_y_10km'])
p=p[p.year<=2018].copy()
h=p[['cell_x_10km','cell_y_10km']].drop_duplicates().reset_index(drop=True)
t=Transformer.from_crs(6933,4326,always_xy=True)
h['longitude'],h['latitude']=t.transform((h.cell_x_10km.to_numpy()+.5)*10000,(h.cell_y_10km.to_numpy()+.5)*10000)
valid=np.isfinite(h.longitude)&np.isfinite(h.latitude)
bad=int((~valid).sum());h=join(h[valid].copy())
h.to_parquet(W/'history_cell_landmass.parquet',index=False)
p=p.merge(h,on=['cell_x_10km','cell_y_10km'],how='left',validate='many_to_one')
p.to_parquet(W/'species_history_landmass.parquet',index=False)
z=c[c.landmass_area_km2<7e6]
receipt={'status':'GEOGRAPHY_JOIN_COMPLETE_NOT_MODEL_RESULT','candidate_cells':len(c),'landmass_assigned_cells':int(c.landmass_id.notna().sum()),'unassigned_cell_centres':int(c.landmass_id.isna().sum()),'noncontinental_cells_area_below_7million_km2':len(z),'noncontinental_landmasses':int(z.landmass_id.nunique()),'port_records':len(ports),'history_cells':len(h),'invalid_history_centres':bad,'shoreline_sha256':hashlib.sha256((W/'gshhg-shp-2.3.7.zip').read_bytes()).hexdigest(),'shoreline_source':'https://github.com/GenericMappingTools/gshhg-gmt/releases','port_source':'https://gis.ngdc.noaa.gov/arcgis/rest/services/nccos/PIRO_DigitalAtlas/MapServer/30','limitations':['GSHHG high resolution level1 landmass polygons; lakes not subtracted; no ocean route inferred.','Centres outside coastlines remain unknown, not absences or islands.','Area7million threshold only descriptive mainland/island partition; model uses continuous landmass area.','2019 ports are retrospective proxies, ineligible as2011/2015 forecast inputs.','Static shoreline not annual historical coast reconstruction.']}
(O/'island_port_source_audit.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2),flush=True)
