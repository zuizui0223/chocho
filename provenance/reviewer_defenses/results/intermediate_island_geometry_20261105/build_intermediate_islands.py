from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,geopandas as gpd,shapely
from pyproj import Geod
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs'
lock=O/'intermediate_island_geometry_lock.json'
assert lock.exists()
mapping=pd.read_parquet(W/'sea_route_mapping.parquet')
metrics=pd.read_parquet(W/'sea_route_metrics.parquet')
routes=mapping.drop_duplicates('route_key').merge(metrics,on='route_key',validate='one_to_one')
assert len(routes)==175625 and metrics.sample_step_km.max()<=.2500001
land=gpd.read_file(next(W.rglob('GSHHS_h_L1.shp')))
geoms=land.geometry.to_numpy();ids=land.sibling_id.astype(str).to_numpy();areas=land.groupby(land.sibling_id.astype(str)).area.max().to_dict()
shapely.prepare(geoms);tree=shapely.STRtree(geoms);geod=Geod(ellps='WGS84')
rows=[]
for r in routes[routes.route_endpoint_unknown.eq(1)].itertuples():
    rows.append(dict(route_key=r.route_key,intermediate_islands=np.nan,intermediate_landmasses=np.nan,intermediate_identity_unknown=1,method='UNKNOWN_ENDPOINT'))
for r in routes[routes.route_endpoint_unknown.eq(0)&routes.sampled_marine_segments.lt(2)].itertuples():
    rows.append(dict(route_key=r.route_key,intermediate_islands=0,intermediate_landmasses=0,intermediate_identity_unknown=0,method='NO_SAMPLED_SEA_LAND_SEA_INTERRUPTION'))
selected=routes[routes.route_endpoint_unknown.eq(0)&routes.sampled_marine_segments.ge(2)]
for i,r in enumerate(selected.itertuples()):
    az,_,metres=geod.inv(r.source_lon,r.source_lat,r.longitude,r.latitude)
    n=max(1,int(np.ceil(metres/250)));step=metres/1000/n
    xs,ys,_=geod.fwd(np.full(n,r.source_lon),np.full(n,r.source_lat),np.full(n,az),(np.arange(n)+.5)*metres/n)
    pts=shapely.points(np.column_stack([xs,ys]));pairs=tree.query(pts)
    hits=pairs[:,shapely.intersects(geoms[pairs[1]],pts[pairs[0]])]
    sea=np.ones(n,dtype=bool);sea[hits[0]]=False
    segments=int((np.diff(np.r_[False,sea,False].astype(int))==1).sum())
    assert segments==r.sampled_marine_segments
    np.testing.assert_allclose(sea.sum()*step,r.sampled_marine_km,atol=1e-8)
    seaidx=np.flatnonzero(sea);assert len(seaidx)>0
    inside=(hits[0]>seaidx[0])&(hits[0]<seaidx[-1])
    h=pd.DataFrame({'point':hits[0,inside],'id':ids[hits[1,inside]]}).drop_duplicates()
    ambiguous=bool((h.groupby('point').id.nunique()>1).any())
    intermediate=set(h.id)-{str(r.source_landmass),str(r.landmass_id)}
    rows.append(dict(route_key=r.route_key,intermediate_islands=np.nan if ambiguous else sum(areas[x]<7e6 for x in intermediate),intermediate_landmasses=np.nan if ambiguous else len(intermediate),intermediate_identity_unknown=int(ambiguous),method='SAMPLED_INTERIOR_LANDMASS_IDS'))
    if (i+1)%1000==0:print('identity routes',i+1,'/',len(selected),flush=True)
out=pd.DataFrame(rows);assert len(out)==len(routes) and out.route_key.is_unique
out.to_parquet(W/'intermediate_island_route_metrics.parquet',index=False)
features=mapping[['speciesKey','period','cell','route_key']].merge(out,on='route_key',validate='many_to_one')
features['log_intermediate_islands']=np.log1p(features.intermediate_islands.fillna(0))
features.to_parquet(W/'intermediate_island_features.parquet',index=False)
audit=dict(status='GEOMETRY_EXTRACTION_COMPLETE_NOT_MODEL_OR_RELAY_RESULT',routes=len(out),detailed_routes=len(selected),mapped_rows=len(features),route_unknown=int(out.intermediate_identity_unknown.sum()),routes_with_intermediate_island=int(out.intermediate_islands.gt(0).sum()),source_hashes={name:hashlib.sha256((W/name).read_bytes()).hexdigest() for name in ['sea_route_mapping.parquet','sea_route_metrics.parquet']},lock_sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),limits='Sampled geographic potential only. All225 candidates preserved; no butterfly response or predictive comparison fitted here.')
(O/'intermediate_island_geometry_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit),flush=True)
