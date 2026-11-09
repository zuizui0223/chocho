from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,geopandas as gpd,shapely
from pyproj import Geod
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs';D=W/'route_terrain_design_parts';D.mkdir(exist_ok=True)
lock=O/'route_terrain_design_lock.json';assert lock.exists()
mapping=W/'sea_route_mapping.parquet';routes=pd.read_parquet(mapping).drop_duplicates('route_key').sort_values('route_key').reset_index(drop=True)
assert len(routes)==175625
routes['route_index']=np.arange(len(routes));routes.to_parquet(W/'route_terrain_routes.parquet',index=False)
land=gpd.read_file(next(W.rglob('GSHHS_h_L1.shp')));geoms=land.geometry.to_numpy();shapely.prepare(geoms);tree=shapely.STRtree(geoms);geod=Geod(ellps='WGS84')
support=[];tiles={}
def tile(x,y):
    lon,lat=int(np.floor(x)),int(np.floor(y))
    return f'Copernicus_DSM_COG_30_{"N" if lat>=0 else "S"}{abs(lat):02d}_00_{"E" if lon>=0 else "W"}{abs(lon):03d}_00_DEM'
for start in range(0,len(routes),1000):
    q=routes.iloc[start:start+1000];xs=[];ys=[];indices=[];steps=[]
    for r in q.itertuples():
        az,_,distance=geod.inv(r.source_lon,r.source_lat,r.longitude,r.latitude)
        n=max(1,int(np.ceil(distance/2000)));step=distance/n
        x,y,_=geod.fwd(np.full(n,r.source_lon),np.full(n,r.source_lat),np.full(n,az),(np.arange(n)+.5)*step)
        xs.extend(x);ys.extend(y);indices.extend([r.route_index]*n);steps.append((r.route_index,n,step/1000))
    pts=shapely.points(np.column_stack([xs,ys]));pairs=tree.query(pts);hits=pairs[:,shapely.intersects(geoms[pairs[1]],pts[pairs[0]])]
    is_land=np.zeros(len(pts),bool);is_land[hits[0]]=True
    points=pd.DataFrame(dict(route_index=indices,longitude=xs,latitude=ys));points=points[is_land].copy()
    points['tile']=[tile(x,y) for x,y in zip(points.longitude,points.latitude)]
    points.to_parquet(D/f'part_{start:06d}.parquet',index=False)
    counts=points.groupby('route_index').size()
    for idx,n,step in steps:support.append(dict(route_index=idx,planned_samples=n,land_samples=int(counts.get(idx,0)),sample_step_km=step))
    for name,count in points.tile.value_counts().items():tiles[name]=tiles.get(name,0)+int(count)
    if start%10000==0:print('route terrain design',start+len(q),'/',len(routes),flush=True)
pd.DataFrame(support).to_parquet(W/'route_terrain_sampling_support.parquet',index=False)
tile_rows=[]
for name,count in sorted(tiles.items()):
    old=W/'dem_samples'/(name+'.json')
    status=json.loads(old.read_text())['status'] if old.exists() else 'NOT_PREVIOUSLY_ACCESSED'
    tile_rows.append(dict(tile=name,route_land_sample_points=count,previous_cell_sample_tile_status=status))
pd.DataFrame(tile_rows).to_csv(O/'route_terrain_required_tiles.csv',index=False)
audit=dict(status='ROUTE_SAMPLE_DESIGN_COMPLETE_ELEVATIONS_NOT_ACQUIRED',routes=len(routes),planned_samples=sum(x['planned_samples'] for x in support),planned_land_samples=sum(x['land_samples'] for x in support),required_tiles=len(tiles),previous_tile_status_counts=pd.Series([x['previous_cell_sample_tile_status'] for x in tile_rows]).value_counts().to_dict(),max_step_km=max(x['sample_step_km'] for x in support),source_mapping_sha256=hashlib.sha256(mapping.read_bytes()).hexdigest(),lock_sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),limits='Tile metadata reuse is not route-elevation acquisition. No route terrain model fitted. Land mask and2km sampling are approximate.')
(O/'route_terrain_design_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit),flush=True)
