from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,geopandas as gpd,shapely
from pyproj import Geod
from sklearn.neighbors import BallTree
R=Path(__file__).resolve().parents[2];W=R/'work/barriers';O=R/'outputs';geod=Geod(ellps='WGS84')
p=pd.read_parquet(W/'coastal_control_design.parquet')
assert not p.duplicated(['port_id','landmass_id','control_number']).any()
ports=json.loads((W/'ports_2019.geojson').read_text())['features']
xy=np.array([v['geometry']['coordinates'][:2] for v in ports]);tree=BallTree(np.radians(xy[:,::-1]),metric='haversine')
nearest=tree.query(np.radians(p[['latitude','longitude']]),k=1)[0][:,0]*6371.0088
assert np.allclose(nearest,p.nearest_any_port_km)
assert (nearest>=5).all()
_,_,offset=geod.inv(p.shoreline_anchor_longitude,p.shoreline_anchor_latitude,p.longitude,p.latitude)
assert np.allclose(offset,250,atol=1e-5)
land=gpd.read_file(next(W.rglob('GSHHS_h_L1.shp')));land['landmass_id']=land.sibling_id.astype(str)
geom=land.geometry.to_numpy();shapely.prepare(geom);ltree=shapely.STRtree(geom)
points=shapely.points(p.longitude,p.latitude);pairs=ltree.query(points)
pairs=pairs[:,shapely.contains(geom[pairs[1]],points[pairs[0]])]
matched=set(zip(pairs[0],land.landmass_id.to_numpy()[pairs[1]]))
assert all((i,str(lid)) in matched for i,lid in enumerate(p.landmass_id))
# Original generator checked the shoreline anchor against prior offset points.
# Enforce spacing on the actual sampled points, before any vegetation is read.
keep=[]
for _,q in p.groupby(['port_id','landmass_id'],sort=False):
    accepted=[]
    for index,row in q.iterrows():
        if accepted:
            prev=p.loc[accepted]
            _,_,distance=geod.inv(np.full(len(prev),row.longitude),np.full(len(prev),row.latitude),prev.longitude.to_numpy(),prev.latitude.to_numpy())
            if np.min(distance)<2000:continue
        accepted.append(index)
    keep.extend(accepted)
q=p.loc[keep].copy().reset_index(drop=True);q['sample_id']=np.arange(len(q))
q.to_parquet(W/'coastal_control_validated_design.parquet',index=False)
receipt=dict(status='GEOMETRY_VALIDATED_NOT_VEGETATION_RESULT',original_points=len(p),retained_points=len(q),
    spacing_removed=len(p)-len(q),ports=int(q.port_id.nunique()),islands=int(q.landmass_id.nunique()),
    pairs_with3plus=int((q.groupby(['port_id','landmass_id']).size()>=3).sum()),
    input_sha256=hashlib.sha256((W/'coastal_control_design.parquet').read_bytes()).hexdigest(),
    checks=['all control points strictly inside claimed original-coordinate landmass',
            'all anchors offset250m geodesically','all points >=5km from catalogued ports (haversine)',
            'retain >=2km geodesic spacing within each port-island group before reading vegetation'],
    limitations=['geometry resolution and position uncertainty remain','not random shoreline sampling',
                 'near-port probe and offshore control shoreline distances are not exactly matched'])
(O/'coastal_control_design_validation.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt),flush=True)
