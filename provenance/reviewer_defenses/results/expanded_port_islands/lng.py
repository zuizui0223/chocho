from pathlib import Path
import hashlib,json
import sys
import requests,pandas as pd,numpy as np,geopandas as gpd, shapely
from pyproj import Geod

out=Path(__file__).resolve().parent
last=int(sys.argv[1]) if len(sys.argv)>1 else 2021
assert last in [2021,2022,2023]
if last!=2021:
    out=out/f"through_{last}"
    out.mkdir(exist_ok=True)
w=Path("work/barriers")
cache=w/"lng_t3"
cache.mkdir(exist_ok=True)
meta=requests.get("https://zenodo.org/api/records/19571058",timeout=30).json()
hashes={}
for f in meta["files"]:
    if f["key"] not in ["LNG_terminal.csv","LNG_tanker_voyage.csv"]:continue
    path=cache/f["key"]
    if not path.exists():
        r=requests.get(f["links"]["self"],timeout=60);r.raise_for_status();path.write_bytes(r.content)
    md5=hashlib.md5(path.read_bytes()).hexdigest()
    assert f["checksum"]=="md5:"+md5
    hashes[path.name]=md5
v=pd.read_csv(cache/"LNG_tanker_voyage.csv").drop_duplicates()
v=v[v.start_date.ge("2020-01-01")&v.end_date.le(f"{last}-12-31")&v.end_date.ge(v.start_date)].copy()
t=pd.read_csv(cache/"LNG_terminal.csv")
t=t[t.name.isin(set(v.from_terminal)|set(v.to_terminal))][["name","lon","lat"]].drop_duplicates()
land=gpd.read_file(w/"GSHHS_shp/h/GSHHS_h_L1.shp")
tree=shapely.STRtree(land.geometry.to_numpy())
ids=land.sibling_id.astype(int).to_numpy()
geod=Geod(ellps="WGS84")
def lookup(lon,lat):
    point=shapely.Point(lon,lat)
    hits=set(ids[tree.query(point,predicate="intersects")])
    if not hits:
        ring=[geod.fwd(lon,lat,a,5000)[:2] for a in range(0,360,10)]
        if np.ptp([x for x,y in ring])>180:return np.nan
        hits=set(ids[tree.query(shapely.Polygon(ring),predicate="intersects")])
    return hits.pop() if len(hits)==1 else np.nan
t["island"]=[lookup(r.lon,r.lat) for r in t.itertuples()]
geo=t.groupby("name").island.agg(lambda x:x.iloc[0] if x.notna().all() and x.nunique()==1 else np.nan)
v["a"]=v.from_terminal.map(geo)
v["b"]=v.to_terminal.map(geo)
v["year"]=v.end_date.str[:4].astype(int)
study=set(pd.read_csv("outputs/island_port_early_inventory_status.csv").landmass_id)-{15,26,456}
# Independent endpoint coverage even when the other end lies outside the study.
valid=v[v.a.notna()&v.b.notna()&v.a.ne(v.b)].copy()
nodes=[]
for direction,col in [('out','a'),('in','b')]:
    for (i,y),g in v[v[col].isin(study)&v.from_terminal.ne(v.to_terminal)].groupby([col,'year']):
        nodes.append(dict(island=int(i),year=y,direction=direction,voyages=len(g)))
pd.DataFrame(nodes).to_csv(out/"endpoint_coverage.csv",index=False)
e=valid[valid.a.isin(study)&valid.b.isin(study)].groupby(["a","b","year"]).size().rename("voyages").reset_index()
e.to_csv(out/"links.csv",index=False)
t[t.island.isin(study)].to_csv(out/"terminal_join.csv",index=False)
receipt=dict(doi=meta['doi'],license=meta['metadata']['license']['id'],checksums=hashes,unique_global_voyages_in_window=len(v),mapped_global_different_landmass=len(valid),study_voyages=int(e.voyages.sum()),directed_study_pairs=len(e[['a','b']].drop_duplicates()),last_year=last,vessel_scope='LNG only; includes export and return')
(out/"source.json").write_text(json.dumps(receipt,indent=2))
print(receipt);print(pd.DataFrame(nodes).groupby('direction').island.nunique())
