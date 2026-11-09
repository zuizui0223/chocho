from pathlib import Path
import json,hashlib
import requests,pandas as pd, numpy as np, geopandas as gpd, shapely
from pyproj import Geod
from functools import lru_cache
out=Path(__file__).resolve().parent
w=Path('work/barriers')
cache=w/'ncfo';cache.mkdir(exist_ok=True)
geo=pd.read_csv('outputs/port_buffer_landmass_audit.csv').query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1")
study=set(geo.single_landmass_id.astype(int))
land=gpd.read_file(w/'GSHHS_shp/h/GSHHS_h_L1.shp')
tree=shapely.STRtree(land.geometry.to_numpy());ids=land.sibling_id.astype(int).to_numpy()
geod=Geod(ellps='WGS84')
@lru_cache(None)
def locate(lon,lat):
    if not np.isfinite(lon+lat) or abs(lon)>180 or abs(lat)>90:return None
    s=set(ids[tree.query(shapely.Point(lon,lat),predicate='intersects')])
    if not s:
        ring=[geod.fwd(lon,lat,a,1000)[:2] for a in range(0,360,15)]
        if np.ptp([x for x,y in ring])>180:return None
        s=set(ids[tree.query(shapely.Polygon(ring),predicate='intersects')])
    return s.pop() if len(s)==1 else None
receipts={};rows=[]
def load(code):
    url='https://data.bts.gov/resource/'+code+'.json?$limit=50000'
    f=cache/(code+'.json')
    if not f.exists():
        res=requests.get(url,timeout=60);res.raise_for_status();f.write_bytes(res.content)
    t=pd.DataFrame(json.loads(f.read_bytes()));assert 0<len(t)<50000
    receipts[code]=dict(url=url,rows=len(t),sha256=hashlib.sha256(f.read_bytes()).hexdigest())
    return t
for y,seg,term,op in [(2019,'gn77-pp24','ca7h-i9yt','ggca-ddee'),(2022,'iyde-f8k7','ke6h-ga46','as3z-m8rd')]:
    t=load(term);s=load(seg);v=load(op)
    assert t.data_year.astype(int).eq(y).all() and v.data_year.astype(int).eq(y).all()
    for c in ['longitude','latitude']:t[c]=pd.to_numeric(t[c],errors='coerce')
    t['island']=[locate(x,y0) for x,y0 in zip(t.longitude,t.latitude)]
    ok=t.groupby('terminal_id').island.agg(lambda x:x.iloc[0] if x.notna().all() and x.nunique()==1 else np.nan)
    good=set(v.loc[pd.to_numeric(v.trips_per_year,errors='coerce').gt(0),'segment_id'])
    s=s[s.segment_id.isin(good)].drop_duplicates('segment_id')
    for r in s.itertuples():
        if len(str(r.segment_name).split(' - '))>2:continue
        a=ok.get(r.seg_terminal1_id);b=ok.get(r.seg_terminal2_id)
        if a==b or (a not in study and b not in study):continue
        rows.append(dict(year=y,segment=r.segment_id,name=r.segment_name,a=a,b=b))
pd.DataFrame(rows).to_csv(out/'ncfo_segments.csv',index=False)
(out/'ncfo_sources.json').write_text(json.dumps(receipts,indent=2))
print(pd.DataFrame(rows).groupby('year').size())
