"""Extract auditable positive traffic evidence from published regional files."""
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import shapefile
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

base=Path(sys.argv[1])
w=base/"work/barriers"
cache=w/"regional_ferry_extension"
calls=w/"multiyear_calls"
out=Path(__file__).resolve().parents[1]/"results/regional_maritime_links"
out.mkdir(parents=True,exist_ok=True)
# Read the frozen URL crosswalk if rebuilding outside the acquisition workspace.
if (cache/"urls.json").exists():
    urls=json.loads((cache/"urls.json").read_text())
else:
    urls={Path(r.file).stem:r.url for r in pd.read_csv(out/"sources.csv").itertuples() if r.file.endswith(".pdf")}
records=[]; receipts=[]
def receipt(path,url):
    receipts.append(dict(file=path.name,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
def leg(source,year,origin,destination,a,b,n,unit,twoway,raw):
    assert n>0
    records.append(dict(source=source,year=year,origin=origin,destination=destination,a=a,b=b,
                        positive_count=n,unit=unit,two_way=twoway,source_file=raw))
def pdf(key):
    receipt(cache/(key+".pdf"),urls[key])
    import pymupdf
    with pymupdf.open(cache/(key+".pdf")) as document:
        return "\n".join(page.get_text(sort=True) for page in document)
# Use positive evidence in a known subperiod; counts are not annualized.
for year in range(2019,2024):
    key=f"bc{year}"
    txt=pdf(key)
    for route,code,a,b,origin,dest in [("06","CROF",50,None,"Crofton","Vesuvius Bay"),
                                       ("11","PR",None,86,"Prince Rupert","Skidegate")]:
        row=re.search(r"(?m)^\s*"+route+r"\s+"+code+r"\s+([\d,]+)\s",txt)
        assert row,(key,route)
        leg("BC Ferries",year,origin,dest,a,b,int(row[1].replace(",","")),"December vehicles at named terminal",True,key+".pdf")
    key=f"ws{year}";txt=pdf(key)
    for origin,dest,b in [("Mukilteo","Clinton",551),("Fauntleroy","Vashon",1531)]:
        row=re.search(origin+r"-\s+([\d,]+)[^\n]*\n\s*"+dest,txt)
        if not row:
            row=re.search(origin+r"-\s+"+dest+r"\s+([\d,]+)",txt)
        assert row,(origin,key)
        leg("Washington State Ferries",year,origin,dest,None,b,int(row[1].replace(",","")),"first reported quarter riders",True,key+".pdf")
txt=pdf("hk")
for year,n in zip(range(2019,2024),[632491,638641,366653,412471,551395]):
    assert f"{n:,}" in txt
    leg("MLIT Hokkaido FY ending March",year,"Hakodate","Aomori",27,13,n,"fiscal-year passengers",True,"hk.pdf")

# Completed inbound calls at all five Balearic ports. Restrict island assignments to
# unique UN/LOCODE coordinates in a study polygon or within 5 km of a same-country WPI port.
ports=pd.read_csv(base/"outputs/port_buffer_landmass_audit.csv").query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1")
f={x["properties"]["OBJECTID"]:x for x in json.loads((w/"ports_2019.geojson").read_text(encoding="utf8"))["features"]}
uk={int(r.single_landmass_id) for r in ports.itertuples() if f[r.port_id]["properties"]["COUNTRY"] in ("GB","IE")}
universe=set(ports.single_landmass_id.astype(int))-uk
assert len(universe)==200
pts=[(f[r.port_id]["properties"]["COUNTRY"],int(r.single_landmass_id),*f[r.port_id]["geometry"]["coordinates"][:2]) for r in ports.itertuples() if r.single_landmass_id in universe]
wk=pd.DataFrame(pts,columns=["country","island","lon","lat"])
reader=shapefile.Reader(str(w/"GSHHS_shp/h/GSHHS_h_L1.shp"))
geoms=[]; ids=[]
for r in reader.iterRecords(fields=["sibling_id"]):
    if int(r[0]) in universe:
        ids.append(int(r[0]));geoms.append(shape(reader.shape(r.oid).__geo_interface__))
tree=STRtree(geoms)
loc=pd.read_csv(calls/"locodes.csv",dtype=str,keep_default_na=False)
loc=loc[loc.Function.str.startswith("1")].copy()
loc["code"]=loc.Country+loc.Location
coords=loc.groupby("code").Coordinates.agg(set).to_dict()
receipt(calls/"locodes.csv","https://raw.githubusercontent.com/datasets/un-locode/0fb7229b64476966c67cb4b90f000d554154591c/data/code-list.csv")
@lru_cache(None)
def resolve(code):
    cs=coords.get(code,set())
    if len(cs)!=1:return None
    m=re.fullmatch(r"(\d{2})(\d{2})([NS])\s+(\d{3})(\d{2})([EW])",next(iter(cs)).strip())
    if not m:return None
    a,b,h,c,d,k=m.groups()
    lat=(int(a)+int(b)/60)*(1 if h=="N" else -1)
    lon=(int(c)+int(d)/60)*(1 if k=="E" else -1)
    hits={ids[i] for i in tree.query(Point(lon,lat),predicate="intersects")}
    p=wk[wk.country.eq(code[:2])]
    arc=np.sin(np.radians(p.lat-lat)/2)**2+np.cos(np.radians(lat))*np.cos(np.radians(p.lat))*np.sin(np.radians(p.lon-lon)/2)**2
    hits.update(p.loc[12742.0176*np.arcsin(np.sqrt(arc.clip(0,1)))<=5,"island"])
    return hits.pop() if len(hits)==1 else None
counters=Counter();seen=set()
for path in sorted(calls.glob("atraques_*.csv")):
    if not "201901"<=path.stem[-6:]<="202401":continue
    receipt(path,"https://posidoniaweb.portsdebalears.com/docs/files/"+path.name)
    with path.open(encoding="utf-8-sig",errors="replace",newline="") as ff:
        for r in csv.DictReader(ff,delimiter=";"):
            origin=r["Locode origen"].strip()
            port=r["Código puerto"]
            if port not in ("A","P","I","F","M"):continue
            yy=re.findall(r"\b20\d\d\b",r["Fecha/Hora llegada"])
            if len(yy)!=1 or not 2019<=int(yy[0])<=2023:continue
            if not r["Fecha/Hora salida"].strip() or not r["Escala"].strip():continue
            key=(port,r["Escala"])
            if key in seen:continue
            seen.add(key)
            a=resolve(origin);b=152 if port in ("A","P") else None
            if a is not None or b is not None:
                counters[int(yy[0]),origin,port,a,b]+=1
for (y,orig,port,a,b),n in counters.items():
    leg("Balearic completed calls",y,orig,port,a,b,n,"distinct completed calls",False,"atraques monthly archive")
pd.DataFrame(records).to_csv(out/"positive_routes.csv",index=False)
pd.DataFrame(receipts).to_csv(out/"sources.csv",index=False)
print(pd.DataFrame(records).groupby("source").agg({"year":"nunique","origin":"size"}))
