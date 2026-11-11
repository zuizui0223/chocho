"""Search the fixed island panel for baseline hosts, immature butterflies and later island records."""
import argparse
import hashlib
import json
import math
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import shapely
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from acquire_route_free_port_flora import ROOT, OUT as SOURCE, GEOD, load_frame

OUT=ROOT/"results/island_breeding_relay"
LOCAL=threading.local()
LAY=re.compile(r"\b(ovipositing|laying (?:an? |her |its )?eggs?|(?:observed|seen|photographed) (?:to )?oviposit(?:ing)?|(?:female|butterfly) (?:was )?(?:observed )?(?:laid eggs?|oviposited))\b",re.I)
NEG=re.compile(r"\b(no|not|possibly|probably|possible|unsure|attempting|appeared|seemed|potential|unsuccessful)\b",re.I)


def frame(workspace):
    ports,a=load_frame(workspace)
    islands=pd.read_csv(SOURCE/"island_registry.csv")
    centers={}
    for i,g in ports[ports.island.isin(islands.island)].groupby("island"):
        lon,lat=np.radians(g[["lon","lat"]].to_numpy()).T
        x,y,z=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]).mean(axis=1)
        centers[i]=np.degrees([np.arctan2(y,x),np.arctan2(z,np.hypot(x,y))])
    migration=pd.read_csv(ROOT/"data/exploratory/public51_author_migrant_20261109/migration_gbif_author2025_crosswalk.csv")
    keys=migration[migration.iloc[:,2].str.startswith("REVIEW_SELECTED")].speciesKey
    assert len(keys)==111 and len(centers)==125
    cells=pd.read_csv(SOURCE/"species_island_hosts.csv")
    cells=cells[cells.key.isin(keys)&(cells.native+cells.introduced).gt(0)]
    a=a[a.island.isin(centers)&a.speciesKey.isin(keys)].copy()
    first=a.groupby(["speciesKey","island"],as_index=False).year.min().rename(columns={"speciesKey":"key","island":"destination","year":"first_year"})
    dist=[]
    for i,xy in centers.items():
        near=sorted((GEOD.inv(*xy,*v)[2]/1000,j) for j,v in centers.items() if i!=j)
        for rank,(km,j) in enumerate(near,1):
            if km<=500:dist.append((i,j,km,rank))
    pairs=cells.merge(pd.DataFrame(dist,columns=["island","destination","km","neighbor_rank"]),on="island").merge(first,on=["key","destination"])
    pairs=pairs[pairs.first_year.ge(2012)]
    # A middle-island record in the original annual series must precede the destination.
    dates=a[a.year.between(2011,2023)].groupby(["speciesKey","island"]).year.min().rename("earliest_middle_candidate")
    pairs=pairs.join(dates,on=["key","island"]).loc[lambda z:z.earliest_middle_candidate<z.first_year]
    pairs.to_csv(OUT/"candidate_links.csv",index=False)
    print('candidate',len(pairs),pairs.island.nunique(),pairs.key.nunique(),flush=True)
    return pairs


def search(job):
    island,geom,species,start,end,kind,cache=job
    if not hasattr(LOCAL,"session"):
        LOCAL.session=requests.Session()
        LOCAL.session.mount("https://",HTTPAdapter(max_retries=Retry(total=5,backoff_factor=.8,status_forcelist=[429,500,502,503,504])))
    xmin,ymin,xmax,ymax=geom.bounds
    params=dict(taxonKey=sorted(map(int,species)),decimalLatitude=f"{math.floor(ymin*1e5)/1e5},{math.ceil(ymax*1e5)/1e5}",
                decimalLongitude=f"{math.floor(xmin*1e5)/1e5},{math.ceil(xmax*1e5)/1e5}",hasCoordinate="true",hasGeospatialIssue="false",
                occurrenceStatus="PRESENT",year=f"{start},{end}",limit=300)
    if kind in ("Larva","Egg"):params["lifeStage"]=kind
    else:params["q"]=kind
    token=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()[:16]
    path=cache/f"{island}_{token}.json"
    if path.exists():return json.loads(path.read_text(encoding="utf8"))
    matches=[];receipts=[];off=0
    while True:
        resp=LOCAL.session.get("https://api.gbif.org/v1/occurrence/search",params=dict(params,offset=off),timeout=75)
        resp.raise_for_status();d=resp.json()
        receipts.append(dict(island=island,kind=kind,url=resp.url,count=d["count"],offset=off,sha256=hashlib.sha256(resp.content).hexdigest()))
        for r in d["results"]:
            if r.get("speciesKey") not in species:continue
            lon,lat=r["decimalLongitude"],r["decimalLatitude"]
            if not geom.covers(shapely.Point(lon,lat)):continue
            stage=r.get("lifeStage","")
            words=" ".join(str(r.get(k,"")) for k in ("behavior","occurrenceRemarks","fieldNotes"))
            textlay=bool(LAY.search(words) and not NEG.search(words))
            if kind not in ("Larva","Egg") and not textlay:continue
            matches.append(dict(island=island,key=r["speciesKey"],year=r["year"],gbif_id=r["key"],stage=stage,query=kind,
                                laying=textlay,evidence_text=words.strip() if textlay else "",lon=lon,lat=lat,
                                dataset_key=r["datasetKey"],uncertainty_m=r.get("coordinateUncertaintyInMeters"),basis=r.get("basisOfRecord")))
        off+=len(d["results"])
        if d["endOfRecords"]:break
        assert 0<off<100000
    result=dict(receipts=receipts,matches=matches)
    path.write_text(json.dumps(result),encoding="utf8")
    return result


def collect(workspace,pairs):
    base=workspace/"work/barriers"
    shape=gpd.read_file(base/"GSHHS_shp/h/GSHHS_h_L1.shp",where="sibling_id IN ("+",".join("'"+str(i)+"'" for i in pairs.island.unique())+")").to_crs(4326)
    cache=base/"relay_stage_cache";cache.mkdir(exist_ok=True)
    jobs=[]
    for i,g in pairs.groupby("island"):
        geom=shapely.union_all(shape.loc[shape.sibling_id.astype(int).eq(i)].geometry.to_numpy())
        assert not geom.is_empty
        for kind in ("Larva","Egg","ovipositing","laying eggs","oviposited"):
            jobs.append((int(i),geom,set(g.key),int(g.earliest_middle_candidate.min()),int(g.first_year.max()-1),kind,cache))
    with ThreadPoolExecutor(max_workers=5) as pool:
        results=list(pool.map(search,jobs))
    rec=pd.DataFrame([r for d in results for r in d["matches"]])
    rec.to_csv(OUT/"stage_search_hits.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    pd.DataFrame([r for d in results for r in d["receipts"]]).to_csv(OUT/"query_receipts.csv",index=False)
    return rec


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("workspace",type=Path)
    workspace=ap.parse_args().workspace
    OUT.mkdir(exist_ok=True,parents=True)
    collect(workspace,frame(workspace))
