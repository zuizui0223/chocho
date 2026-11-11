"""Compile 2001–2010 flora at unambiguous WPI island ports with annual butterfly records."""
import argparse
import hashlib
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from pyproj import Geod
from requests.adapters import HTTPAdapter
from shapely import Polygon, to_wkt
from urllib3.util.retry import Retry

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/route_free_port_hosts"
GEOD=Geod(ellps="WGS84")
LOCAL=threading.local()


def load_frame(workspace):
    w=workspace/"work/barriers"
    pts=json.loads((w/"ports_2019.geojson").read_text())["features"]
    registry={p["properties"]["OBJECTID"]:p for p in pts}
    ports=pd.read_csv(workspace/"outputs/port_buffer_landmass_audit.csv").query("status=='SINGLE_LANDMASS' and islands_intersecting_5km_buffer==1").copy()
    ports["island"]=ports.single_landmass_id.astype(int)
    ports["country"]=ports.port_id.map(lambda i:registry[i]["properties"]["COUNTRY"])
    ports["lon"]=ports.port_id.map(lambda i:registry[i]["geometry"]["coordinates"][0])
    ports["lat"]=ports.port_id.map(lambda i:registry[i]["geometry"]["coordinates"][1])
    annual=pd.read_parquet(w/"global_original_species_island_year.parquet").query("records>0 and 2001<=year<=2024").copy()
    annual["island"]=annual.landmass_id.astype(int)
    ports=ports[ports.island.isin(annual.island)]
    return ports,annual


def inventory(p,w):
    cache=w/"route_free_early_facets"
    cache.mkdir(exist_ok=True)
    stem=f"{p.port_id}.json"
    for directory in ("port_flora_2001_2010","expanded_port_flora_2001_2010","route_free_early_facets","port_flora_inventory"):
        f=w/directory/stem
        if not f.exists():continue
        d=json.loads(f.read_text())
        if directory=="port_flora_2001_2010" and d["status"] in ("COMPLETE_EARLY_FACETS","CACHED_ZERO_EARLY_RECORDS"):
            return d["count"],[int(x["name"]) for x in d["species_counts"]],directory,hashlib.sha256(f.read_bytes()).hexdigest()
        if directory in ("expanded_port_flora_2001_2010","route_free_early_facets"):
            return d["count"],d["keys"],directory,"|".join(d["hashes"])
        if directory=="port_flora_inventory" and d["status"]=="COMPLETE_FACETS" and not any(2001<=int(x["name"])<=2010 for x in d["facets"]["YEAR"]):
            return 0,[],directory,hashlib.sha256(f.read_bytes()).hexdigest()
    polygon=Polygon([GEOD.fwd(p.lon,p.lat,-a,5000)[:2] for a in range(0,360,15)])
    assert polygon.bounds[2]-polygon.bounds[0]<180
    params=dict(taxonKey=7707728,hasCoordinate="true",hasGeospatialIssue="false",occurrenceStatus="PRESENT",year="2001,2010",
                geometry=to_wkt(polygon,rounding_precision=6),limit=0,facet="speciesKey",facetLimit=1000)
    if not hasattr(LOCAL,"client"):
        LOCAL.client=requests.Session()
        LOCAL.client.mount("https://",HTTPAdapter(max_retries=Retry(total=6,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
    keys=[];hashes=[];count=None
    while True:
        response=LOCAL.client.get("https://api.gbif.org/v1/occurrence/search",params=dict(params,facetOffset=len(keys)),timeout=60)
        response.raise_for_status()
        d=response.json()
        assert count is None or count==d["count"]
        count=d["count"];hashes.append(hashlib.sha256(response.content).hexdigest())
        batch=[int(x["name"]) for facet in d["facets"] if facet["field"]=="SPECIES_KEY" for x in facet["counts"]]
        assert not set(keys).intersection(batch)
        keys.extend(batch)
        if len(batch)<1000:break
    (cache/stem).write_text(json.dumps(dict(count=count,keys=keys,hashes=hashes,params=params,retrieved=datetime.now(timezone.utc).isoformat())))
    return count,keys,"route_free_early_facets","|".join(hashes)


def acquire(workspace):
    ports,_=load_frame(workspace)
    ports=ports[~ports.island.isin(ports.loc[ports.country.isin(["DK","DE"]),"island"])]
    rows=[];taxa=[]
    with ThreadPoolExecutor(max_workers=5) as executor:
        tasks=[(p,executor.submit(inventory,p,workspace/"work/barriers")) for p in ports.itertuples()]
        for p,task in tasks:
            n,keys,source,digest=task.result()
            rows.append((p.island,p.port_id,p.country,n,source,digest))
            taxa.extend((p.port_id,k) for k in keys)
    pd.DataFrame(rows,columns=["island","port_id","country","records","source","hashes"]).to_csv(OUT/"port_coverage.csv",index=False)
    pd.DataFrame(taxa,columns=["port_id","species_key"]).to_csv(OUT/"port_taxa.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    print("Coverage complete:",len(rows),len(taxa),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("workspace",type=Path)
    OUT.mkdir(exist_ok=True,parents=True)
    acquire(parser.parse_args().workspace)
