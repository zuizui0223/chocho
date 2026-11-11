"""Audit dated native host occurrences for the ten named, regionally covered cases."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests
import shapely
from pyproj import Geod
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/native_hosts_ten_island_cases"
FIG = ROOT / "figures/integrated_island_port_butterfly_v03/FigS7_local_native_hosts"
SOURCE = ROOT / "results/top_host_island_chronology/events.csv"
API = "https://api.gbif.org/v1/"
GEOD = Geod(ellps="WGS84")


def acquire(workspace):
    from datetime import datetime, timezone
    base = workspace / "work/barriers"
    original = workspace / "work/remote_sensing/host_source/extracted/inputs/original"
    pairs = pd.read_csv(original/"insect_host_accepted.csv",low_memory=False)
    ranges = pd.read_csv(original/"native_extant_nondoubtful_wgsrpd3.csv")
    xwalk = pd.read_parquet(base/"host_flora_name_join_audit.parquet").set_index("host_name")
    cases = pd.read_csv(SOURCE).loc[lambda z:z.global_cell.str.startswith("Already")].copy()
    assert len(cases)==10
    ports = json.loads((base/"ports_2019.geojson").read_text())["features"]
    shore = gpd.read_file(base/"GSHHS_shp/h/GSHHS_h_L1.shp",where="sibling_id IN ('169','276','454','551','821')").to_crs(4326)
    ids = shore.sibling_id.astype(int)
    session = requests.Session()
    session.mount("https://",HTTPAdapter(max_retries=Retry(total=5,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
    def get(endpoint, params):
        r=session.get(API+endpoint,params=params,timeout=75)
        r.raise_for_status()
        return r
    catalogue=[]
    for c in cases.itertuples():
        plant_ids=pairs.loc[pairs.insect_species.eq(c.species),"accepted_plant_name_id"]
        names=ranges[ranges.accepted_plant_name_id.isin(plant_ids)&ranges.area_code_l3.eq(c.region)].accepted_name.unique()
        for name in sorted(names):
            row=xwalk.loc[name]
            if row.join_status=="SPECIES_FACET_COMPATIBLE":
                taxkey=int(row.accepted_gbif_species_key)
                rank="species"
            else:
                assert row.join_status=="EXACT_INFRASPECIFIC_REQUIRES_FINER_RECORDS"
                m=get("species/match",{"name":name,"strict":"true","kingdom":"Plantae"}).json()
                assert m["matchType"]=="EXACT" and m["scientificName"]==name and m["rank"] in ["VARIETY","SUBSPECIES"]
                taxkey=int(m.get("acceptedUsageKey",m["usageKey"]))
                rank="infraspecific"
            catalogue.append((c.island,c.key,c.species,c.year,c.region,name,taxkey,rank))
    cat=pd.DataFrame(catalogue,columns=["island","butterfly_key","butterfly","first_year","region","native_host","plant_key","rank"])
    cat.to_csv(OUT/"candidate_hosts.csv",index=False)
    records=[]
    receipts=[]
    for island, group in cat.groupby("island"):
        poly=shapely.union_all(shore.loc[ids.eq(island)].geometry.values)
        assert not poly.is_empty
        xmin,ymin,xmax,ymax=poly.bounds
        cases_here=cases[cases.island.eq(island)]
        port_names=set(" | ".join(cases_here.ports).split(" | "))
        registry={p["properties"]["PORT_NAME"]:p["geometry"]["coordinates"][:2] for p in ports
                  if p["properties"]["PORT_NAME"] in port_names and
                  xmin-.1<=p["geometry"]["coordinates"][0]<=xmax+.1 and ymin-.1<=p["geometry"]["coordinates"][1]<=ymax+.1}
        assert set(registry)==port_names
        scope=group[["plant_key","rank"]].drop_duplicates()
        end=int(group.first_year.max())-1
        params=dict(taxonKey=sorted(scope.plant_key.astype(int)),year=f"2001,{end}",occurrenceStatus="PRESENT",
                    hasCoordinate="true",hasGeospatialIssue="false",decimalLongitude=f"{math.floor(xmin*1e5)/1e5},{math.ceil(xmax*1e5)/1e5}",
                    decimalLatitude=f"{math.floor(ymin*1e5)/1e5},{math.ceil(ymax*1e5)/1e5}",limit=300)
        offset=0
        seen=set()
        while True:
            resp=get("occurrence/search",dict(params,offset=offset))
            page=resp.json()
            receipts.append(dict(island=int(island),offset=offset,total=page["count"],returned=len(page["results"]),
                                 retrieved_utc=datetime.now(timezone.utc).isoformat(),url=resp.url,sha256=hashlib.sha256(resp.content).hexdigest()))
            for r in page["results"]:
                assert r["key"] not in seen
                seen.add(r["key"])
                lon,lat=r["decimalLongitude"],r["decimalLatitude"]
                if not poly.covers(shapely.Point(lon,lat)):
                    continue
                for k,rank in scope.itertuples(index=False,name=None):
                    if (r.get("speciesKey")==k if rank=="species" else k in (r.get("acceptedTaxonKey"),r.get("taxonKey"))):
                        distances={n:GEOD.inv(*xy,lon,lat)[2] for n,xy in registry.items()}
                        near=min(distances,key=distances.get)
                        records.append((island,k,r["key"],r["year"],r.get("eventDate"),r["datasetKey"],r.get("coordinateUncertaintyInMeters"),near,round(distances[near],2)))
            if page["endOfRecords"]:
                assert len(seen)==page["count"]
                break
            offset+=len(page["results"])
            assert offset and offset<100000
        print(island,len(seen),"bbox records",flush=True)
    pd.DataFrame(records,columns=["island","plant_key","gbif_id","year","event_date","dataset_key","uncertainty_m","nearest_focal_port","distance_m"]).to_csv(OUT/"occurrences.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    (OUT/"query_receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("workspace",type=Path)
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    acquire(args.workspace)
