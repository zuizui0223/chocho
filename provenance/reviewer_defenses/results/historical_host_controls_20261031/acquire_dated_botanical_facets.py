#!/usr/bin/env python3
"""Pilot historical GBIF botanical facets at source-checked 2019 WPI ports.

No original occurrence coordinates or restricted host edges are published.
A plant nonrecord is NOT host absence, nor does a public port observation
demonstrate intentional or ship-borne introduction.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

B=Path("provenance/reviewer_defenses/results")
PLANTS=B/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
PORTS=B/"port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv"
GATE=B/"historical_host_controls_20261031/pair_priority.csv"
PINS={
    PLANTS:"0a01831056c2cd59a4bd41bd1dbd4c84d77cf922",
    PORTS:"fdf73b11ddcb1644f46730123350ae52cc61cae6",
}
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
GBIF="https://api.gbif.org/v1/occurrence/search"
COUNTRIES={
    "US":{"US","USA","UNITED STATES","UNITED STATES OF AMERICA"},
    "DK":{"DK","DENMARK"},
    "DE":{"DE","GERMANY","DEUTSCHLAND"},
}
AGENT="chocho-20261031-historical-botanical-source-audit/1.0"

def read_rows(root,p):
    raw=(root/p).read_bytes()
    if p in PINS:
        sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
        if sha!=PINS[p]:raise ValueError(f"Git source blob changed: {p} {sha}")
    return list(csv.DictReader(raw.decode("utf8").splitlines()))

def get(url, params):
    full=url+"?"+urlencode(params)
    for trial in range(4):
        try:
            req=Request(full,headers={"User-Agent":AGENT,"Accept":"application/json"})
            with urlopen(req,timeout=60) as response:
                data=json.load(response)
            if not isinstance(data,dict) or data.get("error"):
                raise ValueError("Source JSON failed or explicit API error: "+str(data.get("error"))[:300])
            return data
        except (HTTPError,URLError,TimeoutError) as ex:
            if isinstance(ex,HTTPError) and ex.code == 400:
                body=ex.read(1200).decode("utf8",errors="replace")
                raise ValueError(f"Source HTTP 400: {body}; URL={full[:550]}") from ex
            if trial==3:raise
            if isinstance(ex,HTTPError) and ex.code not in (408,425,429,500,502,503,504):
                raise
            time.sleep(min(10,2**trial))
    raise AssertionError("unreachable")

def choose(root):
    plant=read_rows(root,PLANTS)
    allports=read_rows(root,PORTS)
    gate=read_rows(root,GATE)
    # Includes both preselected focal and the five status-unknown priority pairs.
    islands={i for r in gate if r["gate"] in
        ("PRIORITY_UNCLASSIFIED","FOCAL_HOST_OBSERVED") for i in (r["a"],r["b"])}
    portlist=[x for x in plant if x["landmass_id"] in islands]
    if (len(islands),len(portlist))!=(9,23):
        raise ValueError("Expected 9 islands and 23 source-indexed WPI ports")
    countries={r["port_id"]:r["country"] for r in allports}
    selected={}
    for r in portlist:
        island=r["landmass_id"]
        if island not in selected or (int(r["early_plant_records"]),-int(r["port_id"])) > (
             int(selected[island]["early_plant_records"]),-int(selected[island]["port_id"])):
            selected[island]=r
    return portlist,selected,countries

def get_ports(all_port_rows,countries):
    ids=sorted(int(r["port_id"]) for r in all_port_rows)
    data=get(WPI,{"where":"OBJECTID IN ("+",".join(map(str,ids))+")",
         "outFields":"OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE",
         "returnGeometry":"false","f":"json"})
    features=data.get("features",[])
    if len(features)!=len(ids):
        raise ValueError("WPI mirror port id count mismatch: "+str(len(features)))
    out={}
    for x in features:
        a=x["attributes"];p=int(a["OBJECTID"])
        if p in out or p not in ids:raise ValueError("WPI duplicate/unexpected OBJECTID")
        country=str(a.get("COUNTRY") or "").strip().upper()
        expected=countries[str(p)]
        if country not in COUNTRIES.get(expected,{expected}):
            raise ValueError(f"WPI identity mismatches original country for {p}: {expected} vs {country}")
        lat=float(a["LATITUDE"]);lon=float(a["LONGITUDE"])
        if not (-90<=lat<=90 and -180<=lon<=180):raise ValueError("Invalid WPI WGS84 coordinates")
        out[p]=dict(port_id=p,name=a.get("PORT_NAME"),country=country,
                    expected_country=expected,latitude=lat,longitude=lon,index_no=a.get("INDEX_NO"))
    return out

def year_facets(port,year):
    params={"taxonKey":"6","year":str(year),"hasCoordinate":"true",
            "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
            "geoDistance":f"5km,{port['latitude']:.7f},{port['longitude']:.7f}",
            "facet":"speciesKey","facetLimit":"1000","limit":"0"}
    keys={};offset=0;total=None
    while True:
        params["facetOffset"]=str(offset)
        data=get(GBIF,params)
        n=int(data["count"])
        if total is not None and n!=total:raise ValueError("GBIF index changed within faceted pages")
        total=n
        v=next((z.get("counts",[]) for z in data.get("facets",[])
            if z.get("field","").upper()=="SPECIES_KEY"),[])
        for item in v:
            k=str(item["name"])
            if k in keys:raise ValueError("Overlapping speciesKey facets")
            keys[k]=int(item["count"])
        if len(v)<1000:break
        offset+=1000
        if offset>=50000:raise ValueError("GBIF facet limit exceeded")
    if sum(keys.values())>total:raise ValueError("Facets exceed overall indexed records")
    return total,keys,offset//1000+1

def run(root,summary_path,cache_path):
    allports,selected,countries=choose(root)
    source=get_ports(allports,countries)
    annual=[];facets=[];site=[]
    for iid, r in sorted(selected.items(),key=lambda x:int(x[0])):
        p=source[int(r["port_id"])]
        year_sum=0
        for year in range(2001,2011):
            n,keys,pages=year_facets(p,year)
            year_sum+=n
            annual.append(dict(island_id=iid,port_id=p["port_id"],year=year,
                plant_records_indexed=n,species_keys_faceted=len(keys),
                records_without_species_key=n-sum(keys.values()),pages=pages))
            for k,count in sorted(keys.items()):
                facets.append((iid,p["port_id"],year,k,count))
            time.sleep(0.12)
        site.append(dict(island_id=iid,port_id=p["port_id"],port_name=p["name"],
            country=p["expected_country"],original_early_index_records=int(r["early_plant_records"]),
            independently_requeried_early_records=year_sum,
            count_difference=year_sum-int(r["early_plant_records"])))
        print("HISTORICAL_PLANTS",iid,p["port_id"],year_sum,"previous",r["early_plant_records"],flush=True)
    cache_path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(cache_path,"wt",newline="",encoding="utf8") as out:
        writer=csv.writer(out)
        writer.writerow(["island_id","wpi_port_id","year","gbif_speciesKey","indexed_records"])
        writer.writerows(facets)
    hashed=hashlib.sha256(cache_path.read_bytes()).hexdigest()
    summary={
        "schema":"chocho_20261031_plant_species_year_facets_pilot_v1",
        "status":"PUBLIC_BOTANICAL_OBSERVATION_AUDIT_NOT_HISTORIC_INTRODUCTION",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "source_git_blob_pins":{str(p):sha for p,sha in PINS.items()},
        "wpi_mirror_url":WPI,"gbif_occurrence_api":GBIF,
        "wpi_23_port_country_identities_matched":len(source),
        "pilot_islands":9,"one_port_per_island_by_max_original_early_record_index":True,
        "query_params":{"kingdom_taxonKey":6,"year_range":"2001-2010",
            "hasCoordinate":True,"hasGeospatialIssue":False,"occurrenceStatus":"PRESENT",
            "geoDistance":"5km from WPI coordinate","facet":"speciesKey"},
        "source_port_summary":site,"yearly_facets_recording_support":annual,
        "taxon_year_archive":{"filename":cache_path.name,"rows":len(facets),
            "sha256":hashed,"uploaded_to_action_artifacts_only":True},
        "not_source_equivalent_until_5km_shapes_and_original_species_filters_reconciled":True,
        "not_host_absence_or_introductory_pathway":True,
        "GBIF_may_reindex_since_original_collection":True,
        "not_butterfly_breeding_or_inter_island_migration":True,
        "original_butterfly_and_restricted_host_edges_not_exported":True,
    }
    summary_path.parent.mkdir(parents=True,exist_ok=True)
    summary_path.write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print("BOTANICAL_AUDIT",json.dumps(dict(islands=9,port_id_verified=len(source),
        annual_site_year_rows=len(annual),facets=len(facets),yearly_index_sum=sum(x["plant_records_indexed"] for x in annual),
        matched_source_counts=sum(x["count_difference"]==0 for x in site)),sort_keys=True),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path("."))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--artifact",type=Path,required=True)
    args=parser.parse_args()
    run(args.root,args.output,args.artifact)
