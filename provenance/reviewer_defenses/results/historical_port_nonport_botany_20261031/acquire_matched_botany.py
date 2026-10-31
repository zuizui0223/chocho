#!/usr/bin/env python3
"""Historical 2001-2010 GBIF botanical speciesKey comparison at island ports vs coast.

No butterfly future information determines the controls. Uses the same GSHHG
2.3.7 high-resolution L1 landmass IDs and WPI port source as prior chocho work.
These are newly generated fixed-coordinate coastal controls, NOT recovery of
the earlier unpublished coastal_control_design.parquet.

The 10 plants are previously identified global contributors to *potential*
butterfly host opportunity, not confirmed historically exotic, introduced
by ship, or locally used by larvae. Recording is not biological absence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import time
import zipfile
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Point
from shapely.ops import transform
from pyproj import CRS,Transformer

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
DATA=Path("provenance/reviewer_defenses/results")
SOURCES={
  "islands":(DATA/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json",
     "54c4515947d2a3ee7545e4007985ef82a2af6868"),
  "ports":(DATA/"island_dated_host_match_20261030/island_port_early_plant_availability.csv",
     "0a01831056c2cd59a4bd41bd1dbd4c84d77cf922"),
  "top_hosts":(DATA/"geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json",
     "d25a6198d70861d014b156f4f214691b62297634"),
}
WPI="https://services-eu1.arcgis.com/BuS9rtTsYEV5C0xh/ArcGIS/rest/services/World_Port_Index/FeatureServer/0/query"
GBIF="https://api.gbif.org/v1/occurrence/search"
GSHHG_RELEASE="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GSHHG_SHA256="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
UA="chocho-historical-island-port-offport-comparison/2026-10-31 research audit"

def blob(path):
    b=(ROOT/path).read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def read_sources():
    out={}
    for key,(path,expected) in SOURCES.items():
        actual=blob(path)
        if actual!=expected:raise ValueError("source Git blob changed: "+str(path)+" "+actual)
        raw=(ROOT/path).read_text(encoding="utf8")
        out[key]=list(csv.DictReader(io.StringIO(raw))) if path.suffix==".csv" else json.loads(raw)
    islands={str(x["island_id"]) for x in out["islands"]["islands"]}
    ports=[x for x in out["ports"] if x["landmass_id"] in islands]
    hosts=[x for x in out["top_hosts"]["source_contributor_plants"]]
    assert len(islands)==17 and len(ports)==38 and len(hosts)==10
    assert len(set(int(x["port_id"]) for x in ports))==38
    assert len(set(str(x["GBIF_accepted_speciesKey"]) for x in hosts))==10
    return islands,ports,hosts

def get_json(url,params=None,retries=4):
    address=url+("?"+urlencode(params,doseq=True) if params else "")
    last=None
    for trial in range(retries):
        try:
            req=Request(address,headers={"User-Agent":UA,"Accept":"application/json"})
            with urlopen(req,timeout=45) as response:
                if response.status!=200:raise ValueError("source HTTP "+str(response.status))
                d=json.load(response)
            if not isinstance(d,dict) or d.get("error"):
                raise ValueError("response JSON error: "+str(d.get("error"))[:350])
            return d
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
            last=exc
            if isinstance(exc,HTTPError) and exc.code not in (408,425,429,500,502,503,504):
                break
            if trial+1<retries:time.sleep(min(4.0,0.4*2**trial))
    raise RuntimeError("source retrieval failed "+address[:160]+" : "+repr(last))

def get_wpi(needed):
    fields="OBJECTID,INDEX_NO,PORT_NAME,COUNTRY,LATITUDE,LONGITUDE"
    # Direct verification of every source original WPI ID, not fuzzy port names.
    direct=get_json(WPI,{"where":"OBJECTID IN ("+",".join(map(str,sorted(needed)))+")",
        "outFields":fields,"returnGeometry":"false","f":"json"})
    v=direct.get("features",[])
    if len(v)!=len(needed):raise ValueError("Missing one or more source WPI OBJECTIDs: "+str(len(v)))
    out={}
    for f in v:
        row=f["attributes"]
        key=int(row["OBJECTID"])
        lat=float(row["LATITUDE"]);lon=float(row["LONGITUDE"])
        if not(-90<=lat<=90 and -180<=lon<=180):
            raise ValueError("Invalid port coordinates")
        out[key]={"port_id":key,"latitude":lat,"longitude":lon,
                  "name":str(row.get("PORT_NAME") or ""),"country":str(row.get("COUNTRY") or "")}
    # All ports for an actual non-port control distance screen.
    global_rows=[];offset=0;page_size=1000
    while True:
        p=get_json(WPI,{"where":"1=1","outFields":"OBJECTID,LATITUDE,LONGITUDE",
            "returnGeometry":"false","f":"json","orderByFields":"OBJECTID ASC",
            "resultOffset":offset,"resultRecordCount":page_size})
        block=p.get("features",[])
        global_rows.extend(block)
        offset+=len(block)
        if not block:break
        if offset>20000:raise ValueError("Unexpected WPI size")
        if len(block)<page_size and not p.get("exceededTransferLimit",False):break
    points={}
    for feature in global_rows:
        a=feature["attributes"]
        try:
            key=int(a["OBJECTID"]);lat=float(a["LATITUDE"]);lon=float(a["LONGITUDE"])
            if -90<=lat<=90 and -180<=lon<=180:points[key]=(lon,lat)
        except (ValueError,TypeError,KeyError):continue
    if not all(id in points for id in needed):
        raise ValueError("All-WPI locator failed to include all required source ports")
    if len(points)<1500:
        raise ValueError("The global WPI reference has unexpectedly few ports")
    return out,points

def load_land(zip_file,extract_dir,ids):
    data=zip_file.read_bytes()
    actual=hashlib.sha256(data).hexdigest()
    if actual!=GSHHG_SHA256:raise ValueError("GSHHG 2.3.7 ZIP SHA mismatch: "+actual)
    extract_dir.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        selected=[p for p in z.namelist() if p.endswith((".shp",".shx",".dbf",".prj"))
                  and "/h/GSHHS_h_L1." in "/"+p]
        if len(selected)<3:raise ValueError("Missing GSHHG 2.3.7 high-resolution L1 shapefile")
        for path in selected:z.extract(path,extract_dir)
        shp=next(extract_dir.rglob("GSHHS_h_L1.shp"))
    land=gpd.read_file(shp)
    if "sibling_id" not in land.columns:raise ValueError("GSHHG sibling_id missing")
    land["landmass_id"]=land.sibling_id.astype(int).astype(str)
    data={}
    for name,group in land[land.landmass_id.isin(ids)].groupby("landmass_id"):
        polygon=shapely.union_all(group.geometry.to_numpy())
        if polygon.is_empty or not polygon.is_valid:
            raise ValueError("Invalid original GSHHG island polygon "+name)
        data[name]=polygon
    if set(data)!=ids:
        raise ValueError("Unmatched GSHHG island IDs: "+str(sorted(ids-set(data))))
    return data,actual

def make_controls(islands,source_ports,wpi,all_wpi,ncontrols):
    controls=[];audit=[]
    by_island=defaultdict(list)
    for record in source_ports:
        by_island[str(record["landmass_id"])].append(int(record["port_id"]))
    for island_id in sorted(islands,key=int):
        local=by_island[island_id]
        center_lat=sum(wpi[p]["latitude"] for p in local)/len(local)
        center_lon=sum(wpi[p]["longitude"] for p in local)/len(local)
        crs=CRS.from_proj4("+proj=aeqd +lat_0="+str(center_lat)+
                          " +lon_0="+str(center_lon)+" +datum=WGS84 +units=m +no_defs")
        forward=Transformer.from_crs(4326,crs,always_xy=True)
        reverse=Transformer.from_crs(crs,4326,always_xy=True)
        polygon=transform(forward.transform,islands[island_id])
        if polygon.is_empty:raise ValueError("Empty projected island")
        # All WPI sites, including sites not selected for prior botanical surveys,
        # supply the NON-PORT exclusion (>=5 km).
        lons=np.array([x[0] for x in all_wpi.values()])
        lats=np.array([x[1] for x in all_wpi.values()])
        projected=np.column_stack(forward.transform(lons,lats))
        projected=projected[np.isfinite(projected).all(axis=1)]
        for port_id in sorted(local):
            point=Point(*forward.transform(wpi[port_id]["longitude"],wpi[port_id]["latitude"]))
            shoreline_distance=point.distance(polygon)
            p={"landmass_id":island_id,"port_id":port_id,
               "site_id":"PORT_"+str(port_id),"site_type":"port",
               "latitude":wpi[port_id]["latitude"],"longitude":wpi[port_id]["longitude"],
               "shoreline_km":shoreline_distance/1000,"nearest_wpi_km":0.,
               "port_distance_km":0.}
            if shoreline_distance>3000:
                audit.append({"island_id":island_id,"port_id":port_id,
                  "status":"EXCLUDE_SOURCE_PORT_OVER_3KM_FROM_ISLAND_LAND",
                  "port_distance_to_land_m":shoreline_distance,"n_controls":0})
                continue
            # Historical port query radius is fixed at 1km here, regardless of
            # its earlier 5km source; the controls use exactly the same radius.
            selected=[];candidate_counter=0
            for radius in (6000,9000,14000,21000,30000,44000,62000):
                if len(selected)>=ncontrols:break
                for step in range(72):
                    angle=2*math.pi*((step*23+port_id)%72)/72
                    candidate=Point(point.x+radius*math.cos(angle),
                                    point.y+radius*math.sin(angle))
                    if not polygon.contains(candidate):continue
                    coastal_distance=candidate.distance(polygon.boundary)
                    if not 150<=coastal_distance<=1250:continue
                    d=np.sqrt(np.min(np.sum((projected-np.array([candidate.x,candidate.y]))**2,axis=1)))
                    if d<5000:continue
                    if any(candidate.distance(old)<3000 for old in selected):continue
                    candidate_counter+=1
                    selected.append(candidate)
                    if len(selected)>=ncontrols:break
            if len(selected)<ncontrols:
                audit.append({"island_id":island_id,"port_id":port_id,
                  "status":"EXCLUDE_INSUFFICIENT_GEOMETRIC_NONPORT_CONTROLS",
                  "port_distance_to_land_m":shoreline_distance,
                  "n_controls":len(selected),
                  "reason":"same island, within 62 km, 150-1250m from coast, >=5km from all WPI"})
                continue
            controls.append(p)
            for ix,c in enumerate(selected,1):
                lon,lat=reverse.transform(c.x,c.y)
                dist=float(np.sqrt(np.min(np.sum((projected-np.array([c.x,c.y]))**2,axis=1))))
                controls.append({"landmass_id":island_id,"port_id":port_id,
                  "site_id":"OFF_"+str(port_id)+"_"+str(ix),"site_type":"offport",
                  "latitude":lat,"longitude":lon,
                  "shoreline_km":c.distance(polygon.boundary)/1000,
                  "nearest_wpi_km":dist/1000,
                  "port_distance_km":c.distance(point)/1000})
            audit.append({"island_id":island_id,"port_id":port_id,
              "status":"MATCHED_GEOMETRIC_DESIGN_ONLY","n_controls":len(selected),
              "port_distance_to_land_m":shoreline_distance,
              "valid_candidate_count_at_selection":candidate_counter})
    if len({p["site_id"] for p in controls})!=len(controls):
        raise ValueError("Duplicate fixed design site")
    return controls,audit

def plant_facets(site,hosts):
    p={"taxonKey":6,"year":"2001,2010","hasCoordinate":"true",
       "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
       "geoDistance":("%.6f,%.6f,1000m"%(site["latitude"],site["longitude"])),
       "facet":["speciesKey","datasetKey"],"facetLimit":1000,"facetMincount":1,
       "limit":0}
    species={}; datasets=set(); count=None
    for offset in range(0,3001,1000):
        q=dict(p);q["facetOffset"]=offset
        result=get_json(GBIF,q)
        if count is None:count=int(result["count"])
        # GBIF returns canonical facet FIELD ids (e.g. SPECIES_KEY), not
        # necessarily camel-case query names (speciesKey). Never silently
        # convert indexed records into zero botanical species on case mismatch.
        facets={str(x.get("field","")).upper():x.get("counts",[])
                for x in result.get("facets",[])}
        if count and "SPECIES_KEY" not in facets:
            raise ValueError("GBIF returned indexed plant records but no SPECIES_KEY facet: "+
                             repr(sorted(facets)))
        for entry in facets.get("SPECIES_KEY",[]):
            if entry.get("name"):
                species[str(entry["name"])]=int(entry.get("count",0))
        for entry in facets.get("DATASET_KEY",[]):
            if entry.get("name"):datasets.add(str(entry["name"]))
        if len(facets.get("SPECIES_KEY",[]))<1000:break
        if offset==3000:raise ValueError("SpeciesKey facet pagination exceeded safety cap")
    host_key={str(x["GBIF_accepted_speciesKey"]) for x in hosts}
    return dict(site,query_status="COMPLETE",
      plant_occurrence_index_count=count,plant_indexed_species_keys=len(species),
      species_keys=";".join(sorted(species,key=int)),
      target_host_species_keys=";".join(sorted(species.keys() & host_key,key=int)),
      n_target_hosts=len(species.keys() & host_key),
      provider_datasets_minimum=len(datasets))

def query_all(sites,hosts,max_workers):
    out=[]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures={executor.submit(plant_facets,s,hosts):s for s in sites}
        for future in as_completed(futures):
            site=futures[future]
            try:out.append(future.result())
            except Exception as exc:
                out.append(dict(site,query_status="FAILED_"+type(exc).__name__,
                  query_error=str(exc)[:300],
                  plant_occurrence_index_count=None,plant_indexed_species_keys=None,
                  species_keys="",target_host_species_keys="",n_target_hosts=None,
                  provider_datasets_minimum=None))
    return sorted(out,key=lambda x:x["site_id"])

def fmean(values):
    return sum(values)/len(values) if values else None

def evaluate_sites(sites):
    by_pair=defaultdict(list)
    for r in sites:by_pair[r["port_id"]].append(r)
    rows=[]
    for port_id,g in sorted(by_pair.items()):
        ps=[x for x in g if x["site_type"]=="port"]
        off=[x for x in g if x["site_type"]=="offport"]
        assert len(ps)==1
        port=ps[0]
        if port["query_status"]!="COMPLETE" or len(off)<2 or any(x["query_status"]!="COMPLETE" for x in off):
            status="QUERY_FAILED"
        else:status="PAIRED_INDEXED"
        d={"port_id":port_id,"island_id":port["landmass_id"],"status":status,
           "n_controls":len(off),
           "port_plant_species_keys":port["plant_indexed_species_keys"],
           "offport_plant_species_keys_mean":fmean([x["plant_indexed_species_keys"] for x in off])
                   if status=="PAIRED_INDEXED" else None,
           "port_top10_count":port["n_target_hosts"],
           "offport_top10_count_mean":fmean([x["n_target_hosts"] for x in off])
                   if status=="PAIRED_INDEXED" else None}
        for threshold in (1,5,10,20,50):
            valid=status=="PAIRED_INDEXED" and port["plant_indexed_species_keys"]>=threshold and all(
                x["plant_indexed_species_keys"]>=threshold for x in off)
            def fraction(x):return x["n_target_hosts"]/x["plant_indexed_species_keys"]
            d["meets_min_species_"+str(threshold)]=bool(valid)
            d["port_minus_control_top10_index_fraction_"+str(threshold)]=(
                fraction(port)-fmean([fraction(x) for x in off]) if valid else None)
        rows.append(d)
    return rows

def report(sites,pairs,layout,archived,allwpi,hosts,unmatched,maxports):
    summary={}
    for th in (1,5,10,20,50):
        key="port_minus_control_top10_index_fraction_"+str(th)
        eligible=[x for x in pairs if x["meets_min_species_"+str(th)]]
        island_means=defaultdict(list)
        for p in eligible:island_means[p["island_id"]].append(p[key])
        summary[str(th)]={"eligible_port_pairs":len(eligible),
          "eligible_islands":len(island_means),
          "island_equal_mean_port_minus_offport_host_fraction":fmean(
            [fmean(x) for x in island_means.values()]),
          "positive_port_pairs":sum(x[key]>0 for x in eligible),
          "negative_port_pairs":sum(x[key]<0 for x in eligible),
          "equal_port_pairs":sum(x[key]==0 for x in eligible)}
    return {
      "schema":"chocho_pre2011_port_nonport_island_top10_botanical_control_v1",
      "status":"RETROSPECTIVE_BOTANICAL_MATCHED_COASTAL_DIAGNOSTIC_NOT_SHIPPING_OR_BUTTERFLY_ESTABLISHMENT",
      "source_git_blobs":{k:{"path":str(x[0]),"sha":x[1]} for k,x in SOURCES.items()},
      "wpi_source":WPI,"gshhg_release":GSHHG_RELEASE,"gshhg_zip_sha256":archived,
      "gbif_source":GBIF,
      "window":"2001-2010 historical BEFORE 2011-2022 butterfly index",
      "port_buffer_km":1,"offport_buffer_km":1,
      "landmass_original_GSHHG_sibling_ids":17,"original_port_candidates":38,
      "source_ports_in_design":len(set(x["port_id"] for x in layout)),
      "sites_geometrically_matched":len(layout),
      "all_wpi_checked_for_nonport_distance":allwpi,
      "historical_plant_query_complete":sum(s["query_status"]=="COMPLETE" for s in sites),
      "historical_plant_query_failed":sum(s["query_status"]!="COMPLETE" for s in sites),
      "original_GEB_top_host_candidates":len(hosts),
      "excluded_geometric_pairs":unmatched,
      "maximum_ports_for_pilot":maxports,
      "results_by_min_background_plant_species_keys":summary,
      "limits":[
        "Ports are retroactively located from modern 2019 WPI, not historically verified active 2001-2010 harbor facilities.",
        "New controls were sampled on GSHHG coast BEFORE botanical results and are NOT identical to 2021 WorldCover matched points.",
        "Controls exclude all catalogued WPI within 5km but do not exclude unlisted historical harbors or nonport urban infrastructure.",
        "Every plant query is 1km buffer; land and sea sampled fractions may differ for harbor vs inland-offset offport sites.",
        "GEB top10 candidate plants are potential butterfly larval hosts globally, NOT verified introduced plants on a given island in 2001-2010.",
        "GBIF botanical records/datasets are unevenly sampled, and non-indexed species are not true absent species.",
        "Top10 ratio is an indexing proportion, not host biomass, host encounter rate, larval survival or butterfly colonization.",
        "No ship cargo or route movement is established; urbanization and agriculture may drive planting.",
        "The earlier private butterfly-year Parquet is not reconstructed from these botanical data.",
        "The experiment is retrospective and exploratory, NOT registered confirmatory analysis.",
      ]
    }

def save_csv(path,rows,fields):
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,required=True)
    p.add_argument("--work-dir",type=Path,required=True)
    p.add_argument("--max-ports",type=int,default=0)
    p.add_argument("--n-controls",type=int,default=2)
    p.add_argument("--max-workers",type=int,default=6)
    p.add_argument("--design-only",action="store_true")
    args=p.parse_args()
    if args.n_controls<2 or args.n_controls>3:raise ValueError("Require 2 or 3 independent botanical controls per port")
    island_ids,port_rows,hosts=read_sources()
    ids={int(x["port_id"]) for x in port_rows}
    wpi,all_wpi=get_wpi(ids)
    land,digest=load_land(args.zip,args.work_dir/"shape",island_ids)
    if args.max_ports>0:
        # Stratified island-balanced pilot: first source port per island before a
        # second site at any island. Does not select by botanical outcomes.
        selected=[];by_island=defaultdict(list)
        for r in port_rows:by_island[r["landmass_id"]].append(r)
        for z in range(max(len(v) for v in by_island.values())):
            for i in sorted(by_island,key=int):
                if z<len(by_island[i]):selected.append(by_island[i][z])
        port_rows=selected[:args.max_ports]
    site_design,design_status=make_controls(land,port_rows,wpi,all_wpi,args.n_controls)
    valid_port_pairs=len([r for r in design_status if r["status"]=="MATCHED_GEOMETRIC_DESIGN_ONLY"])
    HERE.mkdir(parents=True,exist_ok=True)
    save_csv(HERE/"fixed_geometric_sites.csv",site_design,
         ["landmass_id","port_id","site_id","site_type","latitude","longitude",
          "shoreline_km","nearest_wpi_km","port_distance_km"])
    save_csv(HERE/"geometry_admissions.csv",design_status,
         ["island_id","port_id","status","port_distance_to_land_m","n_controls",
          "valid_candidate_count_at_selection","reason"])
    print(json.dumps({"source_islands":17,"source_ports":38,"pilot_source_ports":len(port_rows),
         "geometric_admitted_port_pairs":valid_port_pairs,
         "fixed_botanical_query_sites":len(site_design),
         "wpi_all_ports":len(all_wpi),"gshhg_sha256":digest}),flush=True)
    if args.design_only:return
    sites=query_all(site_design,hosts,args.max_workers)
    # SpeciesKey strings are public, but never include original GBIF occurrence
    # ID/individual coordinates or restricted host network.
    save_csv(HERE/"historical_plant_site_facets.csv",sites,
         ["landmass_id","port_id","site_id","site_type","latitude","longitude",
          "query_status","plant_occurrence_index_count","plant_indexed_species_keys",
          "provider_datasets_minimum","n_target_hosts","target_host_species_keys",
          "species_keys","query_error"])
    pairs=evaluate_sites(sites)
    save_csv(HERE/"port_offport_historical_pairs.csv",pairs,list(pairs[0]) if pairs else
         ["port_id","island_id","status"])
    receipt=report(sites,pairs,site_design,digest,len(all_wpi),hosts,
         len(design_status)-valid_port_pairs,args.max_ports)
    (HERE/"historical_port_offport_results.json").write_text(
         json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":receipt["status"],"complete":receipt["historical_plant_query_complete"],
         "failed":receipt["historical_plant_query_failed"],
         "support":receipt["results_by_min_background_plant_species_keys"]},
         indent=2,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
