#!/usr/bin/env python3
"""Audit land-area fractions within original 1 km historical GBIF buffers.

Quality sensitivity: whether a port circle mostly over water is being compared
with a nonport circle mostly on land. No botanical records are requested here;
site coordinates/harbor identities and GSHHG 2.3.7 are independently source
locked. The thresholds are exploratory follow-up diagnostics, not preregistered
or proof that maritime plants arrived by vessel.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from pyproj import CRS, Transformer
from shapely.geometry import Point
from shapely.ops import transform

from acquire_matched_botany import load_land, GSHHG_SHA256

HERE=Path(__file__).resolve().parent
ZIP_SHA=GSHHG_SHA256
INPUT_BLOBS={
 "sites":("fixed_geometric_sites.csv","3373121189ac80b237d99e2a1b0a172101cfacb8"),
 "plants":("historical_plant_site_facets.csv","eb09a3f3ec053c00b43e8bc0678d6643ce02697b"),
 "pairs":("port_offport_historical_pairs.csv","ae7fa4986e1a5be73a3c50b89669d05bc5138ac1"),
 "receipt":("historical_port_offport_results.json","53ef6dddbe1a8efceaf518812bf2b3c77cf285d6"),
}

def source(name):
    file,sha=INPUT_BLOBS[name]
    data=(HERE/file).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if actual!=sha:
        raise ValueError("source historical site record changed: "+file+" "+actual)
    if file.endswith(".json"):return json.loads(data)
    return list(csv.DictReader(data.decode("utf8").splitlines()))

def means(x):
    return sum(x)/len(x) if x else None

def quantile(vals,p):
    if not vals:return None
    a=sorted(vals);x=(len(a)-1)*p
    i=math.floor(x);j=math.ceil(x)
    return a[i]+(a[j]-a[i])*(x-i)

def compute(zipfile,work):
    fixed=source("sites")
    plants=source("plants")
    pairs=source("pairs")
    receipt=source("receipt")
    if (len(fixed),len(plants),len(pairs))!=(114,114,38):
        raise ValueError("full38 source coverage changed")
    ids={r["landmass_id"] for r in fixed}
    geoms,digest=load_land(zipfile,work/"GSHHG",ids)
    if digest!=ZIP_SHA:raise ValueError("GSHHG ZIP source mismatch")
    sites_by_id={r["site_id"]:r for r in plants}
    by_island=defaultdict(list)
    for row in fixed:
        if row["site_type"]=="port":
            by_island[row["landmass_id"]].append(row)
    areas={}
    samples=[]
    for island_id in sorted(ids,key=int):
        focal=by_island[island_id]
        lon=sum(float(x["longitude"]) for x in focal)/len(focal)
        lat=sum(float(x["latitude"]) for x in focal)/len(focal)
        local=CRS.from_proj4(f"+proj=aeqd +lon_0={lon} +lat_0={lat} +datum=WGS84 +units=m +no_defs")
        proj=Transformer.from_crs(4326,local,always_xy=True)
        poly=transform(proj.transform,geoms[island_id])
        for s in [x for x in fixed if x["landmass_id"]==island_id]:
            xy=Point(*proj.transform(float(s["longitude"]),float(s["latitude"])))
            circle=xy.buffer(1000,resolution=32)
            frac=circle.intersection(poly).area/circle.area
            if not -1e-6<=frac<=1+1e-6:
                raise ValueError("Invalid 1km terrestrial fraction")
            frac=max(0.,min(1.,frac))
            areas[s["site_id"]]=frac
            samples.append({
                "site_id":s["site_id"],"landmass_id":island_id,
                "port_id":s["port_id"],"site_type":s["site_type"],
                "land_area_fraction_1km":frac,
                "shoreline_km":s["shoreline_km"],
                "gbif_2001_2010_plant_species_keys":sites_by_id[s["site_id"]]["plant_indexed_species_keys"],
                "gbif_2001_2010_top10_host_candidates":sites_by_id[s["site_id"]]["n_target_hosts"],
                "provider_datasets_minimum":sites_by_id[s["site_id"]]["provider_datasets_minimum"],
            })
    by_port=defaultdict(list)
    for s in samples:by_port[s["port_id"]].append(s)
    p_lookup={r["port_id"]:r for r in pairs}
    contrasts=[]
    for port_id,group in sorted(by_port.items(),key=lambda x:int(x[0])):
        p=next(s for s in group if s["site_type"]=="port")
        far=[s for s in group if s["site_type"]=="offport"]
        if len(far)!=2:raise ValueError("Need two offport controls")
        difference=[p["land_area_fraction_1km"]-s["land_area_fraction_1km"] for s in far]
        record=p_lookup[port_id]
        row={"port_id":port_id,"island_id":p["landmass_id"],
             "port_land_fraction":p["land_area_fraction_1km"],
             "mean_nonport_land_fraction":means([x["land_area_fraction_1km"] for x in far]),
             "max_abs_port_control_land_fraction_diff":max(abs(x) for x in difference),
             "mean_port_minus_nonport_land_fraction":means(difference)}
        for th in (1,5,10,20,50):
            row["botanical_support_"+str(th)]=record["meets_min_species_"+str(th)]
            row["host_share_difference_"+str(th)]=record["port_minus_control_top10_index_fraction_"+str(th)]
        contrasts.append(row)
    if len(contrasts)!=38:raise ValueError("38 port land fractions expected")

    qc={}
    for th in (1,5,10,20):
        for tol in (.15,.25,.40):
            available=[x for x in contrasts
               if x["botanical_support_"+str(th)]=="True"
               and x["max_abs_port_control_land_fraction_diff"]<=tol]
            by=defaultdict(list)
            for r in available:by[r["island_id"]].append(float(r["host_share_difference_"+str(th)]))
            means_by=[means(x) for x in by.values()]
            qc[f"min{th}_maxabsland{tol:g}"]={
               "matched_ports":len(available),
               "matched_islands":len(by),
               "island_equal_top10_port_minus_nonport_share":means(means_by),
               "positive_port_count":sum(float(x["host_share_difference_"+str(th)])>0 for x in available),
               "negative_port_count":sum(float(x["host_share_difference_"+str(th)])<0 for x in available)
            }
    return samples,contrasts,{
      "schema":"chocho_original_full38_1km_land_fraction_quality_v1",
      "status":"SOURCE_LOCKED_HISTORICAL_BOTANY_GEOMETRIC_AREA_QA_NOT_CAUSAL_EFFECT",
      "source_git_blobs":{k:{"file":file,"sha":sha} for k,(file,sha) in INPUT_BLOBS.items()},
      "gshhg_release_sha256":digest,
      "source_islands":len(ids),"ports":len(contrasts),"sites":len(samples),
      "land_fraction_mean_port":means([x["land_area_fraction_1km"] for x in samples if x["site_type"]=="port"]),
      "land_fraction_mean_nonport":means([x["land_area_fraction_1km"] for x in samples if x["site_type"]=="offport"]),
      "port_control_maxabs_landfrac_diff_median": sorted(
          x["max_abs_port_control_land_fraction_diff"] for x in contrasts)[len(contrasts)//2],
      "tolerance_sensitivity":qc,
      "notes":[
         "Geometric circle intersection uses 2017 GSHHG 2.3.7 high-resolution L1; not year-specific land change.",
         "Port 2019 WPI retrospective locations may not be historically active in 2001-2010.",
         "Results stratified by terrestrial circle area are exploratory; thresholds selected for quality sensitivity after basic botanical pilot.",
         "Fewer islands and non-independent port pairs reduce precision.",
         "Host plant taxonomy/introduced status and local butterfly feeding remain unverified.",
         "WPI exclusion does not independently map ports or urbanization from the actual years.",
         "No causal effect or ecological rescue inferred."
      ]
    }

def save_csv(file,rows):
    with file.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    args=p.parse_args()
    s,p,r=compute(args.zip,args.work)
    save_csv(HERE/"land_fraction_sites.csv",s)
    save_csv(HERE/"land_fraction_pairs.csv",p)
    (HERE/"land_fraction_quality_results.json").write_text(
        json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(r,ensure_ascii=False,indent=2))

if __name__=="__main__":main()
