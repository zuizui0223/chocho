#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from shapely.geometry import Point, shape

GBIF="https://api.gbif.org/v1"
UA="chocho-full-history-temporal-hindcast/0.1"


def get_json(path:str,params:dict,attempts:int=12,timeout:float=30.0)->dict:
    url=f"{GBIF}/{path}?{urlencode(params)}"
    last=None
    for i in range(attempts):
        try:
            req=Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
            with urlopen(req,timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except HTTPError as exc:
            last=exc
            if i+1>=attempts:
                break
            if exc.code==429:
                raw=exc.headers.get("Retry-After")
                try:
                    delay=float(raw)
                except Exception:
                    delay=15.0+5.0*i
                time.sleep(min(120.0,max(10.0,delay)))
            else:
                time.sleep(min(30.0,0.75*(2**i)))
        except Exception as exc:
            last=exc
            if i+1<attempts:
                time.sleep(min(30.0,0.75*(2**i)))
    raise RuntimeError(f"GBIF request failed: {url} ({last})")


def load_pairs(path:Path):
    out=defaultdict(set)
    with path.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        for row in r:
            s=str(row["insect_species"]).strip(); h=str(row["accepted_plant_name_id"]).strip()
            if s and h: out[s].add(h)
    return out


def load_units(path:Path):
    out=defaultdict(set)
    with path.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        for row in r:
            h=str(row["accepted_plant_name_id"]).strip(); u=str(row["area_code_l3"]).strip()
            if h and u: out[h].add(u)
    return out


def load_geometry(path:Path):
    payload=json.loads(path.read_text(encoding="utf-8"))
    out={}
    for feature in payload.get("features",[]):
        props=feature.get("properties") or {}; raw=feature.get("geometry")
        code=str(props.get("LEVEL3_COD") or "").strip()
        if not code or not raw: continue
        geom=shape(raw); rp=geom.representative_point()
        out[code]={
            "level1":str(props.get("LEVEL1_COD") or "").strip(),
            "name":str(props.get("LEVEL3_NAM") or "").strip(),
            "lat":float(rp.y),"lon":float(rp.x),
            "geometry":geom,
        }
    if len(out)<300: raise RuntimeError("WGSRPD3 geometry drift")
    return out


def haversine_km(lat1,lon1,lat2,lon2):
    radius=6371.0088
    p1,p2=math.radians(lat1),math.radians(lat2)
    dlat=math.radians(lat2-lat1); dlon=math.radians(lon2-lon1)
    a=math.sin(dlat/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dlon/2)**2
    return 2*radius*math.asin(min(1.0,math.sqrt(a)))


def quantile_edges(vals):
    if not vals:return []
    a=np.asarray(vals,float)
    return [float(np.quantile(a,1/3)),float(np.quantile(a,2/3))]


def bin_value(v,edges):
    return int(sum(v>e for e in edges))


def build_snapshot_candidates(mapped_csv:Path,pairs,native,cont,geom,effort_threshold:int):
    universe=set(geom); occ=[]
    with mapped_csv.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        for row in r:
            code=str(row["wgsrpd3_code"]).strip()
            if code not in universe: continue
            try:y=int(float(row["year"]))
            except Exception:continue
            if 2010<=y<=2025:
                occ.append((str(row["species"]).strip(),code,y))
    species=sorted({s for s,_,_ in occ if s in pairs})
    baseline_by_sp=defaultdict(set); test_by_sp=defaultdict(set)
    baseline_counts=defaultdict(lambda:defaultdict(int)); test_counts=defaultdict(lambda:defaultdict(int))
    total_base=defaultdict(int); total_test=defaultdict(int)
    for s,u,y in occ:
        if 2010<=y<=2017:
            baseline_by_sp[s].add(u); baseline_counts[u][s]+=1; total_base[u]+=1
        elif 2018<=y<=2025:
            test_by_sp[s].add(u); test_counts[u][s]+=1; total_test[u]+=1

    resource={}
    for s in species:
        n=set(); c=set()
        for h in pairs[s]:
            n |= native.get(h,set()); c |= cont.get(h,set())
        if not n<=c: raise RuntimeError(f"native not subset contemporary for {s}")
        resource[s]={"native":n&universe,"contemporary":c&universe,"added":(c-n)&universe}

    rows=[]
    for s in species:
        base=baseline_by_sp[s]
        if not base: continue
        base_pts=[geom[u] for u in base]
        nres=resource[s]["native"]; cres=resource[s]["contemporary"]; added=resource[s]["added"]
        cand=[]
        for u in sorted(universe-base-nres):
            treatment=int(u in added)
            if not treatment and u in cres: continue
            eb=total_base[u]-baseline_counts[u].get(s,0)
            et=total_test[u]-test_counts[u].get(s,0)
            if eb<effort_threshold or et<effort_threshold: continue
            g=geom[u]
            dist=min(haversine_km(g["lat"],g["lon"],b["lat"],b["lon"]) for b in base_pts)
            cand.append({
                "species":s,"wgsrpd3_code":u,"wgsrpd3_name":g["name"],"level1":g["level1"],
                "treatment":treatment,"snapshot_outcome":int(u in test_by_sp[s]),
                "baseline_other_effort":eb,"test_other_effort":et,"distance_km":dist,
            })
        if not cand: continue
        eb_edges=quantile_edges([math.log1p(r["baseline_other_effort"]) for r in cand])
        et_edges=quantile_edges([math.log1p(r["test_other_effort"]) for r in cand])
        d_edges=quantile_edges([math.log1p(r["distance_km"]) for r in cand])
        for r in cand:
            r["baseline_effort_bin"]=bin_value(math.log1p(r["baseline_other_effort"]),eb_edges)
            r["test_effort_bin"]=bin_value(math.log1p(r["test_other_effort"]),et_edges)
            r["distance_bin"]=bin_value(math.log1p(r["distance_km"]),d_edges)
            r["stratum_key"]="|".join(map(str,(s,r["level1"],r["baseline_effort_bin"],r["test_effort_bin"],r["distance_bin"])))
            rows.append(r)
    return rows


def resolve_species(name:str):
    payload=get_json("species/match",{"name":name,"strict":"false"})
    key=payload.get("usageKey")
    if not key: raise RuntimeError(f"no GBIF taxon match for {name}: {payload}")
    return int(key),payload


def bbox_wkt(geom):
    minx,miny,maxx,maxy=geom.bounds
    return f"POLYGON(({minx} {miny},{maxx} {miny},{maxx} {maxy},{minx} {maxy},{minx} {miny}))"


def first_exact_record_year(taxon_key:int,geom,start_year:int,end_year:int):
    box=bbox_wkt(geom)
    payload=get_json("occurrence/search",{
        "taxonKey":taxon_key,"hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT","year":f"{start_year},{end_year}",
        "geometry":box,"facet":"year","facetLimit":max(300,end_year-start_year+1),
        "facetMincount":1,"limit":1,
    })
    years=[]
    for facet in payload.get("facets") or []:
        if str(facet.get("field") or "").upper()!="YEAR":continue
        for x in facet.get("counts") or []:
            try:y=int(x.get("name"))
            except Exception:continue
            if start_year<=y<=end_year: years.append(y)
    for year in sorted(set(years)):
        offset=0
        while True:
            page=get_json("occurrence/search",{
                "taxonKey":taxon_key,"hasCoordinate":"true","hasGeospatialIssue":"false",
                "occurrenceStatus":"PRESENT","year":year,"geometry":box,
                "limit":300,"offset":offset,
            })
            results=page.get("results") or []
            exact=[]
            for item in results:
                lat=item.get("decimalLatitude"); lon=item.get("decimalLongitude")
                if lat is None or lon is None:continue
                if geom.covers(Point(float(lon),float(lat))):
                    exact.append(item)
            if exact:
                return {
                    "first_year":year,
                    "first_record_count":len(exact),
                    "first_basis":";".join(sorted({str(x.get("basisOfRecord") or "") for x in exact})),
                    "first_dataset_count":len({str(x.get("datasetKey") or "") for x in exact}),
                }
            if bool(page.get("endOfRecords",True)) or not results:break
            offset += len(results)
            if offset>=100000:break
    return {"first_year":None,"first_record_count":0,"first_basis":"","first_dataset_count":0}


def null_summary(rows,key_fields,permutations,seed):
    if not rows:
        return {"rows":0}
    y=np.asarray([int(r["full_history_outcome"]) for r in rows],dtype=np.int8)
    t=np.asarray([int(r["treatment"]) for r in rows],dtype=np.int8)
    observed=int(np.sum(y*t))
    strata=defaultdict(list)
    for i,r in enumerate(rows):
        strata[tuple(r[k] for k in key_fields)].append(i)
    rng=np.random.default_rng(seed); B=permutations
    null=np.zeros(B,dtype=np.int32)
    for idxs in strata.values():
        idx=np.asarray(idxs,dtype=int); n_t=int(t[idx].sum())
        if n_t==0:continue
        m=len(idx); positives=int(y[idx].sum())
        if n_t==m:
            draws=np.full(B,positives,dtype=np.int16)
        else:
            draws=rng.hypergeometric(positives,m-positives,n_t,size=B).astype(np.int16)
        null += draws
    return {
        "key_fields":list(key_fields),"strata":len(strata),"observed_treatment_hits":observed,
        "null_mean":float(null.mean()),"null_median":float(np.median(null)),
        "null_q025":float(np.quantile(null,.025)),"null_q975":float(np.quantile(null,.975)),
        "observed_minus_null_mean":float(observed-null.mean()),
        "observed_to_null_mean_ratio":None if null.mean()==0 else float(observed/null.mean()),
        "p_high":float((1+np.sum(null>=observed))/(B+1)),
        "permutations":B,"seed":seed,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--occurrence-mapped-csv",type=Path,required=True)
    ap.add_argument("--insect-host-csv",type=Path,required=True)
    ap.add_argument("--native-distribution-csv",type=Path,required=True)
    ap.add_argument("--contemporary-distribution-csv",type=Path,required=True)
    ap.add_argument("--level3-geojson",type=Path,required=True)
    ap.add_argument("--effort-threshold",type=int,default=10)
    ap.add_argument("--history-start",type=int,default=1800)
    ap.add_argument("--history-end",type=int,default=2025)
    ap.add_argument("--permutations",type=int,default=99999)
    ap.add_argument("--seed",type=int,default=20261007)
    ap.add_argument("--workers",type=int,default=8)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    geom=load_geometry(args.level3_geojson)
    pairs=load_pairs(args.insect_host_csv); native=load_units(args.native_distribution_csv); cont=load_units(args.contemporary_distribution_csv)
    rows=build_snapshot_candidates(args.occurrence_mapped_csv,pairs,native,cont,geom,args.effort_threshold)

    snapshot={
        "candidate_cells":len(rows),
        "treatment_cells":sum(int(r["treatment"]) for r in rows),
        "control_cells":sum(1-int(r["treatment"]) for r in rows),
        "treatment_snapshot_hits":sum(int(r["treatment"]) and int(r["snapshot_outcome"]) for r in rows),
        "control_snapshot_hits":sum((1-int(r["treatment"])) and int(r["snapshot_outcome"]) for r in rows),
    }
    expected={"candidate_cells":1477,"treatment_cells":442,"control_cells":1035,"treatment_snapshot_hits":9,"control_snapshot_hits":11}
    if snapshot!=expected:
        raise RuntimeError(f"frozen candidate reconstruction drift: {snapshot} != {expected}")

    species=sorted({r["species"] for r in rows})
    resolved={s:resolve_species(s) for s in species}

    def task(i,r):
        key,_=resolved[r["species"]]
        x=first_exact_record_year(key,geom[r["wgsrpd3_code"]]["geometry"],args.history_start,args.history_end)
        return i,x

    results=[None]*len(rows)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(task,i,r) for i,r in enumerate(rows)]
        for done,fut in enumerate(as_completed(futures),start=1):
            i,x=fut.result(); results[i]=x
            if done%100==0:
                print(json.dumps({"history_cells_completed":done,"total":len(rows)}),flush=True)

    for r,x in zip(rows,results):
        first=x["first_year"]
        r["historical_first_record_year"]="" if first is None else first
        r["historical_record_pre2018"]=int(first is not None and first<=2017)
        r["full_history_eligible"]=int(first is None or first>=2018)
        r["full_history_outcome"]=int(first is not None and 2018<=first<=2025)
        r["historical_first_basis"]=x["first_basis"]
        r["historical_first_dataset_count"]=x["first_dataset_count"]

    clean=[r for r in rows if int(r["full_history_eligible"])]
    treat=[r for r in clean if int(r["treatment"])]
    control=[r for r in clean if not int(r["treatment"])]
    th=sum(int(r["full_history_outcome"]) for r in treat); ch=sum(int(r["full_history_outcome"]) for r in control)
    tr=th/len(treat) if treat else None; cr=ch/len(control) if control else None

    null_species=null_summary(clean,["species"],args.permutations,args.seed)
    null_level1=null_summary(clean,["species","level1"],args.permutations,args.seed)
    null_full=null_summary(clean,["species","level1","baseline_effort_bin","test_effort_bin","distance_bin"],args.permutations,args.seed)

    first_year_t=[int(r["historical_first_record_year"]) for r in treat if int(r["full_history_outcome"])]
    first_year_c=[int(r["historical_first_record_year"]) for r in control if int(r["full_history_outcome"])]
    payload={
        "schema":"chocho_butterfly_temporal_resource_hindcast_full_history_v0.1",
        "status":"POSTHOC_FULL_HISTORY_CORRECTED_TEMPORAL_HINDCAST",
        "snapshot_candidate_reconstruction":snapshot,
        "history_window":[args.history_start,args.history_end],
        "history_audit":{
            "cells_with_pre2018_record_removed":sum(int(r["historical_record_pre2018"]) for r in rows),
            "treatment_pre2018_removed":sum(int(r["historical_record_pre2018"]) and int(r["treatment"]) for r in rows),
            "control_pre2018_removed":sum(int(r["historical_record_pre2018"]) and not int(r["treatment"]) for r in rows),
            "clean_cells":len(clean),"clean_treatment_cells":len(treat),"clean_control_cells":len(control),
        },
        "observed":{
            "treatment_first_records_2018_2025":th,"treatment_rate":tr,
            "control_first_records_2018_2025":ch,"control_rate":cr,
            "raw_risk_difference":None if tr is None or cr is None else tr-cr,
            "raw_risk_ratio":None if tr is None or not cr else tr/cr,
            "treatment_first_years":sorted(first_year_t),
            "control_first_years":sorted(first_year_c),
        },
        "null_ladder":{
            "species_only":null_species,
            "species_plus_level1":null_level1,
            "full_protocol_match":null_full,
        },
        "decision":{
            "placement_signal_beyond_level1":bool(null_level1.get("p_high",1)>0 and null_level1.get("p_high",1)<=0.05),
            "placement_signal_beyond_full_match":bool(null_full.get("p_high",1)>0 and null_full.get("p_high",1)<=0.05),
        },
        "claim_boundary":[
            "The original 53,434-record snapshot is used only to define an outcome-blind candidate/effort-matched cohort; full GBIF history is queried independently for every candidate cell before defining temporal eligibility and outcome.",
            "A first GBIF record is a first documented detection, not a colonization or establishment date.",
            "The full-history correction removes cells with any focal-butterfly record through 2017, but it cannot prove biological absence before 2018.",
            "Contemporary WCVP introduced ranges remain time-invariant, so host presence by 2017 is not established in this stage."
        ]
    }

    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    args.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))


if __name__=="__main__":
    main()
