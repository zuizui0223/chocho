#!/usr/bin/env python3
"""GBIF source-observed butterfly × WGSRPD3 records; capped stratified search.

Exactly two temporal windows and two systematic index pages (<=300 records
each) per butterfly per window are selected WITHOUT inspecting whether a
GBIF butterfly record coincides with an introduced plant host.

This is a BOUNDED NONRANDOM SAMPLE, not a complete GBIF download or a
biodiversity detection/non-detection survey. Output is aggregated only.
See protocol BUTTERFLY_239_GBIF_OCCURRENCE_EFFORT_MATCHED_20261019.md.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from shapely.geometry import Point, shape
from shapely.strtree import STRtree

API="https://api.gbif.org/v1"
WINDOWS=(("2010_2017","2010,2017"),("2018_2025","2018,2025"))
PAGE_SIZE=300
MAX_ORDINAL=100000
PRIMARY_BASIS={"HUMAN_OBSERVATION","PRESERVED_SPECIMEN",
               "OBSERVATION","MACHINE_OBSERVATION","OCCURRENCE"}
CORE_BASIS={"HUMAN_OBSERVATION","PRESERVED_SPECIMEN"}
USER_AGENT="chocho-global-239-butterfly-ecological-test/0.1 (https://github.com/zuizui0223/chocho)"
WAIT_LOCK=threading.Lock()
NEXT_TIME=[0.]
REQUEST_INTERVAL_SECONDS=.23


def polite_wait():
    with WAIT_LOCK:
        now=time.monotonic()
        delay=max(0., NEXT_TIME[0]-now)
        NEXT_TIME[0]=max(now,NEXT_TIME[0]) + REQUEST_INTERVAL_SECONDS
    if delay:
        time.sleep(delay)


def get_json(path,params,*,timeout=18.,tries=5):
    url=API+"/"+path+"?"+urlencode(params)
    last=None
    for attempt in range(tries):
        polite_wait()
        req=Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
        try:
            with urlopen(req,timeout=timeout) as response:
                return json.load(response)
        except HTTPError as exc:
            last=f"HTTP_{exc.code}"
            if exc.code not in (429,500,502,503,504):
                raise RuntimeError(f"{last} on GBIF query {path}") from exc
            retry_after=exc.headers.get("Retry-After")
            try:
                delay=min(30.,max(1.,float(retry_after))) if retry_after else min(24.,2.0**attempt)
            except ValueError:
                delay=min(24.,2.0**attempt)
        except (TimeoutError,URLError,ValueError,OSError) as exc:
            last=f"{type(exc).__name__}: {exc}"
            delay=min(24.,2.0**attempt)
        if attempt<tries-1:
            time.sleep(delay)
    raise RuntimeError(f"GBIF request unsuccessful ({last}) on {path}")


def resolve_exact(name):
    response=get_json("species/match",{"name":name,"strict":"true"})
    canonical=str(response.get("canonicalName") or "").strip()
    if (canonical.lower()!=name.lower() or
        response.get("rank")!="SPECIES" or
        response.get("matchType")!="EXACT" or
        response.get("status") not in ("ACCEPTED","SYNONYM") or
        not response.get("usageKey")):
        return None, {"gbif_result":response,"reason":"NO_EXACT_SPECIES_MATCH"}
    key=int(response.get("acceptedUsageKey") or response["usageKey"])
    return key, {
        "key":key,"canonical":canonical,"match_type":response["matchType"],
        "taxon_status":response["status"],
        "resolver":"GBIF_species_match_strict_canonical_exact",
        "synonym_accepted_key":response.get("acceptedUsageKey")
    }


def chosen_offsets(count,pages):
    """Bounded systematic sample across a retrievable ordinal index."""
    count=int(count)
    if count<=0:return []
    available=min(count,MAX_ORDINAL)
    if available<=PAGE_SIZE:
        return [0]
    number=min(pages,int(math.ceil(available/PAGE_SIZE)))
    if number==1:
        return [min(MAX_ORDINAL-PAGE_SIZE,
                    int((.5*available)//PAGE_SIZE)*PAGE_SIZE)]
    if available<=PAGE_SIZE*number:
        return list(range(0,available,PAGE_SIZE))[:number]
    lo_max=min(MAX_ORDINAL-PAGE_SIZE,available-PAGE_SIZE)
    offsets=[
        min(lo_max,max(0,int((available*(k+.5)/number)//PAGE_SIZE)*PAGE_SIZE))
        for k in range(number)
    ]
    result=sorted(set(offsets))
    if len(result)!=number:
        for i in range(0,lo_max+1,PAGE_SIZE):
            if i not in result:
                result.append(i)
                if len(result)==number:break
    return sorted(result)


def georecords_for_species(species,pages,*,timeout):
    key,taxon=resolve_exact(species)
    if key is None:return {
        "species":species,"status":"TAXON_UNRESOLVED","taxon":taxon,"windows":{},
        "records":[],"errors":[]
    }
    windows={}
    coords=[]
    ids=set()
    basis=Counter()
    errors=[]
    duplicates=0
    year_mismatch=0
    taxon_mismatch=0
    too_uncertain=0
    missing_coords=0
    for tag,years in WINDOWS:
        params={"taxonKey":key,"hasCoordinate":"true",
                "hasGeospatialIssue":"false",
                "occurrenceStatus":"PRESENT",
                "year":years}
        try:
            count=int(get_json("occurrence/search",{**params,"limit":1},
                               timeout=timeout)["count"])
            offsets=chosen_offsets(count,pages)
            ok=True
            received=0
            for off in offsets:
                resp=get_json("occurrence/search",{
                    **params,"limit":PAGE_SIZE,"offset":off},
                    timeout=timeout)
                if "results" not in resp or resp.get("count") is None:
                    raise ValueError("Missing GBIF result page")
                for item in resp["results"]:
                    received+=1
                    try:
                        rec_key=int(item.get("key") or item.get("gbifID"))
                    except (TypeError,ValueError):
                        continue
                    if rec_key in ids:
                        duplicates+=1
                        continue
                    ids.add(rec_key)
                    sy_key=item.get("speciesKey")
                    other=item.get("taxonKey")
                    if ((sy_key is not None and int(sy_key)!=key) or
                        (sy_key is None and (other is None or int(other)!=key))):
                        taxon_mismatch+=1
                        continue
                    yr=item.get("year")
                    if not isinstance(yr,int) or not int(years[:4])<=yr<=int(years[-4:]):
                        year_mismatch+=1
                        continue
                    lat=item.get("decimalLatitude")
                    lon=item.get("decimalLongitude")
                    try:
                        lat=float(lat);lon=float(lon)
                        if (not math.isfinite(lat) or not math.isfinite(lon)
                            or not -90<=lat<=90 or not -180<=lon<=180
                            or (lat==0. and lon==0.)):
                            missing_coords+=1;continue
                    except (TypeError,ValueError):
                        missing_coords+=1;continue
                    uncertainty=item.get("coordinateUncertaintyInMeters")
                    try:
                        if uncertainty is not None and float(uncertainty)>100000:
                            too_uncertain+=1;continue
                    except (TypeError,ValueError):
                        pass
                    original_basis=str(item.get("basisOfRecord") or "").upper()
                    if original_basis not in PRIMARY_BASIS:
                        basis["EXCLUDED_"+original_basis]+=1
                        continue
                    basis[original_basis]+=1
                    coords.append((tag,lat,lon,original_basis))
            windows[tag]={
                "year_filter":years,"reported_GBIF_count":count,
                "index_offsets":offsets,"pages_requested":len(offsets),
                "received_records_before_quality_filter":received,
                "total_search_more_than_100000":count>MAX_ORDINAL,
                "complete":True
            }
        except Exception as exc:
            errors.append({"window":tag,"error":str(exc)[:330]})
            windows[tag]={
                "year_filter":years,"complete":False,
                "reason":"PAGE_RETRIEVAL_FAILED_EXPLICIT_NOT_ABSENCE"
            }
    status=("COMPLETE" if all(z.get("complete") for z in windows.values())
            and len(windows)==len(WINDOWS) else "PARTIAL_FETCH_NOT_ELIGIBLE")
    return {
        "species":species,"status":status,"taxon":taxon,"windows":windows,
        "records":coords,"errors":errors,"raw_taxon_mismatch":taxon_mismatch,
        "raw_year_mismatch":year_mismatch,"raw_duplicates":duplicates,
        "raw_too_uncertain":too_uncertain,"raw_bad_coords":missing_coords,
        "basis_counts":dict(sorted(basis.items()))
    }


def source_polygon_map(geojson):
    d=json.loads(Path(geojson).read_text(encoding="utf-8"))
    shapes=[];ids=[]
    for feature in d["features"]:
        code=str((feature.get("properties") or {}).get("LEVEL3_COD") or "").strip()
        geom=feature.get("geometry")
        if geom and code:
            shapes.append(shape(geom))
            ids.append(code)
    assert 360<=len(shapes)<=450
    assert len(set(ids))==len(ids)
    tree=STRtree(shapes)
    by_id={id(g):i for i,g in enumerate(shapes)}
    return tree,shapes,ids,by_id


def region_for_latlon(bundle,lat,lon):
    tree,shapes,ids,by_id=bundle
    p=Point(lon,lat)
    codes=[]
    for value in tree.query(p):
        idx=(by_id[id(value)] if hasattr(value,"geom_type") else int(value))
        if shapes[idx].covers(p):
            codes.append(ids[idx])
    if not codes:
        return None,0
    return sorted(set(codes))[0],len(set(codes))


def save_json(path,record):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(record,sort_keys=True,indent=2)+"\n",encoding="utf-8")


def run(args):
    with args.species_csv.open(newline="",encoding="utf-8") as file:
        names=sorted({r["species"] for r in csv.DictReader(file)})
    assert len(names)==239, f"Expected original global panel 239, got {len(names)}"
    requested=names[:args.max_species]
    geo=source_polygon_map(args.level3_geojson)
    results=[]
    counts=Counter()
    sample_status=Counter()
    coord_statistics=Counter()
    starting=datetime.now(timezone.utc).isoformat()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        fs={pool.submit(georecords_for_species,s,args.pages_per_window,
                        timeout=args.timeout):s for s in requested}
        for idx,future in enumerate(as_completed(fs),1):
            species=fs[future]
            try:
                received=future.result()
            except Exception as exc:
                received={"species":species,"status":"ERROR_RETRIEVAL",
                          "taxon":None,"windows":{},"records":[],
                          "errors":[{"error":str(exc)[:400]}]}
            sample_status[received["status"]]+=1
            rowcount=Counter()
            for window,lat,lon,basis in received.pop("records",[]):
                code,geometries=region_for_latlon(geo,lat,lon)
                if not code:
                    coord_statistics["unmapped"]+=1
                    continue
                coord_statistics["mapped"]+=1
                if geometries>1:
                    coord_statistics["overlapping_TDWG_polygons"]+=1
                rowcount[window,code,basis]+=1
            if received["status"]=="COMPLETE":
                for (window,code,basis),n in rowcount.items():
                    counts[species,window,code,basis]+=n
            else:
                coord_statistics["records_discarded_incomplete_species"]+=sum(rowcount.values())
            results.append(received)
            if idx%20==0 or idx==len(requested):
                print(json.dumps({"processed":idx,"requested":len(requested),
                    "status":dict(sample_status),"mapped":coord_statistics["mapped"]}),
                    flush=True)
    out_rows=[]
    for (species,window,code,basis),n in sorted(counts.items()):
        out_rows.append({"species":species,"window":window,"wgsrpd3":code,
                         "basis":basis,"sampled_record_count":n})
    args.output_counts.parent.mkdir(parents=True,exist_ok=True)
    with args.output_counts.open("w",encoding="utf-8",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=[
            "species","window","wgsrpd3","basis","sampled_record_count"])
        writer.writeheader()
        writer.writerows(out_rows)
    audit={
        "schema":"chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01",
        "status":"COMPLETE_PAGED_QUERY_AUDIT" if sample_status["ERROR_RETRIEVAL"]==0
          and sample_status["PARTIAL_FETCH_NOT_ELIGIBLE"]==0 else
          "PARTIAL_GBIF_SOURCE_MISSINGNESS_NO_UNSAMPLED_ABSENCE_INFERENCE",
        "started_utc":starting,
        "completed_utc":datetime.now(timezone.utc).isoformat(),
        "gbif_api":API,
        "queries":{"windows":WINDOWS,"limit_per_page":PAGE_SIZE,
                    "max_pages_per_window":args.pages_per_window,
                    "max_ordinal_per_window":MAX_ORDINAL,
                    "query_filter":{
                      "hasCoordinate":True,"hasGeospatialIssue":False,
                      "occurrenceStatus":"PRESENT"
                    }},
        "species_total_panel":len(names),
        "species_requested":len(requested),
        "sampling_is_systematic_and_not_random_or_exhaustive":True,
        "species_taxon_resolved_and_complete":sample_status["COMPLETE"],
        "species_statuses":dict(sample_status),
        "sample_coordinate_diagnostics":dict(coord_statistics),
        "species_receipts":sorted(results,key=lambda x:x["species"]),
        "outputs":{"aggregated_region_basis_rows":len(out_rows),
                   "source_records_mapped":sum(counts.values())},
        "scope_warning":"Not representative occurrence download, not confirmed survey non-detection, cannot estimate biological host acceptance or colonisation; BGIF taxon and GBIF record sampling limits retained."
    }
    save_json(args.output_manifest,audit)
    if args.require_complete_panel and (
        len(requested)!=239 or sample_status["COMPLETE"]!=239):
        raise RuntimeError("INCOMPLETE_239_GBIF_ACQUISITION: retained diagnostics and output, do not interpret unqueried species as absent")
    print(json.dumps({
        "status":audit["status"],
        "species_completed":sample_status["COMPLETE"],
        "species_requested":len(requested),
        "mapped_records":audit["outputs"]["source_records_mapped"],
        "aggregated_rows":len(out_rows),
        "unmapped":coord_statistics["unmapped"],
        "incomplete_species":sample_status["PARTIAL_FETCH_NOT_ELIGIBLE"]
    },sort_keys=True),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--species-csv",type=Path,required=True)
    p.add_argument("--level3-geojson",type=Path,required=True)
    p.add_argument("--output-counts",type=Path,required=True)
    p.add_argument("--output-manifest",type=Path,required=True)
    p.add_argument("--max-species",type=int,default=239)
    p.add_argument("--pages-per-window",type=int,default=2)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--timeout",type=float,default=18.)
    p.add_argument("--require-complete-panel",action="store_true")
    a=p.parse_args()
    if not (1<=a.max_species<=239 and 1<=a.pages_per_window<=4
            and 1<=a.workers<=5):
        raise SystemExit("Disallowed retrospective sampling-budget changes")
    run(a)
