#!/usr/bin/env python3
"""Recover ONLY the original 96 GBIF incompletely sampled butterfly species.

Preserves the original prespecified two ordinal 300-record windows per
butterfly × time band, but makes four 75-record TRANSPORT ONLY subrequests
per 300-index window. A previous GBIF API curl exit28 with HTTP200 is
a RETRYABLE INCOMPLETE BODY, not a completed source or scientific absence.
Names and original accepted GBIF speciesKeys frozen independently from the
239-attempt original artifact. No historical host-plant studies involved.

Only store derived species×WGSRPD3×year-window×basis counts and audit.
No raw GBIF coordinates or WCVP source geography are exported.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode

from acquire_gbif_butterfly_239_region_occurrences_20261019 import (
    API,chosen_offsets,source_polygon_map,region_for_latlon,PRIMARY_BASIS
)

WINDOWS=(("2010_2017","2010,2017"),("2018_2025","2018,2025"))
PAGE=300
TRANSPORT=75
LIMIT=100000
BASES=PRIMARY_BASIS
_LOCK=threading.Lock()
_NEXT=[0.]
UA="chocho-GBIF-partial96-original-ordinal-recovery/0.1 (https://github.com/zuizui0223/chocho)"


def request_json(path,params,*,timeout=28.0,retries=4):
    url=API+"/"+path+"?"+urlencode(params)
    last=""
    for attempt in range(retries):
        with _LOCK:
            now=time.monotonic()
            delay=max(0.,_NEXT[0]-now)
            _NEXT[0]=max(now,_NEXT[0])+.14
        if delay:time.sleep(delay)
        command=[
            "curl","--location","--silent","--show-error","--compressed",
            "--connect-timeout","8","--max-time",str(timeout),
            "--header",f"User-Agent: {UA}","--header","Accept: application/json",
            "--write-out","\n%{http_code}",url]
        try:
            p=subprocess.run(command,capture_output=True,text=True,
                             timeout=timeout+3,check=False)
            body,sep,status=(p.stdout or "").rpartition("\n")
            if not sep:
                status="000"
            try:http=int(status.strip())
            except ValueError:http=0
            if p.returncode==0 and http==200:
                obj=json.loads(body)
                if "results" not in obj or "count" not in obj:
                    raise ValueError("Incomplete GBIF response schema")
                return obj
            if http and http not in (200,429,500,502,503,504):
                raise RuntimeError(f"nonretryable HTTP{http} on {path}")
            # rc28 and http200 is absolutely retryable; prior full-source run
            # erroneously considered this an HTTP200 nonretryable fatality.
            last=f"curl_rc{p.returncode}_http{http}_bytes{len(body)}_stderr{(p.stderr or '')[:65]}"
        except (subprocess.TimeoutExpired,ValueError,json.JSONDecodeError) as ex:
            last=f"{type(ex).__name__}: {str(ex)[:75]}"
        if attempt<retries-1:time.sleep(min(12.,1.6**attempt))
    raise RuntimeError(f"GBIF source exhausted retryable transport: {last}")


def read_source_keys(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==96
    keys={}
    for row in rows:
        n=row["species"].strip()
        k=int(row["accepted_gbif_species_key"])
        if n in keys or k<=0:raise ValueError("Wrong GBIF source-name/key audit")
        keys[n]=k
    return keys


def process_one(species,key,geom):
    audit={"species":species,"accepted_species_key":key,
           "status":"PARTIAL_TRANSPORT_NOT_ABSENCE","errors":[],
           "windows":{},"record_basis_counts":{}}
    included=[]
    ids=set()
    counter=Counter()
    basis=Counter()
    for tag,years in WINDOWS:
        params={
            "taxonKey":key,"hasCoordinate":"true",
            "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
            "year":years}
        try:
            snapshot=request_json("occurrence/search",{**params,"limit":1})
            count=int(snapshot["count"])
            offsets=chosen_offsets(count,2)
            received=[]
            deviations=[]
            for offset in offsets:
                # No changed sampling design: exact same ORDINAL 300-record
                # window originally selected, merely requested 75 at a time.
                for step in range(0,PAGE,TRANSPORT):
                    start=offset+step
                    expected=min(TRANSPORT,max(0,count-start))
                    if expected==0:continue
                    response=request_json("occurrence/search",{
                        **params,"limit":expected,"offset":start})
                    if len(response["results"])!=expected:
                        raise RuntimeError(f"short page at {start}: {len(response['results'])}/{expected}")
                    newer=int(response["count"])
                    if newer!=count:
                        deviations.append(newer-count)
                    received.extend(response["results"])
            if len(deviations)>0 and max(abs(v) for v in deviations)>max(5,int(count*.001)):
                raise RuntimeError("significant GBIF count drift during ordinal sample")
            audit["windows"][tag]={
                "complete":True,"year_filter":years,
                "original_ordinal_page_size":PAGE,
                "transport_subrequest_size":TRANSPORT,
                "reported_GBIF_count":count,
                "selected_original_page_offsets":offsets,
                "sampled_original_ordinal_records":len(received),
                "max_count_drift":max((abs(x) for x in deviations),default=0),
                "source_count_more_than_100000":count>LIMIT
            }
            for obj in received:
                try:
                    rid=int(obj.get("key") or obj.get("gbifID") or 0)
                except (TypeError,ValueError):continue
                if rid<=0 or rid in ids:
                    counter["duplicate_or_missing_key"]+=1
                    continue
                ids.add(rid)
                spkey=obj.get("speciesKey")
                recordkey=obj.get("taxonKey")
                try:
                    belongs=(int(spkey)==key if spkey is not None
                             else recordkey is not None and int(recordkey)==key)
                except (TypeError,ValueError):belongs=False
                if not belongs:counter["other_species_key"]+=1;continue
                yr=obj.get("year")
                if not isinstance(yr,int) or yr<int(years[:4]) or yr>int(years[-4:]):
                    counter["year_mismatch"]+=1;continue
                try:
                    lat=float(obj.get("decimalLatitude"));lon=float(obj.get("decimalLongitude"))
                    if (not math.isfinite(lat) or not math.isfinite(lon)
                        or not -90<=lat<=90 or not -180<=lon<=180
                        or (lat==0 and lon==0)):
                        counter["bad_coords"]+=1;continue
                except (ValueError,TypeError):counter["bad_coords"]+=1;continue
                try:
                    uncert=obj.get("coordinateUncertaintyInMeters")
                    if uncert is not None and float(uncert)>100000:
                        counter["too_uncertain"]+=1;continue
                except (ValueError,TypeError):pass
                recordbasis=str(obj.get("basisOfRecord") or "").upper()
                if recordbasis not in BASES:
                    counter["excluded_basis"]+=1;continue
                code,duplicates=region_for_latlon(geom,lat,lon)
                if code is None:
                    counter["unmapped"]+=1;continue
                if duplicates>1:counter["overlapping_polygons"]+=1
                basis[recordbasis]+=1
                included.append((tag,code,recordbasis))
        except Exception as e:
            audit["windows"][tag]={
                "complete":False,"year_filter":years,
                "reason":"SOURCE_FETCH_INCOMPLETE_NOT_BIOLOGICAL_ABSENCE"}
            audit["errors"].append({"window":tag,"error":str(e)[:240]})
    if len(audit["windows"])==2 and all(d["complete"] for d in audit["windows"].values()):
        audit["status"]="COMPLETE"
        audited=Counter(included)
    else:
        # Critical: original initial GBIF run did NOT count any species
        # that failed one of the 2 windows, even when one year had records.
        audited=Counter()
    audit["n_mapped_complete"]=sum(audited.values())
    audit["record_basis_counts"]=dict(sorted(basis.items()))
    audit["source_validity_counters"]=dict(sorted(counter.items()))
    return audit,audited


def main(args):
    originals=read_source_keys(args.keys)
    ordered=list(originals.items())
    if not (0<=args.start_index<96) or args.start_index+args.max_species>96:
        raise ValueError("Shard must stay within original frozen 96 GBIF taxa")
    target=ordered[args.start_index:args.start_index+args.max_species]
    geom=source_polygon_map(args.tdwg_geojson)
    totals=Counter()
    audits=[]
    status=Counter()
    start=datetime.now(timezone.utc).isoformat()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(process_one,sp,k,geom):sp for sp,k in target}
        for i,task in enumerate(as_completed(pending),1):
            name=pending[task]
            try:audit,group=task.result()
            except Exception as ex:
                audit={"species":name,"accepted_species_key":originals[name],
                       "status":"PARTIAL_TRANSPORT_NOT_ABSENCE",
                       "errors":[{"error":str(ex)[:240]}],
                       "windows":{},"n_mapped_complete":0}
                group=Counter()
            assert name==audit["species"]
            if audit["status"]=="COMPLETE":
                for (window,region,basis),count in group.items():
                    totals[name,window,region,basis]+=count
            audits.append(audit)
            status[audit["status"]]+=1
            if i%10==0 or i==len(target):
                print(json.dumps({"processed":i,"attempted":len(target),
                      "complete":status["COMPLETE"],"partial":status["PARTIAL_TRANSPORT_NOT_ABSENCE"],
                      "mapped_complete":sum(totals.values())}),flush=True)
    args.output_counts.parent.mkdir(parents=True,exist_ok=True)
    with args.output_counts.open("w",encoding="utf-8",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=[
            "species","window","wgsrpd3","basis","sampled_record_count"])
        writer.writeheader()
        for (name,window,region,basis),n in sorted(totals.items()):
            writer.writerow({"species":name,"window":window,"wgsrpd3":region,
                "basis":basis,"sampled_record_count":n})
    report={
        "schema":"chocho_GBIF_96_partial_species_transportreread_20261019_v01",
        "original_source":"Exact key/name original GBIF 239-attempt run 39233336804, 96 PARTIAL_FETCH_NOT_ELIGIBLE",
        "start_utc":start,"end_utc":datetime.now(timezone.utc).isoformat(),
        "species_attempted":len(target),
        "shard_start_index":args.start_index,
        "shard_end_exclusive":args.start_index+len(target),
        "original_frozen_96_species":96,
        "complete_species":status["COMPLETE"],
        "not_completed_species":status["PARTIAL_TRANSPORT_NOT_ABSENCE"],
        "regional_basis_count_rows":len(totals),
        "mapped_complete_records":sum(totals.values()),
        "original_ordinal_pages_per_window":2,"original_ordinal_size":300,
        "transport_only_subrequest_size":75,
        "retained_original_source_sampling_protocol":True,
        "species_receipts":sorted(audits,key=lambda x:x["species"]),
        "caveat":"Nonrandom GBIF source page sample, partial species not absent, no butterfly host use/colonisation/recruitment implied",
    }
    args.output_manifest.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!="species_receipts"},sort_keys=True))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--keys",type=Path,required=True)
    parser.add_argument("--tdwg-geojson",type=Path,required=True)
    parser.add_argument("--output-counts",type=Path,required=True)
    parser.add_argument("--output-manifest",type=Path,required=True)
    parser.add_argument("--max-species",type=int,default=96)
    parser.add_argument("--start-index",type=int,default=0)
    parser.add_argument("--workers",type=int,default=3)
    args=parser.parse_args()
    if not (1<=args.max_species<=96 and 1<=args.workers<=4):
        raise SystemExit("Original missing-species recovery budget changed")
    main(args)
