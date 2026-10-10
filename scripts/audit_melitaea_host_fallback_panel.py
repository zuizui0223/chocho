#!/usr/bin/env python3
"""Source-verified, fail-closed Melitaea host-decline feasibility (no effect fit)."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

URLS=[
 "https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.905qfttrg/download",
 "https://datadryad.org/downloads/file_stream/3349312",
]
NAME="empirical_models/data/RAWDATA/fall_survey_2004_2013.csv"
COLS={"Patch","Year","Network","Area","Occupancy","Nest_count","Pl","Vs"}
MAX=25_000_000

def acquire(opener=urlopen):
    attempts=[]
    for url in URLS:
        try:
            req=Request(url,headers={"User-Agent":"chocho-scientific-source-audit/1.0"})
            with opener(req,timeout=30) as r:
                data=r.read(MAX+1)
                status=getattr(r,"status",None)
            if len(data)>MAX or not zipfile.is_zipfile(io.BytesIO(data)):
                attempts.append({"url":url,"http_status":status,"error":"non-zip or oversized"})
                continue
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                match=[s for s in z.namelist() if s.endswith(NAME)]
                if len(match)!=1:
                    attempts.append({"url":url,"error":"exact source CSV not present exactly once"})
                    continue
                raw=z.read(match[0])
            return raw,{"archive_url":url,"archive_sha256":hashlib.sha256(data).hexdigest(),
                        "file_sha256":hashlib.sha256(raw).hexdigest(),"file_in_zip":match[0]},attempts
        except HTTPError as e:
            attempts.append({"url":url,"http_status":e.code,"error":"HTTP access denied"})
        except (URLError,TimeoutError,OSError,zipfile.BadZipFile) as e:
            attempts.append({"url":url,"error":type(e).__name__})
    return None,None,attempts

def decode(raw):
    text=raw.decode("utf-8-sig")
    reader=csv.DictReader(io.StringIO(text))
    if not COLS.issubset(reader.fieldnames or []):
        raise ValueError("missing original required columns: "+str(sorted(COLS-set(reader.fieldnames or []))))
    records={}
    n=0
    for r in reader:
        n+=1
        patch=r["Patch"].strip()
        if not patch or not r["Year"].strip():
            raise ValueError("empty Patch or Year key")
        try:
            year=int(r["Year"])
            occ=int(r["Occupancy"])
            pl=int(r["Pl"])
            vs=int(r["Vs"])
        except (ValueError,TypeError) as e:
            raise ValueError("missing or noninteger year, occupancy, or host cover") from e
        key=(patch,year)
        if key in records:
            raise ValueError("duplicate Patch x Year key")
        if occ not in (0,1) or pl not in (0,1,2,3) or vs not in (0,1,2,3):
            raise ValueError("unexpected host/occupancy category")
        records[key]={"patch":patch,"year":year,"network":r["Network"].strip(),
                                "occupancy":occ,"pl":pl,"vs":vs}
    return records,n

def count_events(records):
    groups=defaultdict(lambda:{"events":0,"patches":set(),"years":set(),"networks":set(),"extinctions":0})
    risk=0
    for (patch,t),v in records.items():
        before=records.get((patch,t-1))
        after=records.get((patch,t+1))
        if before is None or after is None or v["occupancy"]!=1:
            continue
        risk+=1
        # Only prospective host decline measured by t-1 -> t; next-year outcome is after exposure.
        decline=before["vs"]>v["vs"] and before["vs"]>=1
        g=("decline" if decline else "nondecline", "high" if v["pl"]>=2 else "low")
        x=groups[g]
        x["events"]+=1
        x["patches"].add(patch)
        x["years"].add(t)
        if v["network"]:x["networks"].add(v["network"])
        x["extinctions"]+=int(after["occupancy"]==0)
    output={}
    for k,v in groups.items():
        output["_".join(k)]={"patch_years":v["events"],"unique_patches":len(v["patches"]),
                            "networks":len(v["networks"]),"years":len(v["years"]),
                            "next_year_nonoccupancy":v["extinctions"]}
    a=output.get("decline_high",{})
    b=output.get("decline_low",{})
    checks=[a.get("patch_years",0)>=40,b.get("patch_years",0)>=40,
            a.get("unique_patches",0)>=30,b.get("unique_patches",0)>=30,
            a.get("networks",0)>=5,b.get("networks",0)>=5,
            a.get("years",0)>=3,b.get("years",0)>=3,
            a.get("next_year_nonoccupancy",0)+b.get("next_year_nonoccupancy",0)>=20]
    return {"risk_set_occupied_t_triples":risk,"event_support":output,
            "independent_support_gate_passed":all(checks),
            "gate_checks":checks}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receipt",type=Path,required=True)
    args=ap.parse_args()
    raw, provenance, attempts=acquire()
    receipt={"schema":"chocho_melitaea_host_fallback_source_gate_v01",
             "source_doi":"10.5061/dryad.905qfttrg",
             "attempts":attempts,"source_verified":raw is not None,
             "effect_estimated":False,"ecological_causality_claimed":False,
             "main_GEB_changed":False}
    if raw is not None:
        receipt["source_provenance"]=provenance
        try:
            records,n=decode(raw)
            receipt.update({"status":"SOURCE_VERIFIED","raw_rows":n,
                            "feasibility":count_events(records)})
            if n!=36704:
                receipt["status"]="SOURCE_VERSION_ROW_COUNT_MISMATCH"
                receipt["feasibility"]["independent_support_gate_passed"]=False
        except (UnicodeDecodeError,ValueError) as e:
            receipt.update({"status":"SOURCE_SCHEMA_GATE_FAILED","failure":str(e)})
    else:receipt["status"]="SOURCE_ACCESS_BLOCKED"
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
