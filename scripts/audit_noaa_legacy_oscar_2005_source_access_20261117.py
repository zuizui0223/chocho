#!/usr/bin/env python3
"""Public, historical OSCAR legacy NOAA ERDDAP u/v data-access smoke test.

NOAA mirror is OLD 1/3-degree, ~5-day OSCAR (1992-2018) — distinct from
the protected NASA 0.25-degree daily OSCAR v2.0. Sampling a valid current
cell is NOT a particle trajectory, coastal arrival, or ecological effect.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import math
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DATASET="yearly_336c_0b32_9cd3"
BASE="https://oceanwatch.pifsc.noaa.gov/erddap/griddap/"
INFO="https://oceanwatch.pifsc.noaa.gov/erddap/info/yearly_336c_0b32_9cd3/index.json"
# NOAA's source metadata: 5-day 1992-2018, z=15m, lat 80:-0.333:-80, lon 20:0.333:420.
TEST_SAMPLES=[
    {"label":"Northwest Pacific offshore", "time_index":930, "depth_index":0, "lat_index":150, "lon_index":390},
    {"label":"Central North Pacific offshore", "time_index":930, "depth_index":0, "lat_index":180, "lon_index":555},
]

def fetch(url: str, max_bytes: int=1_500_000) -> bytes:
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-source-velocity-feasibility/2026.11",
        "Accept":"application/json,text/csv,text/plain,*/*"
    })
    with urllib.request.urlopen(req,timeout=35) as rsp:
        if rsp.status!=200:raise RuntimeError("Non-200 data response")
        payload=rsp.read(max_bytes+1)
        if len(payload)>max_bytes:raise RuntimeError("ERDDAP response larger than expected preflight")
        return payload

def sample(var: str, indices: dict) -> dict:
    dims=[indices[k] for k in ("time_index","depth_index","lat_index","lon_index")]
    query=var+"".join(f"[{x}]" for x in dims)
    # urlencode each brackets is supported by ERDDAP (query is a single constraint expression)
    url=BASE+DATASET+".csv?"+urllib.parse.quote(query,safe="[],")
    rows=list(csv.reader(io.StringIO(fetch(url).decode("utf-8-sig"))))
    if len(rows)<3:raise ValueError(f"Unexpected ERDDAP CSV response for {var}: {rows[:2]}")
    hdr=rows[0]
    if not all(x in hdr for x in ("time","depth","latitude","longitude",var)):
        raise ValueError(f"Actual velocity column absent from ERDDAP {var} CSV: {hdr}")
    vals=rows[-1]
    r=dict(zip(hdr,vals))
    value=float(r[var])
    if not math.isfinite(value) or abs(value)>5:
        raise ValueError(f"Velocity unavailable or physically implausible {var}={value}")
    # Source first says 1992-2018; don't force an exact 2005 day without reading returned timestamp.
    year=int(r["time"][:4])
    if year<2004 or year>2007:
        raise ValueError(f"Returned index year was not a historical 2005-era sample: {r['time']}")
    return {"time":r["time"],"depth":float(r["depth"]),"lat":float(r["latitude"]),
            "lon":float(r["longitude"]),"value_m_per_s":value,"variable":var}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    response={
      "source":"NOAA Fisheries ERDDAP public mirrored legacy OSCAR surface velocities",
      "dataset_id":DATASET,"source_metadata":INFO,
      "source_resolution":"1/3 degree","source_temporal":"approximately 5 days, 1992-2018",
      "not_nasa_v2":True,"status":"SOURCE_NOT_YET_TESTED","samples":[],
      "is_physical_seed_drift_prediction":False,
      "notes":["Only the first two offshore u/v points are queried; no ship routes or seed arrival modeled",
               "The original PO.DAAC third-degree legacy product has been retired; NOAA mirror is not daily 0.25deg v2",
               "A successful public u/v sample unlocks only a future physically directional pilot"]
    }
    try:
        # Read actual official mirror metadata before claiming source identity
        metadata=json.loads(fetch(INFO))
        rows=metadata.get("table",{}).get("rows",[])
        mtxt="\n".join(str(r) for r in rows)
        if not ("OSCAR" in mtxt and "Ocean Surface Zonal Currents" in mtxt and "Ocean Surface Meridional Currents" in mtxt):
            raise ValueError("No expected source metadata and u/v variables")
        for q in TEST_SAMPLES:
            u=sample("u",q)
            v=sample("v",q)
            if (u["time"],u["depth"],u["lat"],u["lon"]) != (v["time"],v["depth"],v["lat"],v["lon"]):
                raise ValueError("u and v samples refer to different source cells")
            response["samples"].append({**{k:q[k] for k in q},
                "time":u["time"],"depth_m":u["depth"],"latitude":u["lat"],
                "longitude":u["lon"],"u_m_per_s":u["value_m_per_s"],
                "v_m_per_s":v["value_m_per_s"]})
        response["status"]="ACTUAL_HISTORICAL_UV_FROM_NOAA_MIRROR"
    except urllib.error.HTTPError as e:
        response["status"]="NOAA_HTTP_ACCESS_ERROR"
        response["error"]={"http_code":e.code,"type":type(e).__name__}
    except (urllib.error.URLError,TimeoutError,OSError,ValueError,RuntimeError) as e:
        response["status"]="NOAA_ACCESS_OR_SOURCE_VALIDATION_UNRESOLVED"
        response["error"]={"type":type(e).__name__,"message":str(e)[:300]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(response,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("NOAA_LEGACY_OSCAR_SOURCE_TEST",json.dumps({
         "status":response["status"],"samples":response["samples"],
         "error":response.get("error"),"not_nasa_v2":True
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
