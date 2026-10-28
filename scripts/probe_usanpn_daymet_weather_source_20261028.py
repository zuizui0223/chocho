#!/usr/bin/env python3
"""Check public ORNL Daymet access and one-year CSV calendar/columns before touching original NPN sites."""
import csv,io,json
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from pathlib import Path
import argparse

URL="https://daymet.ornl.gov/single-pixel/api/data"
QS={"lat":"43.1","lon":"-85.3","vars":"tmax,tmin,prcp","years":"2023"}
SOURCE="https://doi.org/10.3334/ORNLDAAC/2361"

def inspect(body:bytes):
    if len(body)>800000 or not body:raise ValueError("Daymet source response size invalid")
    txt=body.decode("utf-8-sig")
    if "<html" in txt.lower() or "login" in txt[:500].lower():raise ValueError("Daymet returned HTML or authorization page")
    lines=txt.splitlines()
    heads=[i for i,l in enumerate(lines) if "year" in l.lower() and "yday" in l.lower() and "tmax" in l.lower() and "prcp" in l.lower()]
    if len(heads)!=1:raise ValueError("Daymet year/yday/tmin/tmax/prcp header not verified")
    reader=csv.DictReader(io.StringIO("\n".join(lines[heads[0]:])))
    cols=list(reader.fieldnames or [])
    if not all(any(z in c.lower() for c in cols) for z in ("year","yday","tmax","tmin","prcp")):
        raise ValueError("Daymet source columns incomplete")
    data=list(reader)
    if len(data)!=365:raise ValueError(f"Daymet source calendar mismatch: {len(data)}")
    return {"calendar_days":len(data),"columns":cols, "format_lines_before_header":heads[0]}

def run(path):
    result={"schema":"chocho_daymet_endpoint_probe_20261028","original_source":SOURCE,"original_sites_requested":0}
    from urllib.parse import urlencode
    req=Request(URL+"?"+urlencode(QS),headers={"User-Agent":"chocho-usanpn-phenology/1.0","Accept":"text/csv,text/plain,*/*"})
    try:
        with urlopen(req,timeout=45) as r:
            result["status_http"]=r.status
            result["content_type"]=r.headers.get("Content-Type","")
            body=r.read(800001)
        result.update(inspect(body))
        result["status"]="DAYMET_CSV_365DAYS_SOURCE_ACCESS_VERIFIED"
    except (OSError,UnicodeError,ValueError,HTTPError,URLError) as e:
        result["status"]="DAYMET_UNAVAILABLE_OR_UNVERIFIED"
        result["blocker_class"]=type(e).__name__
        result["blocker"]=str(e)[:270]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,sort_keys=True))
    return result["status"]=="DAYMET_CSV_365DAYS_SOURCE_ACCESS_VERIFIED"

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    raise SystemExit(0 if run(a.output) else 2)
