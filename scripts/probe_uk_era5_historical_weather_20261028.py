#!/usr/bin/env python3
"""Verify 1976–2021 climate reanalysis availability for UK historical transects.

The example source point is PUBLIC TEST coordinate 52N 1W and NOT any original
LTMN monitoring station. No butterfly data/site coordinates are queried here.
"""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError,URLError

BASE="https://archive-api.open-meteo.com/v1/archive"
PARAMS={
 "latitude":"52.0","longitude":"-1.0",
 "start_date":"1976-01-01","end_date":"2021-12-31",
 "daily":"temperature_2m_mean,precipitation_sum",
 "models":"era5","timezone":"GMT",
 "temperature_unit":"celsius","precipitation_unit":"mm"
}
def check(body):
    if len(body)>5_000_000:raise ValueError("Original reanalysis response too large")
    raw=json.loads(body)
    if "error" in raw:raise ValueError("Public reanalysis says "+str(raw.get("reason","unknown"))[:200])
    daily=raw.get("daily",{})
    dates=daily.get("time",[])
    thermal=daily.get("temperature_2m_mean",[])
    pr=daily.get("precipitation_sum",[])
    expected=(dt.date(2021,12,31)-dt.date(1976,1,1)).days+1
    if not(len(dates)==len(thermal)==len(pr)==expected):
        raise ValueError(f"Original reanalysis calendar incomplete {len(dates)}/{expected}")
    if dates[0]!="1976-01-01" or dates[-1]!="2021-12-31":
        raise ValueError("Source climate dates do not match 46-year requested period")
    valid=[i for i,(t,p) in enumerate(zip(thermal,pr)) if isinstance(t,(float,int)) and
      isinstance(p,(float,int)) and math.isfinite(t) and math.isfinite(p) and -70<t<60 and 0<=p<1000]
    if len(valid)<.95*expected:raise ValueError("Gap/invalid climate exceeds 5%")
    return {"expected_source_days":expected,"valid_original_daily_tmean_and_precip":len(valid),
            "available_weather_variables":list(daily),
            "fixed_ERA5_product":str(raw.get("generationtime_ms",""))[:20] if False else "ERA5 historical reanalysis",
            "no_site_weather_or_original_butterfly_data_exported":True}

def run(path):
    out={"schema":"chocho_original_UK_1976_2021_ERA5_historical_weather_API_gate_v01",
         "original_provider":"Open-Meteo historical archive (ERA5 reanalysis)",
         "source_api":BASE,"test_date_range":"1976-01-01 to 2021-12-31",
         "this_source_point_is_not_LTMN_site":True,
         "source_query_for_reproducibility":BASE+"?"+urlencode(PARAMS)}
    try:
        with urlopen(Request(out["source_query_for_reproducibility"],
                 headers={"User-Agent":"chocho-UK-historical-climate-hindcast/1.0"}),
                 timeout=90) as r:
            body=r.read(5_000_001)
            out["HTTP_status"]=r.status
            out["content_type"]=r.headers.get("Content-Type","")
        out.update(check(body))
        out["original_public_response_sha256"]=hashlib.sha256(body).hexdigest()
        out["status"]="ERA5_1976TO2021_DAILY_SOURCE_COMPLETE"
    except (OSError,ValueError,HTTPError,URLError) as exc:
        out["status"]="SOURCE_MISSING_OR_INCOMPLETE"
        out["failure_type"]=type(exc).__name__
        out["failure_reason"]=str(exc)[:250]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["status"],
           "n_source_days":out.get("expected_source_days"),
           "n_valid":out.get("valid_original_daily_tmean_and_precip"),
           "failure":out.get("failure_reason")},sort_keys=True))
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    r=run(p.parse_args().output)
    if r["status"]!="ERA5_1976TO2021_DAILY_SOURCE_COMPLETE":raise SystemExit(2)
