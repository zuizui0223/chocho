#!/usr/bin/env python3
"""Validate original NASA POWER DAILY 1981–2021 complete point source, not a UKBMS site."""
from __future__ import annotations
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

BASE="https://power.larc.nasa.gov/api/temporal/daily/point"
PARAMS={
 "parameters":"T2M,PRECTOTCORR",
 "community":"AG",
 "latitude":"52.0",
 "longitude":"-1.0",
 "start":"19810101",
 "end":"20211231",
 "time-standard":"UTC",
 "format":"JSON",
}
SOURCE_CALENDAR=(dt.date(2021,12,31)-dt.date(1981,1,1)).days+1

def inspect(payload):
    j=json.loads(payload)
    prop=j.get("properties",{})
    parameter=prop.get("parameter",{})
    if not isinstance(parameter,dict) or not {"T2M","PRECTOTCORR"}<=set(parameter):
        raise ValueError("NASA POWER expected T2M and PRECTOTCORR daily parameters missing")
    timekeys=parameter["T2M"]
    rainkeys=parameter["PRECTOTCORR"]
    if len(timekeys)!=SOURCE_CALENDAR or set(timekeys)!=set(rainkeys):
        raise ValueError(f"NASA POWER daily calendar key set incomplete {len(timekeys)}/{SOURCE_CALENDAR}")
    if min(timekeys)!="19810101" or max(timekeys)!="20211231":
        raise ValueError("NASA POWER daily start/end sources incompatible")
    complete=0
    for k in timekeys:
        t=timekeys[k];p=rainkeys[k]
        if isinstance(t,(int,float)) and isinstance(p,(int,float)) and math.isfinite(t) and math.isfinite(p) and (-80 <t<65) and (0<=p<1500):
            complete+=1
    if complete < .90*SOURCE_CALENDAR:
        raise ValueError("NASA POWER source >10% incomplete or invalid")
    return {
      "original_daily_keys":SOURCE_CALENDAR,
      "valid_original_daily_temperature_and_precipitation":complete,
      "variables":["T2M","PRECTOTCORR"],
      "provider_parameter_headers":sorted(j.get("header",{}).get("sources",[]) if isinstance(j.get("header",{}).get("sources"),list) else []),
      "no_original_UKBMS_site_id_or_weather_series_exported":True,
    }

def run(path):
    r={
      "schema":"chocho_NASA_POWER_1981to2021_original_daily_public_point_source_probe_v01",
      "source":"NASA POWER original daily API",
      "source_url":BASE,
      "source_info":"https://power.larc.nasa.gov/docs/services/api/temporal/daily/",
      "period":"1981-01-01 to 2021-12-31",
      "source_test_point_is_NOT_a_UKBMS_site":True,
      "era5_prior_failure":"ERA5 Open-Meteo 1976-2021 source HTTP429, not biological finding",
      "time_epoch_changed_from_ERA5_proposal":True,
    }
    try:
        url=BASE+"?"+urlencode(PARAMS)
        req=Request(url,headers={"User-Agent":"chocho-UK-original-NASA-POWER-climate-source/1.0",
                                 "Accept":"application/json"})
        with urlopen(req,timeout=95) as f:
            data=f.read(7_000_001)
            r["HTTP_status"]=f.status
            r["content_type"]=f.headers.get("Content-Type","")
        if len(data)>7_000_000:raise ValueError("NASA source payload over bounded download limit")
        r["source_sha256"]=hashlib.sha256(data).hexdigest()
        r.update(inspect(data))
        r["status"]="NASA_POWER_1981_2021_DAILY_SOURCE_VERIFIED"
    except (ValueError,UnicodeError,HTTPError,URLError,OSError) as e:
        r["status"]="NASA_POWER_SOURCE_UNAVAILABLE_OR_INCOMPLETE"
        r["failure_type"]=type(e).__name__
        r["failure_detail"]=str(e)[:300]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({k:r.get(k) for k in ("status","original_daily_keys","valid_original_daily_temperature_and_precipitation","failure_detail")},sort_keys=True))
    return r

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    q=run(p.parse_args().output)
    if q["status"]!="NASA_POWER_1981_2021_DAILY_SOURCE_VERIFIED":raise SystemExit(2)
