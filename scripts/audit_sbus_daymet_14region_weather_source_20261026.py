#!/usr/bin/env python3
"""E1b source-only check: 14 rounded PUBLIC survey position classes × Daymet.

This stage does not claim 1-km climate for original monitored sites, only a
seasonal regional proxy at one pixel per rounded public source class.
No source site IDs, public rounded coordinates, or daily site climate values
are exported. Original source dataset and API response integrity audited.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import time
import urllib.parse
import urllib.request
import urllib.error
import zipfile

from analyze_sbus_illinois_niche_counts_20261025 import frozen_focal,load_source,read_events
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import field_panel

ROOT=Path(__file__).resolve().parents[1]
BASE="https://daymet.ornl.gov/single-pixel/api/data"
YEAR_RANGE=range(2000,2021)
EXPECTED_COORD_CLASSES=14
SOURCE_MD5="3161b1c8802a0d67771b2044e18da928"


def original_coarse_coords(source,eligible):
    siteids={s for s,year in eligible}
    mapping=defaultdict(Counter)
    info=Counter()
    with zipfile.ZipFile(source) as z:
        candidates=[m for m in z.infolist()
                    if m.filename.endswith("/4_res/complete-share.csv")
                    and not m.filename.startswith("__MACOSX/")]
        if len(candidates)!=1:
            raise ValueError("Original source CSV identity cannot be proved")
        with z.open(candidates[0]) as b:
            reader=csv.DictReader(io.TextIOWrapper(b,encoding="utf-8-sig",newline=""))
            for row in reader:
                if row["program"]!="Illinois Butterfly Monitoring Network" or row["site"] not in siteids:
                    continue
                label=str(row.get("lat-lon precision") or "").strip()
                if label!="Nearest 1":
                    raise ValueError("Original coarse public precision label changed")
                info["coarse_precision_source_rows"]+=1
                lat=float(row["lat"])
                lon=float(row["lon"])
                if not (math.isfinite(lat) and math.isfinite(lon) and 35<lat<44 and -93<lon<-86):
                    raise ValueError("Original public rounded survey coordinate invalid")
                mapping[row["site"]][(round(lat,3),round(lon,3))]+=1
    if len(mapping)!=121:
        raise ValueError("Original source Illinois site join changed")
    best={site:hist.most_common(1)[0][0] for site,hist in mapping.items()}
    if any(hist[point]/sum(hist.values())<.95 for site,hist in mapping.items() for point in [best[site]]):
        raise ValueError("Site's rounded source location is not stable")
    unique=sorted(set(best.values()))
    if len(unique)!=EXPECTED_COORD_CLASSES:
        raise ValueError(f"Unexpected public coarse location class count: {len(unique)}")
    return best,unique,info


def fetch_daymet_coarse_point(coord):
    params={
        "lat":str(coord[0]),
        "lon":str(coord[1]),
        "vars":"tmax,tmin,prcp",
        "years":",".join(map(str,YEAR_RANGE)),
    }
    url=BASE+"?"+urllib.parse.urlencode(params)
    request=urllib.request.Request(url,headers={
        "Accept":"text/csv",
        "User-Agent":"chocho-source-audit-14-public-rounded-region-weather/0.1",
    })
    with urllib.request.urlopen(request,timeout=65) as r:
        status=r.status
        mime=r.headers.get("Content-Type","")
        raw=r.read(3_000_001)
    sha=hashlib.sha256(raw).hexdigest()
    if status!=200 or len(raw)>3_000_000 or b"<html" in raw[:500].lower():
        raise ValueError("Independent Daymet API did not return bounded CSV")
    text=raw.decode("utf-8-sig")
    lines=text.splitlines()
    head=next((i for i,x in enumerate(lines) if x.lower().startswith("year,yday,")),None)
    if head is None:
        raise ValueError("Daymet API output lacks year/yday CSV schema")
    reader=csv.DictReader(io.StringIO("\n".join(lines[head:])))
    fields=reader.fieldnames or []
    def getcol(key):
        choices=[x for x in fields if x.lower().startswith(key)]
        if len(choices)!=1:
            raise ValueError(f"Missing/ambiguous original Daymet variable {key}: {fields}")
        return choices[0]
    tmax,tmin,prcp=map(getcol,("tmax","tmin","prcp"))
    byyear=defaultdict(list)
    for z in reader:
        yr=int(z["year"])
        day=int(z["yday"])
        vals=tuple(float(z[k]) for k in (tmax,tmin,prcp))
        if yr not in YEAR_RANGE or not 1<=day<=365:
            raise ValueError("Year/out of range Daymet daily response")
        if not all(math.isfinite(v) for v in vals) or vals[2]<0:
            raise ValueError("Missing/negative source model daily weather")
        byyear[yr].append((day,*vals))
    if set(byyear)!=set(YEAR_RANGE):
        raise ValueError(f"Incomplete Daymet year class response: {sorted(byyear)}")
    result={}
    for yr,values in sorted(byyear.items()):
        if len(values)!=365 or len({x[0] for x in values})!=365:
            raise ValueError("Incomplete daily Daymet source year")
        spring=[x for x in values if 60<=x[0]<=181]
        summer=[x for x in values if 182<=x[0]<=273]
        if len(spring)<110 or len(summer)<83:
            raise ValueError("Daymet source seasonal coverage gate failed")
        # Do NOT output regional locations. These are just alternative regional
        # weather forcings and NOT measured host-plant phenology.
        gdd5=sum(max(0,(tmax+tmin)/2-5) for doy,tmax,tmin,precip in spring)
        summer_precip=sum(precip for doy,tmax,tmin,precip in summer)
        result[yr]={"spring_GDD5":gdd5,"summer_precip_mm":summer_precip}
    return result,{"status":"COMPLETE_21_MODEL_YEARS","rows":sum(len(v) for v in byyear.values()),
                   "sha256":sha,"bytes":len(raw),"content_type":mime,
                   "original_daily_weather_not_exported":True}


def stage(source_path):
    focal=frozen_focal()
    events,_=read_events(source_path,focal)
    good,keys,sp,vector,panel=field_panel(events)
    mapping,classes,labels=original_coarse_coords(source_path,keys)
    siteyears_per_class=Counter()
    class_id={c:i for i,c in enumerate(classes)}
    for st,year in keys:
        siteyears_per_class[class_id[mapping[st]]]+=1
    weather={}
    receipts=[]
    for i,point in enumerate(classes):
        try:
            vals,rec=fetch_daymet_coarse_point(point)
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,
                ValueError,OSError,UnicodeError) as e:
            return {"status":"REGIONAL_WEATHER_SOURCE_INCOMPLETE_NO_EFFECT",
                    "successful_region_classes":len(weather),
                    "source_api_error_class":type(e).__name__,
                    "source_api_http_status":e.code if isinstance(e,urllib.error.HTTPError) else None,
                    "source_original_classes":len(classes),
                    "climate_values_never_interpreted_if_incomplete":True}
        weather[i]=vals
        receipts.append({"nonidentifying_region_index":i,
                         "response_sha256":rec["sha256"],"response_bytes":rec["bytes"],
                         "weather_model_days":rec["rows"],
                         "region_siteyears_source_count":siteyears_per_class[i]})
        if i<len(classes)-1:time.sleep(.2)
    valid=sum((int(yr) in weather[class_id[mapping[st]]]) for st,yr in keys)
    return {"status":"COARSE_REGION_YEAR_DAYMET_SOURCE_COMPLETE_NONLOCAL_PROXY",
            "public_rounded_coordinate_classes":len(classes),
            "original_source_monitoring_sites":len(mapping),
            "eligible_siteyears":len(keys),
            "daymet_years_expected":len(YEAR_RANGE),
            "siteyears_with_daily_model_weather":valid,
            "siteyear_coverage_fraction":valid/len(keys),
            "source_14_region_response_hash_receipts":receipts,
            "weather_stats_not_exported":True,
            "original_position_precision_label_only":"Nearest 1",
            "regional_model_pixel_weather_not_site_measurement":True,
            "actual_local_host_plant_phenology_not_measured":True,
            "original_source_qc":panel}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        source=Path(temp)/"henry_us_butterflies.zip"
        load_source(source)
        report=stage(source)
    report.update({
        "schema":"chocho_sbus_illinois_14coarse_region_daymet_source_v01_20261026",
        "original_US_source_doi":"10.6084/m9.figshare.27934602.v1",
        "weather_source":"NASA Daymet historical daily gridded interpolated weather, 10.3334/ORNLDAAC/2129",
        "source_test_only_no_eco_mechanism_estimate":True,
    })
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
