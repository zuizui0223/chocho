#!/usr/bin/env python3
"""Check NASA POWER original response-body drift vs daily numeric stability.

No site ID, coordinates, original daily values, plant host rows or butterfly
observations are serialized. This is a source diagnosis, not model evaluation.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,time
from pathlib import Path
from pyproj import Transformer
import audit_uk_ltmn_NASA_POWER_historical_baseline_20261028 as nasa
from audit_uk_ltmn_historical_event_join_20261028 import excel_rows,get_source,FIELDS_VISITS
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk

SCHEMA="chocho_original_NASA_19_UKBMS_sites_body_vs_numeric_fingerprint_20261029_v01"
OLD="051c82f9af21941bee3ea2fa9d8fa739990a36ed6e25e43441c25821e257b806"

def source_numeric_sha(weather):
    out=hashlib.sha256()
    for day in sorted(weather):
        t,p=weather[day]
        out.update(f"{day.isoformat()}|{t:.8f}|{p:.8f}\n".encode())
    return out.hexdigest()

def audit(legacy):
    old=json.loads(legacy.read_text())
    if old.get("source",{}).get("n_NASA_POWER_source_transects")!=19 or        old.get("source",{}).get("NASA_POWER_days_per_site")!=14975 or        old.get("source",{}).get("source_original_time_split_early")!=[1981,1995] or        old.get("source",{}).get("source_original_time_split_late")!=[2005,2021]:
        raise ValueError("Original UKBMS 19-site NASA initial source receipt is not valid")
    saved=old["source"]["original_NASA_POWER_source_response_sha256_sorted_unlinked_to_locations"]
    if len(saved)!=19 or hashlib.sha256("\n".join(sorted(saved)).encode()).hexdigest()!=OLD:
        raise ValueError("Original previous raw NASA metadata sha differs")
    cross,qc=read_crosswalk()
    original,sites,site_events=nasa.nasa_site_effort(
        excel_rows(get_source("visits"),FIELDS_VISITS),cross)
    if len(sites)!=19:raise ValueError("Original NASA 1981–95 site sample changed")
    transformer=Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True)
    current=[];numeric=[];repeated={}
    for i,site in enumerate(sorted(sites),1):
        east,north,_=cross[site]
        lon,lat=transformer.transform(east,north)
        weather,sha,invalid=nasa.original_nasa_power_site_weather(lat,lon)
        if invalid or len(weather)!=14975:raise ValueError("Original site NASA daily numeric values incomplete")
        current.append(sha)
        numeric.append(source_numeric_sha(weather))
        if i==1:
            # One matched original site twice, no site coordinates or values in output.
            second,s2,bad=nasa.original_nasa_power_site_weather(lat,lon)
            repeated={
                "one_site_repeat_raw_SHA_identical":sha==s2,
                "one_site_repeat_daily_numeric_SHA_identical":source_numeric_sha(weather)==source_numeric_sha(second),
                "one_site_repeat_same_14975_daily_records":len(second)==len(weather)==14975,
            }
        if i%5==0:print(json.dumps({"NASA_source_original_site_probes":i}),flush=True)
        time.sleep(.4)
    prior=set(saved);new=set(current)
    result={
        "schema":SCHEMA,"status":"CURRENT_SOURCE_NUMERIC_STABILITY_DIAGNOSED",
        "old_19site_raw_source_aggregate_SHA256":OLD,
        "new_19site_raw_source_aggregate_SHA256":hashlib.sha256("\n".join(sorted(current)).encode()).hexdigest(),
        "new_19site_daily_numeric_aggregate_SHA256":hashlib.sha256("\n".join(sorted(numeric)).encode()).hexdigest(),
        "n_identical_original_raw_19_site_responses":len(prior&new),
        "n_original_UKBMS_georeferenced_sites":len(sites),
        "NASA_source_days_per_site":14975,
        "original_historical_time_cutoff":"1981–2021",
        "repeat_within_run":repeated,
        "old_original_daily_numeric_checksum_was_NOT_persisted":True,
        "older_source_daily_value_equivalence_cannot_be_asserted":True,
        "no_original_UKBMS_sites_daily_series_or_observers_exported":True,
    }
    if not (repeated["one_site_repeat_daily_numeric_SHA_identical"]
           and repeated["one_site_repeat_same_14975_daily_records"]):
        result["status"]="SOURCE_UNSTABLE_WITHIN_SAME_RUN"
    return result

def run(old,output):
    try: result=audit(old)
    except (OSError,ValueError,KeyError) as e:
        result={"schema":SCHEMA,"status":"ORIGINAL_SOURCE_AUDIT_FAILED",
                "blocker_class":type(e).__name__,"blocker":str(e)[:350],
                "no_original_sites_exported":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,sort_keys=True))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-result",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    r=run(a.original_result,a.output)
    if r["status"]!="CURRENT_SOURCE_NUMERIC_STABILITY_DIAGNOSED":raise SystemExit(2)
