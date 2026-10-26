#!/usr/bin/env python3
"""Source-only E0 audit: survey site coordinate suitability + Daymet 1-km API.

Never persist survey coordinates, individual adult counts, site IDs, or actual
site time series. One PUBLIC demonstration location (40,-89) for Daymet test.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import tempfile
import urllib.parse
import urllib.error
import urllib.request
import zipfile

from analyze_sbus_illinois_niche_counts_20261025 import load_source,read_events,frozen_focal
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import field_panel

ROOT = Path(__file__).resolve().parents[1]
DAYMET = "https://daymet.ornl.gov/single-pixel/api/data"

def coordinate_qc(src, eligible_site_keys):
    sites = {site for site, year in eligible_site_keys}
    coords = defaultdict(Counter)
    precision = Counter()
    rows = Counter()
    with zipfile.ZipFile(src) as z:
        cands = [m for m in z.infolist() if not m.filename.startswith("__MACOSX/")
                 and m.filename.endswith("/4_res/complete-share.csv")]
        if len(cands) != 1:
            raise ValueError("Original US survey data CSV not unique")
        with z.open(cands[0]) as byte_stream:
            reader = csv.DictReader(io.TextIOWrapper(byte_stream, encoding="utf-8-sig", newline=""))
            required = {"site", "program", "lat", "lon", "lat-lon precision"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("Original Illinois sample coordinate/precision schema changed")
            for row in reader:
                if (row["program"] != "Illinois Butterfly Monitoring Network"
                    or row["site"] not in sites):
                    continue
                rows["original_source_rows_at_eligible_sites"] += 1
                p = str(row.get("lat-lon precision") or "").strip()
                precision[p or "[missing]"] += 1
                try:
                    lat = float(row["lat"] or "nan")
                    lon = float(row["lon"] or "nan")
                except ValueError:
                    rows["invalid_coordinate_rows"] += 1
                    continue
                if (not math.isfinite(lat) or not math.isfinite(lon)
                    or lat<35 or lat>44 or lon < -93 or lon > -86):
                    rows["invalid_coordinate_rows"] += 1
                    continue
                rows["valid_coordinate_rows"] += 1
                coords[row["site"]][(round(lat, 3), round(lon, 3))] += 1
    stable = 0
    medians = []
    for site in sites:
        observed = coords[site]
        if not observed:
            continue
        common,n=observed.most_common(1)[0]
        n_total=sum(observed.values())
        medians.append(n / n_total)
        if n / n_total >= .95:
            stable += 1
    by_coord=Counter()
    for st, hist in coords.items():
        if hist:
            by_coord[hist.most_common(1)[0][0]] += 1
    duplicates=sum(v-1 for v in by_coord.values() if v>1)
    complete=len(coords)
    result={
        "source_eligible_sites":len(sites),
        "sites_with_valid_coordinates":complete,
        "sites_with_stable_rounded_coords_gte95pct":stable,
        "fraction_stable_sites":stable/len(sites),
        "site_coord_mode_frequency_median":statistics.median(medians) if medians else None,
        "sites_sharing_exact_rounded_coordinates_beyond_first":duplicates,
        "coordinate_distinct_mode_count":len(by_coord),
        "coordinate_precision_labels_counts":dict(precision),
        "original_rows":dict(rows),
        "source_site_values_coordinates_not_saved":True,
    }
    result["site_join_gate"] = (
        "SITE_SOURCE_COORDINATES_POTENTIALLY_LINKABLE"
        if stable >= math.ceil(len(sites)*.8) and duplicates <= math.floor(len(sites)*.2)
        else "SITE_WEATHER_JOIN_NOT_IDENTIFIABLE"
    )
    return result


def daymet_test():
    args={"lat":"40.0","lon":"-89.0","vars":"tmax,tmin,prcp,srad","years":"2015"}
    url=DAYMET + "?" + urllib.parse.urlencode(args)
    result={
        "endpoint":DAYMET,
        "test_coordinates_public_demo_only":True,
        "variables":args["vars"].split(","),
        "test_year":2015,
        "http_status":None,
        "status":"NOT_REQUESTED",
    }
    req=urllib.request.Request(url, headers={
        "Accept":"text/csv,text/plain,application/octet-stream",
        "User-Agent":"chocho-independent-weather-source-availability/0.1",
    })
    try:
        with urllib.request.urlopen(req,timeout=65) as resp:
            result["http_status"]=resp.status
            result["content_type"]=resp.headers.get("Content-Type")
            raw=resp.read(650_001)
        result["response_bytes"]=len(raw)
        result["response_sha256"]=hashlib.sha256(raw).hexdigest()
        if len(raw)>650_000:
            result["status"]="DEMO_RESPONSE_EXCEEDS_LIMIT"
            return result
        if resp.status != 200:
            result["status"]="DEMO_HTTP_NOT_200"
            return result
        s=raw.decode("utf-8-sig",errors="replace")
        if "<html" in s[:1024].lower() or "earthdata login" in s[:1024].lower():
            result["status"]="HTML_OR_AUTH_PAGE_NOT_WEATHER_CSV"
            return result
        lines=s.splitlines()
        headeridx=next((i for i,x in enumerate(lines)
                         if x.lower().startswith("year,yday,")),None)
        if headeridx is None:
            result["status"]="CSV_DAYMET_HEADER_MISSING"
            result["first_lines_schema_only"]=[x[:85] for x in lines[:5]]
            return result
        reader=csv.DictReader(io.StringIO("\n".join(lines[headeridx:])))
        header=reader.fieldnames or []
        rows=list(reader)
        result["columns"]=header
        result["daily_rows"]=len(rows)
        if len(rows) < 360:
            result["status"]="PARTIAL_TEST_YEAR_WEATHER_ROWS"
            return result
        years={r.get("year") for r in rows}
        if years != {"2015"}:
            result["status"]="UNEXPECTED_DAYMET_YEAR_VALUES"
            return result
        result["status"]="DAYMET_PUBLIC_DEMO_365_ROWS_PARSED"
    except urllib.error.HTTPError as e:
        result["status"]="HTTP_ERROR_OR_AUTH_REQUIRED"
        result["http_status"]=e.code
    except (urllib.error.URLError, TimeoutError, OSError, UnicodeError) as e:
        result["status"]="SERVICE_TRANSPORT_UNAVAILABLE"
        result["exception_class"]=type(e).__name__
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    focal=frozen_focal()
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"original_US.zip"
        load_source(p)
        events,_=read_events(p,focal)
        _,keys,species,vector,field=field_panel(events)
        locations=coordinate_qc(p,keys)
    weather=daymet_test()
    status = (
        "READY_TO_PLAN_WEATHER_SITE_YEAR_JOIN"
        if locations["site_join_gate"]=="SITE_SOURCE_COORDINATES_POTENTIALLY_LINKABLE"
        and weather["status"]=="DAYMET_PUBLIC_DEMO_365_ROWS_PARSED"
        else "SOURCE_WEATHER_LINK_NOT_READY"
    )
    out={
        "schema":"chocho_illinois_weather_site_coordinate_qa_v01_20261026",
        "status":status,
        "original_sbus_dataset_doi":"10.6084/m9.figshare.27934602.v1",
        "original_illinois_field_panel":field,
        "coordinate_quality":locations,
        "independent_daily_weather":weather,
        "no_individual_site_coord_or_butterfly_count_export":True,
        "no_abiotic_or_host_niche_effect_estimated":True,
        "disclaimer":"Daymet weather are modelled daily grid fields, not measured local larval host plant phenology",
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
