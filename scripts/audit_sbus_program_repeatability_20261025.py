#!/usr/bin/env python3
"""QC programme-level repeatability of the original US butterfly survey counts.

This is a source gate and contains no ecological effects. No site IDs or
underlying survey rows leave this short-lived GitHub Actions computation.
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
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
URL="https://ndownloader.figshare.com/files/50896473"
MD5="3161b1c8802a0d67771b2044e18da928"
ZIP_SIZE=55066096


def get_archive(path):
    md=hashlib.md5(usedforsecurity=False)
    length=0
    request=urllib.request.Request(URL,headers={"User-Agent":"chocho-ecology/1.0", "Accept":"application/zip"})
    with urllib.request.urlopen(request,timeout=85) as res, path.open("wb") as f:
        if res.status!=200:raise ValueError("Unexpected source HTTP status")
        while True:
            c=res.read(1024*1024)
            if not c:break
            length+=len(c)
            if length>60_000_000:raise ValueError("Unexpected source byte budget")
            f.write(c);md.update(c)
    if md.hexdigest()!=MD5 or length!=ZIP_SIZE:
        raise ValueError("Original archive source identity failed")


def focal_names():
    path=ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    with path.open(newline="",encoding="utf-8") as f:
        names={r["species"] for r in csv.DictReader(f)}
    assert len(names)==239
    return names


def audit(path, focal):
    n=defaultdict(Counter)
    event=defaultdict(set)
    site=defaultdict(set)
    siteyear=defaultdict(set)
    focal_names=defaultdict(set)
    years=defaultdict(set)
    zero_species=defaultdict(set)
    site_year_to_events=defaultdict(set)
    site_to_years=defaultdict(set)
    site_focal_to_years=defaultdict(set)
    source_info={}
    with zipfile.ZipFile(path) as z:
        members=[x for x in z.infolist() if x.filename.endswith("/4_res/complete-share.csv")
                 and not x.filename.startswith("__MACOSX/")]
        if len(members)!=1:raise ValueError("Count data path ambiguous")
        with z.open(members[0]) as binary:
            data=io.TextIOWrapper(binary,encoding="utf-8-sig",newline="")
            reader=csv.DictReader(data)
            columns=set(reader.fieldnames or [])
            expected={"program","site","year","Reported species name","event.id","count",
                      "duration","party.minutes","GU code"}
            if not expected.issubset(columns):raise ValueError("Original count columns changed")
            for r in reader:
                pr=str(r["program"]).strip()
                st=str(r["site"]).strip()
                name=str(r["Reported species name"]).strip()
                e=str(r["event.id"]).strip()
                yy=int(float(r["year"]))
                val=float(r["count"])
                if not (2000<=yy<=2020 and math.isfinite(val) and val>=0):
                    raise ValueError("Unexpected monitoring value")
                k=n[pr]
                k["rows"]+=1
                if val>0:k["positive"]+=1
                else:
                    k["explicit_zero"]+=1
                    zero_species[pr].add(name)
                if name in focal:
                    k["rows_original_239"]+=1
                    focal_names[pr].add(name)
                    site_focal_to_years[(pr,st,name)].add(yy)
                if r.get("duration","").strip() not in {"","NA","NaN"}:
                    try:
                        if math.isfinite(float(r["duration"])) and float(r["duration"])>0:k["positive_duration"]+=1
                    except ValueError:pass
                if r.get("party.minutes","").strip() not in {"","NA","NaN"}:
                    try:
                        if math.isfinite(float(r["party.minutes"])) and float(r["party.minutes"])>0:k["positive_party_minutes"]+=1
                    except ValueError:pass
                site[pr].add(st)
                event[pr].add(e)
                years[pr].add(yy)
                siteyear[pr].add((st,yy))
                site_year_to_events[(pr,st,yy)].add(e)
                site_to_years[(pr,st)].add(yy)
    def consecutive3(v):
        return any(y+1 in v and y+2 in v for y in v)
    output=[]
    for pr in sorted(n):
        info=n[pr]
        source_records=n[pr]
        sites=[ys for (name,s),ys in site_to_years.items() if name==pr]
        fsite=[ys for (name,s,sp),ys in site_focal_to_years.items() if name==pr]
        repeats=[len(xs) for (name,s,y),xs in site_year_to_events.items() if name==pr]
        output.append({
            "program":pr,"rows":info["rows"],"n_positive":info["positive"],
            "n_explicit_zero":info["explicit_zero"],
            "zero_fraction":info["explicit_zero"]/info["rows"],
            "n_239_focal_species":len(focal_names[pr]),
            "n_rows_exact_239":info["rows_original_239"],
            "n_species_with_explicit_zero":len(zero_species[pr]),
            "years":[min(years[pr]),max(years[pr])],
            "years_n":len(years[pr]),
            "sites_n":len(site[pr]),
            "events_n":len(event[pr]),
            "site_year_n":len(siteyear[pr]),
            "site_years_with_2plus_events":sum(x>=2 for x in repeats),
            "site_years_with_5plus_events":sum(x>=5 for x in repeats),
            "sites_with_3consecutive_years":sum(consecutive3(x) for x in sites),
            "focal_site_species_with_3consecutive_years":sum(consecutive3(x) for x in fsite),
            "duration_positive_fraction":info["positive_duration"]/info["rows"],
            "party_minutes_positive_fraction":info["positive_party_minutes"]/info["rows"],
            "site_id_values_not_persisted":True,
        })
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    with tempfile.TemporaryDirectory() as d:
        source=Path(d)/"original"
        get_archive(source)
        readout=audit(source,focal_names())
    out={
        "schema":"chocho_sbus_programme_year_count_gate_v01",
        "source_doi":"10.6084/m9.figshare.27934602.v1",
        "source_archive_md5":MD5,
        "source_rights":"CC BY 4.0, Henry et al. (2025), original survey programme contribution attribution required",
        "status":"PROGRAMME_SOURCE_REPEATABILITY_QC_COMPLETE_NOT_ECOLOGICAL_EFFECT",
        "n_programs":len(readout),
        "data":readout,
        "cautions":[
            "Programme sampling protocol and site effort need further source-specific validation",
            "Zero among explicitly reported observations does not imply every omitted row is zero",
            "Adult counts are not larval survival or direct interspecific competition",
            "One programme has many source rows but not necessarily uniform method or effort",
        ]
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["status"],"programme_gate":sorted([
        {"name":d["program"],"year_count":d["years_n"],"n_focal":d["n_239_focal_species"],
         "repeat_sites":d["sites_with_3consecutive_years"],
         "repeat_focal_site_species":d["focal_site_species_with_3consecutive_years"],
         "n_zero":d["n_explicit_zero"]}
        for d in readout],key=lambda z:-z["repeat_focal_site_species"])[:12]},indent=2))


if __name__=="__main__":
    main()
