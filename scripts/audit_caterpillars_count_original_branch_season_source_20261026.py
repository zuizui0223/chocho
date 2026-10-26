#!/usr/bin/env python3
"""Validate original 77.6MB Caterpillars Count branch-level observation source.

Original author code filters and June–July 68,741 survey denominator checked.
No ecological effects or regression fit. No raw survey rows or coordinates
persisted; source-pinned Git SHA and site-level aggregate only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

ORIGINAL_COMMIT="2c5690ac7c6dba7eca506da886c49781c37617bd"
RAW_BLOB="5e9f8c8c98b7e941c53a3ba32baf9be9f6d0692d"
RAW_BYTES=77_610_567
ORIGIN_BLOB="fd14396dd9b83fee12355de115df20caa501ca3d"
EARLY=(91,151)
PEAK=(152,212)
LATE=(213,273)
EXCLUDE={"Coweeta - BB","Coweeta - BS","Coweeta - RK"}
NEEDED={"ID","Name","Region","Year","julianday","Longitude","WetLeaves",
        "PlantFK","Species","sciName","Family","ObservationMethod","Group","Quantity"}


def git_blob(file, expected):
    n=file.stat().st_size
    dig=hashlib.sha1(b"blob "+str(n).encode()+b"\x00")
    with file.open("rb") as f:
        while b:=f.read(1024*1024):
            dig.update(b)
    if dig.hexdigest()!=expected:
        raise ValueError(f"Source blob drifted {file.name}: {dig.hexdigest()} != {expected}")
    return {"name":file.name,"bytes":n,"sha1_git_blob":dig.hexdigest()}


def read_inferred(p):
    store={}
    with p.open(newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            fk=str(row["PlantFK"] or "").strip()
            name=str(row["InferredSciName"] or "").strip()
            try:conf=float(row["NameConfidence"] or "nan")
            except ValueError:conf=float("nan")
            if fk and conf>=2 and name and name.lower()!="nan":
                if fk in store and store[fk]!=name:
                    raise ValueError("Original expert plant inference disagrees at same PlantFK")
                store[fk]=name
    return store


def read_origin(p):
    out={}
    with p.open(newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            name=str(row.get("scientificName") or "").strip()
            origin=str(row.get("plantOrigin") or "").strip().lower()
            if name and origin in {"native","alien"}:
                if name in out and out[name]!=origin:
                    raise ValueError("Conflicting native vs alien scientificName in source")
                out[name]=origin
    return out


def season(day):
    for name,(lo,hi) in (("early",EARLY),("peak",PEAK),("late",LATE)):
        if lo<=day<=hi:
            return name
    return None


def safe_int(s):
    try:
        v=float(s)
        return int(v) if math.isfinite(v) and v.is_integer() else None
    except (TypeError,ValueError):
        return None


def audit(source, inferences, plant_origin):
    surveys={}
    all_source_rows=0
    source_invalid=Counter()
    by_method=Counter()
    year_source=Counter()
    raw_siteids=set()
    source_groups=Counter()
    source_arb=Counter()
    native_ref=Counter()
    for row in csv.DictReader(source.open(newline="",encoding="utf-8-sig")):
        all_source_rows+=1
        if all_source_rows==1:
            if not NEEDED.issubset(row):
                raise ValueError(f"Source 2025-04-18 headers changed: {sorted(NEEDED-set(row))}")
            header=list(row)
        key=str(row.get("ID") or "").strip()
        if not key:
            source_invalid["missing_survey_ID"]+=1
            continue
        site=str(row.get("Name") or "").strip()
        year=safe_int(row.get("Year"))
        doy=safe_int(row.get("julianday"))
        family=str(row.get("Family") or "").strip()
        method=str(row.get("ObservationMethod") or "").strip()
        source_species=str(row.get("sciName") or "").strip()
        plant_user=str(row.get("Species") or "").strip()
        if plant_user=="N/A":
            source_species=inferences.get(str(row.get("PlantFK") or "").strip(),source_species)
        if source_species.lower() in {"na","n/a","nan","none",""}:
            source_species=""
        try:lon=float(row.get("Longitude") or "nan")
        except ValueError:lon=float("nan")
        flag=str(row.get("WetLeaves") or "").strip().lower()
        wet=(flag in {"true","t","1","yes","y"})
        if flag not in {"true","false","t","f","1","0","yes","no","y","n"}:
            source_invalid["unexpected_WetLeaves_flag"]+=1
        valid=(bool(site and family and source_species and method)
               and year is not None and 2010<=year<=2024
               and doy is not None and 1<=doy<=366
               and math.isfinite(lon) and lon>=-100
               and not wet and site not in EXCLUDE)
        if not valid:
            source_invalid["rows_filtered_by_original_source_rules"]+=1
            continue
        original_season=season(doy)
        origin=plant_origin.get(source_species)
        if key not in surveys:
            surveys[key]={
                "site":site,"year":year,"family":family,"method":method,
                "doy":doy,"season":original_season,
                "species":source_species,"origin":origin,
                "positive":False,"source_rows":0
            }
        x=surveys[key]
        if any(x[k]!=value for k,value in
               (("site",site),("year",year),("family",family),
                ("method",method),("doy",doy),("species",source_species),
                ("origin",origin))):
            source_invalid["duplicate_ID_conflicting_metadata"]+=1
        x["source_rows"]+=1
        if str(row.get("Group") or "")=="caterpillar":
            q=row.get("Quantity")
            try:quantity=float(q)
            except (TypeError,ValueError):quantity=float("nan")
            if math.isfinite(quantity) and quantity>0:
                x["positive"]=True
            elif not math.isfinite(quantity) or quantity<0:
                source_invalid["caterpillar_invalid_quantity"]+=1
        source_groups[str(row.get("Group") or "")]+=1
    original_june_july=sum(s["season"]=="peak" for s in surveys.values())
    survey_seasons=Counter()
    season_origin=defaultdict(Counter)
    sample_by_method=Counter()
    matched=defaultdict(lambda:defaultdict(lambda:Counter()))
    for sid,s in surveys.items():
        se=s["season"]
        if se is None:
            continue
        survey_seasons[se]+=1
        sample_by_method[(se,s["method"])]+=1
        if s["origin"] in {"native","alien"}:
            season_origin[se][s["origin"]+"_surveys"]+=1
            season_origin[se][s["origin"]+"_caterpillar_positive"]+=int(s["positive"])
            k=(s["site"],s["year"],s["family"],s["method"],se)
            matched[se][k][s["origin"]+"_n"]+=1
            if s["positive"]:
                matched[se][k][s["origin"]+"_positive"]+=1
    matched_season={}
    ready=True
    for se in ("early","peak","late"):
        strata={k:v for k,v in matched[se].items()
                if v.get("native_n",0)>0 and v.get("alien_n",0)>0}
        st_sites={k[0] for k in strata}
        totals={
            "matched_strata":len(strata),
            "matched_distinct_sites":len(st_sites),
            "matched_botanical_families":len({k[2] for k in strata}),
            "matched_native_surveys":sum(v["native_n"] for v in strata.values()),
            "matched_alien_surveys":sum(v["alien_n"] for v in strata.values()),
            "matched_native_positive":sum(v["native_positive"] for v in strata.values()),
            "matched_alien_positive":sum(v["alien_positive"] for v in strata.values()),
        }
        accept=(totals["matched_strata"]>=30 and totals["matched_distinct_sites"]>=10 and
                totals["matched_native_surveys"]>=200 and totals["matched_alien_surveys"]>=200 and
                totals["matched_native_positive"]+totals["matched_alien_positive"]>=10)
        totals["prespecified_season_source_gate_pass"]=accept
        matched_season[se]=totals
        ready&=accept
    if original_june_july!=68_741:
        ready=False
        source_invalid["published_author_68741_survey_denominator_not_reproduced"]=1
    if source_invalid["duplicate_ID_conflicting_metadata"] or source_invalid["caterpillar_invalid_quantity"]:
        ready=False
    return {
        "status":"ALL_3_SEASONS_SOURCE_GATE_QUALIFIED" if ready else "SOURCE_SEASON_OR_IDENTITY_GATE_NOT_QUALIFIED",
        "input_raw_source_rows":all_source_rows,
        "header":header[:90],
        "original_author_peak_DOY_152_212_survey_count":original_june_july,
        "published_original_author_peak_survey_count_target":68741,
        "source_usable_surveys_all_year":len(surveys),
        "source_survey_season_counts":dict(survey_seasons),
        "source_origin_surveys_by_season":{k:dict(v) for k,v in season_origin.items()},
        "survey_method_by_season":{"|".join(k):v for k,v in sample_by_method.items()},
        "matched_strata_by_season":matched_season,
        "source_exclusion_and_quality_counts":dict(source_invalid),
        "source_plant_names_inferred_high_confidence_n":len(inferences),
        "source_plant_origin_names_n":len(plant_origin),
        "original_survey_ids_sites_coordinates_raw_not_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=args.original_repo
    raw=root/"data/fullDataset_2025-04-18.csv"
    origin=root/"data/plant_origin_status.csv"
    infer=root/"data/inferredPlantNames_2024-08-20.csv"
    receipts=[git_blob(raw,RAW_BLOB),git_blob(origin,ORIGIN_BLOB)]
    if raw.stat().st_size!=RAW_BYTES:
        raise ValueError("Original fullDataset size drift")
    result=audit(raw,read_inferred(infer),read_origin(origin))
    result["schema"]="chocho_2026_caterpillars_count_season_native_alien_source_gate_v01"
    result["source_repo_commit"]=ORIGINAL_COMMIT
    result["source_receipts"]=receipts
    result["independent_evidence_category"]="Observed caterpillar Lepidoptera on branch survey, NOT butterfly species population"
    result["source_previous_study"]="Whitener et al (2026) Ecosphere DOI 10.1002/ecs2.70615"
    result["source_previous_analysis_reuses_JuneJuly"]=True
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in (
        "status","input_raw_source_rows",
        "original_author_peak_DOY_152_212_survey_count",
        "source_usable_surveys_all_year",
        "source_survey_season_counts","source_origin_surveys_by_season",
        "matched_strata_by_season","source_exclusion_and_quality_counts")},indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
