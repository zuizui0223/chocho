#!/usr/bin/env python3
"""Verify full CC-BY original U.S. butterfly survey-count source without publishing rows.

Requires original Figshare exact archive MD5 and dataset paper schema.
No demographic or host-sharing effect computed.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import urllib.request
import urllib.error
import zipfile

ROOT=Path(__file__).resolve().parents[1]
URL="https://ndownloader.figshare.com/files/50896473"
SOURCE_BYTES=55066096
SOURCE_MD5="3161b1c8802a0d67771b2044e18da928"
DATASET_DOI="10.6084/m9.figshare.27934602.v1"
MAX_BYTES=70_000_000
FROZEN_BLOB="9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"


def verify_names():
    path=ROOT/"data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    r=path.read_bytes()
    if hashlib.sha1(b"blob "+str(len(r)).encode()+b"\x00"+r).hexdigest()!=FROZEN_BLOB:
        raise ValueError("Focal 239 changed")
    with path.open(newline="",encoding="utf-8") as file:
        names=[row["species"].strip() for row in csv.DictReader(file)]
    if len(names)!=239 or len(set(names))!=239:
        raise ValueError("Not original 239 distinct focal species")
    return set(names)


def fetch(tmp:Path):
    result={"status":"NOT_RETRIEVED","source_url":URL}
    request=urllib.request.Request(URL,headers={
        "Accept":"application/zip, application/octet-stream",
        "User-Agent":"chocho-butterfly-community-ecology-source-gate/1.0",
    })
    sha=hashlib.sha256()
    md5=hashlib.md5(usedforsecurity=False)
    try:
        with urllib.request.urlopen(request,timeout=90) as resp:
            result["http_status"]=resp.status
            result["content_type"]=resp.headers.get("Content-Type")
            if resp.status!=200:
                result["status"]="HTTP_NOT_OK"
                return result
            total=0
            with tmp.open("wb") as f:
                while True:
                    data=resp.read(1024*1024)
                    if not data:
                        break
                    total+=len(data)
                    if total>MAX_BYTES:
                        result["status"]="OVERSIZE_OR_DRIFTED_ARCHIVE"
                        return result
                    f.write(data)
                    sha.update(data)
                    md5.update(data)
        result.update(
            status="RETRIEVED",bytes=total,sha256=sha.hexdigest(),
            md5=md5.hexdigest(),
        )
        if md5.hexdigest()!=SOURCE_MD5 or total!=SOURCE_BYTES:
            result["status"]="SOURCE_MD5_OR_BYTE_SIZE_MISMATCH"
    except urllib.error.HTTPError as e:
        result.update(status="HTTP_ERROR",http_status=e.code)
    except (urllib.error.URLError,TimeoutError,OSError) as e:
        result.update(status="TRANSPORT_ERROR",type=type(e).__name__)
    return result


def sample_csv(z,member,known):
    stat={"member":member.filename,"uncompressed_bytes":member.file_size,
          "is_raw_site_source_persisted":False}
    try:
        with z.open(member,"r") as f:
            handle=io.TextIOWrapper(f,encoding="utf-8-sig",errors="replace",newline="")
            reader=csv.DictReader(handle)
            stat["header"]=reader.fieldnames or []
            if not reader.fieldnames:
                stat["status"]="NO_HEADER"
                return stat
            # Only read main record file. Metadata on taxonomic crosswalk stays in source ZIP.
            col= {s.lower().replace(" ","").replace("_","").replace(".",""):s
                  for s in reader.fieldnames}
            def lookup(*patterns):
                for p in patterns:
                    if p in col:
                        return col[p]
                return None
            keys={
                "program":lookup("program"),
                "site":lookup("site","siteid","sitecode"),
                "event":lookup("eventid","surveyid"),
                "year":lookup("year"),
                "species":lookup("reportedspeciesname","species","scientificname","speciesname"),
                "GU_code":lookup("gucode"),
                "count":lookup("count","abundance"),
                "duration":lookup("duration"),
                "party_minutes":lookup("partyminutes"),
            }
            stat["inferred_columns"]=keys
            if any(keys[k] is None for k in ("program","site","event","year","species","count")):
                stat["status"]="UNKNOWN_RECORD_SCHEMA_NO_INTERPRETATION"
                return stat
            programmes=collections.Counter()
            years=collections.Counter()
            names=set()
            sites=set()
            events=set()
            siteyear=set()
            names_239=set()
            pos=zero=negative=bad=missing_count=0
            n_rows=0
            valid_count=0
            finite_duration=finite_party=0
            for row in reader:
                n_rows+=1
                prog=str(row[keys["program"]] or "").strip()
                site=str(row[keys["site"]] or "").strip()
                event=str(row[keys["event"]] or "").strip()
                sp=str(row[keys["species"]] or "").strip()
                year=str(row[keys["year"]] or "").strip()
                value=str(row[keys["count"]] or "").strip()
                if prog:programmes[prog]+=1
                if sp:
                    names.add(sp)
                    if sp in known:
                        names_239.add(sp)
                if prog and site:
                    sites.add((prog,site))
                if event and prog:events.add((prog,event))
                if prog and site and year:
                    try: y=int(float(year))
                    except ValueError: y=-1
                    if 1990<=y<=2026:
                        years[y]+=1
                        siteyear.add((prog,site,y))
                if not value or value.lower() in ("nan","na","null","none"):
                    missing_count+=1
                else:
                    try:v=float(value)
                    except ValueError:
                        bad+=1
                        continue
                    if not math.isfinite(v):
                        bad+=1
                        continue
                    valid_count+=1
                    if v>0:pos+=1
                    elif v==0:zero+=1
                    else:negative+=1
                for k,m in (("duration","duration"),("party_minutes","party")):
                    key=keys[k]
                    if key and row.get(key):
                        try:
                            val=float(row[key])
                            if math.isfinite(val) and val>0:
                                if m=="duration":finite_duration+=1
                                else:finite_party+=1
                        except ValueError:
                            pass
            stat.update(
                status="RECORD_SCHEMA_VERIFIED",
                n_rows=n_rows,n_years=len(years),
                years_range=[min(years),max(years)] if years else None,
                program_rows=dict(programmes.most_common()),
                n_programmes=len(programmes),
                n_distinct_site_program=len(sites),
                n_distinct_survey_events=len(events),
                n_distinct_site_year=len(siteyear),
                n_distinct_species_name=len(names),
                n_exact_239_matched_species=len(names_239),
                original_239_exact_species_names=sorted(names_239),
                count_positive=pos,count_zero=zero,count_negative=negative,
                count_missing=missing_count,count_invalid=bad,
                n_rows_with_valid_count=valid_count,
                duration_positive_rows=finite_duration,
                party_minutes_positive_rows=finite_party,
                year_row_counts=dict(sorted(years.items())),
                site_ids_dates_and_counts_not_persisted=True,
                diagnostic_limitation="CSV may omit zero counts, whose interpretation needs survey/site membership analysis",
            )
    except (zipfile.BadZipFile,UnicodeError,csv.Error,OSError) as e:
        stat.update(status="DATA_PARSING_ERROR",exception_class=type(e).__name__)
    return stat


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    names=verify_names()
    result={
        "schema":"chocho_us_sbus_original_survey_count_source_gate_20261025_v01",
        "source_paper_doi":"10.1038/s41597-025-05513-8",
        "source_dataset_doi":DATASET_DOI,
        "source_license":"CC BY 4.0",
        "source_archive_figshare_id":50896473,
        "original_239_focal_blob":FROZEN_BLOB,
        "no_community_or_competition_effect_estimated":True,
        "source_programme_caution":"Circle counts, casual trips and structured transects must not be pooled without protocol/effort controls",
    }
    with tempfile.TemporaryDirectory() as d:
        arch=Path(d)/"sbus_original"
        response=fetch(arch)
        result["transport"]=response
        if response["status"]!="RETRIEVED":
            result["status"]="SOURCE_ACCESS_OR_IDENTITY_FAILED_NO_BIOLOGICAL_INFERENCE"
        elif not zipfile.is_zipfile(arch):
            result["status"]="ORIGINAL_BINARY_NOT_ZIP"
        else:
            with zipfile.ZipFile(arch) as z:
                bad=z.testzip()
                if bad:
                    result["status"]="ORIGINAL_ZIP_CHECKSUM_FAILED"
                else:
                    members=[m for m in z.infolist() if not m.is_dir()]
                    result["members"]=[
                        {"name":m.filename[:180],"uncompressed_bytes":m.file_size}
                        for m in members
                        if m.filename.lower().endswith((".csv",".txt",".zip",".json"))
                    ][:300]
                    candidates=[m for m in members if m.filename.lower().endswith(".csv")
                        and ("complete-share" in m.filename.lower() or
                             "complete share" in m.filename.lower())]
                    if not candidates:
                        result["status"]="NO_COMPLETE_SHARE_SOURCE_CSV_VERIFIED"
                    elif len(candidates)!=1:
                        result["status"]="AMBIGUOUS_PRIMARY_COUNT_FILES"
                    else:
                        doc=sample_csv(z,candidates[0],names)
                        result["count_data_schema"]=doc
                        result["status"]=("COUNT_ROWS_AVAILABLE_TAXONOMY_AND_EFFORT_REVIEW_PENDING"
                            if doc.get("status")=="RECORD_SCHEMA_VERIFIED"
                            else "COUNT_SCHEMA_NOT_VERIFIED")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "transport":result["transport"]["status"],
        "bytes":result["transport"].get("bytes"),
        "count_data":{key:result.get("count_data_schema",{}).get(key)
        for key in ("status","n_rows","n_distinct_survey_events","n_distinct_site_program",
                    "n_exact_239_matched_species","years_range","n_programmes")},
        "possible_members":[x["name"] for x in result.get("members",[])[:30]],
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
