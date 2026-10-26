#!/usr/bin/env python3
"""Independent monarch egg/instar × Asclepias phenology source gate.

Original Dryad DOI 10.5061/dryad.4j0zpc8q5, published 8 Oct 2026.
No survey-level values, transect names, site coordinates or records committed.
No effect fitted. See source-first protocol.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
import zipfile

DOI="10.5061/dryad.4j0zpc8q5"
ENC="doi%3A10.5061%2Fdryad.4j0zpc8q5"
BASE="https://datadryad.org/api/v2"
ZIP_URL=BASE+"/datasets/"+ENC+"/download"
META_URL=BASE+"/datasets/"+ENC
MAX_ZIP=12_000_000
MAIN="MS_western_monarchs_Oikos_2024.csv"
PLANTS=("ASCO","ASCR","ASER","ASER2","ASFA","ASIN","ASSP")
STAGES=("Eggs","1","2","3","4","5")
PHENO=("Vegetative","Flowering","Fruiting","Senescent")


def http_binary(url,limit):
    request=Request(url,headers={
        "Accept":"application/zip,application/octet-stream,application/json,*/*",
        "User-Agent":"chocho-independent-2026-monarch-stage-source-gate/1.0",
    })
    with urlopen(request, timeout=90) as response:
        code=response.status
        mime=response.headers.get("Content-Type","")
        final=response.url.split("?")[0]
        raw=response.read(limit+1)
    if len(raw)>limit:
        raise ValueError(f"Source exceeds {limit} bytes")
    return raw,{"http_status":code,"mime":mime,"resolved_url":final,
                "byte_length":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}


def digest_metadata():
    result={"url":META_URL,"status":"NOT_RETRIEVED"}
    try:
        data,receipt=http_binary(META_URL,1_000_000)
        result.update(receipt)
        obj=json.loads(data.decode("utf-8-sig"))
        result["status"]="DATASET_METADATA_RECEIVED"
        result["identifier"]=obj.get("identifier")
        result["title"]=obj.get("title")
        result["versionNumber"]=obj.get("versionNumber")
        result["curationStatus"]=obj.get("curationStatus")
        result["publicationDate"]=obj.get("publicationDate")
        result["dataDigestNotPublished"]=True
    except (HTTPError,URLError,TimeoutError,ValueError,OSError,UnicodeError,json.JSONDecodeError) as e:
        result["status"]="DATASET_METADATA_ERROR"
        result["error_class"]=type(e).__name__
        result["http_error"]=getattr(e,"code",None)
    return result


def parse_num(row,name,stats,integer=True):
    val=row.get(name)
    if val is None:
        return None
    value=str(val).strip()
    if value=="" or value.lower() in {"na","null","nan","n/a","none"}:
        stats[name]["missing"]+=1
        return None
    try:
        x=float(value)
    except ValueError:
        stats[name]["bad"]+=1
        return None
    if not math.isfinite(x) or x<0 or (integer and not x.is_integer()):
        stats[name]["invalid"]+=1
        return None
    if x==0:
        stats[name]["zero"]+=1
    else:
        stats[name]["positive"]+=1
    return int(x) if integer else x


def inspect_csv(raw):
    txt=raw.decode("utf-8-sig",errors="strict")
    reader=csv.DictReader(io.StringIO(txt,newline=""))
    headers=reader.fieldnames or []
    if not headers or len(headers)!=len(set(headers)):
        return {"status":"SOURCE_HEADER_ABSENT_OR_DUPLICATED","headers":headers[:170]}
    fields={name.casefold():name for name in headers}
    def field(name):
        return fields.get(name.casefold())
    expected=["Transect_ID","Region","Site","Year","Survey_Round","Survey_Type"]
    missing=[x for x in expected if field(x) is None]
    if missing:
        return {"status":"CRITICAL_TRANSECT_SAMPLE_SCHEMA_MISSING","missing":missing,
                "headers":headers[:170]}
    vars={}
    for p in PLANTS:
        vars[p]={}
        for label in list(STAGES)+list(PHENO)+["Total_Available_Stems","Total_Stems","Total"]:
            key=field(p+"_"+label)
            if key:
                vars[p][label]=key
    schema_by_plant={p:sorted(d) for p,d in vars.items()}
    completeness={
        p:{"has_egg":bool("Eggs" in cols),
           "has_5_instars":all(i in cols for i in ("1","2","3","4","5")),
           "has_4_pheno":all(i in cols for i in PHENO),
           "has_available_stems":"Total_Available_Stems" in cols}
        for p,cols in vars.items()
    }
    stats=defaultdict(Counter)
    years=Counter()
    regions=Counter()
    methods=Counter()
    plant_units=Counter()
    plant_positive_late_units=Counter()
    plant_positive_egg_units=Counter()
    n_rows=0
    bytransect={}
    regionyear=set()
    grouped_count=Counter()
    group_plant_seen=set()
    grouped_plant_positives=defaultdict(set)
    observed_total=Counter()
    invalid_row=Counter()
    missing_fields=Counter()
    origin_keys=Counter()
    method_units=defaultdict(set)
    plant_by_method=Counter()
    positives_in_plant_count=Counter()
    stems_positive=Counter()
    multiple_plant_groups=defaultdict(set)
    inconsistent_stem_total=Counter()
    valid_stage_compare=Counter()
    conflicting_stage_compare=Counter()
    for row in reader:
        n_rows+=1
        try:
            year=int(float(row[field("Year")]))
            round_number=int(float(row[field("Survey_Round")]))
        except (ValueError,TypeError):
            invalid_row["invalid_year_or_round"]+=1
            continue
        if year not in {2017,2018,2019,2023,2024} or round_number not in range(1,9):
            invalid_row["out_of_protocol_year_or_round"]+=1
            continue
        reg=str(row[field("Region")] or "").strip()
        st=str(row[field("Site")] or "").strip()
        transect=str(row[field("Transect_ID")] or "").strip()
        method=str(row[field("Survey_Type")] or "").strip()
        if not all((reg,st,transect,method)):
            invalid_row["missing_source_unit_id"]+=1
            continue
        nkey=(transect,reg,year,round_number,method)
        rowkey=(transect,reg,year,round_number,method,
                row.get(field("Meters_along_Transect"),""))
        origin_keys[rowkey]+=1
        bytransect[transect]=bytransect.get(transect,reg)
        regionyear.add((reg,year))
        years[year]+=1
        regions[reg]+=1
        methods[method]+=1
        grouped_count[nkey]+=1
        method_units[method].add(nkey)
        valid_species_present=0
        pl_sum_eggs=0
        pl_sum_larvae=Counter()
        for p,keys in vars.items():
            if not (completeness[p]["has_available_stems"] and
                    completeness[p]["has_egg"] and
                    completeness[p]["has_5_instars"]):
                continue
            stages={k:parse_num(row,keys[k],stats) for k in STAGES}
            pheno={k:parse_num(row,keys[k],stats) for k in PHENO if k in keys}
            available=parse_num(row,keys["Total_Available_Stems"],stats)
            if available is None:
                missing_fields["missing_stem_availability"]+=1
                continue
            if available>0:
                stems_positive[p]+=1
                group_plant_seen.add((nkey,p))
                plant_by_method[(method,p)]+=1
                multiple_plant_groups[nkey].add(p)
                valid_species_present+=1
            if available==0 and any(v and v>0 for v in stages.values()):
                invalid_row["animal_on_zero_available_host"]+=1
            if available>0 and any(v is not None for v in stages.values()):
                plant_units[p]+=1
                early=sum(stages[x] or 0 for x in ("1","2"))
                late=sum(stages[x] or 0 for x in ("3","4","5"))
                eggs=stages["Eggs"] or 0
                observed_total[(p,"Eggs")]+=eggs
                observed_total[(p,"Early_L1_L2")]+=early
                observed_total[(p,"Late_L3_L5")]+=late
                observed_total[(p,"Available_Stems")]+=available
                if eggs>0:
                    plant_positive_egg_units[p]+=1
                    grouped_plant_positives[(nkey,p)].add("egg")
                if late>0:
                    plant_positive_late_units[p]+=1
                    grouped_plant_positives[(nkey,p)].add("late")
            if all(pheno.get(x) is not None for x in PHENO):
                source_all=sum(pheno[x] for x in PHENO)
                # If available stems must be ≤ all stems (including senescent)
                if source_all<available:
                    inconsistent_stem_total["pheno_less_than_available"]+=1
                if "Total" in keys:
                    source_total=parse_num(row,keys["Total"],stats)
                    if source_total is not None:
                        valid_stage_compare["pheno_total_checked"]+=1
                        if abs(source_total-source_all)>0:
                            inconsistent_stem_total["pheno_sum_vs_reported_total"]+=1
            if all(v is not None for v in stages.values()):
                pl_sum_eggs+=stages["Eggs"]
                for k in ("1","2","3","4","5"):
                    pl_sum_larvae[k]+=stages[k]
        if "Sum_Eggs" in fields.values():
            summary_col=field("Sum_Eggs")
            value=parse_num(row,summary_col,stats)
            if value is not None:
                valid_stage_compare["egg_total_checked"]+=1
                if value != pl_sum_eggs:
                    conflicting_stage_compare["egg_total_mismatch"]+=1
        for k in ("1","2","3","4","5"):
            summary_col=field("Sum_L"+k)
            if summary_col:
                value=parse_num(row,summary_col,stats)
                if value is not None:
                    valid_stage_compare["L"+k+"_checked"]+=1
                    if value!=pl_sum_larvae[k]:
                        conflicting_stage_compare["L"+k+"_mismatch"]+=1
    groupcount=len(grouped_count)
    positive_egg_groups={key for key,positive in grouped_plant_positives.items() if "egg" in positive}
    positive_late_groups={key for key,positive in grouped_plant_positives.items() if "late" in positive}
    eligible=np=len(group_plant_seen)
    collision=sum(n-1 for n in origin_keys.values() if n>1)
    site_groups_with_2hosts=sum(len(v)>=2 for v in multiple_plant_groups.values())
    result={
        "status":"ORIGINAL_CSV_PARSED_QC_PENDING",
        "original_csv_rows":n_rows,
        "headers":headers[:210],
        "header_count":len(headers),
        "date_field_present":field("Date") is not None,
        "strict_required_plant_columns":completeness,
        "plant_header_groups":schema_by_plant,
        "regions_row_counts":dict(sorted(regions.items())),
        "years_row_counts":dict(sorted(years.items())),
        "methods_row_counts":dict(sorted(methods.items())),
        "distinct_transects":len(bytransect),
        "distinct_regions_by_year":len(regionyear),
        "transect_round_sampling_units":groupcount,
        "species_x_transect_round_units_with_available_stems":eligible,
        "transect_round_units_at_least_two_plant_species_available":site_groups_with_2hosts,
        "exact_transect_round_slice_row_key_duplicate_extra_rows":collision,
        "exact_transect_round_slice_row_key_missing_marker_handling":"duplicate keys diagnostic only; not stage biological replicates",
        "plant_taxon_source_row_counts_with_positive_available_stems":dict(stems_positive),
        "plant_taxon_units_with_available_stems_and_some_stage_fields":dict(plant_units),
        "plant_taxon_source_row_units_with_positive_eggs":dict(plant_positive_egg_units),
        "plant_taxon_source_row_units_with_positive_L3_L5":dict(plant_positive_late_units),
        "positive_egg_transect_round_plant_units":len(positive_egg_groups),
        "positive_late_transect_round_plant_units":len(positive_late_groups),
        "original_observed_aggregate_totals":{
            p:{k:observed_total[(p,k)] for k in ("Eggs","Early_L1_L2","Late_L3_L5","Available_Stems")}
            for p in PLANTS},
        "source_variable_missing_and_invalid_aggregates":dict(sorted(missing_fields.items())),
        "stages_valid_vs_summary_totals":dict(valid_stage_compare),
        "stage_summary_total_mismatch_counts":dict(conflicting_stage_compare),
        "pheno_stem_source_contradictions":dict(inconsistent_stem_total),
        "source_semantic_failures":dict(invalid_row),
        "core_data_columns_missing_positive_aggregate":{
            k:v.get("missing",0) for k,v in stats.items() if v.get("missing",0)>0},
        "per_method_distinct_sampling_units":{k:len(v) for k,v in method_units.items()},
        "raw_transect_site_geolocation_date_records_NOT_exported":True,
    }
    complete=[p for p,check in completeness.items() if all(check.values())]
    eligible_hosts=sum(plant_units.get(p,0)>=15 for p in complete)
    latemilkweed_units=len(positive_late_groups)
    egg_units=len(positive_egg_groups)
    totalcheck=valid_stage_compare["egg_total_checked"]+sum(valid_stage_compare["L"+k+"_checked"] for k in ("1","2","3","4","5"))
    totalmismatch=sum(conflicting_stage_compare.values())
    if (len(bytransect)>=15 and len(years)>=3 and eligible_hosts>=2 and
        latemilkweed_units>=30 and egg_units>=30 and eligible>=100 and
        totalmismatch<=max(1,totalcheck)*.01 and
        not invalid_row["animal_on_zero_available_host"]):
        result["status"]="SOURCE_STAGE_AND_HOST_PHENOLOGY_FEASIBLE_NEEDS_GRAIN_QC"
    else:
        result["status"]="SOURCE_STAGE_OR_HOST_QUALITY_GATE_INSUFFICIENT"
    result["gate_counts"]={
        "years":len(years),"plant_species_n_with_gte15_source_stem_positive":eligible_hosts,
        "positive_late_plant_round_units":latemilkweed_units,
        "positive_egg_plant_round_units":egg_units,
        "source_total_checked":totalcheck,"source_total_mismatch":totalmismatch,
        "needs_method_deduplication":len(origin_keys)<n_rows
    }
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result={
        "schema":"chocho_western_monarch_stage_and_milkweed_source_gate_20261026_v01",
        "scientific_question":"availability and host phenophase versus stage-specific actual juvenile occupancy, NOT cohort survival",
        "doi":DOI,
        "original_published_date":"2026-10-08",
        "original_host_and_juvenile_stages_species":"Danaus plexippus",
        "raw_survey_transect_coordinates_or_plant_records_not_exported":True,
    }
    result["metadata"]=digest_metadata()
    try:
        raw,receipt=http_binary(ZIP_URL,MAX_ZIP)
        result["download"]=receipt
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            bad=z.testzip()
            if bad:
                raise ValueError("Dryad ZIP CRC failed")
            members=[x for x in z.infolist() if not x.is_dir()]
            result["zip_members"]=[
                {"name":Path(x.filename).name,"uncompressed_bytes":x.file_size,
                 "source_file_sha256":hashlib.sha256(z.read(x)).hexdigest()}
                for x in members]
            choices=[x for x in members if Path(x.filename).name==MAIN]
            if len(choices)!=1:
                raise ValueError(f"Expected exactly one original {MAIN}, found {len(choices)}")
            original=z.read(choices[0])
            if len(original)>6_000_000:
                raise ValueError("Original main CSV unexpectedly changed size")
            result["csv"]=inspect_csv(original)
            result["status"]=result["csv"]["status"]
    except (HTTPError,URLError,TimeoutError,OSError,ValueError,zipfile.BadZipFile,UnicodeError) as e:
        result["status"]="ORIGINAL_SOURCE_ACQUISITION_OR_QC_FAILED"
        result["error_class"]=type(e).__name__
        result["http_status"]=getattr(e,"code",None)
        result["error"]=(str(e)[:200] if not isinstance(e,HTTPError) else "")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_status":result["status"],
        "transport":result.get("download",{}).get("http_status"),
        "metadata":result["metadata"].get("status"),
        "sample":{key:result.get("csv",{}).get(key) for key in (
            "original_csv_rows","distinct_transects","transect_round_sampling_units",
            "species_x_transect_round_units_with_available_stems",
            "transect_round_units_at_least_two_plant_species_available",
            "years_row_counts","methods_row_counts",
            "positive_egg_transect_round_plant_units",
            "positive_late_transect_round_plant_units","stage_summary_total_mismatch_counts",
            "pheno_stem_source_contradictions","gate_counts")}
    },ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
