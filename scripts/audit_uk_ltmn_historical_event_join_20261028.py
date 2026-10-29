#!/usr/bin/env python3
"""Scientific source-complete LTMN visit vs butterfly count identity audit.

Only source-level aggregates are serialized. No site IDs, locations, observer
identities, visit-date combinations, or original source rows leave the process.
Zero records count ONLY explicit numeric zero for a named species on a verified
visit: absence of a species-count row is unknown, NOT observed non-detection.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,re
from pathlib import Path
from audit_uk_ltmn_xls_schema_20261028 import SOURCES,get_source

SCHEMA="chocho_original_LTMN_ukbms_visited_species_counts_source_gate_20261028_v01"
FIELDS_VISITS=("UKBMS_SITENO","DAY","MONTH","YEAR","START_TIME")
FIELDS_COUNTS=FIELDS_VISITS+("WEEKNO","COMMON_NAME","VISIT_TOTAL")
WINDOWS={
    "early_1976_1995":(1976,1995),
    "late_2005_2022":(2005,2022),
    "early_1993_2003":(1993,2003),
    "late_2007_2015":(2007,2015),
    "mid_1996_2006":(1996,2006),
    "recent_2016_2023":(2016,2023),
}

def key_site(v):
    if v is None or v=="":return None
    if isinstance(v,(int,float)) and float(v).is_integer() and v>0:
        return str(int(v))
    s=str(v).strip()
    return s if s and len(s)<=80 else None

def integral(v,lo,hi):
    if v is None or v=="":return None
    try:
        x=float(v)
        if math.isfinite(x) and x.is_integer() and lo<=x<=hi:
            return int(x)
    except (TypeError,ValueError):pass
    return None

def time_key(t):
    """1-minute precision; preserve raw times for nonstandard strings."""
    if t is None or t=="":return "UNKNOWN"
    if isinstance(t,(int,float)):
        v=float(t)
        if 0<=v<1:return str(int(round(v*1440))%1440)
        if v>=1 and v<=24 and v.is_integer():return str(int(v)*60)
        return f"NUM-{round(v,6)}"
    s=str(t).strip()
    if re.fullmatch(r"\d{1,2}:\d{2}(?::\d{2})?",s):
        h,mi=map(int,s.split(":")[:2])
        if h<24 and mi<60:return str(h*60+mi)
    return s[:35]

def excel_rows(data,required):
    import xlrd
    wb=xlrd.open_workbook(file_contents=data,on_demand=True)
    try:
        for sh in wb.sheets():
            if sh.nrows<2:continue
            headers=[str(x).strip().upper() for x in sh.row_values(0)]
            if len(headers)!=len(set(headers)) or not set(required)<=set(headers):
                raise ValueError(f"Source changed column set on sheet {sh.name}")
            idx={key:headers.index(key) for key in required}
            for i in range(1,sh.nrows):
                vals=sh.row_values(i)
                yield {k:vals[j] for k,j in idx.items()}
    finally:
        wb.release_resources()

def parse_event(row):
    site=key_site(row["UKBMS_SITENO"])
    day=integral(row["DAY"],1,31)
    month=integral(row["MONTH"],1,12)
    year=integral(row["YEAR"],1970,2030)
    if not site or None in (day,month,year):return None
    from datetime import date
    try:date(year,month,day)
    except ValueError:return None
    return (site,year,month,day,time_key(row["START_TIME"]))

def joined_source(visits,counts):
    visit_events=set(); date_visits=collections.defaultdict(set); visit_duplicate=0
    invalid_visits=0; year_visit=collections.Counter(); names=collections.Counter()
    for r in visits:
        ev=parse_event(r)
        if ev is None:
            invalid_visits+=1;continue
        if ev in visit_events:
            visit_duplicate+=1;continue
        visit_events.add(ev)
        date_visits[ev[:4]].add(ev)
        year_visit[ev[1]]+=1
    counted_species={};source_duplicate=0;conflicts=0;bad_count=0;unmatched=0
    aligned=0; matched_date_only=0; year_counts=collections.Counter()
    actual_names=collections.Counter(); zero_species=collections.Counter()
    for r in counts:
        ev=parse_event(r)
        species=str(r["COMMON_NAME"] or "").strip()
        try:num=float(r["VISIT_TOTAL"])
        except (ValueError,TypeError):num=float("nan")
        if ev is None or not species or not math.isfinite(num) or num<0:
            bad_count+=1;continue
        key=(ev,species)
        if key in counted_species:
            if num!=counted_species[key]:conflicts+=1
            else:source_duplicate+=1
            # never sum or turn duplicates into biological abundance
            continue
        counted_species[key]=num
        year_counts[ev[1]]+=1
        actual_names[species]+=1
        if num==0:zero_species[species]+=1
        if ev in visit_events:aligned+=1
        elif ev[:4] in date_visits:matched_date_only+=1
        else:unmatched+=1
    valid_counted_events=set(k[0] for k in counted_species)
    count_event_date=set(k[:4] for k in valid_counted_events)
    visit_sites={k[0] for k in visit_events}
    count_sites={k[0] for k in valid_counted_events}
    by_window={}
    for key,(start,end) in WINDOWS.items():
        e={event for event in visit_events if start<=event[1]<=end}
        c={event for event in valid_counted_events if start<=event[1]<=end}
        mapped=sum(1 for event in e if event in c)
        by_window[key]={
            "original_visit_events":len(e),
            "original_count_events":len(c),
            "shared_exact_visit_count_events":mapped,
            "n_original_visit_sites":len({k[0] for k in e}),
            "n_original_count_sites":len({k[0] for k in c}),
            "zero_count_named_species_rows":sum(1 for (event,sp),num in counted_species.items()
                                                if start<=event[1]<=end and num==0),
            "positive_named_species_rows":sum(1 for (event,sp),num in counted_species.items()
                                                if start<=event[1]<=end and num>0),
        }
    early={ev[0] for ev in visit_events if 1976<=ev[1]<=1995}
    late={ev[0] for ev in visit_events if 2005<=ev[1]<=2022}
    e2={ev[0] for ev in visit_events if 1993<=ev[1]<=2003}
    l2={ev[0] for ev in visit_events if 2007<=ev[1]<=2015}
    zero_total=sum(zero_species.values())
    all_total=len(counted_species)
    event_groups=collections.defaultdict(list)
    for (ev,sp),num in counted_species.items():event_groups[ev].append((sp,num))
    count_events_nspecies=[len(v) for v in event_groups.values()]
    from statistics import median
    return {
       "schema":SCHEMA,
       "source_context":{
         "original_visit_xls_sha256":SOURCES["visits"][1],
         "original_weekly_count_xls_sha256":SOURCES["counts"][1],
         "original_ukbms_transect_records_under_ogl":True,
         "source_dates_and_transect_ids_serialized":False,
         "recorder_personal_identifiers_exported":False,
       },
       "source_visit_qc":{
         "unique_original_visit_events":len(visit_events),
         "duplicate_original_visit_rows":visit_duplicate,
         "invalid_original_visit_rows":invalid_visits,
         "n_unique_visit_transect_identifiers":len(visit_sites),
         "earliest_visit_year":min(year_visit) if year_visit else None,
         "latest_visit_year":max(year_visit) if year_visit else None,
         "original_visits_per_year":dict(sorted(year_visit.items())),
       },
       "source_counts_qc":{
         "unique_named_species_count_records":all_total,
         "unique_original_count_events":len(valid_counted_events),
         "n_unique_count_transect_identifiers":len(count_sites),
         "nonidentical_duplicate_count_source_rows":conflicts,
         "duplicate_identical_count_rows":source_duplicate,
         "invalid_or_unscorable_count_rows":bad_count,
         "exact_matched_species_count_records_to_actual_visits":aligned,
         "source_count_rows_at_matching_visit_date_but_not_starttime":matched_date_only,
         "source_count_rows_without_visit_on_same_calendar_day":unmatched,
         "number_named_species":len(actual_names),
         "explicit_species_count_zero_records":zero_total,
         "explicit_species_count_positive_records":all_total-zero_total,
         "fraction_explicit_zero_among_named_count_rows":zero_total/all_total if all_total else None,
         "median_number_named_species_record_rows_per_count_event":median(count_events_nspecies) if count_events_nspecies else None,
         "year_counts":dict(sorted(year_counts.items())),
         "nonzero_common_name_rows_by_original_taxon":dict(sorted(actual_names.items())),
         "explicit_zero_rows_by_named_taxon":dict(sorted(zero_species.items())),
       },
       "source_periods":by_window,
       "longterm_repeated_transects":{
         "n_with_any_visit_1976_1995_and_2005_2022":len(early&late),
         "n_with_any_visit_1993_2003_and_2007_2015":len(e2&l2),
       },
       "inference_firewall":{
         "missing_named_species_rows_are_not_nondetections":True,
         "explicit_zero_rows_if_present_are_only_named_taxon_visit_zero":True,
         "no_population_or_host_mediation_inferred":True,
         "survey_event_can_be_matched_by_site_and_date":bool(aligned+matched_date_only),
         "full_source_protocol_includes_targeted_surveys_unknown":True,
       },
    }

def run(path):
    v=get_source("visits")
    c=get_source("counts")
    out=joined_source(excel_rows(v,FIELDS_VISITS),excel_rows(c,FIELDS_COUNTS))
    out["status"]="ORIGINAL_UK_TRANSECT_VISIT_COUNT_DATA_AUDITED"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    q=out["source_counts_qc"]
    print(json.dumps({
       "status":out["status"],
       "visits":out["source_visit_qc"]["unique_original_visit_events"],
       "count_events":q["unique_original_count_events"],
       "count_rows":q["unique_named_species_count_records"],
       "named_species":q["number_named_species"],
       "explicit_species_zeros":q["explicit_species_count_zero_records"],
       "matched_exact":q["exact_matched_species_count_records_to_actual_visits"],
       "matched_date_only":q["source_count_rows_at_matching_visit_date_but_not_starttime"],
       "no_visit_same_date":q["source_count_rows_without_visit_on_same_calendar_day"],
       "periods":out["source_periods"],
       "repeated":out["longterm_repeated_transects"]},sort_keys=True),flush=True)
    return out

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--output",required=True,type=Path)
    run(a.parse_args().output)
