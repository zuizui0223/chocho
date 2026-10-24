#!/usr/bin/env python3
"""Privacy-safe indexed-host-annotation availability, 39 fixed GBIF immature-stage labels.

Only aggregates leave this process. No occurrence bodies, observers, coordinates,
dates, comments, IDs, photos or links are persisted. Text matching is evidence
of *metadata wording*, never validated host use, offspring or recruitment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API="https://api.gbif.org/v1/occurrence/search"
BASIS=("HUMAN_OBSERVATION","PRESERVED_SPECIMEN","OBSERVATION","MACHINE_OBSERVATION","OCCURRENCE")
FIELDS=("associatedTaxa","associatedOccurrences","occurrenceRemarks",
        "organismRemarks","habitat","dynamicProperties")
KEYWORDS=re.compile(r"\b(?:feed(?:ing|s)?|fed on|oviposit(?:ion|ing|ed|s)?|laid eggs?|"
                    r"host\s+plant|feeding\s+on|larval\s+host|reared\s+on)\b",re.I)
SCHEMA="gbif_original39_juvenile_host_annotation_coverage_20261023_v1"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def utc():
    return datetime.now(timezone.utc).isoformat()

def source_cohort(juvenile_path,positive48_path,five_path):
    juvenile=json.loads(Path(juvenile_path).read_bytes())
    positive=json.loads(Path(positive48_path).read_bytes())
    five=json.loads(Path(five_path).read_bytes())
    if juvenile.get("status")!="COMPLETE_CONDITIONAL_JUVENILE_GADM1" or len(juvenile.get("attempts",[]))!=15:
        raise ValueError("Previous fifteen conditional juvenile sources incomplete")
    if juvenile.get("gadm1_mapped_juvenile_records")!=39 or juvenile.get("resolved_juvenile_records")!=39:
        raise ValueError("Original juvenile stage totals drifted")
    if positive.get("status")!="COMPLETE_48" or len(positive.get("observations",[]))!=48:
        raise ValueError("Original 48 butterfly-country source changed")
    if (five.get("focal_stage_case_count"),five.get("focal_candidate_plant_links"),
        five.get("focal_unique_candidate_plant_ids"))!=(5,20,8):
        raise ValueError("Original WCVP botanical candidate identity changed")
    old={(r["species"],r["country"]):r for r in positive["observations"]}
    plantlist={}
    for row in five["focal_five_cases"]:
        key=(row["species"],row["country"])
        if key in plantlist or not row["selected_five_before_botanical_outcomes"]:
            raise ValueError("Original five plant-case identities changed")
        species=[]
        for plant in row["candidate_plants"]:
            if plant["taxonomic_name_status"]!="EXACT_SOURCE_BINOMIAL_GBIF_RESOLUTION_REQUIRED":
                raise ValueError("Original plant name no longer accepted")
            if len(plant["botanical_names"])!=1:
                raise ValueError("Ambiguous original plant name")
            species.append(plant["botanical_names"][0])
        plantlist[key]=sorted(set(species))
    result=[]
    seen=set()
    for q in juvenile["attempts"]:
        if q["status"]!="COMPLETE_MATCHED_STAGE" or q["expected_stage_count"]!=q["record_count"]:
            raise ValueError("One of 15 stage facets mismatched original source")
        key=(q["species"],q["country"],q["life_stage"])
        if key in seen or (q["species"],q["country"]) not in old:
            raise ValueError("Repeated or unknown species/country/lifeStage")
        seen.add(key)
        source=old[q["species"],q["country"]]
        if source["taxon_key"]!=q["taxon_key"] or source["facets"]["life_stage"]["categories"].get(q["life_stage"])!=q["expected_stage_count"]:
            raise ValueError("Source-locked species/country/stage counts or taxon key changed")
        result.append({"species":q["species"],"country":q["country"],
                       "taxon_key":q["taxon_key"],"stage":q["life_stage"],
                       "expected_count":q["expected_stage_count"],
                       "candidate_names":plantlist.get((q["species"],q["country"]),[])})
    if (len(result),sum(q["expected_count"] for q in result),
        len({(x["species"],x["country"]) for x in result}))!=(15,39,11):
        raise ValueError("Original 15 queries/39 juvenile labels/11 species-country cells changed")
    return sorted(result,key=lambda x:(x["species"],x["country"],x["stage"])), {
        "juvenile_source_sha256":sha(Path(juvenile_path).read_bytes()),
        "all48_source_sha256":sha(Path(positive48_path).read_bytes()),
        "botanical_five_source_sha256":sha(Path(five_path).read_bytes()),
        "frozen_source_queries":15,"frozen_juvenile_label_records":39,"source_juvenile_countries":11,
        "source_plant_listed_country_cases":5,
    }

def query_url(q):
    params=[
        ("taxonKey",str(q["taxon_key"])),("country",q["country"]),
        ("year","2010,2025"),("lifeStage",q["stage"]),
        ("hasCoordinate","true"),("hasGeospatialIssue","false"),
        ("occurrenceStatus","PRESENT"),("limit","100"),("offset","0")
    ]
    params.extend(("basisOfRecord",x) for x in BASIS)
    return API+"?"+urlencode(params)

def normalized_field(value):
    if value is None:return ""
    if isinstance(value,(list,dict)):
        return json.dumps(value,ensure_ascii=False,sort_keys=True)
    return str(value).strip()

def inspect_response(q,raw):
    """Return aggregate flags ONLY. The original records never leave this call."""
    try:
        doc=json.loads(raw)
    except (ValueError,UnicodeDecodeError) as exc:
        raise ValueError("Truncated or malformed source result") from exc
    if not isinstance(doc,dict):
        raise ValueError("GBIF source must be JSON object")
    n=doc.get("count")
    records=doc.get("results")
    if (type(n) is not int or n!=q["expected_count"] or
        not isinstance(records,list) or len(records)!=n or n>=100 or
        doc.get("offset")!=0 or doc.get("endOfRecords") is not True):
        raise ValueError("Indexed source drift or incomplete result page; not a biological zero")
    seen=set()
    coverage=Counter()
    keycounts=Counter()
    names=Counter()
    with_any=mention=mention_with_host_words=0
    missing_id=0
    for x in records:
        if not isinstance(x,dict):raise ValueError("Invalid GBIF occurrence result row")
        key=x.get("key")
        if type(key) is not int:
            missing_id+=1
        elif key in seen:raise ValueError("Duplicate GBIF keys in a stage-query response")
        else:seen.add(key)
        if x.get("countryCode") not in (None,q["country"]):
            raise ValueError("GBIF record country mismatch")
        # DO NOT persist source text. Inspect a fixed set of keys in memory.
        raw_text=[]
        for name in FIELDS:
            text=normalized_field(x.get(name))
            if text:
                coverage[name]+=1
                raw_text.append(text)
        joined="\n".join(raw_text)
        if raw_text:with_any+=1
        local=[]
        for botanical in q["candidate_names"]:
            if re.search(r"(?<![A-Za-z])"+re.escape(botanical)+r"(?![A-Za-z])",joined,re.I):
                names[botanical]+=1
                local.append(botanical)
        if local:
            mention+=1
            if KEYWORDS.search(joined):mention_with_host_words+=1
        if KEYWORDS.search(joined):keycounts["any_host_related_keyword_in_permitted_fields"]+=1
    return {
        "status":"COMPLETE_NO_RECORD_BODY_PERSISTED","frozen_stage_count":n,
        "unique_GBIF_record_keys_checked_not_saved":len(seen),
        "missing_GBIF_record_key_count":missing_id,
        "records_with_any_selected_text_field":with_any,
        "field_nonempty_record_counts":dict(sorted(coverage.items())),
        "records_with_host_related_words_unverified":keycounts["any_host_related_keyword_in_permitted_fields"],
        "records_with_exact_source_candidate_binomial_mention":mention,
        "records_with_candidate_mention_and_host_related_words_unverified":mention_with_host_words,
        "candidate_plant_name_mention_counts":dict(sorted(names.items())),
        "no_raw_observer_location_date_or_text_saved":True
    }

def fetch(q,retries=3):
    url=query_url(q)
    errors=[]
    for attempt in range(retries):
        try:
            req=Request(url,headers={"User-Agent":"chocho-public-GBIF-juvenile-host-text-coverage/0.1",
                                      "Accept":"application/json"})
            with urlopen(req,timeout=40) as response:
                status=response.status
                raw=response.read(3_000_000)
            if status!=200 or len(raw)>=3_000_000:
                raise ValueError("GBIF response body incomplete")
            stats=inspect_response(q,raw)
            return {**{k:v for k,v in q.items() if k!="candidate_names"},
                    "source_candidate_binomial_count":len(q["candidate_names"]),
                    **stats,"query_url":url,"retrieved_utc":utc(),"http_status":status,
                    "raw_response_sha256":sha(raw),"attempts":attempt+1}
        except (HTTPError,URLError,TimeoutError,OSError,ValueError) as exc:
            errors.append({"attempt":attempt+1,"exception_type":type(exc).__name__,"error_message":str(exc)[:150]})
            if attempt<retries-1:time.sleep(2*(attempt+1))
    return {**{k:v for k,v in q.items() if k!="candidate_names"},
            "source_candidate_binomial_count":len(q["candidate_names"]),
            "status":"SOURCE_UNAVAILABLE_OR_INDEX_DRIFT_NOT_ZERO","query_url":url,
            "retrieved_utc":utc(),"attempts":retries,"errors":errors}

def summarise(qs,rows,meta):
    if len(qs)!=15 or len(rows)!=15:return {"status":"INCOMPLETE_LOCKED_QUERY_SET"}
    counters=Counter();by_case=defaultdict(lambda:{"source_stage_records":0,"record_text_annotated":0,
                   "any_host_keyword_records":0,"matched_source_binomial_records":0})
    for row in rows:
        if row["status"]!="COMPLETE_NO_RECORD_BODY_PERSISTED":continue
        n=row["frozen_stage_count"]
        counters["complete_queries"]+=1
        for key in ("records_with_any_selected_text_field","records_with_host_related_words_unverified",
                    "records_with_exact_source_candidate_binomial_mention",
                    "records_with_candidate_mention_and_host_related_words_unverified"):
            counters[key]+=row[key]
        counters["original_juvenile_stage_records_scanned"]+=n
        key=row["species"]+"|"+row["country"]
        g=by_case[key]
        g["source_stage_records"]+=n
        g["record_text_annotated"]+=row["records_with_any_selected_text_field"]
        g["any_host_keyword_records"]+=row["records_with_host_related_words_unverified"]
        g["matched_source_binomial_records"]+=row["records_with_exact_source_candidate_binomial_mention"]
    complete=counters["complete_queries"]==15 and counters["original_juvenile_stage_records_scanned"]==39
    return {"schema":SCHEMA,"status":"COMPLETE_39_METADATA_COVERAGE_NOT_HOST_USE" if complete
            else "INCOMPLETE_SOURCE_METADATA_NOT_EVIDENCE_OF_ABSENCE",
            "executed_utc":utc(),"source":meta,
            "fixed_stage_queries":15,"completed":counters["complete_queries"],
            "original_juvenile_stage_records":39,
            "observed_source_stage_records":counters["original_juvenile_stage_records_scanned"],
            "records_with_nonempty_selected_text_fields":counters["records_with_any_selected_text_field"],
            "records_with_any_host_action_keyword_not_verified":counters["records_with_host_related_words_unverified"],
            "records_with_exact_accepted_candidate_binomial_in_text":counters["records_with_exact_source_candidate_binomial_mention"],
            "records_with_candidate_plant_name_and_host_word_not_verified":counters["records_with_candidate_mention_and_host_related_words_unverified"],
            "species_country_stage_aggregates":rows,
            "species_country_aggregates":dict(sorted(by_case.items())),
            "privacy":"No individual GBIF records or names of observers, record identifiers, source text, exact dates or coordinates are stored.",
            "interpretation":"Metadata availability only. A string mention is not verified larval host use, survival, breeding, introduction date or recruitment."}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--juvenile",type=Path,required=True)
    p.add_argument("--positive48",type=Path,required=True)
    p.add_argument("--host-five",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    cohort,meta=source_cohort(a.juvenile,a.positive48,a.host_five)
    rows=[]
    for q in cohort:
        rows.append(fetch(q))
        time.sleep(1)
    result=summarise(cohort,rows,meta)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GBIF_JUVENILE_HOST_METADATA="+json.dumps(result,separators=(",",":")),flush=True)
    if result["status"]!="COMPLETE_39_METADATA_COVERAGE_NOT_HOST_USE":
        raise SystemExit("Source drift/partial, do not infer negative biological host evidence")

if __name__=="__main__":
    main()
