#!/usr/bin/env python3
"""Taxon-ID / dataset-type / per-record age audits for *candidate* 18–24ka fossils.

The generic Neotoma taxonname filter was proven invalid; numeric taxonid and
age parameters have been checked as potentially functional but NO historical
species occurrence is confirmed without independently authenticated chronology,
site geomorphology, pollen/macro identification and original deposition.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter
from pathlib import Path
import requests
from probe_neotoma_taxonid_sample_age_v12 import source_http

V11="55887306ce8defa7c027fc1d518f2bb798ac2868b6c7f46db288f154e8ea8219"
V12="657d56b7ae45176de48a5f86eab7d32d783a4eaff70cd7a07fe334e57af3c690"
API="https://api.neotomadb.org/v2.0/data/occurrences"
CANDIDATE="SOURCE_18_24KA_METADATA_CANDIDATE_NOT_VERIFIED_FOSSIL_OCCUPANCY"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source_identity(contract,v11,v12):
    if sha(v11)!=V11 or sha(v12)!=V12:raise ValueError("source v11/v12 taxon ID or age metadata hashes drift")
    a=json.loads(Path(v11).read_text())
    b=json.loads(Path(v12).read_text())
    q=json.loads(Path(contract).read_text())
    if (a["status"]!="FOSSIL_TAXONNAME_FILTER_INVALID_SENTINEL_MATCHES_REAL_TAXON"
        or b["status"]!="TAXON_ID_PATH_DIFFERENT_OCCURRENCE_RECORDS_CHRONOLOGY_UNVERIFIED"
        or q["status"]!="BOUNDED_SOURCE_SAMPLE_VERIFICATION_NOT_LGM_OCCUPANCY"):
        raise ValueError("missing source verified taxonomy/filter boundaries")
    chosen={x["name"]:x["neotoma_taxonid"] for x in q["cases"]}
    if chosen!=a["exact_taxon_registry_ids"] or chosen!=b["source_only_target_taxon_ids"]:
        raise ValueError("source Neotoma numeric accepted taxon ID mutated")
    cfg=q["queries"]
    if (cfg["limit"],cfg["offset"],cfg["ageyoung"],cfg["ageold"])!=(30,0,18000,24000):
        raise ValueError("source frozen fossil-query protocol changed")
    return q

def check_sample(record,taxonid,name,datasettype):
    """Never infer glacial occupancy from a missing time/model/taxon source field."""
    src=record.get("sample")
    site=record.get("site")
    age=record.get("age")
    if not isinstance(src,dict) or not isinstance(site,dict) or not isinstance(age,dict):
        return "SOURCE_METADATA_UNVERIFIABLE"
    tid=src.get("taxonid")
    scientific=src.get("taxonname")
    source_proxy=site.get("datasettype")
    if str(tid)!=str(taxonid):return "WRONG_TAXON_ID_FILTER_FAILURE"
    if scientific!=name:return "NOT_EXACT_SPECIES_IN_SOURCE"
    if source_proxy!=datasettype:return "WRONG_DATASETTYPE_FILTER_FAILURE"
    n=age.get("age")
    if not isinstance(n,(int,float)) or not math.isfinite(float(n)):
        return "UNKNOWN_SAMPLE_AGE_DO_NOT_INFER_ABSENCE"
    if not 18000<=float(n)<=24000:return "SAMPLE_CENTER_AGE_NOT_18_24KA"
    return CANDIDATE

def record_source_metadata(record,taxonid,name,proxy):
    state=check_sample(record,taxonid,name,proxy)
    site=record.get("site") or {}
    sample=record.get("sample") or {}
    age=record.get("age") or {}
    def value(key):return age.get(key)
    return {"requested_taxon":name,"requested_Neotoma_taxonid":taxonid,
        "requested_proxy":proxy,
        "occid":str(record.get("occid") or ""),
        "source_sample_taxonid":sample.get("taxonid"),
        "source_sample_taxonname":sample.get("taxonname"),
        "source_datasettype":site.get("datasettype"),
        "source_siteid":site.get("siteid"),
        "source_datasetid":site.get("datasetid"),
        "source_site_location_raw":str(site.get("location") or "")[:250],
        "sample_age_BP":value("age"),"sample_age_younger_BP":value("ageyounger"),
        "sample_age_older_BP":value("ageolder"),
        "source_grade":state,
        "age_model_calibration_independently_verified":False,
        "pollen_species_level_identification_independently_verified":False}

def run(session,spec):
    req=[];occurrences=[]
    cfg=spec["queries"]
    for x in spec["cases"]:
        for proxy in cfg["datasettypes"]:
            params={"taxonid":x["neotoma_taxonid"],"datasettype":proxy,
                "ageyoung":cfg["ageyoung"],"ageold":cfg["ageold"],
                "limit":cfg["limit"],"offset":cfg["offset"]}
            data,meta=source_http(session,API,params)
            req.append({"query_type":"bounded_positive_source","name":x["name"],"taxonid":x["neotoma_taxonid"],
                        "datasettype":proxy,**meta})
            for r in data:
                occurrences.append(record_source_metadata(r,x["neotoma_taxonid"],x["name"],proxy))
    for x in cfg["negative_controls"]:
        params=dict(x,limit=2,offset=0)
        records,meta=source_http(session,API,params)
        req.append({"query_type":"taxonid_or_age_negative_source_control",**params,**meta})
    repeats=Counter(x["occid"] for x in occurrences)
    for x in occurrences:
        x["duplicated_occid_in_bounded_query"]=repeats[x["occid"]]>1
    grades=Counter(x["source_grade"] for x in occurrences)
    if grades["WRONG_TAXON_ID_FILTER_FAILURE"] or grades["WRONG_DATASETTYPE_FILTER_FAILURE"]:
        status="SOURCE_FILTER_INTEGRITY_FAILED_NO_FOSSIL_INFERENCE"
    elif grades[CANDIDATE]:
        status="AGE_WINDOW_METADATA_CANDIDATES_FOUND_BUT_CHRONOLOGY_UNVERIFIED"
    else:
        status="NO_QUALIFIED_AGE_WINDOW_IN_BOUNDED_SOURCE_NO_ABSENCE_INFERENCE"
    by_query=[]
    for r in req:
        k=r.get("name") or str(r.get("taxonid"))
        by_query.append({"source":k,"proxy":r.get("datasettype"),"kind":r["query_type"],
                         "source_status":r["status"],"response_records":r.get("record_count",0)})
    summary={"schema":"chocho_Neotoma_18_24ka_taxonID_sample_proxy_age_candidates_v13",
        "status":status,
        "frozen_taxa":[x["name"] for x in spec["cases"]],
        "sample_records_examined":len(occurrences),
        "distinct_occurrence_ids_examined":len(repeats),
        "source_sample_grade_counts":dict(grades),
        "source_candidate_distinct_occids":len({
           x["occid"] for x in occurrences if x["source_grade"]==CANDIDATE}),
        "response_query_overview":by_query,
        "source_sample_age_calibration_verified":False,
        "species_level_pollen_diagnostic_verified":False,
        "global_21ka_host_plant_distribution_estimated":False,
        "independent_insect_genetic_response_used":False,
        "source_query_not_exhaustive":True,
        "original_negative_host_specificity_result_preserved":True,
        "caveats":[
          "Taxon-ID, source datasettype and sample age numeric fields are source metadata, not independently verified fossil identification or age-depth calibration",
          "Only two temperate WCVP model plants and at most 30 hits per proxy queried; absence cannot be inferred",
          "Pollen species label can be taxonomic assignment from database, not diagnostic pollen morphology",
          "Multiple ages per occurrence and transported pollen can distort true source 21 ka plant occupancy",
          "No new 21ka host range or insect genetic break is inferred; previous TTF and chocho host-identity negative findings remain unchanged"
        ]}
    return req,occurrences,summary

def write_rows(path,rows):
    if not rows:
        Path(path).write_text("status\nNO_SOURCE_ROWS_DO_NOT_INFER_SPECIES_ABSENCE\n")
        return
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with Path(path).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,list) else v for k,v in r.items()})

def main():
    p=argparse.ArgumentParser()
    for field in ("contract","v11","v12","out"):
        p.add_argument("--"+field,required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    q=source_identity(a.contract,a.v11,a.v12)
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho/21ka-source-only-plant-macrofossil-age-identity-audit/0.1"})
    queries,rows,summary=run(session,q)
    write_rows(a.out/"bounded_taxonid_age_source_queries.csv",queries)
    write_rows(a.out/"bounded_18_24ka_sample_taxon_proxy_checks.csv",rows)
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()
