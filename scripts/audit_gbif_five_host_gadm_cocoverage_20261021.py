#!/usr/bin/env python3
"""Source-frozen WCVP accepted host candidates × GBIF botanical GADM1 observations.

ONLY five pre-selected butterfly-country cells with source immature stage
labels and added host-resource opportunity. Never interprets GADM-level
co-detection as larval feeding, establishment or successful recruitment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

TAXON_API="https://api.gbif.org/v1/species/match"
OCC_API="https://api.gbif.org/v1/occurrence/search"
BASIS=("HUMAN_OBSERVATION","PRESERVED_SPECIMEN","OBSERVATION",
       "MACHINE_OBSERVATION","OCCURRENCE")
SCHEMA="chocho_strict_plant_taxon_gadm_cocoverage_five_juvenile_country_cases_20261021_v1"
SOURCE_SCHEMA="chocho_host_plant_candidates_5_juvenile_added_country_cases_20261021_v1"
JUV_SCHEMA="chocho_juvenile_conditioned_GADM1_source_facet_20261021_v1"
ORIGINAL_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"


def digest(b):
    return hashlib.sha256(b).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def load_fixed(source_path,juvenile_path,all_stage_path):
    sourcebytes=Path(source_path).read_bytes()
    src=json.loads(sourcebytes)
    juvenile=json.loads(Path(juvenile_path).read_bytes())
    all_stage=json.loads(Path(all_stage_path).read_bytes())
    if src.get("schema")!=SOURCE_SCHEMA or src.get("status")!="SOURCE_VERIFIED_CANDIDATE_PLANTS_NOT_LOCAL_HOST_USE":
        raise ValueError("Not the original WCVP source-verified host candidate list")
    if src["source_verification"].get("author_inner_zip_sha256")!=ORIGINAL_SHA:
        raise ValueError("Original author-only source checksum mismatched")
    if (src.get("country_positive_cells"),src.get("source_added_flag_cells"),
        src.get("focal_stage_case_count"),src.get("focal_unique_candidate_plant_ids"),
        src.get("focal_candidate_plant_links"))!=(48,26,5,8,20):
        raise ValueError("Source candidate host cohort changed")
    if juvenile.get("schema")!=JUV_SCHEMA or juvenile.get("status")!="COMPLETE_CONDITIONAL_JUVENILE_GADM1":
        raise ValueError("Juvenile-labeled GADM source incomplete")
    if all_stage.get("status")!="COMPLETE_48":
        raise ValueError("Original 48 source-positive butterfly country cells incomplete")
    juv={(r["species"],r["country"]):r for r in juvenile["cell_gadm1_juvenile_summary"]}
    all48={(r["species"],r["country"]):r for r in all_stage["observations"]}
    sourcecases={}
    host_names={}
    pairs=[]
    for cell in src["focal_five_cases"]:
        key=(cell["species"],cell["country"])
        if key in sourcecases or key not in juv or key not in all48:
            raise ValueError("Source case missing juvenile or full butterfly positive country data")
        sourcecases[key]=cell
        if not cell["selected_five_before_botanical_outcomes"] or cell["added_host_country"]!=1:
            raise ValueError("Case not source-selected before plant outcome lookup")
        if cell["candidate_introduced_host_count"]!=len(cell["candidate_plants"]):
            raise ValueError("Source host candidate count mismatch")
        for item in cell["candidate_plants"]:
            if item["taxonomic_name_status"]!="EXACT_SOURCE_BINOMIAL_GBIF_RESOLUTION_REQUIRED":
                raise ValueError("Source botanical candidate not an exact accepted binomial")
            if len(item["botanical_names"])!=1 or not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+",item["botanical_names"][0]):
                raise ValueError("Original host candidate has invalid or ambiguous scientific name")
            name=item["botanical_names"][0]
            if name in host_names and host_names[name]!=item["wcvp_accepted_id"]:
                raise ValueError("Same botanical name maps to conflicting WCVP IDs")
            host_names[name]=item["wcvp_accepted_id"]
            pairs.append({"species":cell["species"],"country":cell["country"],
                "plant_name":name, "wcvp_accepted_id":item["wcvp_accepted_id"],
                "wgsrpd_source_intro_units":item["country_material_introduced_wgsrpd3"]})
    if len(sourcecases)!=5 or len(pairs)!=20 or len(host_names)!=8:
        raise ValueError("Fixed five/20/eight host-source cohort drifted")
    return (sorted(sourcecases),host_names, sorted(pairs,key=lambda x:(x["species"],x["country"],x["plant_name"])),
            juv,all48,{"frozen_botanical_candidate_sha256":digest(sourcebytes),
                "source_focal_case_n":5, "wcvp_species_names":8,"country_host_links":20})


def get(url,tries=3):
    errors=[]
    for attempt in range(tries):
        try:
            req=Request(url,headers={"Accept":"application/json",
                                      "User-Agent":"chocho-host-observed-GADM1-source-locked/0.1"})
            with urlopen(req,timeout=50) as response:
                status=response.status
                content=response.read(1_000_000)
            if status!=200 or len(content)>=1_000_000:
                raise ValueError("GBIF HTTP response invalid")
            return {"status":"SOURCE_HTTP_OK","http_status":status,
                    "retrieved_utc":now(),"url":url,
                    "raw_sha256":digest(content),
                    "attempts":attempt+1,"json":json.loads(content)}
        except (HTTPError,URLError,ValueError,OSError,TimeoutError) as e:
            errors.append({"attempt":attempt+1,"type":type(e).__name__,"message":str(e)[:200]})
            if attempt<tries-1:time.sleep(2*(attempt+1))
    return {"status":"UNAVAILABLE_NOT_ZERO","retrieved_utc":now(),
            "url":url,"attempts":tries,"errors":errors}


def strict_match(name,raw):
    if raw["status"]!="SOURCE_HTTP_OK":
        return {**raw,"status":"TAXON_RESOLUTION_UNAVAILABLE","plant_name":name}
    obj=raw["json"]
    accepted=(obj.get("matchType")=="EXACT" and obj.get("rank")=="SPECIES"
              and obj.get("status")=="ACCEPTED" and obj.get("kingdom")=="Plantae"
              and obj.get("canonicalName")==name
              and type(obj.get("usageKey")) is int)
    return {"plant_name":name,"status":"EXACT_ACCEPTED_GBIF_PLANT" if accepted else "UNRESOLVED_OR_NOT_EXACT_ACCEPTED_SPECIES",
            "accepted_species_key":obj.get("usageKey") if accepted else None,
            "matchType":obj.get("matchType"),"canonicalName":obj.get("canonicalName"),
            "taxonomicStatus":obj.get("status"),"taxonRank":obj.get("rank"),
            "kingdom":obj.get("kingdom"),
            "url":raw["url"],"raw_sha256":raw["raw_sha256"],
            "retrieved_utc":raw["retrieved_utc"]}


def match_url(name):
    return TAXON_API+"?"+urlencode([("name",name),("kingdom","Plantae"),("rank","SPECIES"),("strict","true")])


def occurrence_url(pair,accepted_species_key):
    params=[("taxonKey",str(accepted_species_key)),
        ("country",pair["country"]),("year","2010,2025"),
        ("hasCoordinate","true"),("hasGeospatialIssue","false"),
        ("occurrenceStatus","PRESENT"),("limit","0"),
        ("facet","gadmLevel1Gid"),("facetLimit","300")]
    params += [("basisOfRecord",x) for x in BASIS]
    return OCC_API+"?"+urlencode(params)


def parse_plant_gadm(raw):
    if raw["status"]!="SOURCE_HTTP_OK":
        return {"status":"QUERY_UNAVAILABLE_NOT_PLANT_ABSENCE","url":raw["url"]}
    doc=raw["json"]
    n=doc.get("count")
    if type(n) is not int or n<0 or doc.get("limit")!=0 or doc.get("results")!=[]:
        raise ValueError("Invalid occurrence limit/total and indexed facet response")
    f=doc.get("facets")
    if (not f and n==0):
        f=[{"field":"GADM_LEVEL_1_GID","counts":[]}]
    if not isinstance(f,list) or len(f)!=1 or str(f[0].get("field")).upper().replace("_","")!="GADMLEVEL1GID":
        raise ValueError("GBIF plant GADM1 facet absent")
    rows=f[0].get("counts")
    if not isinstance(rows,list) or len(rows)>=300:
        raise ValueError("Truncated GBIF plant GADM1 facet")
    cats={}
    for row in rows:
        key,k=row.get("name"),row.get("count")
        if not isinstance(key,str) or "." not in key or type(k) is not int or k<0 or key in cats:
            raise ValueError("Malformed plant GADM1 count")
        cats[key]=k
    if sum(cats.values())>n:
        raise ValueError("Plant mapped GADM1 counts exceed GBIF total")
    return {"status":"COMPLETE_INDEXED_PLANT_FACET",
            "indexed_country_plant_records":n,
            "gadm1_count":len(cats),
            "gadm1_mapped_records":sum(cats.values()),
            "gadm1_unmapped_records":n-sum(cats.values()),
            "gadm1_counts":dict(sorted(cats.items())),
            "url":raw["url"],"raw_sha256":raw["raw_sha256"],
            "retrieved_utc":raw["retrieved_utc"]}


def build_report(cases,plant_names,pairs,juv,all48,meta,matching,plant_results):
    bycase=[]
    for sp,ct in cases:
        plantlist=[p for p in pairs if (p["species"],p["country"])==(sp,ct)]
        obs=[r for r in plant_results if (r["species"],r["country"])==(sp,ct)]
        if len(obs)!=len(plantlist):
            raise ValueError("Original botanical subset not queried for every fixed case")
        valid=[p for p in obs if p["status"]=="COMPLETE_INDEXED_PLANT_FACET"]
        pending=[p["plant_name"] for p in obs if p["status"]!="COMPLETE_INDEXED_PLANT_FACET"]
        union=Counter()
        for p in valid:
            for region,n in p["gadm1_counts"].items():union[region]+=n
        immature=juv[sp,ct]["juvenile_gadm1_counts"]
        all_butterfly=all48[sp,ct]["facets"]["gadm1"]["categories"]
        matched={g:n for g,n in immature.items() if union.get(g,0)>0}
        candidate_overlap_count=sum(matched.values())
        bycase.append({"species":sp,"country":ct,
            "source_accepted_host_candidates":len(plantlist),
            "exact_accepted_resolved_plants":sum(matching[p["plant_name"]]["status"]=="EXACT_ACCEPTED_GBIF_PLANT" for p in plantlist),
            "completed_plant_gadm_facets":len(valid),
            "plant_gadm_missing_query_or_taxonomy":pending,
            "indexed_plant_detected_gadm_units":sorted(union),
            "juvenile_labeled_focal_gadm1_counts":immature,
            "juvenile_records_all":sum(immature.values()),
            "juvenile_labeled_gadm_units":len(immature),
            "juvenile_labeled_records_in_plant_detected_gadm_units":candidate_overlap_count,
            "juvenile_labeled_gadm_units_with_plant_detection":len(matched),
            "all_focal_butterfly_gadm_count":len(all_butterfly),
            "all_focal_butterfly_record_count":sum(all_butterfly.values()),
            "candidate_plant_positive_gadm_among_all_focal_gadm_units":len(set(union)&set(all_butterfly)),
            "all_focal_records_in_candidate_plant_gadm":sum(n for g,n in all_butterfly.items() if g in union),
            "gadm1_cross_records_are_not_same_site_host_use":True,
            "metadata_status":"ALL_SOURCE_PLANTS_ASSESSED" if not pending else "INCOMPLETE_PLANT_TAXON_OR_FACET_NOT_ZERO"})
    full=all(x["metadata_status"]=="ALL_SOURCE_PLANTS_ASSESSED" for x in bycase)
    return {"schema":SCHEMA, "source":meta, "reported_utc":now(),
            "status":"ALL_FIVE_CASES_SOURCE_ASSESSED" if full else "PARTIAL_TAXON_OR_PLANT_RECORD_SOURCES_NOT_NEGATIVE",
            "taxon_resolved":sum(x["status"]=="EXACT_ACCEPTED_GBIF_PLANT" for x in matching.values()),
            "total_taxon_names":len(matching),
            "candidate_links":len(pairs),
            "plant_facet_responses":sum(x["status"]=="COMPLETE_INDEXED_PLANT_FACET" for x in plant_results),
            "exact_taxonomy":matching, "per_source_host_country":plant_results,
            "focal_case_results":bycase,
            "interpretation":"Matched source botanical taxa and same GADM1 *indexed records*, not actual local larval feeding, fitness or butterfly colonization."}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-five",type=Path,required=True)
    p.add_argument("--juvenile-gadm",type=Path,required=True)
    p.add_argument("--all-stage",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    cases,plants,pairs,juvenile,all48,meta=load_fixed(
        a.source_five,a.juvenile_gadm,a.all_stage)
    matching={}
    for name in sorted(plants):
        raw=get(match_url(name))
        matching[name]=strict_match(name,raw)
        time.sleep(.7)
    queries=[]
    for pair in pairs:
        info=matching[pair["plant_name"]]
        if info["status"]!="EXACT_ACCEPTED_GBIF_PLANT":
            queries.append({**pair,"status":"SOURCE_PLANT_TAXON_NOT_EXACT_NOT_ZERO"})
            continue
        url=occurrence_url(pair,info["accepted_species_key"])
        try:
            v=parse_plant_gadm(get(url))
        except ValueError as e:
            v={"status":"INVALID_PLANT_FACET_DO_NOT_INTERPRET","exception":str(e)[:250],"url":url}
        queries.append({**pair,"accepted_species_key":info["accepted_species_key"],**v})
        time.sleep(.7)
    report=build_report(cases,plants,pairs,juvenile,all48,meta,matching,queries)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("GBIF_FIVE_HOST_GADM_SUMMARY="+json.dumps({k:v for k,v in report.items()
        if k!="per_source_host_country"},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
