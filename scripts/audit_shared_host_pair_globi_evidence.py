#!/usr/bin/env python3
"""Non-confirmatory GloBI original-event feasibility for exact two-butterfly host use."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
SPECIES=["Junonia coenia","Anartia jatrophae"]
PLANT="Plantago lanceolata"
BASE="https://api.globalbioticinteractions.org/interaction.csv"
FIELDS=[
"source_taxon_name","target_taxon_name","target_taxon_path","interaction_type",
"event_date","latitude","longitude","source_specimen_life_stage","study_source_citation",
"study_citation","study_external_id","source_specimen_occurrence_id","locality"
]
PAGE_SIZE=256
PAGES=2

def exact_species(actual,wanted):
    actual=" ".join(str(actual or "").strip().split())
    return actual==wanted or actual.startswith(wanted+" ")

def stage_larval(x):
    x=str(x or "").lower()
    return "larva" in x or "caterpillar" in x

def valid_geo_and_year(row):
    year=re.search(r"(?<!\d)(19\d{2}|20[0-2]\d)(?!\d)",str(row.get("event_date") or ""))
    try:
        lat=float(row.get("latitude") or "nan")
        lon=float(row.get("longitude") or "nan")
        geo=math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180 and (lat!=0 or lon!=0)
    except (TypeError,ValueError):
        geo=False
    return bool(year and geo)

def count_source_rows(rows,species):
    counts={k:0 for k in ["returned_rows","exact_source","exact_plant","larval","geo_dated","independent_provenance"]}
    ids=[]
    for r in rows:
        counts["returned_rows"]+=1
        if not exact_species(r.get("source_taxon_name"),species):continue
        counts["exact_source"]+=1
        if not exact_species(r.get("target_taxon_name"),PLANT):continue
        counts["exact_plant"]+=1
        if not stage_larval(r.get("source_specimen_life_stage")):continue
        counts["larval"]+=1
        if not valid_geo_and_year(r):continue
        counts["geo_dated"]+=1
        citation=" | ".join(str(r.get(k) or "") for k in
                          ("study_source_citation","study_citation","study_external_id"))
        if not citation.strip(" |"):continue
        if "globalbioticinteractions/hosts" in citation.lower() or "nhm hosts" in citation.lower():
            continue
        counts["independent_provenance"]+=1
        ids.append({"date":r.get("event_date"),"latitude":r.get("latitude"),
                    "longitude":r.get("longitude"),"source":citation[:250],
                    "original_occurrence_id":r.get("source_specimen_occurrence_id")})
    return counts,ids

def fetch(species,offset,opener=urlopen):
    query=urlencode({"sourceTaxon":species,"interactionType":"eats",
                     "includeObservations":"true","limit":PAGE_SIZE,"offset":offset,
                     "fields":",".join(FIELDS)})
    url=BASE+"?"+query
    try:
        request=Request(url,headers={"Accept":"text/csv","User-Agent":"chocho-butterfly-pair-source-gate/1.0"})
        with opener(request,timeout=40) as response:
            data=response.read(12_000_001)
        if len(data)>12_000_000:
            return [],{"url":url,"error":"API returned >12MB"},False
        return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig")))),{
            "url":url,"sha256":hashlib.sha256(data).hexdigest(),"rows_bytes":len(data)},True
    except (HTTPError,URLError,OSError,TimeoutError,UnicodeDecodeError,csv.Error) as e:
        return [],{"url":url,"error":str(e)[:200]},False

def audit():
    output={"schema":"chocho_shared_host_pair_globi_evidence_v01",
            "source":"GloBI observed interactions; this is NOT the frozen HOSTS list",
            "target_plant":PLANT,"source_query_limited":True,
            "same_individual_plant_overlap_demonstrated":False,
            "cross_species_pathogen_transfer_demonstrated":False,
            "ecological_fitness_effect_estimated":False,
            "GEB_PR38_unchanged":True,"species":{}}
    for species in SPECIES:
        allrows=[];pages=[];error=False
        for k in range(PAGES):
            rows,meta,ok=fetch(species,k*PAGE_SIZE)
            pages.append(meta)
            if not ok:
                error=True
                break
            allrows.extend(rows)
            if len(rows)<PAGE_SIZE:
                break
        counts,ids=count_source_rows(allrows,species)
        output["species"][species]={"counts":counts,"pages":pages,
            "last_page_full":len(allrows)>=PAGE_SIZE*PAGES,
            "source_incomplete":bool(error or len(allrows)>=PAGE_SIZE*PAGES),
            "strict_candidates":ids,
            "original_links_not_adjudicated":True}
    output["strict_local_co_use_confirmed"]=False
    output["decision"]="SOURCE_SEARCH_ONLY_NOT_A_VERIFIED_LOCAL_COUSE_RESULT"
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--receipt",type=Path,required=True)
    args=p.parse_args()
    result=audit()
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
