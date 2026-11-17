#!/usr/bin/env python3
"""Reconcile six historical Sri Lanka butterfly names using the 2024-2026
regional nomenclature and GBIF v1 matching + exact candidate search.

This is taxonomic name coverage, NOT a GBIF native-range or endemism
classifier. Neither historical E labels nor contemporary regional endemic
descriptions should be projected to every island or historical date.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1]
REGIONAL="results/butterfly_GBIF_effort_origin_gate_20261117/SRI_LANKA_2006_SIX_UNRESOLVED_NAME_CURRENT_SOURCES_AND_HOSTS_20261117.csv"
PANEL="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
PINS={
 REGIONAL:"1ffa9c008a8ed0244470e45e0f922183f093a79b",
 PANEL:"b6a3c98eb76864052d45d520a55b7d830816126c"
}
BASE="https://api.gbif.org/v1/species"

def load(path):
    f=ROOT/path
    b=f.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    if sha!=PINS[path]:raise ValueError("Frozen Git blob mismatch "+path)
    with f.open(encoding="utf-8-sig",newline="") as h:
        return list(csv.DictReader(h))

def get(url):
    exc=None
    for attempt in range(3):
        try:
            req=Request(url,headers={
                "Accept":"application/json",
                "User-Agent":"chocho-six-srilanka-name-crosswalk/2026-11-17; publication audit"
            })
            with urlopen(req,timeout=35) as res:
                raw=res.read(512001)
            if len(raw)>512000:raise ValueError("Oversized GBIF species query")
            data=json.loads(raw.decode())
            if not isinstance(data,dict):raise ValueError("Unexpected GBIF response schema")
            return data,hashlib.sha256(raw).hexdigest()
        except (HTTPError,URLError,TimeoutError,json.JSONDecodeError) as e:
            exc=e
            if isinstance(e,HTTPError) and e.code in (400,401,403,404):
                raise RuntimeError(f"GBIF permanent HTTP {e.code} {url}") from e
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("GBIF source lookup unavailable "+str(exc))

def accepted_match_key(obj):
    # GBIF name matching is a taxonomic index, not an endemicity database.
    if (obj.get("matchType")!="EXACT" or obj.get("rank")!="SPECIES"
        or obj.get("status") not in ("ACCEPTED","SYNONYM")
        or not isinstance(obj.get("confidence"),(int,float))
        or obj["confidence"]<90):
        return None
    if obj.get("status")=="SYNONYM":
        key=obj.get("acceptedUsageKey")
    else:key=obj.get("usageKey")
    return int(key) if isinstance(key,int) else None

NUB_DATASET_UUID="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"

def search_candidate_keys(obj,name):
    # The species/search endpoint includes usage keys FROM MANY CHECKLISTS.
    # Only a verified backbone nubKey may be compared with frozen GBIF speciesKey.
    records=obj.get("results",[])
    if not isinstance(records,list):raise ValueError("GBIF search results are not a list")
    same=[]
    for record in records:
        if not isinstance(record,dict):continue
        canonical=str(record.get("canonicalName","")).strip().casefold()
        if canonical!=name.casefold() or record.get("rank")!="SPECIES":
            continue
        status=record.get("taxonomicStatus",record.get("status"))
        if status not in ("ACCEPTED","SYNONYM"):
            continue
        key=record.get("acceptedKey") if status=="SYNONYM" else record.get("key")
        if isinstance(key,int):
            same.append({"checklist_acceptedKey":key,"status":status,
                         "scientificName":record.get("scientificName"),
                         "canonicalName":record.get("canonicalName"),
                         "key":record.get("key"),
                         "datasetKey":record.get("datasetKey"),
                         "nubKey":record.get("nubKey"),
                         "taxonKey":record.get("taxonKey")})
    return same

def independently_verify_nub_key(nub):
    u=BASE+"/"+str(nub)
    record,sha=get(u)
    status=record.get("taxonomicStatus",record.get("status"))
    valid=(record.get("datasetKey")==NUB_DATASET_UUID and
            record.get("rank")=="SPECIES" and
            status in ("ACCEPTED","SYNONYM") and
            record.get("key")==nub)
    if valid and status=="SYNONYM":
        accepted=record.get("acceptedKey")
    elif valid:
        accepted=nub
    else:
        accepted=None
    if not isinstance(accepted,int):
        accepted=None
    return accepted,{"url":u,"response_sha256":sha,"data":record,
                     "verified_GBIF_backbone":valid}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,default=Path("_srilanka_six_modern_names_GBIF"))
    args=ap.parse_args()
    regional=load(REGIONAL)
    orig=load(PANEL)
    original_keys={int(z["speciesKey"]) for z in orig}
    original_names={z["species"] for z in orig}
    if (len(regional),len(original_keys),len(original_names),len(orig))!=(6,225,225,11475):
        raise ValueError("Frozen source populations changed")
    output=[]
    raw_receipts={}
    seen=set()
    for item in regional:
        h=item["historical_2006_binomial"];m=item["current_regional_species_name"]
        if h in seen or len(h.split())!=2 or len(m.split())!=2:
            raise ValueError("Bad/missing historical name "+h)
        seen.add(h)
        if not item["endemic_status_evidence"] or not item["endemic_source_url"].startswith("https://"):
            raise ValueError("Modern independently sourced regional endemic evidence missing")
        matchurl=BASE+"/match?"+urlencode({"name":m,"kingdom":"Animalia","rank":"SPECIES","strict":"true","verbose":"true"})
        searchurl=BASE+"/search?"+urlencode({"q":m,"rank":"SPECIES","limit":50})
        ma,masha=get(matchurl)
        se,sesha=get(searchurl)
        direct=accepted_match_key(ma)
        alternatives=search_candidate_keys(se,m)
        unique_search_keys=sorted({r["checklist_acceptedKey"] for r in alternatives})
        nubs=sorted({r["nubKey"] for r in alternatives if isinstance(r.get("nubKey"),int)})
        backbone_receipt=None
        if direct is not None:
            key=direct
            provenance="strict_GBIF_backbone_exact_accepted_or_synonym"
        elif len(nubs)==1:
            key,backbone_receipt=independently_verify_nub_key(nubs[0])
            provenance=("verified_GBIF_backbone_nubKey_from_checklist_search"
                        if key is not None else
                        "checklist_search_nubKey_not_verified_as_species")
        else:
            key=None
            provenance="search_name_only_no_verified_backbone_key"
        relation="unresolved"
        if isinstance(key,int):
            relation="same_accepted_key_as_original225" if key in original_keys else "no_accepted_key_overlap_with_original225"
        out={
            "historical_2006_binomial":h,"regional_modern_species_name":m,
            "historical_to_modern_name_evidence":item["relationship_historical_to_regional"],
            "regional_endemic_status":item["endemic_status_evidence"],
            "regional_endemic_source_url":item["endemic_source_url"],
            "GBIF_strict_matched_name":ma.get("scientificName"),
            "GBIF_match_type":ma.get("matchType"),
            "GBIF_status":ma.get("status"),
            "GBIF_confidence":ma.get("confidence"),
            "GBIF_strict_candidate_key":direct,
            "GBIF_exact_search_checklist_usage_keys_NOT_BACKBONE":unique_search_keys,
            "GBIF_search_nubKeys":nubs,
            "GBIF_backbone_lookup_verified":bool(key is not None),
            "GBIF_accepted_key_if_usable":key,
            "GBIF_resolution_status":provenance,
            "overlap_to_original225_speciesKey":relation,
            "documented_larval_hosts":item["observed_larval_host_names"],
            "recorded_host_dependency":item["resource_dependency_interpretation"],
            "native_range_global_review":"NOT_PERFORMED",
            "same_island_breeding_establishment":"NOT_INFERRED_FROM_REGIONAL_CATALOGUE",
            "GBIF_match_url":matchurl,
            "GBIF_search_url":searchurl,
            "match_response_sha256":masha,
            "search_response_sha256":sesha,
        }
        output.append(out)
        raw_receipts[h]={"modern_name":m,"match":{"url":matchurl,"sha256":masha,"data":ma},
                         "search":{"url":searchurl,"sha256":sesha,"data":se},
                         "verified_backbone_usage":backbone_receipt}
        print("SIX_NAME",h,"=>",m,provenance,"key",key,relation,flush=True)
    categories={}
    for row in output:
        key=row["GBIF_resolution_status"]
        categories[key]=categories.get(key,0)+1
    summary={
        "schema":"chocho_SriLanka2006_six_historical_E_to_modern_2026_species_with_GBIF_backbone_source_gate_v1",
        "status":"MODERN_LOCAL_TAXONOMY_AND_HOST_SOURCE_MATCH_NOT_NATIVE_RANGE",
        "retrieved_utc":datetime.now(timezone.utc).isoformat(),
        "original_panel":{"butterflies":225,"species_island_cells":11475,
                          "source_git_blobs":PINS},
        "historical_unresolved_six":len(output),
        "resolution_categories":categories,
        "taxa":output,
        "limits":[
          "Historic 2006 E-marked names are not interchangeable with GBIF accepted names without outside contemporary taxon resolution.",
          "Peer reviewed/current Sri Lankan endemic claims are source assertions, not range-wide specimen reanalysis.",
          "Six species cover only a preselected historic local endemic shortlist, not the complete modern island fauna.",
          "Host plant records from literature cannot by themselves quantify feeding survival, effective host breadth or colonisation.",
          "Do not infer true absence, butterfly native range or population establishment from original225 selected species names.",
          "GBIF taxon names can reflect snapshot variation and nomenclatural differences.",
      "GBIF species/search checklist usage keys are NOT the original GBIF Backbone speciesKeys: compare only strict species/match or verified NUB dataset IDs."
        ]
    }
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"SRI_LANKA_SIX_MODERN_GBIF_BACKBONE_AUDIT.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    (args.out/"SRI_LANKA_SIX_GBIF_RAW_QUERY_RECEIPTS.json").write_text(
        json.dumps(raw_receipts,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with (args.out/"SRI_LANKA_SIX_MODERN_GBIF_SOURCE_CLASSIFICATION.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]))
        w.writeheader()
        for record in output:
            record=record.copy()
            record["GBIF_exact_search_checklist_usage_keys_NOT_BACKBONE"]=";".join(map(str,record["GBIF_exact_search_checklist_usage_keys_NOT_BACKBONE"]))
            record["GBIF_search_nubKeys"]=";".join(map(str,record["GBIF_search_nubKeys"]))
            w.writerow(record)
    print("BEGIN_SIX_GBIF_REGIONAL_TAXONOMY_COUNTS",json.dumps(categories,sort_keys=True),flush=True)
    print("END_SIX_GBIF_REGIONAL_TAXONOMY_COUNTS",flush=True)

if __name__=="__main__":
    main()
