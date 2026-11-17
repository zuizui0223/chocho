#!/usr/bin/env python3
"""GBIF backbone reconciliation of 18 E-marked Sri Lanka historical butterfly names.

Independent 2006 E labels are NOT modern endemism/origin labels. Only classify
whether accepted GBIF backbone taxonKey is represented in frozen 225 speciesKey
source. Fail closed on fuzzy, doubtful, coarse-rank, low-confidence matches.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1]
HIST="results/butterfly_GBIF_effort_origin_gate_20261117/SRILANKA_2006_HISTORICAL_ENDEMIC_E_18_NAME_MATCH_20261117.csv"
GBIF="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
PINS={HIST:"59e2c707f06e32eee04b45ffe3eda0bbcef798cf",
      GBIF:"b6a3c98eb76864052d45d520a55b7d830816126c"}
API="https://api.gbif.org/v1/species/match"
EXPECTED_SOURCE_N=18

def read(path):
    p=ROOT/path
    b=p.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    if sha!=PINS[path]:raise ValueError("Source pin mismatch "+path)
    with p.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))

def match(name):
    params={"name":name,"kingdom":"Animalia","rank":"SPECIES","strict":"true","verbose":"true"}
    url=API+"?"+urlencode(params)
    last=None
    for attempt in range(3):
        try:
            req=Request(url,headers={"User-Agent":"chocho-public-taxonomy-source-audit/2026-11-17 (research; 18 taxa only)","Accept":"application/json"})
            with urlopen(req,timeout=30) as resp:
                payload=resp.read(512_000)
                if len(payload)>=512_000:raise ValueError("Oversized GBIF species response")
            obj=json.loads(payload.decode("utf-8"))
            if not isinstance(obj,dict):raise ValueError("GBIF API schema expected object")
            return url,obj,hashlib.sha256(payload).hexdigest()
        except (HTTPError,URLError,TimeoutError,UnicodeDecodeError,json.JSONDecodeError) as exc:
            last=exc
            if isinstance(exc,HTTPError) and exc.code in (400,401,403,404):
                raise RuntimeError("GBIF taxon lookup permanent API error "+name+" "+str(exc)) from exc
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("GBIF taxon lookup unavailable "+name+" "+str(last))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output-dir",type=Path,default=Path("_historical_SriLanka_GBIF_taxonomic_reconciliation"))
    args=ap.parse_args()
    hist=read(HIST)
    gbif=read(GBIF)
    if (len(hist),len(gbif),len({x["speciesKey"] for x in gbif}))!=(18,11475,225):
        raise ValueError("Frozen source populations changed")
    frozen_keys={int(z["speciesKey"]) for z in gbif}
    if any(z["current_accepted_name_reviewed"]!="no" or z["current_endemicity_reviewed"]!="no" for z in hist):
        raise ValueError("Historical E labels must remain explicitly unreviewed")
    rows=[]
    raw_response={}
    for i,z in enumerate(hist):
        name=z["historical_2006_binomial"]
        url,response,sha=match(name)
        usage=response.get("usageKey")
        accepted=response.get("acceptedUsageKey")
        status=response.get("status","UNKNOWN")
        typ=response.get("matchType","NONE")
        rank=response.get("rank","UNKNOWN")
        conf=response.get("confidence")
        candidate=(accepted if status=="SYNONYM" else usage)
        ok=(typ=="EXACT" and rank=="SPECIES" and status in ("ACCEPTED","SYNONYM") and
            isinstance(candidate,int) and isinstance(conf,(int,float)) and conf>=90)
        if not ok:category="needs_taxonomic_review"
        elif candidate in frozen_keys:category="included_in_original225_by_current_GBIF_accepted_key"
        else:category="not_in_original225_by_current_GBIF_accepted_key"
        o={
          "historical_2006_binomial":name,
          "2006_page":int(z["printed_page"]),
          "GBIF_legacy_match_type":typ,"GBIF_matched_rank":rank,"GBIF_taxonomic_status":status,
          "GBIF_match_confidence":conf,"GBIF_usageKey":usage,"GBIF_acceptedUsageKey":accepted,
          "reconciled_candidate_accepted_key":candidate if ok else None,
          "GBIF_canonicalName":response.get("canonicalName"),
          "GBIF_scientificName":response.get("scientificName"),
          "taxon_key_comparison_to_original225":category,
          "current_Sri_Lanka_endemicity":"UNRESOLVED_NOT_FROM_GBIF_NAME_MATCH",
          "current_Sri_Lanka_native_origin":"UNRESOLVED_NOT_FROM_GBIF_NAME_MATCH",
          "api_url":url,"api_response_sha256":sha
        }
        rows.append(o)
        raw_response[name]={"url":url,"sha256":sha,"response":response}
        print("GBIF_TAXON",i+1,"/",len(hist),name,category,candidate,flush=True)
    cnt={tag:sum(x["taxon_key_comparison_to_original225"]==tag for x in rows) for tag in
         ("included_in_original225_by_current_GBIF_accepted_key",
          "not_in_original225_by_current_GBIF_accepted_key",
          "needs_taxonomic_review")}
    if sum(cnt.values())!=18:raise ValueError("Unaccounted historical E names")
    result={
      "schema":"chocho_2006_SriLanka_historical_endemic_18_GBIF2026_accepted_taxonKeys_vs_frozen225_v1",
      "status":"GBIF_BACKBONE_TAXON_KEY_SELECTION_SCREEN_NOT_ENDENISM_NATIVE_STATUS",
      "date_retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "frozen_source_git_blobs":PINS,
      "GBIF_api_source":API,"strict_species_Animalia_matching":True,
      "original_225_species_keys":len(frozen_keys),
      "historical_individually_E_marked_names":len(rows),
      "name_resolution_counts":cnt,
      "taxa":rows,
      "limitations":[
        "GBIF accepted backbone name match is not verification of Sri Lanka native origin, current endemicity, range restriction or conservation status.",
        "Historical 2006 source E classification uses taxonomy from 2006, not contemporary species concepts.",
        "Original frozen 225 GBIF speciesKeys from a separate API indexing snapshot may differ after backbone revisions.",
        "This audit tests only the 18 explicitly E-marked names selected from original source, not all 20 endemics in a historical narrative total.",
        "Fuzzy, higher-rank, doubtful or low-confidence matches remain review-needed and are never automatically forced to a source species.",
        "No genus-level name matching, guessed synonym assignments or inferred ecological absences."
      ]}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/"GBIF_BACKBONE_HISTORICAL_18_ENDENIC_SOURCE_AUDIT_20261117.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (args.output_dir/"GBIF_MATCH_RAW18_QUERY_RECEIPTS_20261117.json").write_text(
      json.dumps(raw_response,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    with (args.output_dir/"HISTORIC_ENDENIC_GBIF_TAXON_RECONCILIATION_18.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    print("BEGIN_HISTORICAL_E_TAXON_RECONCILIATION_COUNTS")
    print(json.dumps({"count":len(rows),"counts":cnt},sort_keys=True))
    print("END_HISTORICAL_E_TAXON_RECONCILIATION_COUNTS")

if __name__=="__main__":main()
