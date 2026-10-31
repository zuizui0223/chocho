#!/usr/bin/env python3
"""Recover an author-created GBIF taxonKey IN predicate for migrant butterflies.

Indirect source: Chowdhury et al. 2025, author-code source for GBIF download.
The author supplied a list derived from the migrant-butterfly review, matched
names to GBIF usage keys, and used those keys in a single occurrence download.
This is NOT the original review supplement, a nonmigrant classifier, or an
independent observation of butterfly migration to an island.

No occurrence records or original review spreadsheet are downloaded.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import Counter
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import csv
import hashlib
import io
import json
import re
import time

HERE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[4]
ORIGIN="https://github.com/ShawanChowdhury/MigratoryButterfliesPA/blob/main/01.gbifData.R"
ORIGIN_GIT_BLOB="f972cf78d1a0e7ce3293ca1036025da7a09b74f8"
GBIF_DOWNLOAD="0297462-220831081235567"
API="https://api.gbif.org/v1"
PANEL="provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"
PANEL_BLOB="39fe84981cab2b0e47063399f0a5c0c9731e800d"
GAPS="provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/matched_species_screen.csv"
HEADERS={"Accept":"application/json","User-Agent":"chocho-butterfly-maritime-evidence/1.0"}
BINOMIAL=re.compile(r"^\s*([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)\s*$")

def request_json(url,attempts=3):
    last=None
    for i in range(attempts):
        try:
            req=Request(url,headers=HEADERS)
            with urlopen(req,timeout=25) as r:
                if r.status!=200:raise ValueError(f"Non-200 source HTTP: {r.status}")
                return json.load(r)
        except (URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
            last=exc
            if isinstance(exc,HTTPError) and exc.code not in (408,425,429,500,502,503,504):
                break
            if i+1<attempts:time.sleep(0.4*(2**i))
    raise ValueError(f"Source {url}: {type(last).__name__} {str(last)[:180]}")

def taxon_keys_from_predicate(document):
    d=document.get("request")
    if not isinstance(d,dict):
        raise ValueError("Missing GBIF download request object")
    if isinstance(d.get("predicate"),str):
        pred=json.loads(d["predicate"])
    else:
        pred=d.get("predicate",d)
    if not isinstance(pred,dict):
        raise ValueError("Missing GBIF predicate")
    hits=[]
    def visit(node):
        if isinstance(node,list):
            for x in node:visit(x)
        elif isinstance(node,dict):
            if str(node.get("type","")).upper()=="IN" and str(node.get("key","")).upper() in ("TAXON_KEY","TAXONKEY"):
                values=node.get("values")
                if not isinstance(values,list):
                    raise ValueError("Taxon IN should have explicit values")
                hits.append(values)
            for k,v in node.items():
                if k in ("predicates","predicate","value"):
                    visit(v)
    visit(pred)
    if len(hits)!=1:
        raise ValueError(f"Expected exactly one explicit taxon IN predicate, got {len(hits)}")
    values=[str(x) for x in hits[0]]
    if not (200<=len(values)<=800) or len(set(values))!=len(values) or not all(v.isdigit() for v in values):
        raise ValueError(f"Taxon keys are not a unique 200–800 list: {len(values)}")
    return values

def taxon_record(key):
    try:
        data=request_json(f"{API}/species/{key}")
        name=data.get("canonicalName","")
        if not isinstance(name,str):
            name=""
        m=BINOMIAL.fullmatch(name)
        status="RESOLVED_EXACT_BINOMIAL" if m and str(data.get("rank","")).upper()=="SPECIES" else "UNRESOLVED_NON_SPECIES_OR_NAME"
        return (str(key),name.strip() if m else None,status)
    except Exception as err:
        return (str(key),None,"TAXON_ENDPOINT_FAILED:"+type(err).__name__)

def input225():
    raw=(ROOT/PANEL).read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if blob!=PANEL_BLOB:raise ValueError("Source 225 candidate name file changed")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf8"))))
    if len(rows)!=225:raise ValueError("Butterfly candidates changed")
    return rows

def analyze():
    receipt={"schema":"chocho_migrant_2025_author_GBIF_predicate_crosswalk_v1",
             "status":"NOT_RUN","download_api":f"{API}/occurrence/download/{GBIF_DOWNLOAD}",
             "author_source_url":ORIGIN,"author_source_git_blob":ORIGIN_GIT_BLOB,
             "candidate_225_blob":PANEL_BLOB,
             "claim_level":"INDIRECT_PUBLISHED_AUTHOR_SELECTED_TAXON_KEYS_NOT_ORIGINAL_REVIEW_TABLE",
             "not_listed_must_be_unknown":True,
             "dispersal_on_islands_not_directly_measured":True,
             "no_absence_or_causal_steppingstone_inference":True}
    try:
        metadata=request_json(receipt["download_api"])
        keys=taxon_keys_from_predicate(metadata)
        receipt.update({"GBIF_metadata_doi":metadata.get("doi"),"gbif_download_status":metadata.get("status"),
                        "GBIF_predicate_taxon_keys":len(keys),
                        "GBIF_taxon_keys_sha256":hashlib.sha256((",".join(keys)).encode()).hexdigest()})
        records=[]
        with ThreadPoolExecutor(max_workers=12) as pool:
            tasks=[pool.submit(taxon_record,k) for k in keys]
            for fut in as_completed(tasks):records.append(fut.result())
        records.sort(key=lambda x:int(x[0]))
        good=[t[1] for t in records if t[2]=="RESOLVED_EXACT_BINOMIAL"]
        counters=Counter(t[2] for t in records)
        if counters["RESOLVED_EXACT_BINOMIAL"]<200:
            raise ValueError("Fewer than 200 exactly resolved species; do not classify")
        if len(good)!=len(set(good)):
            raise ValueError("Nonunique GBIF canonical names require taxonomy QA")
        candidate=input225()
        valid=set(good)
        matched=[x for x in candidate if x["species"] in valid]
        rows=[{"species":x["species"],"speciesKey":x["taxon_key"],
               "author_2025_GBIF_requested_taxon_key_proxy":
                   "REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST" if x["species"] in valid
                   else "NOT_LISTED_OR_TAXON_UNRESOLVED_UNKNOWN"}
              for x in candidate]
        table=list(csv.DictReader((ROOT/GAPS).open(encoding="utf8",newline="")))
        if len(table)!=193:raise ValueError("193 species gap panel changed")
        withnames={x["species"]:x for x in table}
        documented=[withnames[x["species"]] for x in matched if x["species"] in withnames]
        uncertain=[x for x in table if x["species"] not in {z["species"] for z in documented}]
        def stat(grp):
            n=sum(int(x["n"]) for x in grp)
            k=sum(int(x["far500"]) for x in grp)
            return {"species":len(grp),"indexed_events":n,"record_gaps_over_500km":k,
                    "record_gap_fraction_event_weighted":k/n if n else None,
                    "mean_species_record_gap_fraction":sum(int(x["far500"])/int(x["n"]) for x in grp)/len(grp) if grp else None}
        receipt.update({"status":"SOURCE_RECOVERED_INDIRECT_MIGRANT_PROXY",
                        "source_taxon_resolution":dict(counters),
                        "unique_exact_binomial_keys":len(good),
                        "candidate_225_exact_matched":len(matched),
                        "candidate_225_not_matched_or_unresolved_UNKNOWN":225-len(matched),
                        "gap_193_comparison":{
                            "matched_author_selected_migrant_names":stat(documented),
                            "unmatched_or_unknown_names_NOT_NONMIGRANTS":stat(uncertain)}
                        })
        return receipt,rows
    except Exception as e:
        receipt.update({"status":"GBIF_MIGRANT_PROXY_NOT_VERIFIED","source_blocker":str(e)[:400]})
        return receipt,[]

def main():
    result,rows=analyze()
    (HERE/"migration_gbif_predicate_audit.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    if rows:
        with (HERE/"migration_gbif_author2025_crosswalk.csv").open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader();w.writerows(rows)
    print(json.dumps({k:result.get(k) for k in
        ("status","GBIF_predicate_taxon_keys","candidate_225_exact_matched",
         "gap_193_comparison","source_blocker")},ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
