#!/usr/bin/env python3
"""Audit and probe Neotoma v2 taxonname filters without accepting fossil presence.

Source-only scientific integrity check:
- audit original 15 taxa x 2 proxy v10 queries for exact duplicated payloads;
- compare same proxy/time-age/limit with two source WCVP modelled names and an
  intentionally nonexistent taxon name (negative control);
- resolve Neotoma accepted binomials in taxon REGISTRY, then use taxon-ID
  occurrence route if and only if IDs resolve exactly.

NO historical occupation, valid age-model, taxon-level pollen identification,
true botanical absence or insect genetic response may be inferred.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import defaultdict,Counter
from pathlib import Path
from urllib.parse import urlsplit,parse_qs
import requests

EXPECTED={
    "prior15":"f53676979098828d9467e261c5d27c60a867a09f72c6e254f6c5c9546b9ee90e",
    "prior_summary":"18e5d38ecaa3951e4433293b6e6ec793b99beb2d43b03600f0f66a0012fec4dc",
    "source_models":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733",
}
SOURCE_ROOT="https://api.neotomadb.org"
BAD="Definitelynotaneotomataxon alpha9999"
NAMES=("Alnus glutinosa","Betula pendula",BAD)
AGE=(18000,24000)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for buf in iter(lambda:fh.read(1<<20),b""):h.update(buf)
    return h.hexdigest()

def read_rows(path):
    with Path(path).open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))

def original_replay(csv_path,summary_path):
    if sha(csv_path)!=EXPECTED["prior15"] or sha(summary_path)!=EXPECTED["prior_summary"]:
        raise ValueError("original Neotoma source metadata checkpoint changed")
    rows=read_rows(csv_path)
    q=json.loads(Path(summary_path).read_text())
    if len(rows)!=30 or (q["attempted_WCVP_taxa"],q["successful_JSON_proxy_queries"],
        q["unavailable_proxy_queries"])!=(15,12,18):
        raise ValueError("source 15x2 API record status changed")
    bodies=defaultdict(list);counts=Counter()
    for r in rows:
        counts[r["status"]]+=1
        if r["status"]=="SOURCE_RESPONSE_UNREVIEWED":
            if r["http_status"]!="200" or not re.fullmatch(r"[0-9a-f]{64}",r["json_sha256"]):
                raise ValueError("successful original response lacks HTTP200/hash")
            params=parse_qs(urlsplit(r["endpoint"]).query)
            if params.get("taxonname")!=[r["name"]] or params.get("datasettype")!=[r["proxy"]]:
                raise ValueError("original Neotoma request provenance mismatch")
            if params.get("ageyoung")!=["18000"] or params.get("ageold")!=["24000"]:
                raise ValueError("original source requested fossil ages changed")
            bodies[r["json_sha256"]].append({"taxon":r["name"],"proxy":r["proxy"]})
    repeat=[{"sha256":k,"query_count":len(rows),"distinct_taxon_names":len({r["taxon"] for r in rows}),
             "proxies":sorted({r["proxy"] for r in rows}),
             "taxa":sorted({r["taxon"] for r in rows})}
             for k,rows in sorted(bodies.items())]
    if (len(repeat)!=2 or sorted(r["query_count"] for r in repeat)!=[5,7]
        or not all(r["distinct_taxon_names"]==r["query_count"] for r in repeat)):
        raise ValueError("original exact duplicate query-body anomaly changed")
    return {"schema":"chocho_neotoma_15_taxa_repeated_content_source_audit_v11",
       "prior_proxy_queries":30,"prior_http200_json":12,
       "prior_request_unavailable":18,"unique_http200_payload_sha256":len(repeat),
       "repeated_payload_clusters":repeat,
       "source_claim":"Nonempty identical returned JSON for distinct accepted plant-name queries prevents confirming taxonname filter integrity",
       "not_a_confirmed_fossil_occurrence":True,"absent_taxon_not_inferred":True,
       "source_sha256":{"prior15":EXPECTED["prior15"],"prior_summary":EXPECTED["prior_summary"]}}

def clean_json_response(response):
    if response.status_code!=200:raise RuntimeError(f"HTTP_{response.status_code}")
    if "json" not in response.headers.get("Content-Type","").lower():
        raise RuntimeError("not_JSON_content_type")
    raw=response.content
    if len(raw)>250000:
        raise RuntimeError("unexpectedly_large_source_metadata")
    obj=response.json()
    if not isinstance(obj,(list,dict)):raise RuntimeError("unexpected_source_response_root")
    if isinstance(obj,dict) and str(obj.get("status","")).lower() in ("error","failed"):
        raise RuntimeError("API_response_declares_error")
    vals=obj if isinstance(obj,list) else obj.get("data")
    if not isinstance(vals,list):raise RuntimeError("missing_source_data_list")
    ids=[]
    for x in vals:
        if not isinstance(x,dict):raise RuntimeError("expected_occurrence_or_taxon_records")
        idval=next((x.get(k) for k in ("occid","occurrenceid","taxonid","id") if x.get(k) is not None),None)
        ids.append(str(idval) if idval is not None else None)
    return obj,{"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),
                "record_count":len(vals),"record_ids":ids,
                "first_keys":sorted(vals[0].keys()) if vals else [],
                "top_keys":sorted(obj.keys()) if isinstance(obj,dict) else [],
                "has_nonempty_data":len(vals)>0}

def live_request(session,url,params):
    try:
        r=session.get(url,params=params,timeout=(7,17),headers={"Accept":"application/json"})
        obj,info=clean_json_response(r)
        return obj,{"status":"HTTP200_SOURCE_CONTENT_UNVERIFIED","actual_url":r.url,**info}
    except (requests.RequestException,RuntimeError,ValueError,KeyError) as exc:
        return None,{"status":"API_UNAVAILABLE_OR_UNEXPECTED_NO_ABSENCE_CLAIM",
                     "error":f"{type(exc).__name__}:{str(exc)[:220]}"}

def strict_registry_taxon_id(payload,name):
    if not isinstance(payload,(dict,list)):return None
    vals=payload.get("data",[]) if isinstance(payload,dict) else payload
    matches=[]
    for x in vals:
        if not isinstance(x,dict):continue
        nm=next((x.get(k) for k in ("taxonname","taxonName","taxon_name","taxon") if isinstance(x.get(k),str)),None)
        tid=next((x.get(k) for k in ("taxonid","taxonId","taxon_id") if isinstance(x.get(k),(int,str))),None)
        if nm==name and tid is not None and str(tid).isdigit():
            matches.append(int(tid))
    return matches[0] if len(set(matches))==1 and matches else None

def audit_live(session,contract):
    x=contract["probe"]
    if (x["positive_WCVP_model_names"]!=["Alnus glutinosa","Betula pendula"]
        or x["sentinel_negative_name"]!=BAD or x["proxy"]!="pollen"
        or (x["ageyoung_BP"],x["ageold_BP"])!=AGE
        or x["query_max_items"]!=2):
        raise ValueError("frozen source integrity probes changed")
    rows=[];raw_lookup={};taxon_ids={}
    for name in NAMES:
        params={"taxonname":name,"limit":2,"offset":0}
        payload,info=live_request(session,SOURCE_ROOT+"/v2.0/data/taxa",params)
        rows.append({"endpoint":"taxon_registry","requested_taxon_name":name,**info})
        tid=strict_registry_taxon_id(payload,name)
        if tid is not None and name!=BAD:taxon_ids[name]=tid
        raw_lookup[(name,"registry")]=payload
    for name in NAMES:
        params={"taxonname":name,"datasettype":"pollen","ageyoung":18000,
                "ageold":24000,"limit":2,"offset":0}
        payload,info=live_request(session,SOURCE_ROOT+"/v2.0/data/occurrences",params)
        rows.append({"endpoint":"filtered_occurrences","requested_taxon_name":name,**info})
        raw_lookup[(name,"occurrence")]=payload
    for name in NAMES[:2]:
        if name not in taxon_ids:
            rows.append({"endpoint":"taxon_id_occurrences","requested_taxon_name":name,
                         "status":"NOT_ATTEMPTED_TAXON_ID_NOT_VERIFIED"})
            continue
        payload,info=live_request(session,SOURCE_ROOT+
            f"/v2.0/data/taxa/{taxon_ids[name]}/occurrences",{"limit":2,"offset":0})
        rows.append({"endpoint":"taxon_id_occurrences","requested_taxon_name":name,
                     "verified_registry_taxonid":taxon_ids[name],**info})
    occult={r["requested_taxon_name"]:r for r in rows if r["endpoint"]=="filtered_occurrences"}
    success={name:occult[name] for name in NAMES
             if occult[name]["status"]=="HTTP200_SOURCE_CONTENT_UNVERIFIED"}
    dup=[]
    for x,aa in success.items():
        for y,bb in success.items():
            if (x<y and aa["has_nonempty_data"] and bb["has_nonempty_data"] and
                    aa["sha256"]==bb["sha256"]):
                dup.append([x,y])
    if any(BAD in pair for pair in dup):
        status="FOSSIL_TAXONNAME_FILTER_INVALID_SENTINEL_MATCHES_REAL_TAXON"
    elif dup:
        status="FOSSIL_TAXONNAME_FILTER_UNVERIFIED_IDENTICAL_DISTINCT_TAXA"
    elif len(success)<3:
        status="SOURCE_TIMEOUT_OR_INCOMPLETE_NO_FILTER_INTEGRITY_CONCLUSION"
    elif success[BAD]["has_nonempty_data"]:
        status="IMPOSSIBLE_TAXON_RETURNS_RECORDS_FILTER_UNVERIFIED"
    else:
        status="PROVISIONAL_DISTINCT_FILTER_RESPONSES_SAMPLES_STILL_UNVERIFIED"
    return rows,{"schema":"chocho_Neotoma_21ka_occurrence_taxon_filter_source_controls_v11",
        "status":status,
        "queried_species":list(NAMES),
        "exact_taxon_registry_ids":taxon_ids,
        "success_occurrence_query_n":len(success),
        "success_records_response_sha_groups":{
            name:{"sha256":row["sha256"],"n":row["record_count"]}
            for name,row in success.items()},
        "identical_distinct_name_response_pairs":dup,
        "fossil_species_21ka_sample_occurrence_verified":False,
        "genetic_response_used":False,
        "sample_age_model_or_pollen_species_identification_verified":False,
        "original_lgm_host_climate_specificity_result_reopened":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--original-status",type=Path,required=True)
    p.add_argument("--original-summary",type=Path,required=True)
    p.add_argument("--models",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if sha(a.models)!=EXPECTED["source_models"]:raise ValueError("original WCVP model source SHA drift")
    c=json.loads(a.contract.read_text())
    if c["status"]!="TECHNICAL_NEGATIVE_CONTROL_BEFORE_NEW_FOSSIL_EVIDENCE":
        raise ValueError("fossil source filter contract not frozen")
    old=original_replay(a.original_status,a.original_summary)
    (a.out/"original_pilot_payload_identity_audit.json").write_text(
        json.dumps(old,indent=2,sort_keys=True)+"\n")
    sess=requests.Session()
    sess.headers.update({"User-Agent":"chocho/strict-source-fossil-species-negative-control/0.1 academic research"})
    rows,result=audit_live(sess,c)
    with (a.out/"live_query_controls.csv").open("w",newline="",encoding="utf8") as f:
        keys=sorted(set().union(*(r.keys() for r in rows)))
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader()
        for row in rows:
            cleaned={k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()}
            w.writerow(cleaned)
    result["original_pilot_audit_status"]="12 HTTP 200 queries collapse to 2 distinct payloads across unrelated requested taxa"
    result["source_openapi"]="https://api.neotomadb.org/swagger.json"
    result["source_original_audit_sha256"]=hashlib.sha256(
        (a.out/"original_pilot_payload_identity_audit.json").read_bytes()).hexdigest()
    result["limitations"]=[
         "No historical species occupancy, absence or fossil calibration conclusion from this source-only endpoint test",
         "A valid Neotoma taxon registry taxonID does not confirm occurrence or per-sample calibrated 21ka age",
         "Neotoma occurrence query age filters can describe dataset envelopes; direct fossil sample age model must be checked",
         "Same returned body for unrelated taxon names under apparent taxon filter invalidates taxon-specific interpretation",
         "No changes to existing TTF/chocho host-vs-alternative-plant negative glacial climate tests"
    ]
    (a.out/"source_filter_integrity_summary.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
