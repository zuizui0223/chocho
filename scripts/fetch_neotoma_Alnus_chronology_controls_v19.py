#!/usr/bin/env python3
"""Retrieve exact official Neotoma chronology objects for fixed 21 source-labelled samples.

This audits API source metadata/chron controls, never declares independently
calibrated pollen ages, occupied 21ka Alnus glutinosa trees or genetics.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,io,re,time
from collections import Counter,defaultdict
from pathlib import Path
import requests

SHA_AGES="85cb546d43e57c8ff195c09d5dba253e6e2cc88a5cbae4e1eddb762e5c67c6e4"
SHA_SUMMARY="934c29c61a9ee8603f54758fcf27f85f1970bd512e39cebf9a44fdebd9a195cf"
BASE="https://api.neotomadb.org/v2.0/data"
CHRONID_KEYS={"chronologyid","chronology_id","chronid"}
CONTROLS_KEYS={"chroncontrols","chronologycontrols","chroncontrol","chronologycontrol"}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def load_csv(path):
    with Path(path).open("r",newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))

def source_gate(contract,ages,summary):
    original=contract["parent_source"]
    if (contract["status"]!="SOURCE_CHRONOLOGY_PREFLIGHT_FROZEN_BEFORE_CONTROL_EVIDENCE"
        or summary["source_original_datasets"]!=6 or
        summary["source_labelled_candidate_age_rows_n"]!=21 or
        summary["accepted_21ka_Alnus_host_occupied_population_sites"]!=0
        or len(ages)!=788):
        raise ValueError("original 6-source chronology and 21 labels drifted")
    expected={str(k):[str(x) for x in v] for k,v in original["source_dataset_to_all_chronology_ids"].items()}
    if len(expected)!=6 or len({x for v in expected.values() for x in v})!=11:
        raise ValueError("expected 11 chronology IDs on six source datasets")
    candidate=[r for r in ages
        if r["source_labelled_18_24ka_cal_BP_candidate"]=="True"
        and r["datum_positive"]=="True"]
    if len(candidate)!=21 or len({(r["datasetid"],r["sampleid"]) for r in candidate})!=21:
        raise ValueError("source positive 21 candidate samples not distinct")
    for r in candidate:
        if (str(r["source_chronologyid"]) not in expected[r["datasetid"]]
            or str(r["source_chronologyid"])!=str(original["source_dataset_to_labelled_candidate_chronology_id"][r["datasetid"]])
            or r["source_agetype_label"].casefold()!="calibrated radiocarbon years bp"):
            raise ValueError("source sample-age candidate chronology ID/type changed")
    assert len(summary["source_sha256"])==6
    return candidate,expected

def response_scan(parsed,cap_nodes=120000):
    """Structural hints only; key name occurrence is not age calibration proof."""
    result={"explicit_chronologyids":[],"observed_paths":[],
            "chroncontrol_container_count":0,
            "chroncontrol_record_count":0,
            "reliability_bound_key_count":0,
            "data_root_type":None}
    payload=parsed.get("data") if isinstance(parsed,dict) else parsed
    result["data_root_type"]=type(payload).__name__
    seen=set();pathcount=Counter();stack=[(payload,"$.data")];steps=0
    while stack:
        v,path=stack.pop()
        steps+=1
        if steps>cap_nodes:raise ValueError("Neotoma chronology JSON source exceeds scan cap")
        if isinstance(v,dict):
            for k,x in v.items():
                kn=k.casefold()
                pathcount[path+"."+k]+=1
                if kn in CHRONID_KEYS and (isinstance(x,(str,int)) and not isinstance(x,bool)):
                    if str(x).strip().isdigit():seen.add(str(int(str(x).strip())))
                if kn in CONTROLS_KEYS and isinstance(x,(list,dict)):
                    result["chroncontrol_container_count"]+=1
                    result["chroncontrol_record_count"]+=len(x) if isinstance(x,list) else 1
                if kn in {"ageboundyounger","ageboundolder","ageboundyoungerbp","ageboundolderbp"}:
                    result["reliability_bound_key_count"]+=1
                if isinstance(x,(dict,list)):stack.append((x,path+"."+k))
        elif isinstance(v,list):
            for w in v:
                if isinstance(w,(dict,list)):stack.append((w,path+"[]"))
    result["explicit_chronologyids"]=sorted(seen,key=int)
    result["observed_paths"]=[{"path":p,"count":n} for p,n in pathcount.most_common(45)]
    result["scanned_nodes"]=steps
    return result

def fetch(session,url,max_bytes,timeout,attempts=2):
    errors=[]
    for i in range(attempts):
        try:
            with session.get(url,headers={"Accept":"application/json"},timeout=timeout,stream=True) as r:
                code=r.status_code
                if code>=400:raise requests.HTTPError(f"HTTP_{code}",response=r)
                if "json" not in r.headers.get("Content-Type","").casefold():
                    raise ValueError("HTTP response not JSON")
                body=bytearray()
                for part in r.iter_content(chunk_size=65536):
                    if part:
                        body.extend(part)
                        if len(body)>max_bytes:raise ValueError("source size limit exceeded")
                payload=json.loads(body)
                if not isinstance(payload,(dict,list)):
                    raise ValueError("source response JSON not mapping/list")
                data=payload.get("data") if isinstance(payload,dict) else payload
                if data is None:raise ValueError("source response missing data")
                return {"status":"SOURCE_RECEIVED_UNCLASSIFIED","http":code,"url":r.url,
                    "bytes":len(body),"sha256":hashlib.sha256(body).hexdigest(),
                    "headers_content_type":r.headers.get("Content-Type",""),
                    "source_structure":response_scan(payload)},bytes(body)
        except (requests.RequestException,ValueError,json.JSONDecodeError) as exc:
            errors.append(type(exc).__name__+": "+str(exc)[:240])
            if i+1<attempts:time.sleep(1)
    return {"status":"SOURCE_UNAVAILABLE_NOT_CHRONOLOGY_ABSENCE",
            "request_url":url,"errors":errors},None

def store_gzip(path,body):
    with Path(path).open("wb") as fd:
        with gzip.GzipFile(fileobj=fd,filename="",mode="wb",mtime=0) as g:
            g.write(body)

def run(session,contract,age_rows,summary,outdir):
    candidates,expected=source_gate(contract,age_rows,summary)
    outdir=Path(outdir);rawdir=outdir/"raw_api_responses";rawdir.mkdir(parents=True,exist_ok=True)
    records=[]
    cap=contract["api"]["max_response_bytes_per_endpoint"]
    timeout=(contract["api"]["connect_timeout_seconds"],contract["api"]["read_timeout_seconds"])
    attempts=contract["api"]["request_retries"]
    for did in sorted(expected,key=int):
        jobs=[("dataset_chronologies",did,f"{BASE}/datasets/{did}/chronologies")]
        jobs.extend(("chronology",chron,f"{BASE}/chronologies/{chron}") for chron in expected[did])
        for kind,id_,url in jobs:
            info,body=fetch(session,url,cap,timeout,attempts=attempts)
            record={"kind":kind,"requested_id":id_,"expected_datasetid":did,**info}
            if body is not None:
                filename=f"{kind}_{id_}.json.gz"
                store_gzip(rawdir/filename,body)
                record["source_gzip_file"]=filename
                record["source_gzip_sha256"]=sha(rawdir/filename)
                explicit=info["source_structure"]["explicit_chronologyids"]
                if kind=="chronology":
                    record["expected_id_explicitly_present"]=id_ in explicit
                    record["status"]="SOURCE_RETURNED_EXPLICIT_CHRONOLOGY_ID" if id_ in explicit else "SOURCE_NOT_IDENTIFIABLE_BY_REQUESTED_CHRONOLOGY_ID"
                else:
                    record["known_dataset_chronology_overlap_n"]=len(set(expected[did])&set(explicit))
            records.append(record)
            print(json.dumps({"kind":kind,"id":id_,"status":record["status"],
                "source_sha256":record.get("sha256"),
                "chronologyids":record.get("source_structure",{}).get("explicit_chronologyids")},sort_keys=True),flush=True)
    fake=str(contract["api"]["negative_control_chronology_id"])
    info,body=fetch(session,f"{BASE}/chronologies/{fake}",cap,timeout,attempts=1)
    no_false_positive=True
    if body is not None:
        explicit=info["source_structure"]["explicit_chronologyids"]
        if fake in explicit:
            # An untrusted route echoes a nonexisting request ID: don't use it.
            no_false_positive=False
        store_gzip(rawdir/"invalid_chronology_negative_control.json.gz",body)
    receipt={"schema":"chocho_neotoma_11_chronology_source_preflight_v19",
      "status":"SOURCE_CHRONOLOGY_METADATA_NOT_INDEPENDENT_AGE_CALIBRATION",
      "candidate_pollen_sample_age_rows":21,
      "candidate_positive_sample_ids":len({(r["datasetid"],r["sampleid"]) for r in candidates}),
      "source_datasets":6,
      "source_chronology_ids":11,
      "requested_official_chronology_endpoints":17,
      "chronology_endpoints_returned_JSON":sum(r["kind"]=="chronology" and "source_structure" in r for r in records),
      "chronology_endpoints_identifiable":sum(r["kind"]=="chronology" and r.get("expected_id_explicitly_present") is True for r in records),
      "chronology_endpoints_unavailable":sum(r["kind"]=="chronology" and "source_structure" not in r for r in records),
      "dataset_chronology_endpoints_returned_JSON":sum(r["kind"]=="dataset_chronologies" and "source_structure" in r for r in records),
      "fake_chronology_negative_control_status":info["status"],
      "fake_chronology_negative_control_is_safe":no_false_positive,
      "original_v18_21_candidate_age_label_rows_unchanged":True,
      "independently_verified_age_depth_model_count":0,
      "independently_verified_LGM_Alnus_occupied_plant_site_count":0,
      "original_genetic_response_used":False,
      "source_urls":["https://api.neotomadb.org/api-docs/",
                    "https://open.neotomadb.org/manual/chronology-age-related-tables-1.html"],
      "limitations":["Source chronology metadata may expose age controls but does not independently validate calibration curves, sediment deposition or Alnus species-diagnostic pollen",
        "Source API may return partial/misleading fields; never count a successful HTTP response as independent validation",
        "Source label candidate samples from 18-24ka do not establish local occupied Alnus glutinosa populations",
        "No original TTF negative host specificity or genetic outcome was reopened"]}
    with (outdir/"all_11_chronology_plus_six_dataset_source_routes.json").open("w",encoding="utf8") as f:
        json.dump(records,f,indent=2,sort_keys=True);f.write("\n")
    (outdir/"summary.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--contract",type=Path,required=True)
    a.add_argument("--ages",type=Path,required=True)
    a.add_argument("--v18-summary",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    if sha(args.ages)!=SHA_AGES or sha(args.v18_summary)!=SHA_SUMMARY:
        raise ValueError("original v18 source sample ages changed")
    contract=json.loads(args.contract.read_text())
    result=run(requests.Session(),contract,load_csv(args.ages),
               json.loads(args.v18_summary.read_text()),args.out)
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
