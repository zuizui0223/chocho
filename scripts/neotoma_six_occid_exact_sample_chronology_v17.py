#!/usr/bin/env python3
"""Six frozen Alnus pollen occurrence-ID to full Neotoma sample chronology probe.

Uses ONE preselected smallest numeric v13 occurrence ID for each source dataset.
Only explicit same-typed sample/analysis-unit IDs permit candidate structural
joins. This does NOT verify radiocarbon calibration, exact pollen morphology,
an occupied 21ka Alnus population, or insect genetics.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import requests

CANDIDATE_SHA="3d47fe9a31cb6a04ce60aa0f1d7fc8dd7fb1508d3966d6f3a7dc6942ddfa9efc"
SITES={"4396":"3392","14556":"9834","22662":"14474",
       "22683":"14485","24197":"15750","41351":"23855"}
API="https://api.neotomadb.org/v2.0/data/occurrences/"
SAMPLE_KEYS={"sampleid","sample_id"}
UNIT_KEYS={"analysisunitid","analysisunit_id"}
CHRON_KEYS={"chronologyid","chronology_id","chronid"}
OCC_KEYS={"occid","occurrenceid","occurrence_id"}

def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def rows(path):
    with Path(path).open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))

def keyid(v):
    if isinstance(v,bool) or v is None:return None
    if isinstance(v,(int,str)) and str(v).strip().isdigit():
        return str(int(str(v).strip()))
    return None

def ids(node,fieldnames,max_nodes=200000):
    answer=set()
    stack=[node];seen=0
    while stack:
        node=stack.pop();seen+=1
        if seen>max_nodes:raise ValueError("metadata record source scan exceeded bound")
        if isinstance(node,dict):
            for k,v in node.items():
                if str(k).lower() in fieldnames:
                    item=keyid(v)
                    if item:answer.add(item)
                elif isinstance(v,(dict,list)):stack.append(v)
        elif isinstance(node,list):
            for v in node:
                if isinstance(v,(dict,list)):stack.append(v)
    return answer

def source_datasets(raw,requested_id,expected_site):
    root=raw.get("data") if isinstance(raw,dict) else raw
    if isinstance(root,list) and len(root)==1:root=root[0]
    if not isinstance(root,dict) or not isinstance(root.get("site"),dict):
        raise ValueError("full download absent unique site structure")
    site=root["site"]
    if keyid(site.get("siteid"))!=expected_site:
        raise ValueError("full download siteid provenance changed")
    unit=site.get("collectionunit")
    if not isinstance(unit,dict) or not isinstance(unit.get("dataset"),dict):
        raise ValueError("no explicit collectionunit.dataset in official download")
    d=unit["dataset"]
    if keyid(d.get("datasetid"))!=requested_id or not isinstance(d.get("samples"),list):
        raise ValueError("dataset ID or sample array mismatch")
    return d

def inspect_download(d):
    samples=d["samples"]
    first=next((s for s in samples if isinstance(s,dict)),None)
    age_sources=[a for sample in samples if isinstance(sample,dict)
                 for a in sample.get("ages",[]) if isinstance(a,dict)]
    features={
      "source_full_samples_n":len(samples),
      "source_sample_keys":sorted(first) if first else [],
      "source_sample_age_records_n":len(age_sources),
      "source_first_age_keys":sorted(age_sources[0]) if age_sources else [],
      "source_unique_sample_ids_n":len(ids(samples,SAMPLE_KEYS)),
      "source_unique_analysisunit_ids_n":len(ids(samples,UNIT_KEYS)),
      "source_distinct_age_chronology_id_n":len(ids(age_sources,CHRON_KEYS)),
      "source_distinct_parent_dataset_chronology_id_n":len(ids(d.get("chronologies",[]),CHRON_KEYS)),
      "source_age_types":sorted({str(a.get("agetype")) for a in age_sources if a.get("agetype")}),
      "source_sample_taxon1299_with_datum_n":sum(any(str(v.get("taxonid"))=="1299"
                  for v in s.get("datum",[]) if isinstance(v,dict))
                  for s in samples if isinstance(s,dict))
    }
    return features

def query_occurrence(session,occid,max_bytes,timeout):
    try:
        with session.get(API+occid,stream=True,timeout=timeout,headers={"Accept":"application/json"}) as r:
            r.raise_for_status()
            if "json" not in r.headers.get("content-type","").lower():
                raise ValueError("non-JSON Neotoma occurrence response")
            body=b""
            for c in r.iter_content(chunk_size=32768):
                body+=c
                if len(body)>max_bytes:raise ValueError("source occurrence response exceeded byte cap")
            j=json.loads(body)
            data=j.get("data") if isinstance(j,dict) else j
            if not isinstance(data,(list,dict)):
                raise ValueError("Neotoma response data unavailable")
            return {"status":"SOURCE_HTTP_JSON_AVAILABLE","source_bytes":len(body),
                "source_sha256":hashlib.sha256(body).hexdigest(),"source_url":r.url},body,data
    except (requests.RequestException,ValueError,json.JSONDecodeError) as exc:
        return {"status":"UNKNOWN_SOURCE_NOT_BIOLOGICAL_ABSENCE",
            "reason":type(exc).__name__+": "+str(exc)[:250]},None,None

def match_explicit_sample(occ_data,d):
    matched=[]
    samples=[s for s in d["samples"] if isinstance(s,dict)]
    occ_sampleids=ids(occ_data,SAMPLE_KEYS)
    occ_unitids=ids(occ_data,UNIT_KEYS)
    for row in samples:
        sampleid=keyid(row.get("sampleid")) or keyid(row.get("sample_id"))
        unitid=keyid(row.get("analysisunitid")) or keyid(row.get("analysisunit_id"))
        modes=[]
        if sampleid and sampleid in occ_sampleids:modes.append("EXACT_SAMPLEID")
        if unitid and unitid in occ_unitids:modes.append("EXACT_ANALYSISUNITID")
        if modes:
            ages=[a for a in row.get("ages",[]) if isinstance(a,dict)]
            matched.append({"sampleid":sampleid,"analysisunitid":unitid,
                "explicit_id_join_modes":modes,
                "sample_chronologyids":sorted(ids(ages,CHRON_KEYS),key=int),
                "sample_age_metadata":[{k:a.get(k) for k in (
                   "age","ageolder","ageyounger","agetype","chronologyid")}
                   for a in ages[:12]],
                "candidate_pollen_taxon_1299_in_download_datum":any(
                    str(t.get("taxonid"))=="1299" for t in row.get("datum",[])
                    if isinstance(t,dict))})
    return matched,{"occurrence_explicit_sample_ids":sorted(occ_sampleids,key=int),
                    "occurrence_explicit_analysisunit_ids":sorted(occ_unitids,key=int)}

def gate_original(contract,candidates):
    if (len(candidates)!=27 or contract["status"]!=
       "SOURCE_ONLY_FROZEN_SIX_FIRST_OCCURRENCE_ID_JOIN_TEST_BEFORE_NEW_HTTP"):
        raise ValueError("frozen six dataset study phase or original 27 candidate rows drift")
    chosen=defaultdict(set)
    for row in candidates:
        did=row["source_datasetid"]
        if (did not in SITES or row["source_siteid"]!=SITES[did]
            or row["source_sample_taxonid"]!="1299"
            or row["source_sample_taxonname"]!="Alnus glutinosa"
            or row["source_datasettype"]!="pollen"):
            raise ValueError("source taxon, proxy or dataset identity drift")
        chosen[did].add(row["occid"])
    if set(chosen)!=set(SITES) or sum(map(len,chosen.values()))!=26:
        raise ValueError("not exactly original six datasets and 26 occurrence identities")
    selected={did:min(values,key=int) for did,values in chosen.items()}
    return dict(sorted(selected.items(),key=lambda kv:int(kv[0])))

def run(session,contract,candidates,rawdir,out):
    selection=gate_original(contract,candidates)
    for k,val in contract["source"]["source_v16_dataset_sha256"].items():
        if k not in selection:raise ValueError("missing previous source dataset")
    timeout=(contract["test"]["max_connection_seconds"],contract["test"]["max_read_seconds"])
    out.mkdir(parents=True,exist_ok=True)
    raw_out=out/"individual_source_responses";raw_out.mkdir(parents=True,exist_ok=True)
    results=[]
    for did,occid in selection.items():
        path=rawdir/f"dataset_{did}.json.gz"
        with gzip.open(path,"rb") as f:body=f.read()
        if hashlib.sha256(body).hexdigest()!=contract["source"]["source_v16_dataset_sha256"][did]:
            raise ValueError("original Neotoma full source dataset bytes changed: "+did)
        d=source_datasets(json.loads(body),did,SITES[did])
        meta=inspect_download(d)
        remote,raw,data=query_occurrence(session,occid,contract["test"]["max_response_bytes"],timeout)
        result={"datasetid":did,"siteid":SITES[did],
            "preselected_occid":occid,"original_full_dataset_sha256":hashlib.sha256(body).hexdigest(),
            **meta,**remote,"source_21ka_occupied_plant_validated":False}
        if raw is not None:
            occurrenceids=ids(data,OCC_KEYS)
            if occid not in occurrenceids:
                result["status"]="SOURCE_OCCURRENCE_ID_NOT_PRESENT_IN_EXPLICIT_ID_FIELDS"
            else:
                matches,fields=match_explicit_sample(data,d)
                result.update(fields)
                result["exact_source_sample_match_count"]=len(matches)
                result["same_type_explicit_sample_matches"]=matches[:10]
                result["status"]=(
                    "ONE_SOURCE_SAMPLE_ID_STRUCTURALLY_JOINED_NOT_CALIBRATED"
                    if len(matches)==1 else
                    "MULTIPLE_EXPLICIT_SAMPLE_CANDIDATES_AMBIGUOUS"
                    if len(matches)>1 else
                    "NO_EXPLICIT_SAMPLE_JOIN_NO_FOSSIL_AGE_ACCEPTED")
            raw_path=raw_out/f"occurrence_{occid}.json.gz"
            with gzip.open(raw_path,"wb") as f:f.write(raw)
            result["source_raw_gzip"]=raw_path.name
        results.append(result)
        print(json.dumps({"datasetid":did,"occid":occid,"status":result["status"],
                          "source_sample_ages":meta["source_sample_age_records_n"],
                          "match_n":result.get("exact_source_sample_match_count")},sort_keys=True),flush=True)
    neg,body,data=query_occurrence(session,"9999999",
        contract["test"]["max_response_bytes"],timeout)
    invalid_returns_sample=False
    if body is not None:
        invalid_returns_sample="9999999" in ids(data,OCC_KEYS)
    q={
      "schema":"chocho_six_Alnus_source_occurrence_to_sample_age_link_v17",
      "status":"SOURCE_ONLY_SAMPLE_ID_CHRONOLOGY_JOIN_AUDIT_NOT_OCCUPIED_GLACIAL_HOST",
      "original_six_datasets":6,"original_Alnus_candidate_rows":27,
      "original_unique_candidate_occurrences":26,
      "preselected_source_occurrence_ids":selection,
      "source_status_counts":dict(Counter(r["status"] for r in results)),
      "sample_age_earlier_than_18ka_and_younger_than_24ka_CALIBRATED_VERIFIED":0,
      "species_diagnostic_Alnus_pollen_morphology_verified":0,
      "LGM_alnus_population_occupancy_accepted":0,
      "invalid_occurrenceid_control_yielded_false_sample":invalid_returns_sample,
      "invalid_occurrence_control_http_source_status":neg["status"],
      "raw_source_data_SHA256":{r["preselected_occid"]:r.get("source_sha256")
                               for r in results if r.get("source_sha256")},
      "previous_source_sha256":{
         "original_27_candidate":contract["source"]["original_candidates_sha256"],
         "v16_full_datasets":contract["source"]["source_v16_dataset_sha256"]},
      "source_only_no_genetic_response_or_21ka_occupancy":True,
      "limitations":[
        "Only minimum original occurrence ID per six fixed source dataset was contacted; this pilot is NOT a survey of all 26 candidate records",
        "Explicit sample-ID join is not a calibrated sample chronology or proof species-diagnostic fossil morphology",
        "Sample ages can be radiocarbon years vs calibrated years BP and multiple age models; no posthoc best age model chosen",
        "No response on a source endpoint is UNKNOWN, not absence of paleoecological occurrence",
        "Previous climate-matched-real-other-plant negative glacial signal remains frozen negative"]}
    with (out/"six_site_first_occurrence_sample_age_link_v17.csv").open("w",newline="",encoding="utf8") as f:
        keys=sorted(set().union(*(r.keys() for r in results)))
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
        w.writeheader()
        for r in results:
            w.writerow({k:json.dumps(v,sort_keys=True,ensure_ascii=False) if isinstance(v,(list,dict))
                        else v for k,v in r.items()})
    (out/"summary.json").write_text(json.dumps(q,indent=2,sort_keys=True)+"\n")
    return q

def main():
    a=argparse.ArgumentParser()
    for k in ("contract","candidates","raw-dir","out"):
        a.add_argument("--"+k,required=True,type=Path)
    args=a.parse_args()
    if file_hash(args.candidates)!=CANDIDATE_SHA:
        raise ValueError("original 27 pollen occurrence source digest mismatch")
    contract=json.loads(args.contract.read_text())
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho-Neotoma-six-site-Alnus-pollen-chronology-join/0.1"})
    result=run(session,contract,rows(args.candidates),args.raw_dir,args.out)
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
