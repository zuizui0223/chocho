#!/usr/bin/env python3
"""Read explicit Neotoma chronology objects for the 21 frozen Alnus age candidates.

This checks SOURCE chronologyid identity and age scale only. It does not choose
a favorable age model, recalibrate radiocarbon measurements, authenticate
species-level pollen morphology, infer occupied glacial hosts or read genetics.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
import requests

LABEL="calibrated radiocarbon years bp"
STATUS="SOURCE_CHRONOLOGY_ID_LINK_QC_NOT_VERIFIED_CALIBRATED_HOST_OCCUPANCY"

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def records(path):
    with Path(path).open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def norm_id(x):
    if x is None or isinstance(x,bool):return None
    x=str(x).strip()
    if not x.isdecimal():return None
    return str(int(x))

def norm_label(s):
    return " ".join(str(s or "").lower().strip().split())

def source_candidates(rows,contract):
    if len(rows)!=788:raise ValueError("original all 788 Alnus-associated age records changed")
    cases=[]
    for row in rows:
        candidate=row["source_labelled_18_24ka_cal_BP_candidate"]=="True"
        if candidate:
            if (row["datum_positive"]!="True" or
                norm_label(row["source_agetype_label"])!=LABEL or
                row["source_age_within_18_24ka"]!="True" or
                not norm_id(row["source_chronologyid"])):
                raise ValueError("original source calibrated-label candidate source identity drift")
            cases.append(row)
    if len(cases)!=21 or len({(v["datasetid"],v["sampleid"]) for v in cases})!=21:
        raise ValueError("original frozen21 distinct source samples changed")
    expected=contract["sources"]["required_source_21_candidate_chron_ids"]
    grouping=defaultdict(set)
    for r in cases:grouping[r["datasetid"]].add(norm_id(r["source_chronologyid"]))
    if {k:sorted(v,key=int) for k,v in grouping.items()}!=expected:
        raise ValueError("original five dataset source assigned chronology IDs drift")
    return sorted(cases,key=lambda r:(int(r["datasetid"]),int(r["sampleid"])))

def parse_chronologies(data,requested_dataset):
    if not isinstance(data,list):raise ValueError("Neotoma chronological data array missing")
    out={}
    for r in data:
        if not isinstance(r,dict):raise ValueError("non-dictionary original chronology")
        chron=r.get("chronology",r)
        if not isinstance(chron,dict):raise ValueError("no typed chronology")
        cid=norm_id(chron.get("chronologyid"))
        if cid is None:raise ValueError("official chronology lacks numeric chronologyid")
        did=norm_id(chron.get("datasetid"))
        if did is not None and did!=requested_dataset:
            raise ValueError("official chronology explicitly belongs to unrelated dataset")
        if cid in out:raise ValueError("duplicate chronologyid in one source dataset response")
        out[cid]=chron
    return out

def source_chron_meta(chron):
    controls=chron.get("controls")
    if controls is None:controls=[]
    if not isinstance(controls,list):raise ValueError("Neotoma chronology controls not a list")
    kinds=Counter()
    geochron=Counter()
    source_controls_with_bounds=0
    for raw in controls:
        if not isinstance(raw,dict):continue
        kinds[str(raw.get("controltype") or "unspecified")]+=1
        geochron[str(raw.get("geochronagetype") or "unspecified")]+=1
        bounds=raw.get("controllimits")
        if isinstance(bounds,dict) and bounds.get("ageolder") is not None and bounds.get("ageyounger") is not None:
            source_controls_with_bounds+=1
    return {
      "chronology_source_agetype":str(chron.get("agetype") or ""),
      "chronology_source_name":str(chron.get("chronologyName") or ""),
      "chronology_source_notes":str(chron.get("chronologynotes") or "")[:750],
      "source_chronology_controls_total":len(controls),
      "source_chronology_controls_with_numerical_bounds":source_controls_with_bounds,
      "source_controls_types":dict(kinds),
      "source_controls_geochronage_types":dict(geochron),
      "original_model_independent_calibration_not_verified":True}

def fetch_chronologies(session,url,cap,timeout):
    try:
        with session.get(url,timeout=timeout,stream=True,headers={"Accept":"application/json"}) as resp:
            resp.raise_for_status()
            if "json" not in resp.headers.get("content-type","").lower():
                raise ValueError("unexpected non JSON chronology source")
            payload=b""
            for chunk in resp.iter_content(chunk_size=32768):
                payload+=chunk
                if len(payload)>cap:raise ValueError("source chronology response exceeded cap")
            j=json.loads(payload)
            data=j.get("data") if isinstance(j,dict) else j
            if not isinstance(data,list):raise ValueError("source chronology data not a list")
            return {"status":"SOURCE_JSON_AVAILABLE",
                    "bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest(),
                    "url":resp.url,"source_chronology_records":len(data)},payload,data
    except (requests.RequestException,ValueError) as exc:
        return {"status":"SOURCE_UNAVAILABLE_NOT_FOSSIL_ABSENCE",
                "error":type(exc).__name__+": "+str(exc)[:280]},None,None

def run(session,contract,ages,out):
    selected=source_candidates(ages,contract)
    settings=contract["acquisition"]
    dids=contract["sources"]["source_dataset_site_map"]
    source_ids=sorted(dids,key=int)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/"official_source_chronologies"
    raw.mkdir(parents=True,exist_ok=True)
    result_by_dataset={};source_receipts=[]
    for did in source_ids:
        url=settings["official_url_template"].replace("{datasetid}",did)
        status,bytes_,data=fetch_chronologies(session,url,settings["max_response_bytes"],
            (settings["connection_timeout_s"],settings["read_timeout_s"]))
        info={"datasetid":did,"siteid":dids[did],**status}
        if bytes_ is not None:
            source=raw/f"chronologies_dataset_{did}.json.gz"
            with gzip.open(source,"wb",compresslevel=6) as f:f.write(bytes_)
            info["archived_source_gzip_file"]=source.name
            info["archived_source_gzip_sha256"]=sha(source)
            chronologies=parse_chronologies(data,did)
            info["explicit_chronologyids"]=sorted(chronologies,key=int)
            result_by_dataset[did]=chronologies
        source_receipts.append(info)
        print(json.dumps({"datasetid":did,"status":info["status"],
              "chronologyids":info.get("explicit_chronologyids",[])},sort_keys=True),flush=True)
    checked=[];admissible=0
    for age in selected:
        did=age["datasetid"];cid=norm_id(age["source_chronologyid"])
        source=result_by_dataset.get(did)
        entry={"datasetid":did,"siteid":age["siteid"],"sampleid":age["sampleid"],
          "source_sample_age_center_BP":age["source_age_BP_numeric"],
          "source_age_type":age["source_agetype_label"],
          "original_explicit_sample_chronologyid":cid,
          "source_age_bounds_both_available":age["has_both_source_age_bounds"]=="True",
          "source_age_younger_BP":age["source_ageyounger_BP"],
          "source_age_older_BP":age["source_ageolder_BP"],
          "independently_accepted_paleo_plant_occupancy":False,
          "independently_verified_calibrated_age_depth_model":False,
          "pollen_species_morphology_independently_confirmed":False}
        if source is None:
            entry["status"]="CHRONOLOGY_SOURCE_UNAVAILABLE"
        elif cid not in source:
            entry["status"]="EXPLICIT_SAMPLE_AGE_CHRONOLOGYID_NOT_IN_DATASET_CHRONOLOGY"
        else:
            detail=source_chron_meta(source[cid])
            entry.update(detail)
            identical=norm_label(entry["source_age_type"])==norm_label(detail["chronology_source_agetype"])
            entry["sample_vs_chronology_source_agetype_match"]=identical
            entry["status"]=("EXACT_CHRONOLOGYID_MATCH_SOURCE_LABEL_CONSISTENT_NOT_CALIBRATED"
                    if identical else "EXACT_CHRONOLOGYID_MATCH_BUT_SOURCE_AGE_SCALE_CONFLICT")
        checked.append(entry)
    flags=Counter(r["status"] for r in checked)
    if len(checked)!=21:raise ValueError("source 21 sample chronology candidates lost")
    src=Counter(r["status"] for r in source_receipts)
    summary={"schema":"chocho_Alnus_six_21_source_sample_calendar_age_chronology_qc_v19",
       "status":STATUS,
       "source_frozen_candidate_samples":21,
       "source_datasets":6,"source_dataset_query_status_counts":dict(src),
       "source_original_21_candidate_chronology_lookup_status_counts":dict(flags),
       "original_dataset_source_SHA":contract["sources"]["original_v18_all_source_age_models_csv_sha256"],
       "source_chronology_response_SHA":{r["datasetid"]:r.get("sha256") for r in source_receipts if r.get("sha256")},
       "sample_age_models_independently_calibrated_n":0,
       "species_diagnostic_pollen_morphology_independently_verified_n":0,
       "source_glacial_occupied_Alnus_populations_accepted":0,
       "insect_genetic_responses_read":False,
       "limits":[
          "Exact Neotoma chronologyid matching verifies only relational database provenance, not independently reviewed radiocarbon calibration or posterior age-depth uncertainty.",
          "A source model agtype label can disagree with a sample age label; do not select a chronology after observing source fossil ages.",
          "Positive Neotoma pollen taxonid1299 datums are not verified species-diagnostic Alnus glutinosa morphology or proof of a local occupied host-plant population.",
          "Insect sampling coordinates in climate graph are modern anchors; no ancient genetics or host-specific LGM effect proven.",
          "Prior negative actual-host vs climate-matched-other-real-plant controls remain frozen and unretuned."
       ]}
    def write_csv(file,items):
        keys=sorted(set().union(*(r.keys() for r in items)))
        with (out/file).open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
            w.writeheader()
            for row in items:
                w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v
                            for k,v in row.items()})
    write_csv("all_21_candidate_direct_sample_chronology_source_crosschecks_v19.csv",checked)
    write_csv("six_source_dataset_chronology_full_fetch_status_v19.csv",source_receipts)
    (out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",required=True,type=Path)
    p.add_argument("--original-sample-age-csv",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    contract=json.loads(args.contract.read_text())
    if sha(args.original_sample_age_csv)!=contract["sources"]["original_v18_all_source_age_models_csv_sha256"]:
        raise ValueError("source v18 direct Alnus datum sample age record changed")
    s=requests.Session()
    s.headers.update({"User-Agent":"chocho-six-Neotoma-Alnus-21-source-sample-chronologyid/0.1"})
    run(s,contract,records(args.original_sample_age_csv),args.out)
if __name__=="__main__":main()
