#!/usr/bin/env python3
"""Audit original Neotoma dataset chronologies and citations for six Alnus pollen sites.

This is verification of metadata availability, NOT acceptance of sample age,
pollen species identification, or actual ancient host plant population range.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
import requests
SOURCE_CSV_SHA="3d47fe9a31cb6a04ce60aa0f1d7fc8dd7fb1508d3966d6f3a7dc6942ddfa9efc"
SOURCE_SUMMARY_SHA="16bc59d362365416ac9144d29bae05ed503341ebb99d109f750ababdaa401bd9"

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rows(path):
    with Path(path).open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def gate_source(contract_path,records_path,summary_path):
    if digest(records_path)!=SOURCE_CSV_SHA or digest(summary_path)!=SOURCE_SUMMARY_SHA:
        raise ValueError("previous bounded 27 Alnus source rows modified")
    q=json.loads(Path(contract_path).read_text())
    s=json.loads(Path(summary_path).read_text())
    rec=rows(records_path)
    if (q["status"]!="EXPLORATORY_DATASET_METADATA_CHECK_NOT_VALIDATED_21KA_OCCUPANCY"
        or len(rec)!=27 or s["source_candidate_distinct_occids"]!=26):
        raise ValueError("previous fossil source identity/candidate count drift")
    if any(r["source_grade"]!="SOURCE_18_24KA_METADATA_CANDIDATE_NOT_VERIFIED_FOSSIL_OCCUPANCY"
           or r["source_sample_taxonid"]!="1299" or r["source_sample_taxonname"]!="Alnus glutinosa"
           or r["source_datasettype"]!="pollen" for r in rec):
        raise ValueError("original Alnus pollen age, taxon and proxy class changed")
    existing={(int(r["source_datasetid"]),int(r["source_siteid"])) for r in rec}
    chosen={(int(x["datasetid"]),int(x["siteid"])) for x in q["source"]["candidate_dataset_sites"]}
    if len(chosen)!=6 or chosen!=existing:
        raise ValueError("six original candidate Neotoma dataset and source site pairs do not match")
    return q,rec

def fetch_metadata(session,url,max_size):
    try:
        r=session.get(url,timeout=(8,20),headers={"Accept":"application/json"})
        r.raise_for_status()
        if "json" not in r.headers.get("Content-Type","").lower() or len(r.content)>max_size:
            raise ValueError("bad Neotoma JSON size or response type")
        obj=r.json()
        data=obj.get("data") if isinstance(obj,dict) else obj if isinstance(obj,list) else None
        if not isinstance(data,(list,dict)):
            raise ValueError("unexpected Neotoma dataset/chronology/publication schema")
        first=data[0] if isinstance(data,list) and data else data if isinstance(data,dict) else None
        return {"status":"SOURCE_METADATA_AVAILABLE_NOT_INDEPENDENTLY_VALIDATED",
            "bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest(),
            "data_type":type(data).__name__,
            "data_count":len(data) if isinstance(data,list) else None,
            "first_source_record_keys":sorted(first) if isinstance(first,dict) else [],
            "first_record_preview":json.dumps(first,sort_keys=True,ensure_ascii=False)[:5500] if first is not None else "",
            "actual_url":r.url}
    except (requests.RequestException,ValueError) as exc:
        return {"status":"UNAVAILABLE_METADATA_NOT_EVIDENCE_OF_MISSING_CHRONOLOGY",
                "error":type(exc).__name__+":"+str(exc)[:250]}

def run(session,contract):
    x=contract["queries"];queries=[];base=x["official_api"]
    for case in sorted(contract["source"]["candidate_dataset_sites"],key=lambda y:y["datasetid"]):
        for suffix in x["per_dataset"]:
            url=base+suffix.replace("{datasetid}",str(case["datasetid"]))
            m=fetch_metadata(session,url,x["max_body_bytes"])
            queries.append({"datasetid":case["datasetid"],"siteid":case["siteid"],
                "source_type":suffix.rsplit("/",1)[-1] if "{" not in suffix.rsplit("/",1)[-1] else "dataset",
                **m})
    # Negative control: a nonexistent dataset must not be an ordinary nonempty chronology.
    url=base+"/v2.0/data/datasets/"+str(x["negative_dataset_id"])+"/chronologies"
    queries.append({"datasetid":x["negative_dataset_id"],"siteid":None,
                    "source_type":"negative_invalid_dataset_chronologies",
                    **fetch_metadata(session,url,x["max_body_bytes"])})
    success=Counter(q["source_type"] for q in queries
                    if q["status"]=="SOURCE_METADATA_AVAILABLE_NOT_INDEPENDENTLY_VALIDATED")
    stale={}
    for field in ("dataset","chronologies","publications"):
        relevant=[r for r in queries if r["source_type"]==field]
        stale[field]={"queried_datasets":len(relevant),
                     "http_available":sum(r["status"]=="SOURCE_METADATA_AVAILABLE_NOT_INDEPENDENTLY_VALIDATED" for r in relevant),
                     "nonempty_source_metadata":sum((r.get("data_count") or 0)>0 for r in relevant)}
    neg=[r for r in queries if r["source_type"]=="negative_invalid_dataset_chronologies"]
    if neg and neg[0].get("data_count",0):
        status="NEGATIVE_DATASET_ID_CONTROL_NONEMPTY_SOURCE_METADATA_INVALID"
    else:
        status="CHRONOLOGY_AND_PUBLICATION_METADATA_AVAILABILITY_AUDITED_NO_FOSSIL_OCCUPANCY"
    return queries,{"schema":"chocho_neotoma_Alnus_six_18_24ka_pollen_datasets_chronology_metadata_v14",
        "status":status,
        "source_Alnus_18_24ka_pollen_candidate_rows":27,
        "source_Alnus_18_24ka_unique_occurrence_ids":26,
        "candidate_source_sites":6,"candidate_source_datasets":6,
        "source_query_count":len(queries),
        "metadata_coverage":stale,
        "sample_ages_independently_chronology_validated":False,
        "pollen_morphology_species_level_validated":False,
        "original_21ka_host_plant_populations_inferred":False,
        "genetic_response_used":False,
        "limits":["Neotoma metadata being available does not imply the exact sample's age/depth model is independently validated",
          "Dataset citations do not prove Alnus species-diagnostic pollen or lack of remote transported pollen",
          "Missing dataset chronology metadata due HTTP source errors never means ancient fossil absence",
          "No historical host occupancy probabilities or genetic barriers have been inferred"]}
def write(path,entries):
    keys=sorted(set().union(*(v.keys() for v in entries)))
    with Path(path).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader()
        for row in entries:
            w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()})

def main():
    p=argparse.ArgumentParser()
    for k in ("contract","prior-records","prior-summary","out"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    contract,original=gate_source(a.contract,a.prior_records,a.prior_summary)
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho/Neotoma-six-dated-Alnus-pollen-sites-independent-chronology-source-review/0.1"})
    requests_,report=run(session,contract)
    write(a.out/"all_six_dataset_chronology_publication_source_queries.csv",requests_)
    (a.out/"summary.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()
