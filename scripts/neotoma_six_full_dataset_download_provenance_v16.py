#!/usr/bin/env python3
"""Source-only Neotoma full-download age provenance on six pre-fixed Alnus datasets.

EXACT source candidate fossil records and original WCVP taxon ID only. This
does not assert 21ka calibrated sample ages, species-diagnostic pollen,
occupied 21ka larval host plants or insect genetic differentiation.

The occurrence-search endpoint's taxonname filter previously failed an
impossible-name control. Full downloads MUST have explicitly matching source
dataset IDs, and no parent chronology may be assigned to a sample by guessing.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,time
from collections import Counter,defaultdict
from pathlib import Path
import requests

SOURCES={
 "candidate":"3d47fe9a31cb6a04ce60aa0f1d7fc8dd7fb1508d3966d6f3a7dc6942ddfa9efc",
 "ageflag":"914270f10a0874ea1c000dac395ab6589c95485c1db59d63a6aef3903751a518"
}
EXPECTED={
 "4396":"3392","14556":"9834","22662":"14474",
 "22683":"14485","24197":"15750","41351":"23855"}
BASE="https://api.neotomadb.org/v2.0/data/downloads/"
ID_KEYS={"datasetid","dataset_id"}
OCC_KEYS={"occid","occurrenceid","occurrence_id"}
CHRON_KEYS={"chronologyid","chronology_id","chronid"}
SAMPLE_KEYS={"sampleid","sample_id","analysisunitid","analysisunit_id"}

def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for x in iter(lambda:stream.read(1<<20),b""):h.update(x)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def _to_id(value):
    if value is None or isinstance(value,bool):return None
    if isinstance(value,(str,int)) and str(value).strip().isdigit():
        return str(int(str(value).strip()))
    return None

def tree_field_summary(payload,*,max_nodes=250000):
    """Bounded scan for explicit typed IDs, never infer chronology associations."""
    ids=defaultdict(set);paths=Counter()
    stack=[(payload,"$")];visited=0
    observed_sample_links=[]
    while stack:
        obj,path=stack.pop()
        visited+=1
        if visited>max_nodes:raise ValueError("nested Neotoma source exceeds predeclared scan budget")
        if isinstance(obj,dict):
            immediate={}
            for key,value in obj.items():
                field=str(key).lower()
                paths[path+"."+str(key)]+=1
                group=("dataset" if field in ID_KEYS else
                       "occurrence" if field in OCC_KEYS else
                       "chronology" if field in CHRON_KEYS else
                       "sample" if field in SAMPLE_KEYS else None)
                if group:
                    norm=_to_id(value)
                    if norm:
                        ids[group].add(norm);immediate[group]=norm
                if isinstance(value,(dict,list)):
                    stack.append((value,path+"."+str(key)))
            # Raw co-location of ID keys is a STRUCTURAL hint only, NOT proof
            # that a calibrated chronology was used to date the pollen sample.
            if "chronology" in immediate and ("sample" in immediate or "occurrence" in immediate):
                observed_sample_links.append(immediate)
        elif isinstance(obj,list):
            for value in obj:
                if isinstance(value,(list,dict)):
                    stack.append((value,path+"[]"))
    return {"all_explicit_datasetids":sorted(ids["dataset"],key=int),
            "all_explicit_occurrenceids":sorted(ids["occurrence"],key=int),
            "all_explicit_chronologyids":sorted(ids["chronology"],key=int),
            "all_explicit_sampleids":sorted(ids["sample"],key=int),
            "raw_node_count":visited,
            "source_key_paths":[{"path":k,"count":v} for k,v in paths.most_common(75)],
            "colocated_sample_chronology_ids_first20":observed_sample_links[:20],
            "colocated_sample_chronology_count":len(observed_sample_links)}

def query_source(session,api_url,max_size,timeout):
    try:
        with session.get(api_url,timeout=timeout,stream=True,
                         headers={"Accept":"application/json"}) as resp:
            resp.raise_for_status()
            content_type=resp.headers.get("content-type","").lower()
            if "json" not in content_type:raise ValueError("non-JSON Neotoma response")
            length=resp.headers.get("content-length")
            if length and int(length)>max_size:
                raise ValueError("size limit exceeded before download")
            parts=[];seen=0
            for part in resp.iter_content(chunk_size=65536):
                if part:
                    seen+=len(part)
                    if seen>max_size:raise ValueError("Neotoma source byte cap exceeded")
                    parts.append(part)
            body=b"".join(parts)
            parsed=json.loads(body)
            data=parsed.get("data") if isinstance(parsed,dict) else parsed
            if not isinstance(data,(dict,list)):
                raise ValueError("Neotoma JSON missing data object")
            return {"status":"HTTP_JSON_SOURCE_RECEIVED",
                "http_status":resp.status_code,"source_sha256":hashlib.sha256(body).hexdigest(),
                "source_bytes":len(body),"returned_url":resp.url,
                "data_payload_type":type(data).__name__,
                "data_top_level_keys":sorted(data)[:40] if isinstance(data,dict)
                 else sorted(data[0])[:40] if data and isinstance(data[0],dict) else [],
                "data_items_n":len(data) if isinstance(data,list) else None},body,data
    except (requests.RequestException,ValueError,json.JSONDecodeError) as exc:
        return {"status":"SOURCE_UNAVAILABLE_OR_INVALID_NOT_FOSSIL_ABSENCE",
           "error":type(exc).__name__+": "+str(exc)[:350],
           "requested_url":api_url},None,None

def source_gate(contract,candidates,ageflags):
    if (len(candidates)!=27 or len(ageflags)!=27
        or contract["status"]!="SOURCE_ONLY_DETERMINISTIC_FROZEN_BEFORE_FULL_DATASET_DOWNLOADS"):
        raise ValueError("sample-date source cohort or source preflight contract changed")
    c=contract["source"]
    pairs={(r["source_datasetid"],r["source_siteid"]) for r in candidates}
    if pairs!=set(EXPECTED.items()) or sorted(c["frozen_datasetids"])!=sorted(map(int,EXPECTED)):
        raise ValueError("original fixed Neotoma six sites were changed")
    occ=[r["occid"] for r in candidates]
    if len(set(occ))!=26 or any(
       r["source_sample_taxonid"]!="1299" or
       r["source_sample_taxonname"]!="Alnus glutinosa" or
       r["source_datasettype"]!="pollen"
       for r in candidates):
        raise ValueError("source taxon identity/proxy/26 occurrence IDs drift")
    if {r["occid"] for r in ageflags}!={r["occid"] for r in candidates}:
        raise ValueError("unreconciled original 27 candidate source ages")
    return {did:{r["occid"] for r in candidates if r["source_datasetid"]==did}
            for did in EXPECTED}

def run(session,contract,candidates,ageflags,out_dir):
    chosen=source_gate(contract,candidates,ageflags)
    cap=contract["source_diagnostic"]["max_response_bytes_per_dataset"]
    timeout=(contract["source_diagnostic"]["timeout_connect_seconds"],
             contract["source_diagnostic"]["timeout_read_seconds"])
    out_dir.mkdir(parents=True,exist_ok=True)
    results=[]
    rawdir=out_dir/"source_original_downloads";rawdir.mkdir(parents=True,exist_ok=True)
    for did in sorted(chosen,key=int):
        response,body,data=query_source(session,BASE+did,cap,timeout)
        row={"datasetid":did,"expected_siteid":EXPECTED[did],
             "original_Alnus_candidate_occurrence_count":len(chosen[did]),**response}
        if body is not None:
            scanned=tree_field_summary(data)
            if did not in scanned["all_explicit_datasetids"]:
                row["status"]="UNVERIFIED_DATASET_ID_IN_DOWNLOAD_BODY"
            else:
                row["status"]="SOURCE_DATASET_ID_VERIFIED_CHRONOLOGY_LINK_NOT_YET_VERIFIED"
            row["original_candidate_occurrenceid_found_as_explicit_field_n"]=len(
                set(scanned["all_explicit_occurrenceids"])&chosen[did])
            row["candidate_occurrenceid_missing_from_explicit_download_fields_n"]=len(chosen[did])-row[
                "original_candidate_occurrenceid_found_as_explicit_field_n"]
            row["raw_record_scan_n"]=scanned["raw_node_count"]
            row["chronology_ID_count"]=len(scanned["all_explicit_chronologyids"])
            row["sample_ID_count"]=len(scanned["all_explicit_sampleids"])
            row["dataset_id_count"]=len(scanned["all_explicit_datasetids"])
            row["colocated_sample_chronology_count_not_age_validation"]=scanned["colocated_sample_chronology_count"]
            row["nested_source_structure"]=scanned["source_key_paths"]
            gzfile=rawdir/f"dataset_{did}.json.gz"
            with gzip.open(gzfile,"wb",compresslevel=6) as f:f.write(body)
            row["source_gzip_file"]=gzfile.name
            row["compressed_sha256"]=file_sha(gzfile)
        results.append(row)
        print(json.dumps({"datasetid":did,"status":row["status"],
           "source_bytes":row.get("source_bytes"),"candidate_explicit_occurrence_count":row.get(
               "original_candidate_occurrenceid_found_as_explicit_field_n")},sort_keys=True),flush=True)
    # The invalid dataset is a test of source endpoint behavior, never a real record.
    badid=str(contract["source_diagnostic"]["negative_control_dataset_id"])
    neg,body,data=query_source(session,BASE+badid,cap,timeout)
    has_fake_dataset=False
    if body is not None:
        links=tree_field_summary(data)
        has_fake_dataset=badid in links["all_explicit_datasetids"]
    summary={
        "schema":"chocho_6_Alnus_Neotoma_full_download_source_schema_v16",
        "status":"FULL_DATASET_SCHEMA_ACQUISITION_UNRESOLVED_SAMPLE_CALIBRATION_AND_TAXONOMY",
        "original_candidate_rows":27,"original_candidate_occurrence_ids":26,
        "source_taxon":"Alnus glutinosa","source_taxonid":1299,
        "source_datasets_n":6,
        "source_success_frozen_dataset_id_n":sum(r["status"]=="SOURCE_DATASET_ID_VERIFIED_CHRONOLOGY_LINK_NOT_YET_VERIFIED" for r in results),
        "source_invalid_dataset_false_positive":has_fake_dataset,
        "invalid_dataset_request":neg,
        "source_query_records":len(results)+1,
        "source_response_availability":dict(Counter(x["status"] for x in results)),
        "source_original_data_sha256":{r["datasetid"]:r.get("source_sha256")
                                       for r in results if r.get("source_sha256")},
        "original_source_used":{"candidate_sha256":SOURCES["candidate"],"age_QC_sha256":SOURCES["ageflag"]},
        "sample_assigned_age_depth_models_verified":0,
        "independent_species_diagnostic_pollen_morphology_verified":0,
        "Alnus_21ka_fossil_site_presence_accepted":0,
        "historical_host_occurrence_or_genetic_response_used":False,
        "limits":[
            "A full Neotoma download may include multiple age chronologies and sample ages; field presence alone is not an exact assigned sample age and not a calibrated age posterior",
            "Even an explicit source occurrence ID is not independent morphology confirmation of Alnus glutinosa in source pollen",
            "Failed or missing full Neotoma datasets do NOT mean a species or any sample was absent in 21ka",
            "Some original source numeric ages were outside coarse published dataset envelope; do not resolve by selecting whichever chronology best fits",
            "No DNA response, LGM occupied resource corridors or observed insect genetic barriers are examined"]}
    with (out_dir/"six_dataset_full_source_schema_v16.csv").open("w",newline="",encoding="utf8") as f:
        fields=sorted(set().union(*(r.keys() for r in results)))
        writer=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        writer.writeheader()
        for row in results:
            writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    (out_dir/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("contract","candidate","ageflags","out"):
        p.add_argument("--"+name,required=True,type=Path)
    args=p.parse_args()
    for k,path in (("candidate",args.candidate),("ageflag",args.ageflags)):
        if file_sha(path)!=SOURCES[k]:
            raise ValueError("original source audit SHA changed: "+k)
    contract=json.loads(args.contract.read_text())
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho-Neotoma-6-datasets-source-chronology-audit/0.1"})
    report=run(session,contract,read_csv(args.candidate),read_csv(args.ageflags),args.out)
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()
