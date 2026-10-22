#!/usr/bin/env python3
"""Inspect sample-level age fields on exact Neotoma taxon-ID occurrence API routes.

The Neotoma generic occurrence taxonname filter was confirmed invalid by an
impossible-name sentinel. This script uses separately verified 1299/3461
taxon IDs and never reports actual 18–24ka host occupation or fossil absence.
The raw JSON response is limited (20/page) and not a comprehensive dataset.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import requests

EXPECTED_V11="55887306ce8defa7c027fc1d518f2bb798ac2868b6c7f46db288f154e8ea8219"
ROOT="https://api.neotomadb.org"

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source_http(session,url,params):
    try:
        r=session.get(url,params=params,timeout=(8,18),headers={"Accept":"application/json"})
        r.raise_for_status()
        if "json" not in r.headers.get("Content-Type","").lower():
            raise ValueError("Neotoma did not return JSON")
        raw=r.content
        if len(raw)>400_000:raise ValueError("bounded metadata response unexpectedly too large")
        x=r.json()
        data=x.get("data") if isinstance(x,dict) else x if isinstance(x,list) else None
        if not isinstance(data,list):raise ValueError("Neotoma metadata has no data list")
        if any(not isinstance(z,dict) for z in data):
            raise ValueError("Neotoma occurrence metadata is not object records")
        return data,{"status":"SOURCE_HTTP200_NOT_FOSSIL_VALIDATION",
                     "sha256":hashlib.sha256(raw).hexdigest(),
                     "actual_url":r.url,"record_count":len(data),"bytes":len(raw),
                     "record_ids":[str(z.get("occid") or z.get("occurrenceid") or "") for z in data]}
    except (requests.RequestException,ValueError) as e:
        return [],{"status":"UNAVAILABLE_OR_INVALID_SOURCE_DO_NOT_INFER_ABSENCE",
                   "error":type(e).__name__+":"+str(e)[:240]}

def value_shape(value,depth=2):
    if value is None:return {"kind":"null"}
    if isinstance(value,list):
        return {"kind":"list","n":len(value),"item_shape":value_shape(value[0],depth-1) if value and depth else None}
    if isinstance(value,dict):
        return {"kind":"object","keys":sorted(value),
                "subshapes":{k:value_shape(v,depth-1) for k,v in value.items()} if depth else None}
    return {"kind":type(value).__name__}

def fossil_age_site_preview(record):
    age=record.get("age")
    sample=record.get("sample")
    site=record.get("site")
    # Public, source metadata only; no WCVP taxon label, no sample-age calibration assertion.
    return {"occurrence_id":str(record.get("occid") or record.get("occurrenceid") or ""),
        "age_shape":value_shape(age),
        "sample_shape":value_shape(sample,1),
        "site_shape":value_shape(site,1),
        "age_metadata_uninterpreted":json.dumps(age,ensure_ascii=False,sort_keys=True)[:2500],
        "sample_id":sample.get("sampleid") if isinstance(sample,dict) else None,
        "site_id":site.get("siteid") if isinstance(site,dict) else None}

def frozen_contract(path,v11_summary):
    q=json.loads(Path(path).read_text())
    x=json.loads(Path(v11_summary).read_text())
    if (digest(v11_summary)!=EXPECTED_V11 or
        x["status"]!="FOSSIL_TAXONNAME_FILTER_INVALID_SENTINEL_MATCHES_REAL_TAXON" or
        x["exact_taxon_registry_ids"]!={"Alnus glutinosa":1299,"Betula pendula":3461}):
        raise ValueError("original source-bound taxon-ID registry and sentinel v11 source drift")
    if (q["status"]!="SOURCE_SAMPLE_LEVEL_CHRONOLOGY_PRECHECK_NOT_LGM_PRESENCE" or
        q["age_target_BP"]!=[18000,24000] or
        q["source_pilot_limits"]["offsets"]!=[0,50] or
        q["source_pilot_limits"]["items_per_page"]!=20):
        raise ValueError("pre-source age/taxon-ID pilot contract changed")
    expected={r["accepted"]:r["neotoma_taxonid"] for r in q["exact_taxon_records"]}
    if expected!=x["exact_taxon_registry_ids"]:
        raise ValueError("model taxon IDs disagree with verified source registry")
    return q

def run_probe(session,q):
    outputs=[];preview=[]
    ages={}
    ids_by_species={}
    for taxon in q["exact_taxon_records"]:
        name=taxon["accepted"];taxonid=taxon["neotoma_taxonid"]
        ids_by_species[name]=set()
        for offset in q["source_pilot_limits"]["offsets"]:
            url=f"{ROOT}/v2.0/data/taxa/{taxonid}/occurrences"
            params={"limit":20,"offset":offset}
            data,meta=source_http(session,url,params)
            outputs.append({"probe":"exact_taxon_id_route","taxon":name,
                            "taxonid":taxonid,"offset":offset,**meta})
            if data:
                for x in data:
                    v=fossil_age_site_preview(x)
                    if v["occurrence_id"] in ids_by_species[name]:
                        # Pagination can overlap without implying independent fossil specimens.
                        v["duplicate_record_across_pages"]=True
                    ids_by_species[name].add(v["occurrence_id"])
                    preview.append({"requested_taxon":name,"requested_taxonid":taxonid,
                                    "offset":offset,**v})
                ages[name]=value_shape(data[0].get("age"))
    for example in q["generic_endpoint_negative_controls"]:
        params=dict(part.split("=",1) for part in example["query"].split("&"))
        data,meta=source_http(session,ROOT+"/v2.0/data/occurrences",params)
        outputs.append({"probe":"generic_occurrence_filter_negative_control",
                        "query":example["query"],**meta})
    generic=[x for x in outputs if x["probe"]=="generic_occurrence_filter_negative_control"
             and x["status"]=="SOURCE_HTTP200_NOT_FOSSIL_VALIDATION"]
    duplicate_generic=[]
    for i,a in enumerate(generic):
        for b in generic[i+1:]:
            if a["sha256"]==b["sha256"] and a["record_count"]>0 and b["record_count"]>0:
                duplicate_generic.append([a["query"],b["query"]])
    direct=[x for x in outputs if x["probe"]=="exact_taxon_id_route"
            and x["status"]=="SOURCE_HTTP200_NOT_FOSSIL_VALIDATION"]
    positive={r["taxon"] for r in direct if r["record_count"]>0}
    if len(positive)==2:
        first=next((r for r in direct if r["taxon"]=="Alnus glutinosa" and r["record_count"]>0),None)
        second=next((r for r in direct if r["taxon"]=="Betula pendula" and r["record_count"]>0),None)
        if first and second and first["record_ids"]==second["record_ids"]:
            status="TAXON_ID_OCCURRENCE_PATH_NOT_DISTINCT_FILTER_UNTRUSTWORTHY"
        else:status="TAXON_ID_PATH_DIFFERENT_OCCURRENCE_RECORDS_CHRONOLOGY_UNVERIFIED"
    else:status="TAXON_ID_PATH_INCOMPLETE_NO_AGE_OR_ABSENCE_INFERENCE"
    summary={"schema":"chocho_neotoma_exact_taxonID_sample_age_metadata_v12",
        "status":status,
        "official_source":"https://api.neotomadb.org/swagger.json",
        "source_only_target_taxon_ids":{"Alnus glutinosa":1299,"Betula pendula":3461},
        "source_taxa_with_nonempty_direct_path":sorted(positive),
        "source_direct_success_response_pages":len(direct),
        "source_direct_occurrence_ids_covered_by_bounded_pages":
           {k:len(v) for k,v in ids_by_species.items()},
        "source_direct_age_data_shapes":ages,
        "generic_negative_controls_identical_payload_pairs":duplicate_generic,
        "source_any_original_species_18_24ka_fossil_occupancy_verified":False,
        "source_sample_level_chronology_verified":False,
        "source_pollen_species_vs_genus_identity_verified":False,
        "source_pagination_exhausted":False,
        "source_missing_response_does_not_mean_species_absence":True,
        "original_TTF_host_climate_negative_control_unchanged":True,
        "limitations":["Source-only limited direct-ID occurrence metadata pages are not a complete fossil survey",
            "No sample-level chronometric calibration/age-depth provenance or species-level pollen taxonomic proof",
            "The global 416 WCVP native plants have not been matched to validated 21ka fossil geographic observations",
            "No insect ancient population or genetic response has been used"]}
    return outputs,preview,summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--v11-summary",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    q=frozen_contract(a.contract,a.v11_summary)
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho/source-only-Neotoma-exact-taxonID-chronology-validation/0.1"})
    records,samples,summary=run_probe(session,q)
    def dump_csv(file,rows):
        if not rows:
            (a.out/file).write_text("source_status\\nNO_DATA_RETURNED_DO_NOT_INFER_ABSENCE\\n")
            return
        keys=sorted(set().union(*(v.keys() for v in rows)))
        with (a.out/file).open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow({k:json.dumps(v,ensure_ascii=False,sort_keys=True) if isinstance(v,(dict,list)) else v
                        for k,v in row.items()})
    dump_csv("taxonid_route_and_generic_controls.csv",records)
    dump_csv("sample_age_site_raw_metadata_shapes.csv",samples)
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()
