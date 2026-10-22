#!/usr/bin/env python3
"""Neotoma v2 source-access preflight for independently dated plant fossil evidence.

This script DOES NOT output LGM plant occupied ranges, confirmed 21ka fossil
occurrences, genus=>species pollen imputation, or any insect genetic response.
Even an API 200 response / apparent 18-24ka query hit is ONLY source metadata
until calibrated per-sample ages, proxy and taxonomic precision are checked.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from pathlib import Path
import requests

ORIGINAL_MODELS_SHA="153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733"
DOMAIN="https://api.neotomadb.org"
PROXIES=("plant macrofossil","pollen")

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def species_of_source(path):
    with Path(path).open(newline="",encoding="utf8") as f:
        q=list(csv.DictReader(f))
    assert len(q)==518
    r={}
    for x in q:
        key=str(int(x["wcvp_id"]))
        if key in r:raise ValueError("duplicate source WCVP ID")
        r[key]=x
    return r

def preflight_spec(openapi,contract):
    if openapi.get("openapi","").split(".")[0]!="3":
        raise ValueError("Neotoma source OpenAPI v3 schema not returned")
    endpoints=contract["Neotoma"]["paths"]
    paths=openapi.get("paths",{})
    for endpoint in endpoints:
        if endpoint not in paths or "get" not in paths[endpoint]:
            raise ValueError("Neotoma API path not found in published schema: "+endpoint)
    components=openapi.get("components",{}).get("parameters",{})
    occurrence_refs=[x.get("$ref","").split("/")[-1] for x in
            paths["/v2.0/data/occurrences"]["get"]["parameters"]]
    resolved={components[k].get("name") for k in occurrence_refs if k in components}
    need={"taxonname","ageold","ageyoung","datasettype","limit"}
    if not need<=resolved:
        raise ValueError(f"v2 occurrence query parameters unexpectedly changed: {sorted(need-resolved)}")
    if (contract["Neotoma"]["per_query_limit"]!=5
        or contract["Neotoma"]["ageyoung_calibrated_BP"]!=18000
        or contract["Neotoma"]["ageold_calibrated_BP"]!=24000):
        raise ValueError("frozen fossil source age/limit contract changed")
    return {"source_paths_checked":endpoints,"occurrence_filter_names":sorted(resolved)}

def describe_json_shape(payload):
    """Do not infer matched fossil taxa, ages or occupancy from a result shape."""
    if isinstance(payload,list):
        return {"payload_type":"list","top_level_elements":len(payload),"keys":sorted(payload[0]) if payload and isinstance(payload[0],dict) else []}
    if isinstance(payload,dict):
        out={"payload_type":"dict","top_level_keys":sorted(payload)}
        for key in ("data","result","results","taxa","occurrences"):
            val=payload.get(key)
            if isinstance(val,list):
                out[key+"_list_length"]=len(val)
                out[key+"_first_keys"]=sorted(val[0]) if val and isinstance(val[0],dict) else []
            elif isinstance(val,dict):
                out[key+"_dict_keys"]=sorted(val)
        return out
    return {"payload_type":type(payload).__name__}

def query_metadata(session,url,params):
    last=None
    for retry in range(2):
        try:
            resp=session.get(url,params=params,timeout=(10,25),headers={"Accept":"application/json"})
            resp.raise_for_status()
            if "json" not in resp.headers.get("Content-Type","").lower():
                raise ValueError("expected Neotoma JSON content type")
            obj=resp.json()
            if isinstance(obj,dict) and obj.get("success") is False:
                raise ValueError("Neotoma explicitly returned success=false")
            return {"state":"SOURCE_RESPONSE_UNREVIEWED",
                "http_status":resp.status_code,"url":resp.url,
                "response_sha256":hashlib.sha256(resp.content).hexdigest(),
                "response_bytes":len(resp.content),
                "response_structure":describe_json_shape(obj)}
        except (requests.RequestException,ValueError) as e:
            last=f"{type(e).__name__}: {str(e)[:180]}"
            if retry==0:time.sleep(.6)
    return {"state":"API_QUERY_FAILED_DO_NOT_INFER_FOSSIL_ABSENCE",
            "http_status":None,"url":url,"reason":last}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",type=Path,required=True)
    p.add_argument("--source-models",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    spec=json.loads(a.contract.read_text())
    if spec["status"]!="SOURCE_ONLY_EXPLORATORY_FEASIBILITY_NO_ANCIENT_OCCUPANCY_ESTIMATES":
        raise ValueError("Neotoma pre-outcome pilot contract changed")
    if sha(a.source_models)!=ORIGINAL_MODELS_SHA:
        raise ValueError("WCVP original 518 name/model source hash changed")
    names=species_of_source(a.source_models)
    selected=spec["selected_source_taxa"]
    if len(selected)!=15 or len({z["wcvp_id"] for z in selected})!=15:
        raise ValueError("frozen 15 pilot name selection drift")
    for x in selected:
        original=names.get(x["wcvp_id"])
        if (original is None or original["accepted_host_name"]!=x["name"]
            or original["model_status"]!="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"):
            raise ValueError("pilot is no longer inside original 416 modeled source")
    api_root=spec["Neotoma"]["api_root"].rstrip("/")
    if api_root!=DOMAIN:raise ValueError("source changed to an unreviewed fossil service")
    s=requests.Session()
    s.headers.update({"User-Agent":"chocho/source-only-glacial-plant-fossil-audit/0.1 (academic data feasibility)"})
    try:
        v=query_metadata(s,spec["Neotoma"]["openapi"],{})
        if v["state"]!="SOURCE_RESPONSE_UNREVIEWED":
            raise ValueError("official Neotoma swagger not fetched: "+str(v.get("reason")))
        swagger=s.get(spec["Neotoma"]["openapi"],timeout=(10,25)).json()
        api_shape=preflight_spec(swagger,spec)
        api_shape["openapi_response_sha256"]=v["response_sha256"]
    except (requests.RequestException,ValueError) as e:
        api_shape={"status":"API_SCHEMA_UNAVAILABLE_FAIL_CLOSED",
            "reason":type(e).__name__+": "+str(e)[:220]}
        (a.out/"source_schema.json").write_text(json.dumps(api_shape,indent=2,sort_keys=True)+"\n")
        raise SystemExit("Official Neotoma source schema cannot be verified; no fossil absence inference")
    (a.out/"source_schema.json").write_text(json.dumps(api_shape,indent=2,sort_keys=True)+"\n")
    outputs=[]
    for seq,x in enumerate(selected):
        for proxy in PROXIES:
            params={"taxonname":x["name"],"datasettype":proxy,
                    "ageyoung":spec["Neotoma"]["ageyoung_calibrated_BP"],
                    "ageold":spec["Neotoma"]["ageold_calibrated_BP"],
                    "limit":spec["Neotoma"]["per_query_limit"],"offset":0}
            query=query_metadata(s,api_root+"/v2.0/data/occurrences",params)
            outputs.append({"wcvp_id":x["wcvp_id"],"name":x["name"],
                "selection_role":x["design"],"proxy":proxy,
                "ageold_query_BP":params["ageold"],"ageyoung_query_BP":params["ageyoung"],
                "requested_max_return":5,"source_query_is_exhaustive":False,
                "status":query["state"],"http_status":query.get("http_status"),
                "endpoint":query["url"],"json_bytes":query.get("response_bytes"),
                "json_sha256":query.get("response_sha256"),
                "structure":query.get("response_structure"),
                "failure":query.get("reason")})
            if (seq+1)%5==0:print(json.dumps({"taxa_attempted":seq+1,"query_responses":len(outputs)}),flush=True)
            time.sleep(.15)
    report={"schema":"chocho_Neotoma_15_source_plant_fossil_API_PILOT_v01",
        "status":"SOURCE_ACCESS_QUERY_RESPONSES_ONLY_NO_21KA_SPECIES_OCCUPANCY_CLAIMS",
        "original_modelled_plants":416,"attempted_WCVP_taxa":15,
        "requested_proxy_queries":len(outputs),
        "Neotoma_query_age_BP":[18000,24000],
        "successful_JSON_proxy_queries":sum(x["status"]=="SOURCE_RESPONSE_UNREVIEWED" for x in outputs),
        "unavailable_proxy_queries":sum(x["status"]!="SOURCE_RESPONSE_UNREVIEWED" for x in outputs),
        "technical_positive_control_n":sum(x["design"]=="API-positive-technical-control" for x in selected),
        "source_models_sha256":ORIGINAL_MODELS_SHA,
        "Neotoma_openapi_sha256":api_shape["openapi_response_sha256"],
        "fossil_species_occurrence_support_verified_n":None,
        "ancient_plant_geography_estimated":False,"no_genetic_response_opened":True,
        "sample_chronology_or_taxon_accuracy_not_yet_verified":True,
        "limitations":["HTTP 200 and a source taxon-name query do not establish exact species-level fossil identification or 21ka sample age",
        "Age-range query can filter dataset ages rather than actual calibrated sample posterior and is not a species occupancy datum",
        "Pollen often identified only at genus level and MUST NOT be assigned to source binomial species",
        "A maximum 5 records per query is a source API feasibility test, not complete retrieval or an absence survey",
        "Twelve WCVP-ID spaced taxa and three deliberately selected temperate tree technical controls are not a representative fossil coverage survey",
        "Previous negative TTF/chocho host-identity glacial comparisons are unchanged and not reopened"]}
    (a.out/"source_pilot_summary.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    with (a.out/"all_15_taxa_source_query_status.csv").open("w",newline="",encoding="utf8") as f:
        fields=["wcvp_id","name","selection_role","proxy","ageold_query_BP","ageyoung_query_BP",
           "requested_max_return","source_query_is_exhaustive","status","http_status",
           "endpoint","json_bytes","json_sha256","structure","failure"]
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader()
        for r in outputs:
            r=dict(r,structure=json.dumps(r["structure"],sort_keys=True) if r["structure"] else "")
            w.writerow(r)
    print(json.dumps(report,sort_keys=True))
if __name__=="__main__":main()
