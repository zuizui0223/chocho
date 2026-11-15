#!/usr/bin/env python3
"""Source-verified 2018–2022 individual GBIF occurrences for three uncertain cases.

The exact three candidate cases were fixed by the 1700–2010 time-consistent
port-host/sea model, before inspecting raw occurrence identities. The original
GSHHG 2.3.7 polygons and original speciesKeys must be retained. A GBIF record
does not establish immigration, breeding, or a host-plant-mediated pathway.
"""
from __future__ import annotations
import argparse, collections, csv, hashlib, json, time
from datetime import datetime, timezone
from pathlib import Path

from build_public51_butterfly_risk_20261109 import GBIF, GSHHG_SHA256, load_polygons, get_json

FIXED = {
    ("2018_2022", "7762841", "26"): "Ireland / Nymphalis polychloros",
    ("2018_2022", "1929270", "391"): "Madeira / Lycaena dispar",
    ("2018_2022", "8049830", "456"): "Orkney Mainland / Pararge aegeria",
}
MAX_RESULTS_PER_CASE = 500
PRESERVE = (
    "gbifID","key","speciesKey","taxonKey","species","scientificName",
    "acceptedScientificName","taxonRank","taxonomicStatus","eventDate",
    "year","basisOfRecord","lifeStage","country","countryCode","stateProvince",
    "locality","verbatimLocality","datasetKey","occurrenceID","references",
    "coordinateUncertaintyInMeters","issues","geospatialIssue",
    "identificationRemarks","occurrenceRemarks","license",
    "decimalLatitude","decimalLongitude",
)

def sha(x):
    return hashlib.sha256(x).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--model-json",type=Path,required=True)
    p.add_argument("--case-crosswalk",type=Path,required=True)
    p.add_argument("--gshhg-zip",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    args=p.parse_args()
    model_b=args.model_json.read_bytes()
    cross_b=args.case_crosswalk.read_bytes()
    model=json.loads(model_b)
    cross=json.loads(cross_b)
    eight=model.get("revised_historical_candidate_first_index_intro_host_cases",{})
    indices={(period,str(x["species"]),str(x["island"])):x
             for period,rows in eight.items() for x in rows}
    checked={(x["period"],x["speciesKey"],x["landmass_id"]):x for x in cross["cases"]}
    if len(indices)!=8 or len(checked)!=8 or set(indices)!=set(checked):
        raise ValueError("Original model and external historical crosswalk disagree on fixed eight cases")
    if not set(FIXED)<=set(indices):
        raise ValueError("Prespecified three targeted cases not in historical model")
    if any(checked[k]["class"] not in (
        "SOURCE_CONFLICT_NEEDS_RECORD_VALIDATION","NO_ESTABLISHMENT_VERIFICATION"
    ) for k in FIXED):
        raise ValueError("External source triage status changed unexpectedly")
    geometries,geoaudit,gshhg=load_polygons(args.gshhg_zip,args.work,
                                          {k[2] for k in FIXED},2300)
    if gshhg!=GSHHG_SHA256:
        raise ValueError("GSHHG checksum mismatch")
    out=[]
    for period,key,island in sorted(FIXED):
        values=geometries[island]
        rows=[];receipts=[];failures=[];ids=set()
        observed_count=0
        for frag,wkt in enumerate(values,1):
            offset=0
            for page in range(10):
                params={
                    "taxonKey":key,"geometry":wkt,"year":"2018,2022",
                    "hasCoordinate":"true","hasGeospatialIssue":"false",
                    "occurrenceStatus":"PRESENT","limit":100,"offset":offset,
                }
                try:
                    response=get_json(GBIF+"/occurrence/search",params)
                    found=response.get("results",[])
                    n=int(response.get("count",-1))
                    if n<0 or not isinstance(found,list):raise ValueError("Invalid GBIF occurrence response")
                    original=sha(json.dumps(response,sort_keys=True,ensure_ascii=False).encode("utf8"))
                    receipts.append({
                        "fragment":frag,"page":page,"offset":offset,"gbif_total_count":n,
                        "response_sha256":original,
                        "WKT_sha256":sha(wkt.encode("utf8")),
                        "source_query":"same original 51 landmass geometry / species taxonKey / indexed years",
                    })
                    observed_count+=len(found)
                    for record in found:
                        rid=str(record.get("gbifID") or record.get("key") or "")
                        if not rid:continue
                        if rid in ids:continue
                        ids.add(rid)
                        rows.append({z:record.get(z) for z in PRESERVE if record.get(z) is not None})
                    offset+=len(found)
                    if offset>=n or not found:break
                    if offset>=MAX_RESULTS_PER_CASE:break
                except Exception as e:
                    failures.append({"fragment":frag,"offset":offset,"error":str(e)[:240]})
                    break
        exact=[x for x in rows if str(x.get("speciesKey"))==key]
        other=[x for x in rows if str(x.get("speciesKey"))!=key]
        status="COMPLETE_WITH_RECORDS" if not failures and exact else (
            "COMPLETE_NO_EXACT_SPECIES_RECORD_POSSIBLE_SOURCE_DRIFT" if not failures
            else "INCOMPLETE_SOURCE_ERROR")
        counter=lambda field:dict(collections.Counter(str(x.get(field,"UNKNOWN")) for x in exact))
        out.append({
            "species":indices[(period,key,island)]["species_name"],
            "speciesKey":key,"island_id":island,"island_source_description":FIXED[(period,key,island)],
            "source_history_class":checked[(period,key,island)]["class"],
            "period":period,"status":status,
            "fragments":len(values),"query_pages":len(receipts),"query_errors":failures,
            "n_gbif_unique_records":len(rows),"n_exact_speciesKey":len(exact),
            "n_nonexact_speciesKey":len(other),
            "first_year":min((x["year"] for x in exact if isinstance(x.get("year"),int)),default=None),
            "last_year":max((x["year"] for x in exact if isinstance(x.get("year"),int)),default=None),
            "basisOfRecord_counts":counter("basisOfRecord"),
            "datasetKey_counts":counter("datasetKey"),
            "lifeStage_counts":counter("lifeStage"),
            "original_model_details":indices[(period,key,island)],
            "source_records":exact[:MAX_RESULTS_PER_CASE],
            "nonexact_records":other[:20],
            "query_receipts":receipts,
        })
        print(json.dumps({
            "case":FIXED[(period,key,island)],"status":status,
            "indexed_exact_records":len(exact),
            "years":counter("year"),"basis":counter("basisOfRecord"),
            "errors":len(failures)
        },ensure_ascii=False),flush=True)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    obj={
        "schema":"chocho_fixed_three_original_GBIF_record_identity_audit_v1",
        "audited_at_utc":datetime.now(timezone.utc).isoformat(),
        "historical_model_sha256":sha(model_b),
        "external_historical_crosswalk_sha256":sha(cross_b),
        "GSHHG_source_zip_sha256":gshhg,
        "GSHHG_fragments":geoaudit,
        "scope":"Three model-preselected source-conflicted or unverified 2018-2022 first GBIF indexes",
        "cases":out,
        "limitations":[
            "Original 51 landmass indexing cohort is selected on preexisting port/butterfly data",
            "GBIF search can drift when indexing is updated",
            "No indexed record does not mean ecological absence",
            "taxonKey searches can include other species/children and exact speciesKey rows must be separated",
            "Published regional checklists and original GBIF records measure different domains",
            "No single adult or larval GBIF observation proves settlement, offspring production or host usage"
        ]
    }
    path=args.out_dir/"fixed_three_original_GBIF_record_audit.json"
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("SUMMARY",json.dumps({x["island_source_description"]: {
        "exact_records":x["n_exact_speciesKey"],"status":x["status"],
        "GBIF_basis":x["basisOfRecord_counts"]
    } for x in out},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
