#!/usr/bin/env python3
"""Small official GBIF backbone/occurrence audit for two unresolved Mauritius cases.

Only source QA. Never infer arrival, native origin, breeding or endemism.
No bulk GBIF record downloads. Manually triggered, four taxon queries and one
capped country-level occurrence query.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

BASE="https://api.gbif.org/v1/"
URLS={
 "current_taxon_7346444":"species/7346444",
 "zizina_antanossa_name_match":"species/match?"+urllib.parse.urlencode({"name":"Zizina antanossa"}),
 "zizina_otis_name_match":"species/match?"+urllib.parse.urlencode({"name":"Zizina otis"}),
 "zerynthia_polyxena_name_match":"species/match?"+urllib.parse.urlencode({"name":"Zerynthia polyxena"}),
 "zerynthia_MU_2018_2022":"occurrence/search?"+urllib.parse.urlencode({
        "taxon_key":5141262,"country":"MU","year":"2018,2022","limit":10})
}
def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"chocho-GBIF-historical-island-record-audit/1.0 (research source validation)","Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=45) as stream:
        raw=stream.read(250000)
    return json.loads(raw),hashlib.sha256(raw).hexdigest(),len(raw)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",type=Path,default=Path("_mauritius_GBIF_taxon_source"))
    args=p.parse_args()
    receipt={"schema":"chocho_mauritius2011_GBIF_taxon_key_and_one_record_source_diagnostic_v1",
      "status":"OFFICIAL_GBIF_TAXON_OBSERVATION_ONLY_NOT_COLONIZATION_OR_NATIVE_ORIGIN",
      "retrieved_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
      "endpoints":{}}
    for id,relative in URLS.items():
        url=BASE+relative
        x,sha,bytes_=get(url)
        if id=="current_taxon_7346444":
            fields=["key","rank","nameType","scientificName","canonicalName","taxonomicStatus","status",
                    "acceptedKey","speciesKey","species","genus","family","authorship"]
            data={k:x.get(k) for k in fields}
        elif id=="zerynthia_MU_2018_2022":
            attrs=["key","speciesKey","acceptedScientificName","scientificName","taxonKey",
                   "decimalLatitude","decimalLongitude","countryCode","year",
                   "basisOfRecord","occurrenceStatus","geodeticDatum","datasetKey",
                   "identifiedBy","recordedBy","issue","issues","coordinateUncertaintyInMeters"]
            data={"count":x.get("count"),"endOfRecords":x.get("endOfRecords"),
                  "records":[{k:z.get(k) for k in attrs} for z in x.get("results",[])[:10]]}
        else:
            fields=["usageKey","scientificName","canonicalName","rank","status","matchType","confidence",
                    "note","synonym","acceptedUsageKey","speciesKey","kingdom","family"]
            data={k:x.get(k) for k in fields}
            # Alternative matches useful for source ambiguity; never automatically call them native.
            data["alternatives"]=[{k:a.get(k) for k in fields} for a in x.get("alternatives",[])[:5]]
        receipt["endpoints"][id]={"url":url,"sha256_source_bytes":sha,"source_bytes":bytes_,"fields":data}
        print(json.dumps({id:data},ensure_ascii=False),flush=True)
    receipt["interpretation"]=[
      "Davis and Barnes 1993 used Zizina antanossa; Yago et al. 2008 proposed including it within Z. otis.",
      "A GBIF name match alone does not confirm the 1993 historical specimen concept for original speciesKey 7346444; inspect authority and accepted taxonomy.",
      "One current country occurrence search can disagree with the archived 2026-11-09 51-island period index due to geography, source or API revision; do not silently delete or relabel original records.",
      "Zerynthia polyxena in Mauritius is a geographically anomalous candidate for raw occurrence QA, not a proved error.",
      "No change to original 317 first-index labels or butterfly native/endemic/breeding status."
    ]
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"GBIF_MAURITIUS_TAXON_AND_OCCURRENCE_SOURCE_AUDIT_20261117.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("MAURITIUS_TARGETED_GBIF_SOURCE_AUDIT_SUCCESS")
if __name__=="__main__": main()
