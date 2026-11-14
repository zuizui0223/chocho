#!/usr/bin/env python3
"""Trace GBIF-indexed Pieris rapae larvae to original observation and plant evidence.

Strict evidence tiers:
- Indexed larva, 2020 and 2022: previously source-verified facts.
- Explicit original 'feeding on' or named host/plant observation field: potential
  source text evidence, retained verbatim and flagged for manual verification.
- An image of a larva on vegetation is NOT feeding or host identity proof.
- Historical port candidate plants from 2001-2010 are area-wide opportunities,
  not spatially/temporally verified local larval food.
No nonindex is treated as true absence. No raw photo files are reproduced.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

RECORD_IDS = ("2814126623", "3772636962")
BASE = "https://api.gbif.org/v1/occurrence/"
INAT_API = "https://api.inaturalist.org/v1/observations/"
INAT = re.compile(r"(?:www\.)?inaturalist\.org/observations/(\d+)")
RELATED = re.compile(r"\b(?:host plant|food plant|feeding|feeding on|cabbage|"
                     r"broccoli|mustard|brassica|lepidium|sisymbrium|"
                     r"nasturtium|plant species|eating|oviposit)\b", re.I)
TEXT_KEYS = (
    "associatedTaxa","associatedOccurrences","occurrenceRemarks",
    "identificationRemarks","eventRemarks","fieldNotes","dynamicProperties",
    "verbatimLocality","locality","habitat","eventDate","lifeStage",
    "references","occurrenceID","organismRemarks","samplingProtocol","behavior",
)
GBIF_WHITELIST = [
    "gbifID","key","species","speciesKey","taxonKey","scientificName",
    "eventDate","year","month","day","lifeStage","basisOfRecord",
    "country","countryCode","stateProvince","locality","verbatimLocality",
    "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
    "geospatialIssue","issues","datasetKey","publishingOrgKey",
    "occurrenceID","references","associatedTaxa","associatedOccurrences",
    "occurrenceRemarks","identificationRemarks","eventRemarks","fieldNotes",
    "dynamicProperties","habitat","organismRemarks","samplingProtocol",
    "behavior","recordedBy","lastInterpreted","license","taxonRank",
]
PORTS = {
    # Frozen WPI 2019 OBJECTID, validated against original 51-island graph.
    "912": {"name": "ISTHMUS COVE", "lat": 33.449999999590645, "lon": -118.4999999999738},
    "1012": {"name": "AVALON", "lat": 33.350000000174745, "lon": -118.31666666674916},
}


def get_json(url: str, timeout=28):
    last = None
    for attempt in range(2):
        try:
            with urlopen(Request(url, headers={
                "User-Agent":"chocho-larval-host-evidence-audit/20261114",
                "Accept":"application/json",
            }), timeout=timeout) as req:
                b = req.read()
                if req.status != 200:
                    raise ValueError(f"HTTP {req.status}")
            return json.loads(b), hashlib.sha256(b).hexdigest(), None
        except (URLError, HTTPError, TimeoutError, OSError, ValueError, json.JSONDecodeError) as exc:
            last = f"{type(exc).__name__}: {str(exc)[:180]}"
            time.sleep(0.5)
    return None, None, last


def observation_id(gbif):
    candidates = []
    for key in ("references", "occurrenceID", "identifier", "associatedMedia"):
        value = gbif.get(key)
        if value is not None:
            candidates.extend(INAT.findall(str(value)))
    return candidates[0] if candidates else None


def haversine_km(a, b, c, d):
    r = 6371.0088
    a,b,c,d = map(math.radians,[a,b,c,d])
    dl, dp = d-b, c-a
    x = math.sin(dp/2)**2 + math.cos(a)*math.cos(c)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1,math.sqrt(x)))


def potential_text_indicators(obj):
    found=[]
    for key in TEXT_KEYS:
        v=obj.get(key)
        if not v:
            continue
        t=v if isinstance(v,str) else json.dumps(v, ensure_ascii=False)
        if RELATED.search(t):
            found.append({"field":key,"text":t[:2000]})
    return found


def extract_inaturalist(j):
    if not isinstance(j,dict):
        return None
    arr=j.get("results")
    if not isinstance(arr,list) or len(arr)!=1:
        return {"status":"unexpected_results", "n_results":len(arr) if isinstance(arr,list) else None}
    z=arr[0]
    fields=[]
    for f in z.get("ofvs",[]) or z.get("observation_fields",[]) or []:
        ff=f.get("field") or f.get("observation_field") or {}
        fields.append({
            "field":ff.get("name") if isinstance(ff,dict) else None,
            "value":f.get("value"),
        })
    photos=[]
    for ph in z.get("photos",[]) or []:
        photos.append({
            "id":ph.get("id"),"url":ph.get("url"),"original_url":ph.get("original_url"),
            "attribution":ph.get("attribution"),"license_code":ph.get("license_code"),
        })
    return {
        "id":z.get("id"),"uri":z.get("uri"),"quality_grade":z.get("quality_grade"),
        "description":z.get("description"),"observed_on":z.get("observed_on"),
        "observed_on_string":z.get("observed_on_string"),
        "place_guess":z.get("place_guess"),
        "taxon": {"id":(z.get("taxon") or {}).get("id"),
                  "name":(z.get("taxon") or {}).get("name")},
        "annotations":[{"controlled_attribute_id":a.get("controlled_attribute_id"),
                        "controlled_value_id":a.get("controlled_value_id")} for a in z.get("annotations",[]) or []],
        "observation_fields":fields,"tags":z.get("tags"),
        "photos":photos,
        "obscured":z.get("obscured"),"geoprivacy":z.get("geoprivacy"),
        "positional_accuracy":z.get("positional_accuracy"),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=Path("results/port_host_original_observation_20261114"))
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    args=ap.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for rid in RECORD_IDS:
        gbif, gbif_hash, gbif_error=get_json(BASE+rid)
        if gbif is None:
            rows.append({"gbifID":rid,"status":"GBIF_UNAVAILABLE","error":gbif_error})
            continue
        if str(gbif.get("speciesKey"))!="1920496" or str(gbif.get("lifeStage","")).lower()!="larva":
            rows.append({"gbifID":rid,"status":"INCONSISTENT_WITH_FACET",
                         "gbif_sha256":gbif_hash,"actual_speciesKey":gbif.get("speciesKey"),
                         "actual_lifeStage":gbif.get("lifeStage")})
            continue
        mid=observation_id(gbif)
        inat=None
        inat_hash=None
        inat_error="NO_EXPLICIT_INAT_OBSERVATION_ID" if mid is None else None
        if mid:
            inat_source, inat_hash, inat_error=get_json(INAT_API+quote(mid))
            inat=extract_inaturalist(inat_source)
        raw, raw_hash, raw_error = get_json(BASE+rid+"/verbatim")
        # Verbatim fields are source properties, potentially carrying plant names
        # that were not retained by interpreted GBIF. Audit their labels only.
        rawfields = raw if isinstance(raw,dict) else {}
        rawmatches=[]
        for key,val in rawfields.items():
            if val is not None and (RELATED.search(key) or RELATED.search(str(val))):
                rawmatches.append({"field":key,"text":str(val)[:2000]})
        lat, lon = gbif.get("decimalLatitude"),gbif.get("decimalLongitude")
        port_km={}
        if isinstance(lat,(int,float)) and isinstance(lon,(int,float)):
            for pid,port in PORTS.items():
                port_km[pid]=round(haversine_km(float(lat),float(lon),
                                             port["lat"],port["lon"]),3)
        iNatText={"description":inat.get("description"),"observation_fields":inat.get("observation_fields"),
                  "tags":inat.get("tags")} if inat else {}
        found={
            "gbif_interpreted":potential_text_indicators(gbif),
            "gbif_verbatim":rawmatches,
            "inaturalist":potential_text_indicators(iNatText),
        }
        media=[{"identifier":m.get("identifier"),"title":m.get("title"),
                "description":m.get("description"),"license":m.get("license"),
                "rightsHolder":m.get("rightsHolder")} for m in gbif.get("media",[]) or []][:12]
        rows.append({
            "gbifID":rid,"status":"GBIF_RETRIEVED",
            "gbif_response_sha256":gbif_hash,
            "gbif_url":BASE+rid,
            "gbif_fields":{k:gbif[k] for k in GBIF_WHITELIST if k in gbif and k!="recordedBy"},
            "media":media,
            "nearest_original_WPI_port_distance_km":port_km,
            "iNaturalist_id":mid,
            "iNaturalist_response_sha256":inat_hash,
            "iNaturalist_metadata":inat,
            "iNaturalist_error":inat_error,
            "gbif_verbatim_sha256":raw_hash,
            "gbif_verbatim_error":raw_error,
            "explicit_plant_or_feeding_text_indicators":found,
            "classification":(
                "EXPLICIT_TEXT_REQUIRES_MANUAL_HOST_VERIFICATION" if any(found.values())
                else "NO_EXPLICIT_HOST_FEEDING_TEXT_IN_RETRIEVED_FIELDS"
            ),
            "important_limitation":"The observed larva is not automatically feeding on a named plant; "
                                   "nearby pre-2011 port plants are not the observed larval host.",
        })
        print(json.dumps({
            "gbifID":rid,"status":rows[-1]["status"],"inat_id":mid,
            "inat_fetch_status":"SUCCESS" if inat else inat_error,
            "gbif_fields_present":[k for k in GBIF_WHITELIST if k in gbif],
            "media_n":len(media),"min_port_distance_km":min(port_km.values()) if port_km else None,
            "host_text_status":rows[-1]["classification"]
        },ensure_ascii=False),flush=True)

    result={
        "study":"Original two-year Pieris rapae larval observation audit",
        "date_utc":datetime.now(timezone.utc).isoformat(),
        "records":rows,"frozen_record_ids":RECORD_IDS,
        "original_port_points":PORTS,
        "cohort":"Selected original island 951 after previous 15-case GBIF life-stage facet",
        "limitations":[
            "No independently witnessed feeding or oviposition unless a specific source explicitly reports it",
            "Camera photo alone may show a leaf without verifying which plant species it is",
            "GBIF/iNaturalist public source metadata may have coordinate masking",
            "Port distance is descriptive; plant records were taken within 5 km of port, not necessarily at the larval observation",
            "Source location and historical plant-origin labels cannot identify a causal food-rescue mechanism"
        ],
    }
    path=args.out_dir/"original_observation_audit.json"
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2,default=str)+"\n",encoding="utf8")
    print("SUMMARY",json.dumps({
        "source_records":len(rows),
        "GBIF_retrieved":sum(x["status"]=="GBIF_RETRIEVED" for x in rows),
        "explicit_text_candidates":sum(x.get("classification","")=="EXPLICIT_TEXT_REQUIRES_MANUAL_HOST_VERIFICATION" for x in rows),
        "output":str(path)},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
