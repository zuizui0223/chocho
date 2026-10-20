#!/usr/bin/env python3
"""Exact source-locked GBIF GADM1 and developmental-stage metadata audit.

48 already-positive species×country cells (20 species, 23 countries) only.
No individual occurrence record bodies or precise coordinate access.
No host-use, establishment, colonization, or independent-population inference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.gbif.org/v1/occurrence/search"
BASIS = ("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
         "MACHINE_OBSERVATION", "OCCURRENCE")
FIELDS = ("gadmLevel1Gid", "lifeStage", "establishmentMeans")
FACET_MAP = {"GADMLEVEL1GID": "gadm1", "LIFESTAGE": "life_stage",
             "ESTABLISHMENTMEANS": "establishment_means"}
ORIGINAL_SHA = "18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2"
COUNTRY_SHA = "d9f5207cd611e737279b7b579fe8f2347132d68bab36fc78a56e727f498c0359"
STRATA_SHA = "6eaa91001b282a8bb8bb3f78847386ef5c4fbefcf273f4b09d173d6a691100ec"
SCHEMA = "chocho_subcountry_48_positive_GBIF_GADM1_lifestage_facet_v1"
PILOT = {("Aguna asander", "DO"): 9, ("Papilio demoleus", "IL"): 86}


def digest(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def key_field(field):
    return re.sub(r"[^A-Z0-9]", "", str(field).upper())


def load_locked_positive_cells(strata_path, taxon_manifest, counts_path):
    strata_b = Path(strata_path).read_bytes()
    names_b = Path(taxon_manifest).read_bytes()
    counts_b = Path(counts_path).read_bytes()
    if (digest(strata_b), digest(names_b), digest(counts_b)) != (
            STRATA_SHA, ORIGINAL_SHA, COUNTRY_SHA):
        raise ValueError("Changed source bytes: original 48-cell or taxonomy/country archives")
    strata = json.loads(strata_b)
    manifest = json.loads(names_b)
    if len(strata) != 255 or len(manifest["species_receipts"]) != 239:
        raise ValueError("Frozen 255 strata / 239 species not available")
    taxa = {
        row["species"]: row["taxon"]["key"]
        for row in manifest["species_receipts"]
        if row["status"] != "TAXON_UNRESOLVED"
    }
    source_counts = Counter()
    with Path(counts_path).open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != ["species", "window", "country", "basis_scope", "record_count"]:
            raise ValueError("Frozen facet CSV header changed")
        for row in reader:
            if row["basis_scope"] == "all_accepted_bases":
                source_counts[(row["species"], row["country"], row["window"])] += int(row["record_count"])
    cells = []
    seen = set()
    for st in strata:
        for cell in st["cells"]:
            if not cell["positive"]:
                continue
            key = (st["species"], cell["country"])
            if key in seen or key[0] not in taxa:
                raise ValueError("Duplicate or unknown positive source species-country pair")
            seen.add(key)
            total = sum(source_counts[(key[0], key[1], w)] for w in ("2010_2017", "2018_2025"))
            if total <= 0:
                raise ValueError("A positive primary source cell has no indexed records")
            cells.append({"species": key[0], "country": key[1],
                          "taxon_key": taxa[key[0]], "added_resource_country": int(cell["added"]),
                          "source_pooled_count": total,
                          "source_late_count": source_counts[(key[0],key[1],"2018_2025")]})
    if len(cells) != 48 or len(seen) != 48:
        raise ValueError("Frozen 48 positive country/source cells missing")
    if (len({r["species"] for r in cells}), len({r["country"] for r in cells}),
            sum(c["added_resource_country"] for c in cells)) != (20, 23, 26):
        raise ValueError("Frozen 20 species / 23 countries / 26 added positives changed")
    for pair, n in PILOT.items():
        record = next((r for r in cells if (r["species"],r["country"]) == pair), None)
        if record is None or record["source_late_count"] != n:
            raise ValueError("Pilot is not the previously audited late-record pair")
    return sorted(cells, key=lambda x:(x["species"],x["country"])), {
        "source_strata_sha256": STRATA_SHA,
        "source_taxon_manifest_sha256": ORIGINAL_SHA,
        "source_country_csv_sha256": COUNTRY_SHA,
        "source_positive_cells": 48,
        "source_distinct_species": 20,
        "source_distinct_countries": 23,
        "source_added_resource_positive_cells": 26,
    }


def url_for(r, *, pilot=False, gadm_only=False):
    params = [
        ("taxonKey",str(r["taxon_key"])), ("country",r["country"]),
        ("year","2018,2025" if pilot else "2010,2025"),
        ("hasCoordinate","true"), ("hasGeospatialIssue","false"),
        ("occurrenceStatus","PRESENT"),("limit","0"), ("facetLimit","300")
    ]
    params += [("basisOfRecord",v) for v in BASIS]
    params += [("facet","gadmLevel1Gid")] if gadm_only else [("facet",v) for v in FIELDS]
    return API+"?"+urlencode(params)


def parse_facets(raw, *, require_gadm=True):
    try:
        body=json.loads(raw)
    except (UnicodeDecodeError, ValueError) as e:
        raise ValueError("Invalid/truncated GBIF facet JSON") from e
    if not isinstance(body,dict) or type(body.get("count")) is not int or body["count"] < 0:
        raise ValueError("GBIF count must be nonnegative integer")
    if body.get("limit") != 0 or body.get("results") != []:
        raise ValueError("Unexpected individual occurrence bodies")
    facets=body.get("facets")
    if not isinstance(facets,list):
        raise ValueError("GBIF omitted facet list")
    bins={}
    for f in facets:
        if not isinstance(f,dict):
            raise ValueError("Malformed facet field")
        kind=FACET_MAP.get(key_field(f.get("field")))
        if not kind or kind in bins:
            raise ValueError("Unknown/duplicate requested GBIF facet")
        parts=f.get("counts")
        if not isinstance(parts,list) or len(parts)>=300:
            raise ValueError("Missing/overflowed category list")
        categories={}
        for p in parts:
            if not isinstance(p,dict):
                raise ValueError("Malformed category")
            code=p.get("name")
            n=p.get("count")
            if not isinstance(code,str) or not code.strip() or len(code)>160 or type(n) is not int or n<0 or code in categories:
                raise ValueError("Invalid/duplicate facet name or count")
            categories[code]=n
        summed=sum(categories.values())
        if summed>body["count"]:
            raise ValueError("Facet exceeds total GBIF records")
        bins[kind]={"categories":dict(sorted(categories.items())),
                    "assigned_records":summed,"unassigned_records":body["count"]-summed}
    if require_gadm and "gadm1" not in bins:
        raise ValueError("No GADM1 facet returned")
    return {"record_count":body["count"],"facets":bins,
            "missing_optional_facets":sorted({"life_stage","establishment_means"}-set(bins))}


def fetch_one(row, *, pilot=False, max_attempts=3):
    attempts=[]
    for gadm_only in (False,True):
        url=url_for(row,pilot=pilot,gadm_only=gadm_only)
        for attempt in range(max_attempts):
            try:
                rq=Request(url,headers={"User-Agent":"chocho-GADM1-source-audit/1.0",
                                        "Accept":"application/json"})
                with urlopen(rq,timeout=35) as res:
                    status=res.status
                    raw=res.read(1_000_000)
                if status!=200 or len(raw)>=1_000_000:
                    raise ValueError(f"HTTP/response capacity: {status}/{len(raw)}")
                parsed=parse_facets(raw)
                return {**row,**parsed,"status":"COMPLETE_GADM1",
                        "source_pooled_delta": (parsed["record_count"]-(
                            row["source_late_count"] if pilot else row["source_pooled_count"])),
                        "source_window":"2018_2025" if pilot else "2010_2025",
                        "requested_gadm_only_fallback":gadm_only,
                        "query_url":url,"retrieved_utc":now(),
                        "raw_sha256":digest(raw),"http_status":status,
                        "attempts":attempts+[{"attempt":attempt+1,"gadm_only":gadm_only,"status":200}]}
            except (HTTPError,URLError,TimeoutError,OSError,ValueError) as e:
                attempts.append({"attempt":attempt+1,"gadm_only":gadm_only,
                                 "error_type":type(e).__name__,"message":str(e)[:350]})
                time.sleep(1.0 if attempt < max_attempts-1 else 0)
        # Try a smaller compulsory GADM-only facet if an optional facet was
        # unsupported; never silently claim the optional metadata complete.
    return {**row,"status":"GADM1_UNAVAILABLE_NOT_ZERO",
            "source_window":"2018_2025" if pilot else "2010_2025",
            "query_url":url_for(row,pilot=pilot),
            "attempts":attempts,"retrieved_utc":now()}


def summarize(cells, pilot):
    complete=[x for x in cells if x["status"]=="COMPLETE_GADM1"]
    exact=[x for x in complete if x["source_pooled_delta"]==0]
    observed=[]
    for x in complete:
        n=x["record_count"]
        bins=x["facets"]["gadm1"]
        mapped=bins["assigned_records"]
        counts=list(bins["categories"].values())
        dominant=(max(counts)/n) if counts and n else 0
        observed.append({**x,"mapped_gadm_units":len(counts),
                         "mapped_records":mapped,
                         "gadm_coverage_fraction_of_all_records":mapped/n if n else None,
                         "dominant_gadm_fraction_of_all_records":dominant,
                         "one_gadm_unit":len(counts)==1,
                         "life_stage_annotated":x["facets"].get("life_stage",{}).get("assigned_records",0),
                         "establishment_annotated":x["facets"].get("establishment_means",{}).get("assigned_records",0)})
    total_records=sum(x["record_count"] for x in complete)
    by_added={}
    for status in (0,1):
        group=[x for x in observed if x["added_resource_country"]==status]
        by_added[str(status)]={
            "positive_cells":len(group),
            "single_gadm_units":sum(x["one_gadm_unit"] for x in group),
            "gadm_unassigned_cells":sum(x["mapped_records"]==0 for x in group),
            "sum_records":sum(x["record_count"] for x in group),
            "sum_mapped_records":sum(x["mapped_records"] for x in group)
        }
    return {
        "schema":SCHEMA,
        "status":"COMPLETE_48" if not pilot and len(exact)==48 else (
            "COMPLETE_PILOT" if pilot and len(exact)==2 else "PARTIAL_OR_INDEX_DRIFT_NOT_ZERO"),
        "pilot":pilot,"eligible_cells":2 if pilot else 48,
        "gadm1_facets_complete":len(complete),
        "source_count_exact_cells":len(exact),
        "source_count_drifted_cells":len(complete)-len(exact),
        "indexed_records_across_completed_queries":total_records,
        "gadm1_assigned_records":sum(x["mapped_records"] for x in observed),
        "life_stage_metadata_records":sum(x["life_stage_annotated"] for x in observed),
        "establishment_metadata_records":sum(x["establishment_annotated"] for x in observed),
        "cells_with_one_mapped_gadm_unit":sum(x["one_gadm_unit"] for x in observed),
        "cells_with_no_mapped_gadm_unit":sum(x["mapped_records"]==0 for x in observed),
        "by_added_resource_flag":by_added,
        "observations":observed,
        "inference":"GADM1 reporting-unit mapping and provider life-stage strings only; not native population, host usage, recruitment or butterfly range expansion"
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--strata",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--counts",type=Path,required=True)
    ap.add_argument("--pilot",action="store_true")
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    cells,source=load_locked_positive_cells(args.strata,args.manifest,args.counts)
    selected=([x for x in cells if (x["species"],x["country"]) in PILOT]
              if args.pilot else cells)
    fetched=[]
    for row in selected:
        fetched.append(fetch_one(row,pilot=args.pilot))
        time.sleep(1)
    result=summarize(fetched,pilot=args.pilot)
    result["source"]=source
    result["executed_utc"]=now()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    compact={k:v for k,v in result.items() if k!="observations"}
    print("GBIF_SUBCOUNTRY_SUMMARY="+json.dumps(compact,separators=(",",":")),flush=True)
    if result["status"] not in ("COMPLETE_48","COMPLETE_PILOT"):
        raise SystemExit("GADM1/POSITIVE48_INCOMPLETE_OR_INDEX_DRIFT; do not infer absence/host use")


if __name__=="__main__":
    main()
