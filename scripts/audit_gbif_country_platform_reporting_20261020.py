#!/usr/bin/env python3
"""All-taxon and Lepidoptera recording denominators for fixed GBIF candidate cells.

Two sources of reporting activity, NOT biological abundance or independent
population evidence. No individual records are downloaded or redistributed.
The old GBIF Lepidoptera order key 797 is conditionally verified for nesting
after the 2026 taxonomy migration rather than assumed transferable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
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
YRS = tuple(str(y) for y in range(2018, 2026))
FACETS = ("all_taxa", "lepidoptera_legacy_797")
SCHEMA = "gbif_butterfly_country_dataset_year_reporting_denominators_v1"


def sha(b):
    return hashlib.sha256(b).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def read_locked_dataset_year(path: Path):
    raw = path.read_bytes()
    p = json.loads(raw)
    if p.get("schema") != "gbif_country_seven_candidate_dataset_by_year_qc_v1":
        raise ValueError("Unexpected source year matrix")
    if p.get("status") != "ALL_SEVEN_SOURCE_CELLS_RECONCILED":
        raise ValueError("Missing year-by-dataset source")
    if p.get("complete_country_cells") != 7 or p.get("complete_dataset_queries") != 11:
        raise ValueError("Frozen seven/eleven cohort changed")
    if not p.get("outcomes_cannot_establish_butterfly_colonization"):
        raise ValueError("Unsafe interpretation flag missing")
    sources = []
    seen = set()
    for cell in p["country_cells"]:
        if cell["status"] != "COMPLETE_RECONCILED":
            raise ValueError("Source cell incomplete")
        for d in cell["datasets"]:
            key = (d["dataset_key"], cell["country"])
            if key in seen:
                raise ValueError("Repeated dataset/country group; separate taxa need explicit pooling")
            seen.add(key)
            if d["status"] != "COMPLETE_MATCHED_SOURCE_COUNT":
                raise ValueError("Incomplete sub-source")
            if d["country"] != cell["country"] or d["species"] != cell["species"]:
                raise ValueError("Dataset identity changed")
            if sum(d["year_counts"].values()) != d["source_dataset_count"]:
                raise ValueError("Source year count differs from frozen total")
            sources.append({"species": cell["species"], "country": cell["country"],
                            "dataset_key": d["dataset_key"], "dataset_title": d["dataset_title"],
                            "source_years": d["year_counts"],
                            "source_total": d["source_dataset_count"],
                            "added_host_country": cell["added_host_country"]})
    if len(seen) != 11 or len(sources) != 11:
        raise ValueError("Not the 11 locked dataset/country reporting units")
    return sorted(sources, key=lambda z: (z["species"], z["country"], z["dataset_key"])), {
        "dataset_year_source_sha256": sha(raw), "source_group_count": len(seen),
        "seven_country_species_cases": len(p["country_cells"]),
        "original_record_years": "2018-2025", "source_query_type": "country/dataset/year facets"
    }


def url_for(row, kind):
    if kind not in FACETS:
        raise ValueError("Unknown reporting denominator")
    params = [
        ("datasetKey", row["dataset_key"]), ("country", row["country"]),
        ("year", "2018,2025"), ("hasCoordinate", "true"),
        ("hasGeospatialIssue", "false"), ("occurrenceStatus", "PRESENT"),
        ("limit", "0"), ("facet", "year"), ("facetLimit", "300")
    ] + [("basisOfRecord", basis) for basis in BASIS]
    if kind == "lepidoptera_legacy_797":
        params.append(("orderKey", "797"))
    return API + "?" + urlencode(params)


def parse_years(raw):
    try:
        o = json.loads(raw)
    except (UnicodeDecodeError, ValueError) as exc:
        raise ValueError("GBIF invalid JSON body") from exc
    if type(o.get("count")) is not int or o["count"] < 0:
        raise ValueError("Invalid GBIF total count")
    if o.get("limit") != 0 or o.get("results") != []:
        raise ValueError("Occurrence bodies must not be returned")
    facet = o.get("facets")
    if not isinstance(facet, list) or len(facet) != 1 or not isinstance(facet[0], dict):
        raise ValueError("GBIF year facet missing")
    if str(facet[0].get("field")).upper().replace("_","") != "YEAR":
        raise ValueError("Unexpected facet")
    counts = facet[0].get("counts")
    if not isinstance(counts, list) or len(counts) >= 300:
        raise ValueError("Invalid facet count rows")
    out = {}
    for c in counts:
        if not isinstance(c, dict):
            raise ValueError("Invalid facet entry")
        y, n = c.get("name"), c.get("count")
        if y not in YRS or type(n) is not int or n <= 0 or y in out:
            raise ValueError("Unacceptable year or duplicated facet entry")
        out[y] = n
    if sum(out.values()) != o["count"]:
        raise ValueError("GBIF year facet totals do not match")
    return {"index_count": o["count"], "year_counts": dict(sorted(out.items()))}


def fetch(row, kind, *, tries=3):
    url = url_for(row, kind)
    failures = []
    for attempt in range(tries):
        try:
            req = Request(url, headers={"User-Agent": "chocho-gbif-platform-denominator/0.1",
                                        "Accept": "application/json"})
            with urlopen(req, timeout=45) as r:
                status = r.status
                raw = r.read(1000000)
            if status != 200 or len(raw) >= 1000000:
                raise ValueError("Incomplete GBIF HTTP response")
            item = parse_years(raw)
            return {**item, "status": "COMPLETE", "url": url,
                    "http_status": status, "raw_response_sha256": sha(raw),
                    "retrieved_utc": now(), "attempts": attempt + 1}
        except (HTTPError, URLError, ValueError, TimeoutError, OSError) as exc:
            failures.append({"attempt": attempt + 1, "kind": type(exc).__name__,
                             "details": str(exc)[:250]})
            if attempt < tries-1:
                time.sleep(2*(attempt+1))
    return {"status": "UNAVAILABLE_NOT_ZERO", "url": url,
            "retrieved_utc": now(), "attempts": tries, "errors": failures}


def assess_nested(focal_years, denominator):
    if denominator["status"] != "COMPLETE":
        return {"status": "UNAVAILABLE_NOT_ZERO"}
    yrs = denominator["year_counts"]
    mismatches = {y: {"focal": k, "denominator": yrs.get(y,0)}
                  for y, k in focal_years.items() if yrs.get(y,0) < k}
    if mismatches:
        return {"status": "DENOMINATOR_NOT_NESTED_DO_NOT_INTERPRET",
                "mismatches": mismatches}
    rates = {y: round(10000 * k/yrs[y], 8) for y,k in focal_years.items() if k and yrs.get(y,0)}
    return {"status": "NESTED_VALID", "focal_rate_per_10000": rates,
            "all_year_total": denominator["index_count"]}


def assemble(rows, fetched, source):
    if len(rows) != 11 or len(fetched) != 22:
        raise ValueError("Missing exact eleven×two denominator attempts")
    results = []
    for p in rows:
        key = (p["dataset_key"], p["country"])
        this = {kind: fetched[(key[0], key[1], kind)] for kind in FACETS}
        grades = {kind: assess_nested(p["source_years"], this[kind]) for kind in FACETS}
        results.append({**p, "denominators": this, "nested_quality": grades})
    primary = next(r for r in results if (r["species"],r["country"],r["dataset_title"]) ==
                   ("Papilio demoleus","IL","iNaturalist Research-grade Observations"))
    inpa = next(r for r in results if (r["species"],r["country"],r["dataset_title"]) ==
                ("Papilio demoleus","IL","INPA - Observations"))
    def compare_years(item, kind):
        f = item["source_years"]
        d = item["denominators"][kind].get("year_counts",{})
        if item["nested_quality"][kind]["status"] != "NESTED_VALID" or not d.get("2024") or not d.get("2025"):
            return {"status": "NOT_COMPARABLE"}
        x24,x25 = int(f.get("2024",0)),int(f.get("2025",0))
        if not x24 or not x25:
            return {"status": "FOCAL_YEAR_ZERO_NOT_RATIO_COMPARABLE"}
        rate24=x24/d["2024"]
        rate25=x25/d["2025"]
        return {"status":"DESCRIPTIVE_RATES_ONLY","focal_2024":x24,"focal_2025":x25,
                "reporting_2024":d["2024"],"reporting_2025":d["2025"],
                "relative_focal_rate_2025_to_2024":round(rate25/rate24,8),
                "focal_2024_per_10000":round(rate24*10000,8),
                "focal_2025_per_10000":round(rate25*10000,8),
                "reporting_volume_ratio_2025_to_2024":round(d["2025"]/d["2024"],8)}
    summary = {
        "israel_papilio": {provider: {kind: compare_years(item,kind) for kind in FACETS}
                          for provider, item in (("iNaturalist",primary),("INPA",inpa))},
        "dominican_aguna": next({ "country":r["country"],
              "source_years":r["source_years"],
              "reporting_years": {kind:r["denominators"][kind].get("year_counts",{})
                                  for kind in FACETS},
              "nested_status": {kind:r["nested_quality"][kind]["status"] for kind in FACETS}}
              for r in results if (r["species"],r["country"]) == ("Aguna asander","DO"))
    }
    all_taxa_ok = sum(r["nested_quality"]["all_taxa"]["status"]=="NESTED_VALID" for r in results)
    Lep_ok = sum(r["nested_quality"]["lepidoptera_legacy_797"]["status"]=="NESTED_VALID" for r in results)
    return {
        "schema": SCHEMA,"source":source,"executed_utc":now(),
        "status": "ALL_BROAD_DENOMINATORS_VALID" if all_taxa_ok==11 else "BROAD_DENOMINATOR_UNAVAILABLE_DO_NOT_INTERPRET",
        "source_conditioned_queries":22,"all_taxa_compatible_groups":all_taxa_ok,
        "lepidoptera_legacy_key_compatible_groups":Lep_ok,
        "cannot_infer_butterfly_abundance_or_breeding":True,
        "cases":results,"focal_summary":summary
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dataset-year-source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    rows,src=read_locked_dataset_year(args.dataset_year_source)
    attempts={}
    for row in rows:
        for kind in FACETS:
            attempts[(row["dataset_key"],row["country"],kind)]=fetch(row,kind)
            time.sleep(1.1)
    result=assemble(rows,attempts,src)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GBIF_PLATFORM_REPORTING_RESULT="+json.dumps(result,separators=(",",":")),flush=True)
    if result["status"]!="ALL_BROAD_DENOMINATORS_VALID":
        raise SystemExit("Platform denominator incomplete, saved source-failure receipts.")


if __name__=="__main__":
    main()
