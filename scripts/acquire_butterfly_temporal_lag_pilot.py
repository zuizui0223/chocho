#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from shapely.geometry import Point, shape

GBIF = "https://api.gbif.org/v1"
UA = "chocho-temporal-lag-pilot/0.1"


def get_json(path: str, params: dict, attempts: int = 5, timeout: float = 30.0) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    last = None
    for i in range(attempts):
        try:
            req = Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            if i + 1 < attempts:
                time.sleep(min(8.0, 0.5 * (2 ** i)))
    raise RuntimeError(f"GBIF request failed: {url} ({last})")


def load_regions(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        geom_raw = feature.get("geometry")
        if code and geom_raw:
            out[code] = {
                "name": str(props.get("LEVEL3_NAM") or "").strip(),
                "geometry": shape(geom_raw),
            }
    if not out:
        raise RuntimeError("no WGSRPD3 regions")
    return out


def facet_years(taxon_key: int, geom) -> list[int]:
    minx, miny, maxx, maxy = geom.bounds
    bbox = (
        f"POLYGON(({minx} {miny},{maxx} {miny},{maxx} {maxy},"
        f"{minx} {maxy},{minx} {miny}))"
    )
    payload = get_json(
        "occurrence/search",
        {
            "taxonKey": taxon_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT",
            "year": "1900,2025",
            "geometry": bbox,
            "facet": "year",
            "facetLimit": 200,
            "facetMincount": 1,
            "limit": 1,
        },
    )
    years = []
    for facet in payload.get("facets") or []:
        if str(facet.get("field") or "").upper() != "YEAR":
            continue
        for row in facet.get("counts") or []:
            try:
                years.append(int(row.get("name")))
            except Exception:
                pass
    return sorted(set(y for y in years if 1900 <= y <= 2025))


def first_exact_record_year(taxon_key: int, geom):
    minx, miny, maxx, maxy = geom.bounds
    bbox = (
        f"POLYGON(({minx} {miny},{maxx} {miny},{maxx} {maxy},"
        f"{minx} {maxy},{minx} {miny}))"
    )
    candidate_years = facet_years(taxon_key, geom)
    for year in candidate_years:
        offset = 0
        found = []
        while True:
            payload = get_json(
                "occurrence/search",
                {
                    "taxonKey": taxon_key,
                    "hasCoordinate": "true",
                    "hasGeospatialIssue": "false",
                    "occurrenceStatus": "PRESENT",
                    "year": year,
                    "geometry": bbox,
                    "limit": 300,
                    "offset": offset,
                },
            )
            results = payload.get("results") or []
            for item in results:
                lat = item.get("decimalLatitude")
                lon = item.get("decimalLongitude")
                if lat is None or lon is None:
                    continue
                if geom.covers(Point(float(lon), float(lat))):
                    found.append(
                        {
                            "gbif_key": item.get("key"),
                            "basis_of_record": item.get("basisOfRecord"),
                            "dataset_key": item.get("datasetKey"),
                            "lat": float(lat),
                            "lon": float(lon),
                        }
                    )
            if found:
                return year, found
            if bool(payload.get("endOfRecords", True)) or not results:
                break
            offset += len(results)
            if offset >= 100000:
                break
    return None, []


def resolve(name: str):
    p = get_json("species/match", {"name": name, "strict": "false"})
    key = p.get("usageKey")
    if not key:
        return None, p
    return int(key), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    regions = load_regions(args.level3_geojson)
    with args.pairs_csv.open(newline="", encoding="utf-8") as f:
        pairs = list(csv.DictReader(f))

    names = sorted({r["plant_name"].strip() for r in pairs})
    resolved = {}
    for name in names:
        key, meta = resolve(name)
        resolved[name] = (key, meta)

    cache = {}
    rows = []
    for row in pairs:
        region = row["region"].strip()
        name = row["plant_name"].strip()
        if region not in regions:
            raise RuntimeError(f"missing WGSRPD3 region {region}")
        key, meta = resolved[name]
        out = dict(row)
        out["gbif_taxon_key"] = "" if key is None else key
        out["gbif_match_type"] = str(meta.get("matchType") or "")
        out["gbif_rank"] = str(meta.get("rank") or "")
        out["gbif_canonical_name"] = str(meta.get("canonicalName") or "")
        if key is None:
            out.update(
                {
                    "first_host_record_year": "",
                    "recorded_resource_lead_years": "",
                    "host_record_by_2017": 0,
                    "first_year_record_count_inside_region": 0,
                    "first_year_basis_of_record": "",
                    "first_year_dataset_count": 0,
                    "status": "NO_GBIF_TAXON_MATCH",
                }
            )
            rows.append(out)
            continue

        ck = (key, region)
        if ck not in cache:
            cache[ck] = first_exact_record_year(key, regions[region]["geometry"])
        first_year, found = cache[ck]
        butterfly_year = int(row["butterfly_first_test_year"])
        if first_year is None:
            out.update(
                {
                    "first_host_record_year": "",
                    "recorded_resource_lead_years": "",
                    "host_record_by_2017": 0,
                    "first_year_record_count_inside_region": 0,
                    "first_year_basis_of_record": "",
                    "first_year_dataset_count": 0,
                    "status": "NO_COORDINATE_HOST_RECORD_1900_2025",
                }
            )
        else:
            out.update(
                {
                    "first_host_record_year": first_year,
                    "recorded_resource_lead_years": butterfly_year - first_year,
                    "host_record_by_2017": int(first_year <= 2017),
                    "first_year_record_count_inside_region": len(found),
                    "first_year_basis_of_record": ";".join(
                        sorted({str(x.get("basis_of_record") or "") for x in found})
                    ),
                    "first_year_dataset_count": len(
                        {str(x.get("dataset_key") or "") for x in found}
                    ),
                    "status": "DATED",
                }
            )
        rows.append(out)

    by_cell = defaultdict(list)
    for r in rows:
        by_cell[(r["species"], r["region"], int(r["butterfly_first_test_year"]))].append(r)

    cell_rows = []
    for (species, region, butterfly_year), vals in sorted(by_cell.items()):
        years = [
            int(v["first_host_record_year"])
            for v in vals
            if str(v["first_host_record_year"]).strip()
        ]
        earliest = min(years) if years else None
        cell_rows.append(
            {
                "species": species,
                "region": region,
                "butterfly_first_test_year": butterfly_year,
                "contributing_hosts": len(vals),
                "hosts_with_dated_records": len(years),
                "earliest_documented_host_year": earliest,
                "recorded_resource_lead_years": None if earliest is None else butterfly_year - earliest,
                "any_host_documented_by_2017": bool(earliest is not None and earliest <= 2017),
            }
        )

    dated_cells = [r for r in cell_rows if r["earliest_documented_host_year"] is not None]
    lags = [int(r["recorded_resource_lead_years"]) for r in dated_cells]
    payload = {
        "schema": "chocho_butterfly_temporal_lag_pilot_v0.1",
        "status": "POSTHOC_DATED_HOST_RECORD_PILOT",
        "host_region_pairs": len(rows),
        "unique_butterfly_region_cells": len(cell_rows),
        "dated_host_region_pairs": sum(r["status"] == "DATED" for r in rows),
        "dated_cells": len(dated_cells),
        "cells_with_any_host_documented_by_2017": sum(bool(r["any_host_documented_by_2017"]) for r in cell_rows),
        "positive_host_first_lag_cells": sum(x > 0 for x in lags),
        "zero_lag_cells": sum(x == 0 for x in lags),
        "negative_host_lag_cells": sum(x < 0 for x in lags),
        "median_recorded_resource_lead_years": None if not lags else sorted(lags)[len(lags)//2],
        "cell_results": cell_rows,
        "claim_boundary": [
            "First GBIF record is a detection date, not an establishment date.",
            "Plant and butterfly recording intensity differ strongly through time.",
            "This pilot is a feasibility and chronology diagnostic on the nine introduced-resource cells with subsequent butterfly detections, not an unbiased estimate of colonization lag."
        ],
    }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
