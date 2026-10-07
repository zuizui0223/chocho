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
UA = "chocho-temporal-lag-symmetric-pilot/0.2"


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
            out[code] = shape(geom_raw)
    if not out:
        raise RuntimeError("no WGSRPD3 regions")
    return out


def resolve(name: str):
    p = get_json("species/match", {"name": name, "strict": "false"})
    key = p.get("usageKey")
    return (None if not key else int(key)), p


def bbox_wkt(geom) -> str:
    minx, miny, maxx, maxy = geom.bounds
    return (
        f"POLYGON(({minx} {miny},{maxx} {miny},{maxx} {maxy},"
        f"{minx} {maxy},{minx} {miny}))"
    )


def candidate_years(taxon_key: int, geom, start_year: int, end_year: int) -> list[int]:
    payload = get_json(
        "occurrence/search",
        {
            "taxonKey": taxon_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT",
            "year": f"{start_year},{end_year}",
            "geometry": bbox_wkt(geom),
            "facet": "year",
            "facetLimit": max(250, end_year - start_year + 1),
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
                y = int(row.get("name"))
            except Exception:
                continue
            if start_year <= y <= end_year:
                years.append(y)
    return sorted(set(years))


def first_exact_record(taxon_key: int, geom, start_year: int, end_year: int):
    for year in candidate_years(taxon_key, geom, start_year, end_year):
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
                    "geometry": bbox_wkt(geom),
                    "limit": 300,
                    "offset": offset,
                },
            )
            results = payload.get("results") or []
            for item in results:
                lat, lon = item.get("decimalLatitude"), item.get("decimalLongitude")
                if lat is None or lon is None:
                    continue
                if geom.covers(Point(float(lon), float(lat))):
                    found.append(
                        {
                            "key": item.get("key"),
                            "basisOfRecord": item.get("basisOfRecord"),
                            "datasetKey": item.get("datasetKey"),
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--start-year", type=int, default=1800)
    ap.add_argument("--end-year", type=int, default=2025)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    regions = load_regions(args.level3_geojson)
    with args.pairs_csv.open(newline="", encoding="utf-8") as f:
        pairs = list(csv.DictReader(f))

    taxa = sorted(
        {r["plant_name"].strip() for r in pairs}
        | {r["species"].strip() for r in pairs}
    )
    resolved = {name: resolve(name) for name in taxa}
    first_cache = {}

    def first_for(name: str, region: str):
        key, meta = resolved[name]
        if key is None:
            return None, [], meta
        ck = (key, region)
        if ck not in first_cache:
            first_cache[ck] = first_exact_record(
                key, regions[region], args.start_year, args.end_year
            )
        y, rec = first_cache[ck]
        return y, rec, meta

    host_rows = []
    by_cell = defaultdict(list)
    for row in pairs:
        species = row["species"].strip()
        region = row["region"].strip()
        plant = row["plant_name"].strip()
        host_year, host_rec, host_meta = first_for(plant, region)
        out = dict(row)
        out.update(
            {
                "host_gbif_key": resolved[plant][0] or "",
                "host_match_type": str(host_meta.get("matchType") or ""),
                "host_first_record_year": "" if host_year is None else host_year,
                "host_first_basis": ";".join(
                    sorted({str(x.get("basisOfRecord") or "") for x in host_rec})
                ),
                "host_first_dataset_count": len(
                    {str(x.get("datasetKey") or "") for x in host_rec}
                ),
            }
        )
        host_rows.append(out)
        by_cell[(species, region)].append(out)

    cell_rows = []
    for (species, region), vals in sorted(by_cell.items()):
        butterfly_year, butterfly_rec, butterfly_meta = first_for(species, region)
        host_years = [
            int(v["host_first_record_year"])
            for v in vals
            if str(v["host_first_record_year"]).strip()
        ]
        earliest_host = min(host_years) if host_years else None
        input_test_year = int(vals[0]["butterfly_first_test_year"])
        lag = (
            None
            if earliest_host is None or butterfly_year is None
            else butterfly_year - earliest_host
        )
        cell_rows.append(
            {
                "species": species,
                "region": region,
                "contributing_hosts": len(vals),
                "hosts_with_historical_records": len(host_years),
                "earliest_documented_host_year": earliest_host,
                "historical_butterfly_first_record_year": butterfly_year,
                "recorded_host_to_butterfly_lag_years": lag,
                "input_2018_2025_first_year": input_test_year,
                "butterfly_window_truncation_years": (
                    None if butterfly_year is None else input_test_year - butterfly_year
                ),
                "butterfly_match_type": str(butterfly_meta.get("matchType") or ""),
                "butterfly_first_basis": ";".join(
                    sorted({str(x.get("basisOfRecord") or "") for x in butterfly_rec})
                ),
                "butterfly_first_dataset_count": len(
                    {str(x.get("datasetKey") or "") for x in butterfly_rec}
                ),
            }
        )

    dated = [
        r for r in cell_rows
        if r["earliest_documented_host_year"] is not None
        and r["historical_butterfly_first_record_year"] is not None
    ]
    lags = [int(r["recorded_host_to_butterfly_lag_years"]) for r in dated]
    trunc = [int(r["butterfly_window_truncation_years"]) for r in dated]
    payload = {
        "schema": "chocho_butterfly_temporal_lag_symmetric_pilot_v0.2",
        "status": "POSTHOC_SYMMETRIC_HISTORICAL_RECORD_WINDOW_PILOT",
        "record_window": [args.start_year, args.end_year],
        "host_region_pairs": len(host_rows),
        "butterfly_region_cells": len(cell_rows),
        "dated_cells": len(dated),
        "positive_host_first_lag_cells": sum(x > 0 for x in lags),
        "zero_lag_cells": sum(x == 0 for x in lags),
        "negative_host_lag_cells": sum(x < 0 for x in lags),
        "median_recorded_host_to_butterfly_lag_years": (
            None if not lags else float(__import__("statistics").median(lags))
        ),
        "median_butterfly_window_truncation_years": (
            None if not trunc else float(__import__("statistics").median(trunc))
        ),
        "cell_results": cell_rows,
        "claim_boundary": [
            "Both plant and butterfly first-record dates are now searched over the same 1800-2025 window.",
            "First GBIF record remains a detection/digitization datum, not an establishment date.",
            "Plant and butterfly collection effort differ through time and among regions.",
            "The nine cells were selected because they had post-2017 butterfly detections, so this pilot diagnoses chronology feasibility rather than estimating an unbiased lag distribution."
        ],
    }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(host_rows[0].keys()))
        writer.writeheader()
        writer.writerows(host_rows)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
