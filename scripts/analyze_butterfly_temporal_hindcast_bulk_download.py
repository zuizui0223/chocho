#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path

from shapely.geometry import Point, shape
from shapely.strtree import STRtree

from audit_butterfly_temporal_hindcast_full_history import null_summary


def load_geometry(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    geoms, codes = [], []
    for feature in payload.get("features", []):
        props = feature.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        raw = feature.get("geometry")
        if code and raw:
            geoms.append(shape(raw))
            codes.append(code)
    if not geoms:
        raise RuntimeError("no WGSRPD3 geometry")
    return geoms, codes


def tree_candidates(tree, geoms, point):
    index_by_id = {id(g): i for i, g in enumerate(geoms)}
    out = []
    for candidate in tree.query(point):
        if hasattr(candidate, "geom_type"):
            out.append(index_by_id[id(candidate)])
        else:
            out.append(int(candidate))
    return out


def map_point(tree, geoms, codes, lon, lat):
    p = Point(float(lon), float(lat))
    matches = []
    for i in tree_candidates(tree, geoms, p):
        if geoms[i].covers(p):
            matches.append(codes[i])
    return None if not matches else sorted(set(matches))[0]


def pick_occurrence_member(zf: zipfile.ZipFile):
    names = [n for n in zf.namelist() if not n.endswith("/")]
    preferred = [n for n in names if Path(n).name.lower() in {"occurrence.txt", "occurrence.csv"}]
    candidates = preferred or [n for n in names if Path(n).suffix.lower() in {".txt", ".csv", ".tsv"}]
    if not candidates:
        raise RuntimeError(f"no tabular occurrence member in GBIF ZIP: {names[:20]}")
    return max(candidates, key=lambda n: zf.getinfo(n).file_size)


def lookup_field(fieldnames, options):
    by_lower = {str(x).lower(): x for x in fieldnames or []}
    for option in options:
        if option.lower() in by_lower:
            return by_lower[option.lower()]
    raise RuntimeError(f"missing one of fields {options}; got {fieldnames}")


def read_history(download_zip: Path, taxon_map: dict[str, str], candidate_keys: set[tuple[str, str]], level3_geojson: Path):
    geoms, codes = load_geometry(level3_geojson)
    tree = STRtree(geoms)
    species_set = {s for s, _ in candidate_keys}
    earliest = {}
    all_earliest = {}
    first_basis = defaultdict(set)
    first_datasets = defaultdict(set)
    scanned = mapped = relevant = 0

    with zipfile.ZipFile(download_zip) as zf:
        member = pick_occurrence_member(zf)
        with zf.open(member) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
            header = text.readline()
            if not header:
                raise RuntimeError("empty GBIF occurrence file")
            delimiter = "\t" if header.count("\t") >= header.count(",") else ","
            reader = csv.DictReader(io.StringIO(header + text.read()), delimiter=delimiter)
            fields = reader.fieldnames or []
            key_field = lookup_field(fields, ["speciesKey", "taxonKey", "acceptedTaxonKey"])
            lat_field = lookup_field(fields, ["decimalLatitude"])
            lon_field = lookup_field(fields, ["decimalLongitude"])
            year_field = lookup_field(fields, ["year"])
            basis_field = next((f for f in fields if f.lower() == "basisofrecord"), None)
            dataset_field = next((f for f in fields if f.lower() == "datasetkey"), None)

            for row in reader:
                scanned += 1
                raw_key = str(row.get(key_field) or "").strip()
                species = taxon_map.get(raw_key)
                if species is None or species not in species_set:
                    continue
                try:
                    year = int(float(str(row.get(year_field) or "").strip()))
                    lat = float(row[lat_field])
                    lon = float(row[lon_field])
                except Exception:
                    continue
                if not (1800 <= year <= 2025):
                    continue
                code = map_point(tree, geoms, codes, lon, lat)
                if code is None:
                    continue
                mapped += 1
                key = (species, code)
                old_all = all_earliest.get(key)
                if old_all is None or year < old_all:
                    all_earliest[key] = year
                if key not in candidate_keys:
                    continue
                relevant += 1
                current = earliest.get(key)
                if current is None or year < current:
                    earliest[key] = year
                    first_basis[key].clear()
                    first_datasets[key].clear()
                if earliest.get(key) == year:
                    if basis_field:
                        first_basis[key].add(str(row.get(basis_field) or ""))
                    if dataset_field:
                        first_datasets[key].add(str(row.get(dataset_field) or ""))
    return earliest, all_earliest, first_basis, first_datasets, {"scanned_records": scanned, "mapped_records": mapped, "candidate_region_records": relevant, "all_species_region_first_records": len(all_earliest)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download-zip", type=Path, required=True)
    ap.add_argument("--taxon-map-json", type=Path, required=True)
    ap.add_argument("--candidates-csv", type=Path, required=True)
    ap.add_argument("--level3-geojson", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=99999)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-all-first-csv", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    taxon_payload = json.loads(args.taxon_map_json.read_text(encoding="utf-8"))
    taxon_map = {str(k): str(v) for k, v in taxon_payload["taxon_key_to_species"].items()}
    with args.candidates_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1477:
        raise RuntimeError(f"expected 1477 frozen candidates; got {len(rows)}")
    candidate_keys = {(r["species"], r["wgsrpd3_code"]) for r in rows}
    if len(candidate_keys) != len(rows):
        raise RuntimeError("candidate key duplication")

    earliest, all_earliest, first_basis, first_datasets, ingest = read_history(
        args.download_zip, taxon_map, candidate_keys, args.level3_geojson
    )
    for r in rows:
        key = (r["species"], r["wgsrpd3_code"])
        first = earliest.get(key)
        r["historical_first_record_year"] = "" if first is None else first
        r["historical_first_basis"] = ";".join(sorted(x for x in first_basis.get(key, set()) if x))
        r["historical_first_dataset_count"] = len({x for x in first_datasets.get(key, set()) if x})
        r["historical_record_pre2018"] = int(first is not None and first <= 2017)
        r["full_history_eligible"] = int(first is None or first >= 2018)
        r["full_history_outcome"] = int(first is not None and 2018 <= first <= 2025)

    clean = [r for r in rows if int(r["full_history_eligible"])]
    treat = [r for r in clean if int(r["treatment"])]
    control = [r for r in clean if not int(r["treatment"])]
    th = sum(int(r["full_history_outcome"]) for r in treat)
    ch = sum(int(r["full_history_outcome"]) for r in control)
    tr = th / len(treat) if treat else None
    cr = ch / len(control) if control else None

    null_species = null_summary(clean, ["species"], args.permutations, args.seed)
    null_level1 = null_summary(clean, ["species", "level1"], args.permutations, args.seed)
    null_full = null_summary(
        clean,
        ["species", "level1", "baseline_effort_bin", "test_effort_bin", "distance_bin"],
        args.permutations,
        args.seed,
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_all_first_csv.open("w", newline="", encoding="utf-8") as handle:
        w_all = csv.DictWriter(handle, fieldnames=["species", "wgsrpd3_code", "historical_first_record_year"])
        w_all.writeheader()
        for (species, code), year in sorted(all_earliest.items()):
            w_all.writerow({
                "species": species,
                "wgsrpd3_code": code,
                "historical_first_record_year": year,
            })

    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    payload = {
        "schema": "chocho_butterfly_temporal_resource_hindcast_bulk_gbif_v0.1",
        "status": "FULL_HISTORY_BULK_DOWNLOAD_TEMPORAL_HINDCAST",
        "history_window": [1800, 2025],
        "download": taxon_payload.get("download", {}),
        "ingest": ingest,
        "history_audit": {
            "candidate_cells": len(rows),
            "cells_with_pre2018_record_removed": sum(int(r["historical_record_pre2018"]) for r in rows),
            "treatment_pre2018_removed": sum(int(r["historical_record_pre2018"]) and int(r["treatment"]) for r in rows),
            "control_pre2018_removed": sum(int(r["historical_record_pre2018"]) and not int(r["treatment"]) for r in rows),
            "clean_cells": len(clean),
            "clean_treatment_cells": len(treat),
            "clean_control_cells": len(control),
        },
        "observed": {
            "treatment_first_records_2018_2025": th,
            "treatment_rate": tr,
            "control_first_records_2018_2025": ch,
            "control_rate": cr,
            "raw_risk_difference": None if tr is None or cr is None else tr - cr,
            "raw_risk_ratio": None if tr is None or not cr else tr / cr,
        },
        "null_ladder": {
            "species_only": null_species,
            "species_plus_level1": null_level1,
            "full_protocol_match": null_full,
        },
        "decision": {
            "placement_signal_beyond_level1": bool(null_level1.get("p_high", 1) <= 0.05),
            "placement_signal_beyond_full_match": bool(null_full.get("p_high", 1) <= 0.05),
        },
        "claim_boundary": [
            "The capped 53,434-record snapshot defines the outcome-blind candidate and effort strata only.",
            "Temporal eligibility and outcome are rebuilt from one bulk GBIF download covering all focal taxa over 1800-2025.",
            "A first GBIF record is a first documented detection, not a colonization or establishment date.",
            "Contemporary WCVP introduced ranges remain time-invariant and do not establish host presence by 2017."
        ],
    }
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
