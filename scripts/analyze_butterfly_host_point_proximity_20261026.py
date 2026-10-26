#!/usr/bin/env python3
"""Source-locked GBIF butterfly / WCVP candidate-host point co-recording audit.

26 pre-existing positive butterfly-country cases, 44 independently fixed
plant-country source taxa. No inferred host use or climate niche entry.
Unredacted coordinates stay in GitHub Actions memory ONLY.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qsl, urlencode, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "chocho_exact26_butterfly_plant_point_20261026_v01"
BASIS = {"HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
         "MACHINE_OBSERVATION", "OCCURRENCE"}
MIN_WAIT_SECONDS = 1.05
PER_QUERY_CAP = 5000
TOTAL_RECORD_CAP = 12000
PAGE_LIMIT = 300
KM_EARTH = 6371.0088
JUVENILE = ("larva", "caterpillar", "egg", "pupa", "nymph")
USE_WORDS = ("feeding", "feeds on", "feeding on", "host plant",
             "larval host", "ovipos", "eating", "on leaves of")


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def source_inputs(inventory_path, butterfly_path, plant_path):
    invraw = Path(inventory_path).read_bytes()
    braw = Path(butterfly_path).read_bytes()
    praw = Path(plant_path).read_bytes()
    inv, butterfly, plant = [json.loads(b) for b in (invraw, braw, praw)]
    if inv.get("schema") != "chocho_all26_positive_country_source_hosts_v1":
        raise ValueError("Wrong original 26-case accepted-host roster")
    if butterfly.get("schema") != "chocho_subcountry_48_positive_GBIF_GADM1_lifestage_facet_v1":
        raise ValueError("Wrong 48-case GBIF butterfly source")
    if plant.get("schema") != "chocho_all26_gbif_plant_gadm1_feasibility_20261023_v1":
        raise ValueError("Wrong plant/source GBIF facet snapshot")
    if (inv["frozen_added_host_cells"], inv["candidate_butterfly_country_plant_links"],
            inv["distinct_source_plant_ids"], inv["distinct_source_plant_country_queries"]
            if "distinct_source_plant_country_queries" in inv
            else inv["distinct_plant_country_source_queries"]) != (26, 46, 24, 44):
        raise ValueError("Frozen cohort margins changed")
    if sha(invraw) != plant["source"]["source_inventory_sha256"] or \
       sha(braw) != plant["source"]["source_positive48_sha256"]:
        raise ValueError("Frozen initial source-body checksums changed")
    if plant["source"]["source_original_zip_sha256"] != \
            "408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323":
        raise ValueError("Original accepted-host source ZIP identity changed")
    if plant.get("plant_country_GADM_facets_complete") != 44 or \
       butterfly.get("gadm1_facets_complete") != 48:
        raise ValueError("Original source query cohort incomplete")
    bs = {(r["species"], r["country"]): r for r in butterfly["observations"]}
    if len(bs) != 48:
        raise ValueError("Duplicated original GBIF butterfly source records")
    plants = {(r["wcvp_accepted_id"], r["country"]): r
              for r in plant["plant_country_results"]}
    if len(plants) != 44:
        raise ValueError("Missing unique source accepted plant country queries")
    cases = {}
    for r in inv["all_48_positive_country_cells"]:
        if not r["added_host_resource_country"]:
            continue
        key = (r["butterfly"], r["country"])
        if key not in bs or bs[key]["added_resource_country"] != 1:
            raise ValueError("Butterfly accepted host opportunity mismatch")
        source_host = r["candidate_accepted_plants"]
        if not source_host or len(source_host) != r["candidate_source_plant_count"]:
            raise ValueError("Source host candidates changed")
        hosts = []
        for host in source_host:
            pk = (host["wcvp_accepted_id"], r["country"])
            if pk not in plants:
                raise ValueError("Missing exact accepted source plant record")
            p = plants[pk]
            if p["plant_name"] not in host["accepted_names"]:
                raise ValueError("Source plant accepted name changed")
            if p["gbif_taxonomy"]["status"] != "EXACT_ACCEPTED_GBIF_PLANT":
                raise ValueError("Unresolved source plant taxonomy")
            hosts.append(pk)
        if key in cases:
            raise ValueError("Duplicate butterfly-country case")
        cases[key] = {"row": bs[key], "hosts": hosts}
    if len(cases) != 26 or len({s for s, c in cases}) != 16 or \
       len({c for s, c in cases}) != 11:
        raise ValueError("Original 26-cases/16 species/11 countries changed")
    if sum(len(r["hosts"]) for r in cases.values()) != 46:
        raise ValueError("Source butterfly-host edges changed")
    keys = {k for r in cases.values() for k in r["hosts"]}
    if keys != set(plants):
        raise ValueError("Original 44 accepted plant-country pairs changed")
    return cases, plants, {
        "frozen_case_count": len(cases),
        "frozen_plant_country_queries": len(plants),
        "frozen_unique_host_species": 24,
        "frozen_butterfly_source_count": sum(r["row"]["record_count"] for r in cases.values()),
        "frozen_plant_source_count": sum(r["indexed_country_plant_records"] for r in plants.values()),
        "original_inventory_sha256": sha(invraw),
        "original_butterfly_GADM_sha256": sha(braw),
        "original_plant_GADM_sha256": sha(praw),
    }


def records_url(source_url, wanted_taxon, wanted_country, offset=0):
    """Only replace facet lookup with documented paginated record request."""
    parsed = urlparse(source_url)
    if parsed.scheme != "https" or parsed.netloc != "api.gbif.org" or \
       parsed.path != "/v1/occurrence/search":
        raise ValueError("Changed original GBIF source endpoint")
    pairs = parse_qsl(parsed.query, keep_blank_values=True)
    by = defaultdict(list)
    for k, v in pairs:
        by[k].append(v)
    expected = {
        "taxonKey", "country", "year", "hasCoordinate",
        "hasGeospatialIssue", "occurrenceStatus", "basisOfRecord", "limit"
    }
    if not expected <= set(by) or set(by) - expected - {"facet", "facetLimit"}:
        raise ValueError("Non-frozen source query parameters")
    if by["taxonKey"] != [str(wanted_taxon)] or \
       by["country"] != [wanted_country] or by["year"] != ["2010,2025"] or \
       by["hasCoordinate"] != ["true"] or by["hasGeospatialIssue"] != ["false"] or \
       by["occurrenceStatus"] != ["PRESENT"] or \
       set(by["basisOfRecord"]) != BASIS or len(by["basisOfRecord"]) != len(BASIS):
        raise ValueError("Original query taxon/filter identity has changed")
    if by["limit"] != ["0"] or offset < 0 or offset + PAGE_LIMIT > 100000:
        raise ValueError("Unexpected original facet or GBIF pagination")
    args = [
        ("taxonKey", str(wanted_taxon)), ("country", wanted_country),
        ("year", "2010,2025"), ("hasCoordinate", "true"),
        ("hasGeospatialIssue", "false"), ("occurrenceStatus", "PRESENT"),
    ]
    args += [("basisOfRecord", basis) for basis in sorted(BASIS)]
    args += [("limit", str(PAGE_LIMIT)), ("offset", str(offset))]
    return "https://api.gbif.org/v1/occurrence/search?" + urlencode(args)


class Throttle:
    def __init__(self):
        self.last = None
        self.requests = 0
        self.retries = 0

    def fetch(self, url):
        for attempt in range(4):
            if self.last is not None:
                sleep = MIN_WAIT_SECONDS - (time.monotonic() - self.last)
                if sleep > 0:
                    time.sleep(sleep)
            self.last = time.monotonic()
            self.requests += 1
            try:
                req = Request(url, headers={
                    "User-Agent": "chocho-butterfly-point-science/0.1 (+https://github.com/zuizui0223/chocho)",
                    "Accept": "application/json",
                })
                with urlopen(req, timeout=40) as response:
                    raw = response.read()
                    if response.status != 200:
                        raise RuntimeError("Non-200 GBIF record response")
                data = json.loads(raw)
                if not isinstance(data, dict) or "count" not in data or \
                   "results" not in data:
                    raise ValueError("Bad GBIF occurrence response schema")
                return data, sha(raw)
            except (HTTPError, URLError, TimeoutError) as exc:
                if isinstance(exc, HTTPError) and exc.code not in \
                   (429, 500, 502, 503, 504):
                    raise
                if attempt == 3:
                    raise RuntimeError("Repeated GBIF transport failure") from exc
                self.retries += 1
                time.sleep(min(2 ** (attempt + 1), 12))
        raise RuntimeError("GBIF fetch exhausted")


def valid_point(record, expected_species_key, country):
    """None = discarded index; missing uncertainty stays missing, not 0."""
    try:
        lat, lon = float(record["decimalLatitude"]), float(record["decimalLongitude"])
        year = int(record["year"])
    except (KeyError, ValueError, TypeError):
        return None, "missing_coordinate_or_year"
    if not (math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90
            and -180 <= lon <= 180) or (lat == 0 and lon == 0) or \
            not 2010 <= year <= 2025:
        return None, "invalid_coordinate_or_year"
    if record.get("countryCode") not in (country, None):
        return None, "country_mismatch"
    key = record.get("speciesKey")
    if key is None:
        key = record.get("acceptedTaxonKey")
    if key is None or int(key) != expected_species_key:
        return None, "accepted_species_key_mismatch"
    uncertainty = record.get("coordinateUncertaintyInMeters")
    try:
        uncertainty = float(uncertainty) if uncertainty is not None else None
    except (TypeError, ValueError):
        uncertainty = None
    if uncertainty is not None and (not math.isfinite(uncertainty) or uncertainty < 0):
        uncertainty = None
    stage = str(record.get("lifeStage") or "").strip().lower()
    juvenile = any(k in stage for k in JUVENILE)
    text = " ".join(str(record.get(f) or "") for f in (
        "associatedTaxa", "occurrenceRemarks", "dynamicProperties", "behavior",
    )).lower()
    return {
        "lat": lat, "lon": lon, "year": year, "uncertainty": uncertainty,
        "strict": uncertainty is not None and uncertainty <= 1000,
        "stage_juvenile": juvenile, "annotation_text": text if juvenile else "",
        "dataset": record.get("datasetKey"),
    }, "valid"


def fetch_pair(q, throttle, remaining):
    expected_key, country = q["gbif_species_key"], q["country"]
    source_count = q["original_faceted_index_count"]
    url = records_url(q["source_url"], expected_key, country, 0)
    receipts = []
    discarded = Counter()
    seen = set()
    out = []
    try:
        first, digest = throttle.fetch(url)
        count = int(first["count"])
        if count < 0:
            raise ValueError("Negative GBIF source count")
        if count > PER_QUERY_CAP or count > remaining:
            return {"status": "OVER_BOUNDED_ACQUISITION_CAP",
                    "count": count, "old_count": source_count,
                    "reason": "source volume exceeds precommitted cap"}, []
        offsets = list(range(0, count, PAGE_LIMIT)) or [0]
        for offset in offsets:
            if offset == 0:
                page, body_sha = first, digest
            else:
                page, body_sha = throttle.fetch(
                    records_url(q["source_url"], expected_key, country, offset)
                )
            if int(page["count"]) != count:
                raise ValueError("GBIF live count drifted during pagination")
            returned = page["results"]
            if not isinstance(returned, list):
                raise ValueError("GBIF source results not a list")
            if len(returned) != min(PAGE_LIMIT, count - offset):
                raise ValueError("GBIF result page short or changed during pagination")
            receipts.append({"offset": offset, "rows": len(returned),
                             "sha256": body_sha})
            for rec in returned:
                identity = rec.get("key")
                if identity is None or str(identity) in seen:
                    raise ValueError("Duplicate or missing GBIF occurrence key")
                seen.add(str(identity))
                p, status = valid_point(rec, expected_key, country)
                discarded[status] += 1
                if p:
                    out.append(p)
        if len(seen) != count:
            raise ValueError("Indexed country count vs retrieved source mismatch")
        return {
            "status": "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT",
            "count": count,
            "old_count": source_count,
            "delta_vs_previous_facet_count": count - source_count,
            "valid_points": len(out), "excluded": dict(discarded),
            "pages": receipts, "source_query_with_public_taxon_and_country": url,
        }, out
    except (ValueError, HTTPError, URLError, RuntimeError) as exc:
        return {"status": "INCOMPLETE_GBIF_SOURCE",
                "old_count": source_count,
                "reason": str(exc)[:280],
                "pages_retrieved_before_error": receipts}, []


def distance_km(a, b):
    a1, a2 = math.radians(a["lat"]), math.radians(b["lat"])
    phi = a2 - a1
    theta = math.radians(b["lon"] - a["lon"])
    h = math.sin(phi / 2) ** 2 + math.cos(a1) * math.cos(a2) * \
        math.sin(theta / 2) ** 2
    return KM_EARTH * 2 * math.asin(min(1., math.sqrt(max(0., h))))


def evaluate(butter, plants, host_names):
    """Complete-candidate source only; all counts are indexed point proximity."""
    spans = {"three_years": 3, "same_year": 0, "all_years": 100}
    radii = {"1km": 1., "5km": 5., "25km": 25.}
    out = {"butterfly_points": len(butter),
           "candidate_plant_points_per_case": len(plants),
           "butterfly_juvenile_indices": sum(x["stage_juvenile"] for x in butter),
           "butterfly_strict_uncertainty_1km": sum(x["strict"] for x in butter),
           "plant_strict_uncertainty_1km": sum(x["strict"] for x in plants)}
    for window, years in spans.items():
        out[window] = {distance: 0 for distance in radii}
    out["strict_known_uncertainty_1km"] = {
        "denominator": sum(x["strict"] for x in butter),
        "5km_3years": 0,
    }
    out["literal_juvenile_host_action_annotations"] = 0
    out["butterfly_index_bins_0p01degree"] = len(
        {(round(x["lat"], 2), round(x["lon"], 2)) for x in butter}
    )
    out["plant_index_bins_0p01degree"] = len(
        {(round(x["lat"], 2), round(x["lon"], 2)) for x in plants}
    )
    out["butterfly_with_any_year_matched_plant"] = 0
    for b in butter:
        if b["stage_juvenile"] and any(
            h.lower() in b["annotation_text"] for h in host_names
        ) and any(w in b["annotation_text"] for w in USE_WORDS):
            out["literal_juvenile_host_action_annotations"] += 1
        minima = {}
        for window, years in spans.items():
            pool = [p for p in plants if abs(b["year"] - p["year"]) <= years]
            if pool:
                minima[window] = min(distance_km(b, p) for p in pool)
            else:
                minima[window] = None
            for rad, km in radii.items():
                if minima[window] is not None and minima[window] <= km:
                    out[window][rad] += 1
        if minima["three_years"] is not None:
            out["butterfly_with_any_year_matched_plant"] += 1
        if b["strict"]:
            pool = [p for p in plants if p["strict"]
                    and abs(b["year"] - p["year"]) <= 3]
            if pool and min(distance_km(b, p) for p in pool) <= 5.:
                out["strict_known_uncertainty_1km"]["5km_3years"] += 1
    return out


def build_queries(cases, plants):
    definitions = {}
    for (name, country), data in sorted(cases.items()):
        row = data["row"]
        key = ("butterfly", name, country)
        definitions[key] = {
            "kind": "butterfly", "label": name, "country": country,
            "gbif_species_key": int(row["taxon_key"]),
            "original_faceted_index_count": int(row["record_count"]),
            "source_url": row["query_url"],
        }
    for (pid, country), row in sorted(plants.items()):
        key = ("plant", pid, country)
        definitions[key] = {
            "kind": "plant", "label": row["plant_name"], "country": country,
            "gbif_species_key": int(row["gbif_taxonomy"]["accepted_species_key"]),
            "original_faceted_index_count": int(row["indexed_country_plant_records"]),
            "source_url": row["url"],
        }
    if len(definitions) != 70:
        raise ValueError("Not exactly 26 butterfly and 44 botanical queries")
    for x in definitions.values():
        records_url(x["source_url"], x["gbif_species_key"], x["country"])
    return definitions


def analyse(cases, plants, source_info, query_results, records, throttle):
    case_rows = []
    for (sp, country), case in sorted(cases.items()):
        butterfly_key = ("butterfly", sp, country)
        host_keys = [("plant", pid, country) for pid, ct in case["hosts"]]
        qs = [query_results[k] for k in [butterfly_key, *host_keys]]
        if any(q["status"] != "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT" for q in qs):
            case_rows.append({"butterfly": sp, "country": country,
                              "status": "INCOMPLETE_REQUIRED_SOURCE",
                              "butterfly_point_count": None,
                              "candidate_plant_index_positive": None})
            continue
        bfly = records[butterfly_key]
        host_records = [p for k in host_keys for p in records[k]]
        host_names = [plants[(pid, country)]["plant_name"]
                      for _, pid, _ in host_keys]
        measurements = evaluate(bfly, host_records, host_names)
        case_rows.append({
            "butterfly": sp, "country": country,
            "source_candidate_plant_species": len(case["hosts"]),
            "source_plants_with_any_valid_point": sum(
                bool(records[k]) for k in host_keys
            ),
            "status": "COMPLETE_OBSERVATIONAL_INDEX_ONLY",
            **measurements,
        })
    valid = [r for r in case_rows if r["status"] == "COMPLETE_OBSERVATIONAL_INDEX_ONLY"]
    index_records = sum(r["butterfly_points"] for r in valid)
    all_pairs = sum(r["candidate_plant_points_per_case"] for r in valid)
    spatial = {window: {
        k: sum(r[window][k] for r in valid) for k in ("1km", "5km", "25km")
    } for window in ("three_years", "same_year", "all_years")}
    scored = [r for r in valid if r["butterfly_points"] >= 5
              and r["source_plants_with_any_valid_point"] > 0]
    source_complete = len(valid) == 26 and \
                      all(q["status"] == "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT"
                          for q in query_results.values())
    source_status = (
        "COMPLETE_LIVE_GBIF_SOURCE_SNAPSHOT" if source_complete
        else "INCOMPLETE_SOURCE_SNAPSHOT"
    )
    science_gate = (
        "DESCRIPTIVE_POINT_INDEX_FEASIBILITY" if source_complete and
        len(scored) >= 5 else "INSUFFICIENT_VALIDATED_POINT_COHORT"
    )
    primary = spatial["three_years"]["5km"]
    out = {
        "schema": SCHEMA,
        "status": source_status, "scientific_information_gate": science_gate,
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "source": source_info,
        "acquisition": {
            "complete_taxon_country_queries": sum(
                q["status"] == "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT"
                for q in query_results.values()),
            "total_queries": len(query_results), "GBIF_network_calls": throttle.requests,
            "rate_limit_retries": throttle.retries,
            "source_rows_total": sum(q.get("count", 0)
                                     for q in query_results.values()
                                     if q["status"] == "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT"),
            "query_receipts": [
                {"kind": k[0], "label": q["label"], "country": q["country"],
                 "gbif_species_key": q["gbif_species_key"], **query_results[k]}
                for k, q in sorted(build_queries_cached.items())
            ] if False else [
                {"kind": key[0], "taxon_label": key[1], "country": key[2],
                 **q} for key, q in sorted(query_results.items())
            ],
            "no_raw_coordinates_preserved": True,
        },
        "analysis": {
            "complete_original_source_positive_cases": len(valid),
            "incomplete_cases": 26-len(valid),
            "butterfly_valid_point_indices": index_records,
            "candidate_plant_valid_point_indices_reused_between_cases": all_pairs,
            "cases_min5_butterfly_points_and_any_candidate_plant": len(scored),
            "cases_with_any_plant_index_record": sum(
                r["candidate_plant_points_per_case"] > 0 for r in valid),
            "butterfly_with_any_3year_plant": sum(
                r["butterfly_with_any_year_matched_plant"] for r in valid),
            "matched_within_1_5_25km_by_source_index_year": spatial,
            "primary_within5km_and3years_indices": primary,
            "primary_fraction_of_coordinate_valid_butterfly_indices": (
                primary / index_records if source_complete and index_records else None
            ),
            "strict_known_uncertainty_1km": {
                "butterfly_denominator": sum(
                    r["strict_known_uncertainty_1km"]["denominator"] for r in valid
                ),
                "within5km_and3years": sum(
                    r["strict_known_uncertainty_1km"]["5km_3years"] for r in valid
                ),
            },
            "literal_juvenile_host_action_annotations": sum(
                r["literal_juvenile_host_action_annotations"] for r in valid
            ),
            "cases": case_rows,
        },
        "red_lines": [
            "Source-positive 26 butterfly-country cases are selected on observed butterfly presence.",
            "GBIF-indexed candidate botanical observations may be cultivated, incomplete or not contemporary foliage.",
            "Point co-recording and same-year proximity do not demonstrate butterfly feeding or establishment.",
            "Source query records are a new live snapshot, not original immutable GBIF occurrence DOI data.",
            "All-valid index coordinates include missing source coordinate uncertainty.",
            "This analysis cannot test actual butterfly climatic niche entry: parent country-level niche source cohort differs.",
        ],
    }
    return out


def run(args):
    cases, plants, meta = source_inputs(
        args.inventory, args.butterflies, args.plants)
    queries = build_queries(cases, plants)
    throttle = Throttle()
    remaining = TOTAL_RECORD_CAP
    result, bodies = {}, {}
    for key, q in sorted(queries.items()):
        r, points = fetch_pair(q, throttle, remaining)
        result[key] = r
        bodies[key] = points
        if r["status"] == "COMPLETE_NEW_LIVE_SOURCE_SNAPSHOT":
            remaining -= r["count"]
    out = analyse(cases, plants, meta, result, bodies, throttle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n")
    x = out["analysis"]
    print(json.dumps({
        "status": out["status"],
        "scientific_information_gate": out["scientific_information_gate"],
        "complete_original_source_positive_cases": x["complete_original_source_positive_cases"],
        "butterfly_valid_point_indices": x["butterfly_valid_point_indices"],
        "butterfly_with_any_3year_plant": x["butterfly_with_any_3year_plant"],
        "primary_within5km_and3years_indices": x["primary_within5km_and3years_indices"],
        "primary_fraction_of_coordinate_valid_butterfly_indices":
            x["primary_fraction_of_coordinate_valid_butterfly_indices"],
        "strict": x["strict_known_uncertainty_1km"],
        "juvenile_host_action_annotations": x["literal_juvenile_host_action_annotations"],
        "source_queries_complete": out["acquisition"]["complete_taxon_country_queries"],
    }, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--butterflies", required=True, type=Path)
    parser.add_argument("--plants", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    run(parser.parse_args())
