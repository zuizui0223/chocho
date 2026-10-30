#!/usr/bin/env python3
"""Baseline-only audit of all120 unseen near-continent / far-WPI candidates.

No 2018–2025 query is constructed. Baseline recording coverage is neither survey
effort nor butterfly absence. The observed96 panel is never silently enlarged.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import acquire_island_butterfly_panel_20261029 as pilot
import acquire_expanded_island_panel_20261029 as expansion

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/"data/exploratory/island_baseline_20261029"
PLAN = ROOT/"provenance/archive/host_hindcast/ISLAND_BASELINE_HOST_PLAN_20261029.md"
PLAN_SHA = "67c25d1f28ff9f98def1122fc0fa0b4782062e08cf3dae5497cd3806795a82e6"
WINDOW = (2010, 2017)


def select_candidates(rows, known_keys):
    if len({r["island_key"] for r in rows}) != len(rows):
        raise ValueError("Duplicate original island identity")
    return sorted([r for r in rows if r["island_key"] not in known_keys
                   and 100 <= float(r["area_km2"]) < 20000
                   and -55 <= float(r["representative_lat"]) <= 55
                   and not int(r["antarctica_could_be_nearer"])
                   and float(r["mainland_coast_gap_upper_km"]) < 100
                   and float(r["nearest_wpi_coast_distance_lower_km"]) > 100],
                  key=lambda r: r["island_key"])


def baseline_payload(wkt, offset=0):
    return pilot.payload(wkt, WINDOW, offset)


def baseline_quality(result):
    if result["status"] != "COMPLETE":
        return {"query_status": result["status"], "records": None,
                "species_assigned_records": None, "recorded_species_bins": None,
                "species_assignment_fraction": None, "positive_years": None,
                "eligible": None, "exclusion_reasons": "query_not_complete"}
    total = result["count"]
    assigned = sum(result["facets"]["SPECIES_KEY"].values())
    years = len(result["facets"]["YEAR"])
    fraction = assigned/total if total else None
    reasons = []
    if assigned < 20:
        reasons.append("fewer_than20_assigned_records")
    if fraction is None or fraction < .9:
        reasons.append("species_assignment_below90pct_or_no_records")
    if years < 3:
        reasons.append("fewer_than3_recorded_years")
    return {"query_status": "COMPLETE", "records": total,
            "species_assigned_records": assigned,
            "recorded_species_bins": len(result["facets"]["SPECIES_KEY"]),
            "species_assignment_fraction": fraction, "positive_years": years,
            "eligible": not reasons, "exclusion_reasons": "|".join(reasons)}


def decision(rows, requested):
    complete = sum(r["query_status"] == "COMPLETE" for r in rows)
    eligible = sum(r["eligible"] is True for r in rows)
    if len(rows) != requested or complete != requested:
        status = "INCOMPLETE_BASELINE_AUDIT"
    elif eligible < 3:
        status = "NO_LATE_ACQUISITION_BASELINE_CELL_INSUFFICIENT"
    else:
        status = "BASELINE_CELL_SCREEN_CLEARED_COHORT_NOT_YET_FIXED"
    return {"status": status, "requested_islands": requested, "completed_islands": complete,
            "eligible_islands": eligible, "late_outcomes_requested": 0,
            "interpretation": "Baseline documentation only; not survey effort or absence. Clearing3 does not satisfy the full model or validate independent replication."}


class BaselineClient(expansion.SourceClient):
    def request(self, data, folder, prefix):
        # The shared client makes at most2 attempts per request invocation.
        if self.new_requests + 2 > 360:
            raise RuntimeError("Bounded360-attempt baseline budget reached; no zero inferred")
        return super().request(data, folder, prefix)


def acquire(row, wkt, cache, client):
    folder = cache/row["island_key"]/"2010_2017"
    folder.mkdir(parents=True, exist_ok=True)
    receipts = []
    started = pilot.utc()
    def fetch(offset, recheck=False):
        prefix = f"page_{offset:05d}"+("_recheck" if recheck else "")
        request_bytes = pilot.canonical(baseline_payload(wkt, offset))
        if (folder/(prefix+"_response.json")).exists():
            data, receipt = expansion.read_cached(folder, prefix, request_bytes)
        else:
            (folder/(prefix+"_request.json")).write_bytes(request_bytes)
            body = client.request(request_bytes, folder, prefix)
            receipt = {"request_sha256": pilot.sha(request_bytes), "response_sha256": pilot.sha(body),
                       "retrieved_utc": pilot.utc(), "facet_offset": offset,
                       "is_first_page_recheck": recheck}
            (folder/(prefix+"_response.json")).write_bytes(body)
            (folder/(prefix+"_receipt.json")).write_bytes(pilot.canonical(receipt))
            data = json.loads(body)
        receipts.append(receipt)
        return data
    result = pilot.query_safely(lambda: pilot.collect_facets(fetch))
    return {"island_key": row["island_key"], "window": "2010_2017",
            "geometry_sha256": row["geometry_sha256"], "started_utc": started,
            "finished_utc": pilot.utc(), "source_pages": receipts, **result}


def publish(output, lock, results, client, finished):
    selected = {r["island_key"]: r for r in lock["selection"]["selected"]}
    rows = [{"island_key": r["island_key"], "window": "2010_2017",
             "area_km2": selected[r["island_key"]]["area_km2"],
             "cluster_500": selected[r["island_key"]]["cluster_500"], **baseline_quality(r)}
            for r in sorted(results, key=lambda r: r["island_key"])]
    decision_result = decision(rows, len(selected))
    complete = [r for r in rows if r["query_status"] == "COMPLETE"]
    summary = {**decision_result, "complete_zero_record_islands": sum(r["records"] == 0 for r in complete),
               "complete_positive_record_islands": sum(r["records"] > 0 for r in complete),
               "total_records": sum(r["records"] for r in complete),
               "atleast3_species_islands": sum(r["recorded_species_bins"] >= 3 for r in complete),
               "eligible_geographic_clusters": len({r["cluster_500"] for r in rows if r["eligible"]}),
               "source_requests_this_invocation": client.new_requests}
    report = {"schema": "chocho_unseen_deficient_cell_baseline_20261029_v1",
              "acquisition_status": "FINISHED" if finished else "RUNNING", "updated_utc": pilot.utc(),
              "selection_sha256": lock["selection_sha256"], "selection_locked_utc": lock["locked_before_queries_utc"],
              "plan_sha256": PLAN_SHA, "taxonomy": {"checklistKey": pilot.CHECKLIST, "taxonKey": "5G9"},
              "endpoint": pilot.ENDPOINT, "window": list(WINDOW),
              "query_contract": baseline_payload("<exact original oriented WKT, cached only>"),
              "interpretation": "New baseline-only recording feasibility audit; all120 candidates were unseen in the prior96. No late data requested; no outcome-based replacement; no biological absences or port effects inferred.",
              "summary": summary, "results": results}
    (output/"unseen_island_baseline_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    if rows:
        pilot.write_csv(output/"unseen_island_baseline_counts.csv", rows, list(rows[0]))
    for field in pilot.FIELDS:
        values = [{"island_key": r["island_key"], "facet_key": key, "record_count": n}
                  for r in results for key,n in sorted(r.get("facets", {}).get(field, {}).items())]
        pilot.write_csv(output/("unseen_island_baseline_"+field.lower()+"_facets.csv"), values,
                        ["island_key", "facet_key", "record_count"])
    return report


def run(source_cache, cache, output, prepare_only=False):
    from shapely import orient_polygons, to_wkt, get_num_coordinates
    from measure_island_port_geography_20261029 import source_landmasses, derived_key
    from analyze_expanded_island_recordings_20261029 import fixed_clusters
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    if pilot.sha(PLAN.read_bytes()) != PLAN_SHA:
        raise ValueError("Frozen continuation plan changed")
    known_path = ROOT/"data/exploratory/island_expansion_20261029/expanded_island_panel_selection.json"
    known = json.loads(known_path.read_bytes())
    if pilot.sha(known_path.read_bytes()) != "5ccc868b89973b8142a7fdb4b3014c4394f134596ae2ef58e1179eaab68f5243":
        raise ValueError("Known96 selection changed")
    exposure = expansion.BASE/"island_geography_exposures.csv"
    if pilot.sha(exposure.read_bytes()) != known["selection"]["geography_exposures_sha256"]:
        raise ValueError("Geographic exposure source changed")
    with exposure.open() as f:
        all_rows = list(csv.DictReader(f))
    land = source_landmasses(source_cache)
    geometries = {derived_key(r.parent_source_ID): r.geometry for _,r in land[land.area_km2 >= 5].iterrows()}
    for row in all_rows:
        point = geometries[row["island_key"]].representative_point()
        row.update(representative_lon=point.x, representative_lat=point.y)
    selected = select_candidates(all_rows, {r["island_key"] for r in known["selection"]["selected"]})
    if len(selected) != 120 or sum(float(r["area_km2"]) >= 1000 for r in selected) != 3:
        raise ValueError("Frozen120 candidate census changed")
    clusters, cluster_details = fixed_clusters(selected, 500)
    wkts, public = {}, []
    for row in selected:
        geom = orient_polygons(geometries[row["island_key"]], exterior_cw=False)
        if not geom.is_valid or geom.is_empty or geom.geom_type not in ("Polygon", "MultiPolygon"):
            raise ValueError("Invalid exact source geometry; no substitute")
        wkt = to_wkt(geom, rounding_precision=-1, trim=True)
        wkts[row["island_key"]] = wkt
        public.append({**row, "geometry_sha256": pilot.sha(wkt.encode()),
                       "geometry_vertices": int(get_num_coordinates(geom)), "cluster_500": clusters[row["island_key"]]})
    selection = {"schema": "chocho_unseen_baseline_fixed_roster_20261029_v1",
                 "plan_sha256": PLAN_SHA, "known96_selection_file_sha256": pilot.sha(known_path.read_bytes()),
                 "geography_exposures_sha256": pilot.sha(exposure.read_bytes()),
                 "gshhg_source": known["selection"]["gshhg_source"],
                 "wpi_source": known["selection"]["wpi_source"],
                 "rule": "All unseen islands100<=area<20000km2, representative latitude within55degrees, noAntarcticambiguity, mainlandupper<100km, WPIlower>100km; excludeknown96; no outcomes used.",
                 "window": list(WINDOW), "selected": public, "fixed500km_clusters": cluster_details}
    lock_path = cache/"baseline_selection_lock.json"
    if lock_path.exists():
        lock = json.loads(lock_path.read_bytes())
        if lock["selection"] != selection or lock["selection_sha256"] != pilot.sha(pilot.canonical(selection)):
            raise ValueError("Baseline selection changed after lock")
    else:
        lock = {"locked_before_queries_utc": pilot.utc(), "selection_sha256": pilot.sha(pilot.canonical(selection)),
                "selection": selection}
        lock_path.write_bytes(pilot.canonical(lock))
    (output/"unseen_island_baseline_selection.json").write_text(json.dumps(lock, indent=2, sort_keys=True)+"\n")
    pilot.write_csv(output/"unseen_island_baseline_selection.csv", public, list(public[0]))
    print(json.dumps({"stage": "BASELINE_ROSTER_LOCKED", "islands": len(public),
                      "clusters500": len(cluster_details), "selection_sha256": lock["selection_sha256"],
                      "locked_utc": lock["locked_before_queries_utc"]}), flush=True)
    if prepare_only:
        return lock
    client = BaselineClient()
    results = []
    for row in public:
        result = acquire(row, wkts[row["island_key"]], cache, client)
        results.append(result)
        publish(output, lock, results, client, False)
        print(json.dumps({"completed": len(results), "total": len(public), "status": result["status"],
                          "records": result["count"], "new_requests": client.new_requests}), flush=True)
    final = publish(output, lock, results, client, True)
    print(json.dumps(final["summary"], indent=2), flush=True)
    return final


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-cache", type=Path, default=Path("/tmp/chocho-island-sources"))
    p.add_argument("--cache-dir", type=Path, default=Path("/tmp/chocho-island-baseline-cache"))
    p.add_argument("--output-dir", type=Path, default=OUTPUT)
    p.add_argument("--prepare-only", action="store_true")
    args = p.parse_args()
    run(args.source_cache, args.cache_dir, args.output_dir, args.prepare_only)
