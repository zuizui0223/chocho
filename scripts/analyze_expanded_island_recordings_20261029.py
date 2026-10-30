#!/usr/bin/env python3
"""Pilot-informed, frozen analysis of exact-island butterfly recording pools.

Finite-record rarefaction does not estimate true community richness, abundance,
survey effort or absence. Geographic associations are fitted only after the
prespecified support gates; the whole selected roster remains visible.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/exploratory/island_expansion_20261029"
PLAN = ROOT / "provenance/archive/host_hindcast/ISLAND_EXPANSION_PLAN_20261029.md"
PLAN_SHA = "6374bdf4f93b83173ec8c9c94a447fb96bd376963d14ddf8de7827643255319a"
WINDOWS = ("2010_2017", "2018_2025")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rarefied_richness(counts, sample):
    """Exact expected bins in a uniform without-replacement record subsample."""
    counts = list(counts)
    if type(sample) is not int or sample < 1 or any(type(n) is not int or n < 0 for n in counts):
        raise ValueError("Integer nonnegative counts and positive sample size required")
    total = sum(counts)
    if total < sample:
        raise ValueError("Insufficient species-assigned records; never extrapolate")
    terms = []
    for n in counts:
        if not n:
            continue
        if total - n < sample:
            terms.append(1.0)
        else:
            log_absence = math.fsum(math.log1p(-n / (total-j)) for j in range(sample))
            terms.append(-math.expm1(log_absence))
    return math.fsum(terms)


def spherical_distance(a, b):
    lon1, lat1, lon2, lat2 = map(math.radians, (
        float(a["representative_lon"]), float(a["representative_lat"]),
        float(b["representative_lon"]), float(b["representative_lat"])))
    h = math.sin((lat2-lat1)/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 2*6371.0088*math.asin(math.sqrt(min(1, max(0, h))))


def fixed_clusters(roster, radius_km):
    """Freeze connected components over ALL selected units, before filtering."""
    rows = sorted(roster, key=lambda r: r["island_key"])
    if radius_km <= 0 or len({r["island_key"] for r in rows}) != len(rows):
        raise ValueError("Invalid cluster radius or duplicated island identity")
    for r in rows:
        lon, lat = float(r["representative_lon"]), float(r["representative_lat"])
        if not math.isfinite(lon+lat) or not -180 <= lon <= 180 or not -90 <= lat <= 90:
            raise ValueError("Invalid representative point")
    parent = list(range(len(rows)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    distances = {}
    for i in range(len(rows)):
        for j in range(i):
            distance = spherical_distance(rows[i], rows[j])
            distances[i, j] = distance
            if distance <= radius_km:
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(len(rows)):
        groups[find(i)].append(i)
    mapping, details = {}, []
    for indices in groups.values():
        keys = sorted(rows[i]["island_key"] for i in indices)
        key = f"r{int(radius_km)}_{keys[0]}"
        mapping.update({k: key for k in keys})
        diameter = max((distances[max(i,j), min(i,j)] for i in indices for j in indices if i != j), default=0)
        details.append({"cluster_id": key, "islands": len(keys), "diameter_km": diameter,
                        "island_keys": keys})
    return mapping, sorted(details, key=lambda r: r["cluster_id"])


def eligibility_reasons(row, sample, min_years=3, min_human=None, max_dataset=None):
    if row["query_status"] != "COMPLETE":
        return ["query_not_complete"]
    reasons = []
    if row["species_assigned_records"] < sample:
        reasons.append("too_few_species_assigned_records")
    if row["species_assignment_fraction"] is None or row["species_assignment_fraction"] < .9:
        reasons.append("species_assignment_below_90pct")
    if row["positive_years"] < min_years:
        reasons.append("fewer_than_required_recorded_years")
    if min_human is not None and (row["human_observation_fraction"] is None or row["human_observation_fraction"] < min_human):
        reasons.append("human_observation_fraction_below_threshold")
    if max_dataset is not None and (row["dominant_dataset_share"] is None or row["dominant_dataset_share"] > max_dataset):
        reasons.append("dominant_dataset_share_above_threshold")
    return reasons


def design_matrix(rows, with_port=False):
    blocks = sorted({r["longitude_block"] for r in rows})
    names = ["intercept", "log10_area", "absolute_latitude_per10", "log10_one_plus_gap"]
    names += ["block_"+b for b in blocks[1:]]
    if with_port:
        names.append("wpi_within50_vs_beyond100")
    values = []
    for r in rows:
        x = [1, math.log10(float(r["area_km2"])), abs(float(r["representative_lat"]))/10,
             math.log10(1+float(r["mainland_coast_gap_upper_km"]))]
        x += [int(r["longitude_block"] == b) for b in blocks[1:]]
        if with_port:
            x.append(int(r["port_class"] == "within50"))
        values.append(x)
    X = np.asarray(values, dtype=float).reshape(len(rows), len(names))
    if not len(rows):
        return X, names, 0, None
    scaled = X.copy()
    sd = X[:, 1:].std(axis=0)
    if np.any(sd == 0):
        return X, names, int(np.linalg.matrix_rank(X)), None
    scaled[:, 1:] = (X[:, 1:] - X[:, 1:].mean(axis=0)) / sd
    return X, names, int(np.linalg.matrix_rank(scaled)), float(np.linalg.cond(scaled))


def gated_model(rows, sample, with_port=False, cluster_column="cluster_500", influence_only=False):
    gap = Counter(r["gap_class"] for r in rows)
    cells = Counter((r["gap_class"], r["port_class"]) for r in rows)
    regions = defaultdict(set)
    for r in rows:
        regions[r["longitude_block"]].add(r["gap_class"])
    clusters = sorted({r[cluster_column] for r in rows})
    gap_clusters = {g: len({r[cluster_column] for r in rows if r["gap_class"] == g}) for g in ("under100", "atleast500")}
    X, names, rank, condition = design_matrix(rows, with_port)
    n, p, G = len(rows), len(names), len(clusters)
    failures = []
    if n < 30:
        failures.append("fewer_than_30_eligible_islands")
    if any(gap[g] < 10 for g in ("under100", "atleast500")):
        failures.append("fewer_than_10_islands_per_gap_class")
    if sum(len(v) == 2 for v in regions.values()) < 2:
        failures.append("fewer_than_2_longitude_blocks_with_both_gap_classes")
    if G < 12:
        failures.append("fewer_than_12_geographic_clusters")
    if any(v < 5 for v in gap_clusters.values()):
        failures.append("fewer_than_5_clusters_per_gap_class")
    if rank < p or n <= p:
        failures.append("rank_deficient_or_no_residual_degrees_of_freedom")
    if condition is None or condition > 100:
        failures.append("standardized_design_condition_exceeds_100_or_undefined")
    if with_port:
        if any(cells[g, q] < 3 for g in ("under100", "atleast500") for q in ("within50", "beyond100")):
            failures.append("fewer_than_3_islands_per_gap_port_cell")
        if any(len({r[cluster_column] for r in rows if r["port_class"] == q}) < 4 for q in ("within50", "beyond100")):
            failures.append("fewer_than_4_clusters_per_port_class")
    result = {"status": "NOT_FIT_SUPPORT_GATE", "islands": n, "gap_class_counts": dict(gap),
              "gap_port_counts": {g+":"+q: cells[g,q] for g in ("under100", "atleast500") for q in ("within50", "beyond100")},
              "longitude_blocks_with_both_gap_classes": sum(len(v) == 2 for v in regions.values()),
              "geographic_clusters": G, "clusters_per_gap_class": gap_clusters,
              "design_rank": rank, "parameters": p, "standardized_condition_number": condition,
              "gate_failures": failures, "sample_records": sample, "cluster_column": cluster_column,
              "with_port": with_port}
    if failures and not influence_only:
        return result
    if influence_only and (rank < p or n <= p or G < 2 or condition is None or condition > 100 or len(gap) != 2):
        result["status"] = "INFLUENCE_NOT_IDENTIFIABLE"
        return result
    y = np.log([r[f"rarefied_species_{sample}"] for r in rows])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    errors = y - X @ beta
    inverse = np.linalg.inv(X.T @ X)
    meat = np.zeros((p, p))
    for group in clusters:
        ix = np.array([r[cluster_column] == group for r in rows])
        score = X[ix].T @ errors[ix]
        meat += np.outer(score, score)
    correction = G/(G-1) * (n-1)/(n-p)
    covariance = inverse @ meat @ inverse * correction
    se = np.sqrt(np.maximum(0, np.diag(covariance)))
    critical = float(t.ppf(.975, G-1))
    effects = {}
    for i, name in enumerate(names):
        effects[name] = {"log_coefficient": float(beta[i]), "cluster_standard_error": float(se[i]),
                         "ratio": math.exp(float(beta[i])),
                         "ratio_ci95_lower": math.exp(float(beta[i]-critical*se[i])),
                         "ratio_ci95_upper": math.exp(float(beta[i]+critical*se[i]))}
    result.update(status="INFLUENCE_ONLY_NOT_PRIMARY" if influence_only else "FIT_EXPLORATORY_ASSOCIATION",
                  cluster_df=G-1, cr1_correction=correction, effects=effects,
                  r_squared=float(1-errors@errors/np.sum((y-y.mean())**2)) if np.var(y) else None)
    return result


def build_metrics(lock, manifest):
    taxonomy = manifest.get("taxonomy", {})
    if taxonomy.get("checklistKey") != "7ddf754f-d193-4cc9-b351-99906754a03b" or taxonomy.get("taxonKey") != "5G9":
        raise ValueError("Source taxonomy is not the frozen COL XR Papilionoidea concept")
    if (manifest.get("schema") != "chocho_geography_only_true_island_expansion_20261029_v1"
            or manifest.get("endpoint") != "https://api.gbif.org/v1/occurrence/search/predicate"
            or manifest.get("windows") != [[2010, 2017], [2018, 2025]]):
        raise ValueError("Source acquisition contract changed")
    contract = manifest.get("query_contract_geometry_redacted", {})
    filters = contract.get("predicate", {}).get("predicates", [])
    expected_filters = [
        {"type": "equals", "key": "TAXON_KEY", "value": "5G9", "checklistKey": taxonomy["checklistKey"]},
        {"type": "greaterThanOrEquals", "key": "YEAR", "value": "2010"},
        {"type": "lessThanOrEquals", "key": "YEAR", "value": "2017"},
        {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
        {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
        {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "PRESENT"},
        {"type": "within", "geometry": "<exact original oriented WKT, cached only>"},
    ]
    if (filters != expected_filters or contract.get("predicate", {}).get("type") != "and"
            or contract.get("parameters") != {"CHECKLIST_KEY": [taxonomy["checklistKey"]]}
            or set(contract.get("facets", [])) != {"SPECIES_KEY", "BASIS_OF_RECORD", "DATASET_KEY", "YEAR"}
            or contract.get("limit") != 0 or contract.get("facetMultiSelect") is not False):
        raise ValueError("Source occurrence query filters changed")
    selected = lock["selection"]["selected"]
    expected_sha = hashlib.sha256(json.dumps(lock["selection"], sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()
    if lock["selection_sha256"] != expected_sha or manifest["selection_sha256"] != expected_sha:
        raise ValueError("Selection identity changed")
    if manifest["status"] != "FINISHED":
        raise ValueError("Acquisition is still running; do not analyze partial selected roster")
    lookup = {(r["island_key"], r["window"]): r for r in manifest["results"]}
    expected = {(r["island_key"], w) for r in selected for w in WINDOWS}
    if len(lookup) != len(manifest["results"]) or set(lookup) != expected:
        raise ValueError("Island-window roster incomplete or duplicated")
    cluster500, details500 = fixed_clusters(selected, 500)
    cluster1000, details1000 = fixed_clusters(selected, 1000)
    rows = []
    for r in selected:
        area_bin, gap, port = r["selection_stratum"].split(":")
        for window in WINDOWS:
            query = lookup[r["island_key"], window]
            if query["geometry_sha256"] != r["geometry_sha256"]:
                raise ValueError("Query does not match the selected exact geometry")
            row = {k: r[k] for k in ("island_key", "selection_stratum", "selection_rank", "provenance_origin",
                                    "area_km2", "representative_lon", "representative_lat",
                                    "mainland_coast_gap_lower_km", "mainland_coast_gap_upper_km",
                                    "nearest_wpi_coast_distance_lower_km", "nearest_wpi_coast_distance_upper_km")}
            quadrant = min(3, int((float(r["representative_lon"])+180)//90))
            row.update(window=window, query_status=query["status"], area_class=area_bin,
                       gap_class=gap, port_class=port, longitude_block=f"longitude_{quadrant}",
                       cluster_500=cluster500[r["island_key"]], cluster_1000=cluster1000[r["island_key"]])
            for field in ("records", "species_assigned_records", "recorded_species_bins", "species_assignment_fraction",
                          "positive_years", "datasets", "dominant_dataset_share", "human_observation_fraction",
                          "rarefied_species_10", "rarefied_species_20", "rarefied_species_50"):
                row[field] = None
            if query["status"] == "COMPLETE":
                total = query["count"]
                facets = query["facets"]
                if type(total) is not int or total < 0:
                    raise ValueError("Invalid complete recording count")
                for field in ("SPECIES_KEY", "BASIS_OF_RECORD", "DATASET_KEY", "YEAR"):
                    if field not in facets or any(type(n) is not int or n <= 0 for n in facets[field].values()):
                        raise ValueError("Malformed positive facet bins")
                    residual = total-sum(facets[field].values())
                    if residual < 0 or residual != query["unassigned_facet_records"][field]:
                        raise ValueError("Facet residual does not reconcile")
                    if field in ("DATASET_KEY", "YEAR") and residual:
                        raise ValueError("Incomplete source margin")
                assigned = sum(facets["SPECIES_KEY"].values())
                start, end = map(int, window.split("_"))
                if any(not start <= int(year) <= end for year in facets["YEAR"]):
                    raise ValueError("Year bins fall outside requested window")
                row.update(records=total, species_assigned_records=assigned,
                           recorded_species_bins=len(facets["SPECIES_KEY"]),
                           species_assignment_fraction=assigned/total if total else None,
                           positive_years=len(facets["YEAR"]), datasets=len(facets["DATASET_KEY"]),
                           dominant_dataset_share=max(facets["DATASET_KEY"].values(), default=0)/total if total else None,
                           human_observation_fraction=facets["BASIS_OF_RECORD"].get("HUMAN_OBSERVATION", 0)/total if total else None)
                for size in (10, 20, 50):
                    if assigned >= size:
                        row[f"rarefied_species_{size}"] = rarefied_richness(facets["SPECIES_KEY"].values(), size)
            reasons = eligibility_reasons(row, 20)
            row.update(primary_eligible=int(not reasons), primary_exclusion_reasons="|".join(reasons))
            rows.append(row)
    return rows, {"500_km": details500, "1000_km": details1000}


def analyze(rows):
    specifications = [
        ("primary_late_m20", "2018_2025", 20, {}),
        ("late_m10", "2018_2025", 10, {}),
        ("late_m50", "2018_2025", 50, {}),
        ("early_m20", "2010_2017", 20, {}),
        ("late_atleast5_years", "2018_2025", 20, {"min_years": 5}),
        ("late_human_atleast90pct", "2018_2025", 20, {"min_human": .9}),
        ("late_dominant_dataset_atmost90pct", "2018_2025", 20, {"max_dataset": .9}),
    ]
    models, support, cohorts = {}, [], {}
    for name, window, size, filters in specifications:
        cohort = [r for r in rows if r["window"] == window and not eligibility_reasons(r, size, **filters)]
        cohorts[name] = cohort
        models[name] = gated_model(cohort, size)
        models[name]["definition"] = {"window": window, "sample_records": size, **filters}
        for stratum in sorted({r["selection_stratum"] for r in rows}):
            group = [r for r in rows if r["window"] == window and r["selection_stratum"] == stratum]
            eligible = [r for r in cohort if r["selection_stratum"] == stratum]
            support.append({"analysis": name, "window": window, "stratum": stratum, "selected_islands": len(group),
                            "complete_windows": sum(r["query_status"] == "COMPLETE" for r in group),
                            "zero_record_windows": sum(r["records"] == 0 for r in group),
                            "atleast3_recorded_species": sum((r["recorded_species_bins"] or 0) >= 3 for r in group),
                            "eligible_islands": len(eligible),
                            "median_rarefied_species": float(np.median([r[f"rarefied_species_{size}"] for r in eligible])) if eligible else None})
    primary = cohorts["primary_late_m20"]
    models["secondary_port_late_m20"] = gated_model(primary, 20, with_port=True)
    models["primary_1000km_clusters"] = gated_model(primary, 20, cluster_column="cluster_1000")
    influence = []
    if models["primary_late_m20"]["status"] == "FIT_EXPLORATORY_ASSOCIATION":
        for column in ("longitude_block", "cluster_500"):
            for value in sorted({r[column] for r in primary}):
                fit = gated_model([r for r in primary if r[column] != value], 20, influence_only=True)
                influence.append({"omit_column": column, "omit_value": value, **fit})
    coverage = {}
    for window in WINDOWS:
        group = [r for r in rows if r["window"] == window]
        complete = [r for r in group if r["query_status"] == "COMPLETE"]
        coverage[window] = {"selected_islands": len(group), "complete_windows": len(complete),
                            "records": sum(r["records"] for r in complete),
                            "species_assigned_records": sum(r["species_assigned_records"] for r in complete),
                            "zero_record_islands": sum(r["records"] == 0 for r in complete),
                            "atleast3_species_islands": sum(r["recorded_species_bins"] >= 3 for r in complete),
                            "primary_eligible_islands": sum(r["primary_eligible"] for r in group),
                            "exclusion_reason_counts": dict(Counter(reason for r in group for reason in eligibility_reasons(r, 20)))}
    paired = defaultdict(list)
    for r in rows:
        paired[r["island_key"]].append(r)
    coverage["paired"] = {"both_windows_atleast3_species": sum(all((r["recorded_species_bins"] or 0) >= 3 for r in v) for v in paired.values()),
                           "both_windows_primary_eligible": sum(all(r["primary_eligible"] for r in v) for v in paired.values())}
    return {"coverage": coverage, "models": models, "influence": influence}, support


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def run(folder):
    if digest(PLAN) != PLAN_SHA:
        raise ValueError("Frozen analysis plan changed; preserve it and record amendments separately")
    selection_path = folder/"expanded_island_panel_selection.json"
    manifest_path = folder/"expanded_island_panel_manifest.json"
    lock = json.loads(selection_path.read_bytes())
    manifest = json.loads(manifest_path.read_bytes())
    metrics, clusters = build_metrics(lock, manifest)
    result, support = analyze(metrics)
    write_csv(folder/"expanded_island_recording_metrics.csv", metrics)
    write_csv(folder/"expanded_island_analysis_support.csv", support)
    result.update(schema="chocho_frozen_expanded_island_recorded_richness_20261029_v1",
                  plan_sha256=PLAN_SHA, source_selection_sha256=digest(selection_path),
                  source_manifest_sha256=digest(manifest_path),
                  selected_geographic_clusters=clusters,
                  interpretation="Pilot-informed exploratory finite-record assemblage association only. Records are not individuals or visits; empty records are not absences. Rarefaction does not remove selective recording, source bias, spatial habitat/climate differences, or make the purposive panel representative. Longitude blocks are not biological realms; 500/1000km connected components are dependence approximations, not validated independent populations.",
                  output_sha256={name: digest(folder/name) for name in ("expanded_island_recording_metrics.csv", "expanded_island_analysis_support.csv")})
    (folder/"expanded_island_recording_analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print(json.dumps({"coverage": result["coverage"], "models": {k: {f:v[f] for f in ("status", "islands", "geographic_clusters", "gate_failures")} for k,v in result["models"].items()}}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT)
    run(parser.parse_args().data_dir)
