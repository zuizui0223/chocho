#!/usr/bin/env python3
"""Source-consistent observed-index bridge from earlier indexed host flora.

This is NOT a test of butterfly immigration, feeding, population growth,
colonization, local origin of plants, or true ecological absence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

PANEL = "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
HOST = "results/route_free_port_hosts/species_island_hosts.csv"
EFFORT = "results/route_free_port_hosts/port_coverage.csv"
EXPOSURES = ("introduced_only", "native_only", "both", "unresolved_only")
CATEGORIES = ("no_indexed_host",) + EXPOSURES
WINDOWS = ("new_index_2011_2017", "new_index_2018_2022", "previously_indexed_reindex_2011_2022")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cat(h: dict[str, str]) -> str:
    n = int(h["native"]) > 0
    i = int(h["introduced"]) > 0
    u = int(h["unresolved"]) > 0
    if n and i:
        return "both"
    if n:
        return "native_only"
    if i:
        return "introduced_only"
    if u:
        return "unresolved_only"
    return "no_indexed_host"


def make_panel(root: Path, history_csv: Path, summary_json: Path):
    original = read(root / PANEL)
    if len(original) != 11475 or len({(x["speciesKey"], x["landmass_id"]) for x in original}) != 11475:
        raise ValueError("Original 225x51 butterfly index has changed")
    if len({x["speciesKey"] for x in original}) != 225 or len({x["landmass_id"] for x in original}) != 51:
        raise ValueError("Original 225 butterflies or 51 landmasses changed")
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    if summary["status"] != "SOURCE_COMPLETE" or summary["original_panel_sha256"] != sha(root / PANEL):
        raise ValueError("Mismatch with uniform historical source")
    history = read(history_csv)
    hm = {(x["speciesKey"], x["landmass_id"]): x for x in history}
    if len(history) != 11475 or len(hm) != 11475:
        raise ValueError("Historical 225x51 identities missing or duplicated")
    if any(x["pre2001_source_complete"] != "1" for x in history):
        raise ValueError("Historical source fragment was not queried successfully")
    if sum(int(x["indexed_first_2011_2017_longer_history"]) for x in history) != 99:
        raise ValueError("Earlier corrected record-index source changed")
    if sum(int(x["indexed_first_2018_2022_longer_history"]) for x in history) != 73:
        raise ValueError("Later corrected record-index source changed")
    if sum(int(x["pre2001_any_index"]) for x in history) != 894:
        raise ValueError("Historical 1700–2000 record-index source changed")

    hosts = read(root / HOST)
    hostmap = {(x["key"], x["island"]): x for x in hosts}
    if len(hostmap) != len(hosts):
        raise ValueError("Duplicate host species-island labels")
    effort = defaultdict(int)
    for p in read(root / EFFORT):
        effort[p["island"]] += int(p["records"])

    result = []
    for p in original:
        key = (p["speciesKey"], p["landmass_id"])
        h = hm[key]
        status = hostmap.get(key)
        bot = effort[p["landmass_id"]]
        if status is None or bot <= 0:
            continue
        prior = int(h["pre2001_indexed_records"]) > 0 or int(p["pre2001_2010_indexed_records"]) > 0
        post = int(p["future2011_2017_indexed_records"]) > 0 or int(p["future2018_2022_indexed_records"]) > 0
        result.append({
            "species": key[0],
            "island": key[1],
            "category": cat(status),
            "effort": bot,
            "prior": int(prior),
            "post": int(post),
            "prior_1700_2000_count": int(h["pre2001_indexed_records"]),
            "prior_2001_2010_count": int(p["pre2001_2010_indexed_records"]),
            "eligible_early": int(h["eligible_2011_2017_longer_history"]),
            "event_early": int(h["indexed_first_2011_2017_longer_history"]),
            "eligible_late": int(h["eligible_2018_2022_longer_history"]),
            "event_late": int(h["indexed_first_2018_2022_longer_history"]),
        })
    return result, {
        "source_paths": [PANEL, HOST, EFFORT],
        "source_sha256": {x: sha(root / x) for x in (PANEL, HOST, EFFORT)},
        "historical_risk_csv_sha256": sha(history_csv),
        "historical_risk_source_summary_sha256": sha(summary_json),
        "original_cells": len(original),
        "history_cells": len(history),
        "history_pre2001_positive": 894,
        "history_early_first_index": 99,
        "history_late_first_index": 73,
        "indexed_botanical_covered_cells": len(result),
    }


def groups(rows: list[dict], minimum: int, window: str):
    subset = [x for x in rows if x["effort"] >= minimum]
    if window == "new_index_2011_2017":
        subset = [x for x in subset if x["eligible_early"]]
        response = "event_early"
    elif window == "new_index_2018_2022":
        subset = [x for x in subset if x["eligible_late"]]
        response = "event_late"
    elif window == "previously_indexed_reindex_2011_2022":
        subset = [x for x in subset if x["prior"]]
        response = "post"
    else:
        raise ValueError(window)
    g = {k: {"risk": 0, "record_index": 0, "islands": set(), "species": set()}
         for k in CATEGORIES}
    for x in subset:
        z = g[x["category"]]
        z["risk"] += 1
        z["record_index"] += x[response]
        z["islands"].add(x["island"])
        z["species"].add(x["species"])
    tab = {}
    for key, v in g.items():
        tab[key] = {
            "risk": v["risk"], "record_index": v["record_index"],
            "index_fraction": v["record_index"] / v["risk"] if v["risk"] else None,
            "island_count": len(v["islands"]), "species_count": len(v["species"]),
        }
    return subset, response, tab


def fe_residualize(values: np.ndarray, a: np.ndarray, b: np.ndarray):
    """Alternating two-way within projection for unbalanced binary risk tables."""
    z = values.astype(float).copy()
    groups = []
    for labels in (a, b):
        _, inverse = np.unique(labels, return_inverse=True)
        count = np.bincount(inverse)
        groups.append((inverse, count))
    converged = False
    for iteration in range(2000):
        max_mean = 0.0
        for inverse, count in groups:
            mu = np.bincount(inverse, weights=z) / count
            max_mean = max(max_mean, float(np.max(np.abs(mu))))
            z -= mu[inverse]
        if max_mean < 1e-10:
            converged = True
            break
    if not converged:
        raise ValueError("Two-way fixed-effect projection did not converge")
    return z


def model(subset: list[dict], response: str, adjust_prior_intensity: bool = False):
    y = np.array([x[response] for x in subset], dtype=float)
    island = np.array([x["island"] for x in subset])
    species = np.array([x["species"] for x in subset])
    X = np.array([[int(x["category"] == key) for key in EXPOSURES]
                  for x in subset], dtype=float)
    if adjust_prior_intensity:
        previous = np.array([
            [np.log1p(x["prior_1700_2000_count"]),
             np.log1p(x["prior_2001_2010_count"])] for x in subset
        ], dtype=float)
        X = np.column_stack((X, previous))
    Xw = np.column_stack([fe_residualize(X[:, j], island, species)
                          for j in range(X.shape[1])])
    yw = fe_residualize(y, island, species)
    rank = np.linalg.matrix_rank(Xw)
    if rank != X.shape[1]:
        return {"status": "NO_FULL_FE_RANK", "rank": int(rank), "n": len(subset),
                "events": int(y.sum()), "covariates": int(X.shape[1])}
    beta = np.linalg.lstsq(Xw, yw, rcond=None)[0]
    resid = yw - Xw @ beta
    bread = np.linalg.pinv(Xw.T @ Xw)

    def clustered_meat(labels: np.ndarray):
        # group scores sum (x_i * residual_i) before outer product
        _, inverse = np.unique(labels, return_inverse=True)
        scores = np.zeros((int(inverse.max()) + 1, Xw.shape[1]))
        np.add.at(scores, inverse, Xw * resid[:, None])
        return scores.T @ scores

    # Approximate two-way clustered covariance with intersection correction;
    # not a decisive statistical test with sparse exposed first indexes.
    cross = np.array([f"{a}|{b}" for a, b in zip(island, species)])
    cov = bread @ (
        clustered_meat(island) + clustered_meat(species) - clustered_meat(cross)
    ) @ bread
    islands = len(set(island))
    if islands > 1:
        cov *= islands / (islands - 1)
    def contrast(v):
        v = np.array(v, dtype=float)
        b = float(v @ beta)
        variance = float(v @ cov @ v)
        # A two-way robust estimate can have nonpositive variance at small support.
        se = float(np.sqrt(variance)) if variance > 0 else None
        return {"percentage_points": 100 * b, "two_way_cluster_se_pp":
                100 * se if se is not None else None,
                "approximate_95pct_pp": [100*(b-1.96*se), 100*(b+1.96*se)]
                if se is not None else None}
    def host_contrast(four):
        return contrast(list(four) + ([0., 0.] if adjust_prior_intensity else []))

    results = {key: host_contrast(np.eye(len(EXPOSURES))[i])
               for i, key in enumerate(EXPOSURES)}
    results["both_minus_native_only"] = host_contrast([0, -1, 1, 0])
    results["both_minus_introduced_only"] = host_contrast([-1, 0, 1, 0])
    results["native_minus_introduced_only"] = host_contrast([-1, 1, 0, 0])
    if adjust_prior_intensity:
        results["ln1p_prior_1700_2000_record_count"] = contrast([0, 0, 0, 0, 1, 0])
        results["ln1p_prior_2001_2010_record_count"] = contrast([0, 0, 0, 0, 0, 1])
    return {
        "status": "EXPLORATORY_TWOWAY_FE_ASSOCIATION",
        "n": len(subset), "events": int(y.sum()), "islands": islands,
        "species": len(set(species)), "rank": int(rank),
        "baseline": "no indexed compatible host in source port inventory; not real host absence",
        "adjust_previous_butterfly_index_counts": adjust_prior_intensity,
        "coefficients": results,
        "model": "LPM with butterfly species/destination island FE, four exclusive host classes"
                 + (", plus two separate log prior source-count controls" if adjust_prior_intensity else ""),
        "uncertainty": "approximate two-way cluster sandwich; no test of ecological arrival or breeding",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--history-csv", type=Path, required=True)
    p.add_argument("--history-summary", type=Path, required=True)
    p.add_argument("--out", type=Path,
                   default=Path("results/native_introduced_bridge_20261116/source_index_bridge.json"))
    a = p.parse_args()
    rows, receipts = make_panel(a.repo_root, a.history_csv, a.history_summary)
    reports = {}
    for minimum in (10, 30):
        for window in WINDOWS:
            subset, response, categories = groups(rows, minimum, window)
            name = f"plants_ge{minimum}__{window}"
            reports[name] = {
                "source_period": window, "min_2001_2010_plant_records": minimum,
                "eligible_source_cells": len(subset),
                "source_index_events": sum(x[response] for x in subset),
                "category_counts": categories,
                "model": model(subset, response),
            }
            if window == "previously_indexed_reindex_2011_2022":
                reports[name]["model_prior_intensity_adjusted"] = model(
                    subset, response, adjust_prior_intensity=True
                )
                # Explicitly diagnostic restrictions: same frozen cohort, not
                # independent external replications or new island sampling.
                restrictions = {
                    "exclude_great_britain_15_and_ireland_26": [
                        x for x in subset if x["island"] not in ("15", "26")
                    ],
                    "prior_2001_2010_index_positive": [
                        x for x in subset if x["prior_2001_2010_count"] > 0
                    ],
                    "prior_pre2001_index_positive": [
                        x for x in subset if x["prior_1700_2000_count"] > 0
                    ],
                }
                reports[name]["diagnostic_reindex_scopes"] = {}
                for scope_name, subscope in restrictions.items():
                    reports[name]["diagnostic_reindex_scopes"][scope_name] = {
                        "risk": len(subscope),
                        "events": sum(x[response] for x in subscope),
                        "categories": {
                            c: [sum(x["category"] == c for x in subscope),
                                sum(x["category"] == c and x[response] for x in subscope)]
                            for c in CATEGORIES
                        },
                        "adjusted": model(
                            subscope, response, adjust_prior_intensity=True
                        ) if len(subscope) > 50 else {
                            "status": "TOO_FEW_CELLS", "n": len(subscope)
                        },
                    }
    result = {
        "schema": "chocho_source_temporal_index_host_class_20261116_v2",
        "status": "SOURCE_DESCRIPTIVE_NOT_COLONIZATION",
        "inputs": receipts,
        "results": reports,
        "interpretation_gate": [
            "Historical GBIF zero records do not imply biological absence.",
            "2001–2010 port botanical host record and contemporary regional WCVP origin do not establish a local introduction date or larval host utilization.",
            "No-indexed-host is an unrecorded/unknown exposure; this is not a known biological zero.",
            "Butterfly recording effort is not fully measured; species and island fixed effects cannot remove species×island collection heterogeneity.",
            "For previously recorded butterflies we additionally adjust for the two historical record counts, but these are only imperfect observation-effort/intensity proxies.",
            "The 2018–2022 first-index risk population excludes recorded butterflies in all earlier windows; the two first-index periods are different risk sets.",
            "The previously indexed reindexing endpoint is GBIF record recurrence, not confirmed residence, population persistence, reproduction or a new colonization.",
            "First-index host-positive events may be too sparse for reliable conditional effects, regardless of a positive or negative model coefficient.",
            "No plant abundance, climate suitability, crossing route, ship transport or independent population establishment was measured.",
            "Geographic and historical subscope sensitivities were chosen after examining the main result; they are descriptive diagnostics, not confirmatory replications.",
        ],
        "preserved_results": ["239-butterfly global candidate-resource expansion +54.9%",
                              "Island relative resource expansion +118.9%; mainland +40.9%",
                              "155-species known-butterfly recording prediction MSE reduction +4.12%",
                              "51-island retrospective conditional host and marine associations"],
    }
    output = a.out if a.out.is_absolute() else a.repo_root / a.out
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("SOURCE_INDEX_BRIDGE", json.dumps({
        k: {
            "risk": v["eligible_source_cells"],
            "events": v["source_index_events"],
            "indexed": {c: [z["risk"], z["record_index"]]
                        for c, z in v["category_counts"].items()},
            "model": {"status": v["model"]["status"],
                      "coefficients": v["model"].get("coefficients", {})},
            "historical_record_count_adjusted": (
                v["model_prior_intensity_adjusted"].get("coefficients")
                if "model_prior_intensity_adjusted" in v else None
            ),
            "geography_and_prior_period_sensitivities": (
                {scope: {
                    "risk": val["risk"], "events": val["events"],
                    "categories": val["categories"],
                    "model_status": val["adjusted"]["status"],
                    "adjusted_coefficients": val["adjusted"].get("coefficients", {})
                } for scope, val in v["diagnostic_reindex_scopes"].items()}
                if "diagnostic_reindex_scopes" in v else None
            )
        }
        for k, v in reports.items()
    }, ensure_ascii=False), flush=True)
    print("SOURCE_INDEX_BRIDGE_OUTPUT", str(output), flush=True)


if __name__ == "__main__":
    main()
