#!/usr/bin/env python3
"""Audit later butterfly GBIF re-indexing versus prior port-indexed known host categories.

Crucial boundary: GBIF source indexing is neither true absence/presence nor
island colonisation. The host list is catalogue compatibility, not observed
larval feeding. Historical risk is corrected using original 1700-2000 screen.

Source data are frozen in chocho plus the source-complete previous Actions artifact.
No GBIF calls, retrospective redefinition of island/species, or manuscript edits.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

INPUTS = {
    "butterflies": "data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv",
    "hosts": "results/route_free_port_hosts/species_island_hosts.csv",
    "ports": "results/route_free_port_hosts/port_coverage.csv",
}
GIT_BLOBS = {
    "butterflies": "b6a3c98eb76864052d45d520a55b7d830816126c",
    "hosts": "8ab83a1ebb5bc58e815e19c45ff338025f7b2bd7",
    "ports": "755ef83fb8ad534a82f20981c02b7b263bd210a4",
}
CATEGORIES = ("introduced_only", "native_only", "both", "unresolved")
ALL_CATEGORIES = ("no_indexed_catalog_host",) + CATEGORIES
CRIT = 2.03  # approximate t(35) near fixed 36-island primary
GB_IRELAND_GSHHG_IDS = frozenset(("15", "26"))


def table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def git_blob(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def two_way_demean(v: np.ndarray, group1: np.ndarray, group2: np.ndarray) -> np.ndarray:
    out = v.astype(float, copy=True)
    if out.ndim == 1:
        out = out[:, None]
    for _ in range(150):
        max_group_mean = 0.
        for group in (group1, group2):
            counts = np.bincount(group)
            means = np.vstack([
                np.bincount(group, weights=out[:, j], minlength=len(counts)) / counts
                for j in range(out.shape[1])
            ]).T
            out -= means[group]
            max_group_mean = max(max_group_mean, float(np.max(np.abs(means))))
        if max_group_mean < 1e-10:
            break
    if max_group_mean > 1e-7:
        raise ValueError("Two-way within estimator did not converge")
    return out


def cluster_diagonal(influence: np.ndarray, groups: np.ndarray) -> np.ndarray:
    n_group = int(np.max(groups)) + 1
    totals = np.column_stack([
        np.bincount(groups, weights=influence[:, k], minlength=n_group)
        for k in range(influence.shape[1])
    ])
    return (totals**2).sum(axis=0)


def fit(rows: list[dict], label: str, adjust_early: bool, later_risk_only: bool) -> dict:
    if later_risk_only:
        rows = [r for r in rows if not r["pre2011"]]
    n = len(rows)
    if n < 100:
        raise ValueError(f"too few source records {label}")
    labels = list(CATEGORIES) + (["pre2011_recorded"] if adjust_early else [])
    species = np.asarray([r["species"] for r in rows])
    islands = np.asarray([r["island"] for r in rows])
    _, sp = np.unique(species, return_inverse=True)
    _, il = np.unique(islands, return_inverse=True)
    _, joint = np.unique(species + "|" + islands, return_inverse=True)
    y = np.asarray([r["later_recorded"] for r in rows], dtype=float)
    x = np.column_stack(
        [np.asarray([int(r["host_category"] == category) for r in rows])
         for category in CATEGORIES]
        + ([np.asarray([r["pre2011"] for r in rows])] if adjust_early else [])
    )
    yd = two_way_demean(y, il, sp)[:, 0]
    xd = two_way_demean(x, il, sp)
    rank = int(np.linalg.matrix_rank(xd))
    if rank != len(labels):
        raise ValueError(f"insufficient group/category rank {label}: {rank}/{len(labels)}")
    bread = np.linalg.inv(xd.T @ xd)
    beta = bread @ (xd.T @ yd)
    residual = yd - xd @ beta
    influence = residual[:, None] * (xd @ bread.T)
    var = (cluster_diagonal(influence, il) + cluster_diagonal(influence, sp)
           - cluster_diagonal(influence, joint))
    n_islands = len(set(islands.tolist()))
    n_species = len(set(species.tolist()))
    var *= n_islands / (n_islands - 1)
    se = np.sqrt(np.maximum(var, 0))
    categories = {}
    for category in ALL_CATEGORIES:
        sub = [r for r in rows if r["host_category"] == category]
        categories[category] = {
            "cells": len(sub),
            "later_GBIF_indexed": sum(r["later_recorded"] for r in sub),
            "earlier_1700_2010_GBIF_indexed": sum(r["pre2011"] for r in sub),
            "later_source_index_fraction": sum(r["later_recorded"] for r in sub) / len(sub)
            if sub else None,
        }
    return {
        "panel": label,
        "n_species_island_cells": n,
        "islands": n_islands,
        "species": n_species,
        "later_indexed_cells": int(y.sum()),
        "earlier_indexed_cells": sum(r["pre2011"] for r in rows),
        "reference": "no_indexed_catalog_host: unknown, not confirmed biological absence",
        "model": "Two-way island/species FE linear probability model; one row per species x island",
        "early_record_control": adjust_early,
        "later_risk_only": later_risk_only,
        "category_counts": categories,
        "estimates": {
            name: {
                "beta_pp": float(100 * beta[i]),
                "approx_two_way_cluster_se_pp": float(100 * se[i]),
                "approx_interval_95_pp": [float(100*(beta[i]-CRIT*se[i])),
                                          float(100*(beta[i]+CRIT*se[i]))],
            }
            for i, name in enumerate(labels)
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--history-csv", type=Path, required=True)
    parser.add_argument("--history-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path,
                        default=Path("results/host_marine_bridge_audit_20261116/port_host_record_crossover.json"))
    args = parser.parse_args()

    files = {k: args.repo_root / p for k, p in INPUTS.items()}
    hashes = {}
    for key, path in files.items():
        sha = git_blob(path)
        if sha != GIT_BLOBS[key]:
            raise ValueError(f"Changed input Git blob: {key}: {sha}")
        hashes[key] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "git_blob": sha}

    butterflies = table(files["butterflies"])
    hosts = table(files["hosts"])
    port_rows = table(files["ports"])
    history = table(args.history_csv)
    summary = json.loads(args.history_summary.read_text(encoding="utf-8"))
    if (len(butterflies) != 225*51 or len(history) != 225*51
            or summary["status"] != "SOURCE_COMPLETE"
            or summary["original_panel_sha256"] != hashes["butterflies"]["sha256"]):
        raise ValueError("Original complete 225 x 51 historical source does not match")
    hist = {(r["speciesKey"], r["landmass_id"]): r for r in history}
    hm = {(r["key"], r["island"]): r for r in hosts}
    if len(hist) != 225 * 51 or len(hm) != len(hosts):
        raise ValueError("Duplicate original host or historical index key")
    if any(row["pre2001_source_complete"] != "1" for row in history):
        raise ValueError("An original island has incomplete historical source")
    effort: dict[str, int] = defaultdict(int)
    for p in port_rows:
        effort[p["island"]] += int(p["records"])

    combined = []
    for b in butterflies:
        h = hm.get((b["speciesKey"], b["landmass_id"]))
        if h is None or effort[b["landmass_id"]] < 1:
            continue
        native, intro, unres = (int(h[z]) > 0 for z in ("native", "introduced", "unresolved"))
        if native and intro:
            category = "both"
        elif native:
            category = "native_only"
        elif intro:
            category = "introduced_only"
        elif unres:
            category = "unresolved"
        else:
            category = "no_indexed_catalog_host"
        hb = hist[b["speciesKey"], b["landmass_id"]]
        early = bool(int(hb["pre2001_indexed_records"]) or
                     int(b["pre2001_2010_indexed_records"]))
        later = bool(int(b["future2011_2017_indexed_records"]) or
                     int(b["future2018_2022_indexed_records"]))
        combined.append({
            "species": b["speciesKey"],
            "island": b["landmass_id"],
            "port_plant_records": effort[b["landmass_id"]],
            "host_category": category,
            "pre2011": int(early),
            "later_recorded": int(later),
        })

    plans = [
        ("primary_ge10_history_adjusted", 10, False, True, False),
        ("ge30_history_adjusted", 30, False, True, False),
        ("ge100_history_adjusted", 100, False, True, False),
        ("ge10_no_GB_Ireland", 10, True, True, False),
        ("ge10_early_source_unindexed_only", 10, False, False, True),
        ("ge10_no_history_control", 10, False, False, False),
    ]
    panels = {}
    for name, minimum, exclude_gb, control_early, risk_only in plans:
        use = [r for r in combined if r["port_plant_records"] >= minimum
               and (not exclude_gb or r["island"] not in GB_IRELAND_GSHHG_IDS)]
        panels[name] = fit(use, name, control_early, risk_only)

    result = {
        "schema": "butterfly_known-host-plant_GBIF_records_descriptive_bridge_20261116",
        "source_git_blobs": hashes,
        "history": {
            "full_11475_1700_2000_SHA256": hashlib.sha256(args.history_csv.read_bytes()).hexdigest(),
            "original_summary_SHA256": hashlib.sha256(args.history_summary.read_bytes()).hexdigest(),
            "screen_complete_cells": len(history),
            "sum_pre2001_indexed": sum(int(r["pre2001_indexed_records"]) > 0 for r in history),
        },
        "panels": panels,
        "scientific_interpretation": {
            "outcome": "at least one later (2011-2022) GBIF index in source species x island; not first arrival",
            "earlier_control": "known 1700-2010 GBIF index from uniform island screen, including 2001-2010",
            "host_categories": "known butterfly larval-host candidate taxa GBIF indexed near original WPI port 2001-2010; WCVP regional origin",
            "all_nonindexed_species_or_hosts": "UNKNOWN not biological absence",
            "excluded_claims": [
                "host use on that island", "larval survival", "true butterfly occupancy",
                "plant introduction dates", "ship transport", "island colonization",
                "population persistence", "causal mediation of global resource expansion"
            ],
            "limitations": [
                "same species x island observational source with incomplete GBIF sampling",
                "geographical/covariate effort differences and host detection vary by taxon",
                "island/species fixed effects do not ensure species-specific climate niche match",
                "host catalogue plant presences may be spatially separated from butterfly records",
                "retrospective association; historical GBIF index not comprehensive island fauna",
                "approximate cluster intervals, same dataset used in multiple diagnostics"
            ],
        },
    }
    target = args.output if args.output.is_absolute() else args.repo_root / args.output
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    simple = {name: {"n": panel["n_species_island_cells"], "events": panel["later_indexed_cells"],
                     "beta_pp": {k: round(v["beta_pp"], 4)
                                 for k, v in panel["estimates"].items()}}
              for name, panel in panels.items()}
    print(json.dumps({"status": "SOURCE_LOCKED", "models": simple, "output": str(target)},
                     ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
