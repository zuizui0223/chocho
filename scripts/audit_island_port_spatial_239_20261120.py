#!/usr/bin/env python3
"""Prespecified spatial autocorrelation sensitivity for fixed-239 island/port study.

Input is the immutable-source independent 126-island result CSV from the parent
audit, plus the same pinned TDWG Level3 GeoJSON. This does not reproduce the
unavailable temporary extra.py or the original 125-island / 2517-species panel.

Centroids describe botanical *regions* (not necessarily single islands). Spatial
Moran weights: binary <=2000 km great-circle separation, no self-links; 1,999
permutations, both unrestricted and within nine TDWG Level1 groups.
Moran eigenvector filter: top 10 positive eigenvectors of H W H, fixed a priori
from the spatial weights alone (not chosen using outcomes or p-values).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from pyproj import Geod

from audit_island_port_isolation_239_20261120 import (
    load_geojson, design, y_standardized,
)

RADIUS_KM = 2000.0
MEM_COUNT = 10
PERMUTATIONS = 1999
SEED = 20261120


def make_weights(frame: pd.DataFrame, geojson: Path):
    regions = load_geojson(geojson, set(frame.code))
    if len(regions) != 126:
        raise RuntimeError("Expected 126 source-pinned botanical region geometries")
    # Shapely centroid of each entire botanical area: a reproducible proximity
    # proxy, potentially poor for multi-island, dateline-spanning regions.
    longitudes = np.array([regions[c].centroid.x for c in frame.code], dtype=float)
    latitudes = np.array([regions[c].centroid.y for c in frame.code], dtype=float)
    n = len(frame)
    geod = Geod(ellps="WGS84")
    _, _, metres = geod.inv(
        np.repeat(longitudes, n), np.repeat(latitudes, n),
        np.tile(longitudes, n), np.tile(latitudes, n)
    )
    dist = np.asarray(metres, dtype=float).reshape(n, n) / 1000
    np.fill_diagonal(dist, 0.0)
    if not np.isfinite(dist).all():
        raise RuntimeError("Non-finite centroid distance")
    weights = (dist <= RADIUS_KM).astype(float)
    np.fill_diagonal(weights, 0.0)
    if not np.array_equal(weights, weights.T):
        raise RuntimeError("Spatial adjacency must be symmetric")
    if weights.sum() < n:
        raise RuntimeError("Not enough near-neighbour edges for a spatial test")
    return weights, {"centroid_distance_km_median_offdiagonal":
                     float(np.median(dist[np.triu_indices(n, k=1)])),
                     "undirected_neighbor_pairs": int(weights.sum()/2),
                     "unconnected_regions": int((weights.sum(axis=1)==0).sum()),
                     "neighbor_count_range": [
                        int(weights.sum(axis=1).min()),
                        int(weights.sum(axis=1).max())]}


def preselected_mem(weights):
    n = weights.shape[0]
    h = np.eye(n) - np.ones((n,n))/n
    eigenvalues, vectors = np.linalg.eigh(h @ weights @ h)
    valid = np.flatnonzero(eigenvalues > 1e-8)[::-1]
    if len(valid) < MEM_COUNT:
        raise RuntimeError(f"Only {len(valid)} positive Moran eigenvectors")
    return vectors[:, valid[:MEM_COUNT]], eigenvalues[valid[:MEM_COUNT]]


def moran(residual, weights, groups, permutations=PERMUTATIONS):
    e = np.asarray(residual, dtype=float)
    e = e - e.mean()
    n = len(e)
    s0 = weights.sum()
    if s0 <= 0 or e @ e <= 0:
        raise RuntimeError("Degenerate Moran calculation")
    def score(v):
        return float((n/s0)*(v @ weights @ v)/(v @ v))
    observed = score(e)
    rng = np.random.default_rng(SEED)
    global_ge = local_ge = 0
    indices = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    for _ in range(permutations):
        perm = rng.permutation(e)
        global_ge += score(perm) >= observed
        perm = e.copy()
        for ind in indices:
            perm[ind] = e[rng.permutation(ind)]
        local_ge += score(perm) >= observed
    return {
        "I": observed,
        "p_positive_global_permutation": float((global_ge+1)/(permutations+1)),
        "p_positive_within_level1_permutation": float((local_ge+1)/(permutations+1))
    }


def term_summary(model, term):
    a,b = model.conf_int().loc[term]
    return {"beta":float(model.params[term]),"HC3_p":float(model.pvalues[term]),
            "HC3_95":[float(a),float(b)]}


def evaluate(frame, weights, mem):
    base = design(frame)
    added = pd.DataFrame(mem, columns=[f"spatial_MEM_{i+1:02d}" for i in range(mem.shape[1])],
                         index=base.index)
    extended = pd.concat([base, added], axis=1)
    if np.linalg.matrix_rank(extended.to_numpy()) != extended.shape[1]:
        raise RuntimeError("Spatial augmented model is rank deficient")
    ys = {
        "plant_introduced": y_standardized(frame, "introduced_plant"),
        "host_native_239": y_standardized(frame, "native_host_239"),
        "host_current_239": y_standardized(frame, "current_host_239"),
        "host_current_minus_native_z":
            y_standardized(frame, "current_host_239")-
            y_standardized(frame, "native_host_239"),
    }
    out = {}
    for name, y in ys.items():
        raw = sm.OLS(y, base).fit(cov_type="HC3")
        filtered = sm.OLS(y, extended).fit(cov_type="HC3")
        out[name] = {
            "plain":{"ports":term_summary(raw,"ports"),
                     "distance":term_summary(raw,"distance"),
                     "residual_moran_2000km":moran(raw.resid,weights,frame.level1.to_numpy()),
                     "r2":float(raw.rsquared)},
            "spatial_MEM10":{"ports":term_summary(filtered,"ports"),
                     "distance":term_summary(filtered,"distance"),
                     "residual_moran_2000km":moran(filtered.resid,weights,frame.level1.to_numpy()),
                     "r2":float(filtered.rsquared)},
        }
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--island-csv", type=Path, required=True)
    p.add_argument("--tdwg-geojson", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a=p.parse_args()
    frame=pd.read_csv(a.island_csv,dtype={"code":str,"level1":str})
    if len(frame)!=126 or not frame.code.is_unique or frame.level1.nunique()!=9:
        raise RuntimeError("Unexpected source panel or Level1 geography")
    weights, graph=make_weights(frame,a.tdwg_geojson)
    mem, eig=preselected_mem(weights)
    results=evaluate(frame,weights,mem)
    report={
      "schema":"chocho_fixed239_2000km_moran_MEM10_20261120_v1",
      "scope":"126 WGSRPD3 botanical regions / fixed 239 known-host butterflies; not original temporary 125-island / 2517 panel",
      "weight_definition":"W_ij=1 when WGS84 geodesic of full-region shapely centroids <= 2000 km, zero otherwise; diagonal=0",
      "moran":"Raw residual I after area, ports, distance, absolute latitude and TDWG Level1; positive-tail 1999 permutations global / within Level1; residual permutation is descriptive, not a spatial-error likelihood",
      "spatial_filter":"10 largest strictly positive eigenvectors of H W H, chosen without outcomes, jointly added to same OLS controls; HC3 conditional intervals",
      "important_caveats":["multi-island botanical region centroids can mislocate geographical neighbours",
                          "MEM10 and nine-cluster CR1/t8 are distinct sensitivity designs",
                          "plain permutation p-values need not control residual confounding",
                          "no causal shipping or butterfly establishment inference"],
      "graph":graph,
      "top_10_positive_moran_eigenvalues":[float(v) for v in eig],
      "results":results,
    }
    a.out_dir.mkdir(parents=True, exist_ok=True)
    dest=a.out_dir/"SOURCE239_MORAN_MEM10_SPATIAL_SENSITIVITY.json"
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[]
    for endpoint,d in results.items():
        for variant,v in d.items():
            for term in ("ports","distance"):
                lines.append({"endpoint":endpoint,"model":variant,"term":term,
                    **v[term],"residual_I":v["residual_moran_2000km"]["I"],
                    "moran_global_p":v["residual_moran_2000km"]["p_positive_global_permutation"],
                    "moran_within_level1_p":v["residual_moran_2000km"]["p_positive_within_level1_permutation"]})
    pd.DataFrame(lines).to_csv(a.out_dir/"SOURCE239_MORAN_MEM10_EFFECTS.csv", index=False)
    print("SPATIAL_239_RECALCULATED",json.dumps({
        "graph":graph,
        "positive_eigenvalues_count":len(eig),
        "endpoints":{k:{"raw_Moran_I":v["plain"]["residual_moran_2000km"]["I"],
                        "mem10_Moran_I":v["spatial_MEM10"]["residual_moran_2000km"]["I"],
                        "mem10_beta_ports":v["spatial_MEM10"]["ports"]["beta"],
                        "mem10_p_ports":v["spatial_MEM10"]["ports"]["HC3_p"],
                        "mem10_beta_distance":v["spatial_MEM10"]["distance"]["beta"],
                        "mem10_p_distance":v["spatial_MEM10"]["distance"]["HC3_p"]}
                     for k,v in results.items()}},ensure_ascii=False))

if __name__=="__main__":
    main()
