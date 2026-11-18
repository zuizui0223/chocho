#!/usr/bin/env python3
"""Source-locked 104-island strict double-unseen transfer of *uniform electrical* sea circuits.

Reuses original WCVP/GSHHG plant pairs, response, ridge(alpha=20), source-only
5 hashed geographic folds from original 2026-11-17 shortest-water-path control.
Only marine predictor changes to authenticated undirected finite-volume
Circuitscape resistance; no ship routes or ocean-current direction measured.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np

import audit_WCVP126_crossed_native_alien_geometry_port_predictions_20261117 as src
import audit_WCVP126_native_alien_double_unseen_islands_20261117 as strict

SOURCE_CIRCUIT_CI = 44013658955
SOURCE_ORIGINAL_PAIR_CI = 44010307551
SOURCE_COAST_CI = 44010673580
SOURCE_PARENT_SINGLE_UNSEEN = 44013825017
SOURCE_STRICT_SHORTEST = 44011411479
PAIR_N = 5355
ISLAND_N = 104
FOLD_N = 5
BOOT_N = 4999
SEED = 20261118

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def contrasts(fold_mat):
    # Identical to the frozen original double-unseen implementation.
    B, C, P, PC = fold_mat.mean(axis=0)
    if B<=0 or P<=0 or C<=0:raise RuntimeError("Nonpositive reference MSE")
    return {
        "circuit_beyond_geography_pct":float(100*(B-C)/B),
        "port_beyond_geography_pct":float(100*(B-P)/B),
        "circuit_and_port_beyond_geography_pct":float(100*(B-PC)/B),
        "circuit_given_port_pct":float(100*(P-PC)/P),
        "port_given_circuit_pct":float(100*(C-PC)/C),
    }

def matrix(report):
    rows = report["folds"]
    if len(rows)!=FOLD_N:raise RuntimeError("Expected same five original source folds")
    return np.array([[r["MSE_"+m] for m in ("B","G","P","GP")] for r in rows], dtype=float)

def main():
    ap=argparse.ArgumentParser()
    for name in ("plants","geometry","uniform-circuit","circuit-receipt","island-area","out-dir"):
        ap.add_argument("--"+name, required=True, type=Path)
    a=ap.parse_args()
    cr=json.loads(a.circuit_receipt.read_text(encoding="utf-8"))
    if (cr.get("status")!="SOURCE_VERIFIED_ORIGINAL_UNIFORM_MARINE_CIRCUIT_104_REGIONS_NOT_OCEAN_CURRENTS" or
        cr.get("source_original_104_unordered_pair_resistances")!=5356 or
        cr.get("source_component_marine_nodes")!=189550):
        raise RuntimeError("Original authentic 104-coastal-region circuit source lock not verified")
    inp=SimpleNamespace(plants=a.plants,geometry=a.geometry)
    entries,undefined=src.input_pair(inp)
    if len(entries)!=PAIR_N or len(undefined)!=1:raise RuntimeError("WCVP native/alien pair source drift")
    areas={r["code"]:float(r["area_km2"]) for r in read_csv(a.island_area)}
    if len(areas)!=126:raise RuntimeError("Original 126 WCVP island area roster drift")
    islands={e["island1"] for e in entries}|{e["island2"] for e in entries}
    if len(islands)!=ISLAND_N:raise RuntimeError("Original 104 accepted coast regions drift")
    xy={}
    for r in read_csv(a.uniform_circuit):
        key=tuple(sorted((r["source_region"],r["destination_region"])))
        if key in xy:raise RuntimeError("Duplicate marine resistance pair")
        v=float(r["log1p_normalized_uniform_circuit"])
        if not (math.isfinite(v) and v>0):
            raise RuntimeError("Invalid original electrical resistance")
        xy[key]=v
    if len(xy)!=5356:raise RuntimeError("Original 5356 electrical marine dyads drift")
    shortest={
        "native":strict.evaluate(entries,"native",areas),
        "alien":strict.evaluate(entries,"alien",areas),
    }
    fixed={"native":1.73510722916521,"alien":1.8576534605954167}
    for kind,expected in fixed.items():
        fit=shortest[kind]
        got=fit["macro_equal_fold_relative_mse_improvement_percent"]["G_vs_B_pct"]
        if abs(got-expected)>1e-8:raise RuntimeError(f"Original shortest-water double-unseen score changed ({kind}): {got}")
        if fit["n_double_unseen_test_pairs"]!=1089:
            raise RuntimeError("Strict original double-unseen pair count changed")
    for e in entries:
        key=tuple(sorted((e["island1"],e["island2"])))
        if key not in xy:raise RuntimeError("Missing source-resolved original circuit dyad")
        e["sea_detour"]=xy[key]
    # Same exact target, standardization, ridge alpha, five hashed folds and test
    # indices as previously source-verified geometry-only double-unseen test.
    circuit={
        "native":strict.evaluate(entries,"native",areas),
        "alien":strict.evaluate(entries,"alien",areas),
    }
    for kind in ("native","alien"):
        ref,aug=shortest[kind],circuit[kind]
        if ref["n_double_unseen_test_pairs"]!=aug["n_double_unseen_test_pairs"]:
            raise RuntimeError("Source pair population not identical")
        if [q["held_island_codes"] for q in ref["folds"]] != [q["held_island_codes"] for q in aug["folds"]]:
            raise RuntimeError("Held island fold membership changed")
        for r,z in zip(ref["folds"],aug["folds"]):
            if r["n_train_pairs"]!=z["n_train_pairs"] or r["n_test_double_unseen_pairs"]!=z["n_test_double_unseen_pairs"]:
                raise RuntimeError("Fold train/test split changed")
            for m in ("B","P"):
                if abs(r["MSE_"+m]-z["MSE_"+m])>1e-12:
                    raise RuntimeError(f"Frozen geographical baseline/port score changed ({kind}, {m})")
        if not aug["test_islands_absent_from_train_endpoints"]:
            raise RuntimeError("Geographical island endpoint leakage detected")
    nmat=matrix(circuit["native"])
    amat=matrix(circuit["alien"])
    actual_native=contrasts(nmat)
    actual_alien=contrasts(amat)
    paired_diff=actual_alien["circuit_beyond_geography_pct"]-actual_native["circuit_beyond_geography_pct"]
    n_test= circuit["native"]["n_double_unseen_test_pairs"]
    rng=np.random.default_rng(SEED)
    qn={k:[] for k in actual_native}
    qa={k:[] for k in actual_alien}
    qdelta=[]
    for _ in range(BOOT_N):
        idx=rng.integers(0,FOLD_N,FOLD_N)
        nstat=contrasts(nmat[idx])
        astat=contrasts(amat[idx])
        for k,v in nstat.items():qn[k].append(v)
        for k,v in astat.items():qa[k].append(v)
        qdelta.append(astat["circuit_beyond_geography_pct"]-nstat["circuit_beyond_geography_pct"])
    ci=lambda v:[float(np.quantile(v,.025)),float(np.quantile(v,.975))]
    out={
        "schema":"chocho_104_WCVP_islands_5355_pairs_uniform_electrical_circuit_double_unseen_origin_cross_v1",
        "status":"SOURCE_REPLAYED_BOTH_ENDPOINTS_GEOGRAPHICALLY_UNSEEN_NO_CAUSAL_CURRENT_OR_SHIP_MECHANISM",
        "protocol":"results/WCVP126_circuit_strict_geographic_transfer_20261118/PROTOCOL.md",
        "source_runs":{
            "plant_geography":SOURCE_ORIGINAL_PAIR_CI,
            "source_coast":SOURCE_COAST_CI,
            "uniform_GSHHG_circuit":SOURCE_CIRCUIT_CI,
            "original_single_unseen_circuit":SOURCE_PARENT_SINGLE_UNSEEN,
            "original_double_unseen_shortest_water":SOURCE_STRICT_SHORTEST,
        },
        "source_common_regions":len(islands),"source_common_pairs":len(entries),
        "source_n_test_double_unseen_pairs":n_test,
        "n_source_folds":FOLD_N,
        "uses_original_hashed_source_only_fold_allocation":True,
        "source_model_ridge_alpha":src.RIDGE_ALPHA,
        "original_shortest_double_unseen_verified_before_electrical_fit":{
            k:v["macro_equal_fold_relative_mse_improvement_percent"] for k,v in shortest.items()
        },
        "native_uniform_circuit":circuit["native"],
        "alien_uniform_circuit":circuit["alien"],
        "origin_paired_comparison":{
            "same_folds_same_pairs":True,
            "native_circuit_gain_pct":actual_native["circuit_beyond_geography_pct"],
            "alien_circuit_gain_pct":actual_alien["circuit_beyond_geography_pct"],
            "alien_minus_native_relative_mse_gain_points":paired_diff,
            "descriptive_shared_5fold_resample_4999_interval":ci(qdelta),
            "bootstrap_seed":SEED,
            "bootstrap_draws":BOOT_N,
            "native_MSE_gain_descriptive_4999_interval":ci(qn["circuit_beyond_geography_pct"]),
            "alien_MSE_gain_descriptive_4999_interval":ci(qa["circuit_beyond_geography_pct"]),
        },
        "main_claim_boundary":[
            "This is an exploratory source-reused stricter geographically unseen holdout, specified after original source gains were known.",
            "The electrical GSHHG circuit is UNDIRECTED sea geometry, NOT observed directional ocean currents, propagule survival or ship voyage/cargo transport.",
            "Although native and alien share island pairs, their Jaccard source margins differ; paired relative error gain is not a causal comparative effect.",
            "The five outcome-blind held-island bins have overlapping training and plant source data, so bootstrap intervals are descriptive, not independent effect confidence intervals.",
            "Marginal botanical richness of each unseen island is known to the model; there is no historically prospective test or butterfly colonization observation.",
        ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    target=a.out_dir/"WCVP126_104_true_uniform_circuit_strict_double_unseen.json"
    target.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("WCVP104_TRUE_UNIFORM_CIRCUIT_STRICT_DOUBLE_UNSEEN",json.dumps({
        "status":out["status"],"n_test":n_test,
        "native":actual_native,"alien":actual_alien,
        "alien_minus_native_gain_points":paired_diff,
        "paired_descriptive_interval":out["origin_paired_comparison"]["descriptive_shared_5fold_resample_4999_interval"],
        "native_descriptive_interval":out["origin_paired_comparison"]["native_MSE_gain_descriptive_4999_interval"],
        "alien_descriptive_interval":out["origin_paired_comparison"]["alien_MSE_gain_descriptive_4999_interval"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
