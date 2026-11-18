#!/usr/bin/env python3
"""Source-locked Jaccard floor diagnostics in WCVP104 island origin × sea circuits.

Exploratory follow-up after positive electrical-circuit gains were known.
Use original source-only five hash island bins, 5,355 matched Jaccard dyads,
original GSHHG undirected uniform-electrical circuit, ridge(alpha=20), geography,
port and known same-origin richness; no historically observed dispersal.
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

N_PAIRS = 5355
N_ISLANDS = 104
N_STRICT = 1089
N_BINS = 5
BOOT_N = 4999
BOOT_SEED = 20261118
MODELS = ("B", "G", "P", "GP")
TARGETS = ("native", "alien")

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))

def numeric(value):
    if value is None or str(value).strip() == "":
        return None
    z = float(value)
    return z if math.isfinite(z) else None

def positive_quantiles(series):
    pos = series[series > 0]
    return {
        "n":int(series.size),
        "positive_n":int(pos.size),
        "zero_n":int((series == 0).sum()),
        "zero_fraction":float(np.mean(series == 0)),
        "positive_fraction":float(np.mean(series > 0)),
        "median_all":float(np.median(series)),
        "median_positive":float(np.median(pos)) if len(pos) else None,
        "q25_positive":float(np.quantile(pos, 0.25)) if len(pos) else None,
        "q75_positive":float(np.quantile(pos, 0.75)) if len(pos) else None,
        "mean_all":float(np.mean(series)),
    }

def comparison(matrix):
    values=matrix.mean(axis=0)
    B,G,P,GP=map(float,values)
    if min(B,G,P,GP)<=0:
        raise RuntimeError("Source Brier/MSE loss nonpositive")
    return {
        "circuit_vs_base_gain_pct":100*(B-G)/B,
        "port_vs_base_gain_pct":100*(B-P)/B,
        "circuit_port_vs_base_gain_pct":100*(B-GP)/B,
        "circuit_after_port_gain_pct":100*(P-GP)/P,
        "port_after_circuit_gain_pct":100*(G-GP)/G,
    }

def evaluate(entries, outcome, y, area_by, *, train_positive_only=False, score_mask=None, binary=False):
    X=src.make_X(entries, outcome, area_by)
    groups={}
    ids=sorted({r["island1"] for r in entries}|{r["island2"] for r in entries})
    for code in ids:
        groups.setdefault(strict.label(code),set()).add(code)
    if len(groups)!=N_BINS:
        raise RuntimeError("Original source 5-fold geographic bins changed")
    folds=[]
    for k in range(N_BINS):
        held=groups[k]
        tr=np.array([i for i,e in enumerate(entries) if
            e["island1"] not in held and e["island2"] not in held],dtype=int)
        te=np.array([i for i,e in enumerate(entries) if
            e["island1"] in held and e["island2"] in held],dtype=int)
        if len(tr)<2500 or len(te)<50:
            raise RuntimeError(f"Source geographic fold {k} does not have sufficient support")
        if train_positive_only:
            tr=tr[y[tr]>0]
            te=te[y[te]>0]
            if len(tr)<30 or len(te)<5:
                return {"status":"NOT_ESTIMABLE_POSITIVE_ONLY_SUPPORT",
                        "failed_fold":k,"positive_train":int(len(tr)),"positive_test":int(len(te))}
        if score_mask is not None:
            te=te[score_mask[te]]
            if len(te)<5:
                return {"status":"NOT_ESTIMABLE_MATCHED_POSITIVE_SUPPORT",
                        "failed_fold":k,"positive_pair_test":int(len(te))}
        losses={}
        for m in MODELS:
            pred=src.ridge(X[m][tr],y[tr],X[m][te])
            losses["MSE_"+m]=float(np.mean((y[te]-pred)**2))
        record={
            "fold":k,
            "n_held_islands":len(held),
            "n_train_pairs":len(tr),
            "n_test_pairs":len(te),
            "n_positive_test":int(np.sum(y[te]>0)),
            **losses,
        }
        if binary:
            base_rate=float(np.mean(y[tr]))
            record["train_shared_taxon_prevalence"]=base_rate
            record["held_test_shared_taxon_prevalence"]=float(np.mean(y[te]))
            record["null_constant_training_prevalence_MSE"]=float(np.mean((y[te]-base_rate)**2))
        folds.append(record)
    mat=np.array([[f["MSE_"+m] for m in MODELS] for f in folds])
    point=comparison(mat)
    weighted=np.average(mat,weights=np.array([f["n_test_pairs"] for f in folds]),axis=0)
    return {
        "status":"ESTIMATED_EXPLORATORY_SAME_SOURCE",
        "n_test_pairs":int(sum(f["n_test_pairs"] for f in folds)),
        "n_held_bins":N_BINS,
        "source_hash_folds_unchanged":True,
        "folds":folds,
        "equal_fold_mean_mse":dict(zip(MODELS,map(float,mat.mean(axis=0)))),
        "equal_fold_relative_gain_pct":point,
        "pair_weighted_relative_gain_pct":{
            "circuit_vs_base_gain_pct":float(100*(weighted[0]-weighted[1])/weighted[0]),
            "port_vs_base_gain_pct":float(100*(weighted[0]-weighted[2])/weighted[0])
        },
        "n_folds_circuit_better":sum(f["MSE_G"]<f["MSE_B"] for f in folds),
        "n_folds_port_better":sum(f["MSE_P"]<f["MSE_B"] for f in folds),
        "equal_fold_constant_null_mse":float(np.mean([f["null_constant_training_prevalence_MSE"] for f in folds])) if binary else None,
        "predictions_are_historic_arrival":False,
    }

def paired_interval(native, alien):
    if native["status"]!="ESTIMATED_EXPLORATORY_SAME_SOURCE" or alien["status"]!="ESTIMATED_EXPLORATORY_SAME_SOURCE":
        return {"status":"NOT_ESTIMABLE"}
    if [r["fold"] for r in native["folds"]] != [r["fold"] for r in alien["folds"]]:
        raise RuntimeError("Origin folds lost exact identity")
    nmat=np.array([[f["MSE_"+m] for m in MODELS] for f in native["folds"]])
    amat=np.array([[f["MSE_"+m] for m in MODELS] for f in alien["folds"]])
    rng=np.random.default_rng(BOOT_SEED)
    native_draw=[];alien_draw=[];difference=[]
    for _ in range(BOOT_N):
        ix=rng.integers(0,N_BINS,N_BINS)
        n=comparison(nmat[ix])["circuit_vs_base_gain_pct"]
        a=comparison(amat[ix])["circuit_vs_base_gain_pct"]
        native_draw.append(n);alien_draw.append(a);difference.append(a-n)
    interval=lambda z:[float(np.quantile(z,.025)),float(np.quantile(z,.975))]
    g_native=comparison(nmat)["circuit_vs_base_gain_pct"]
    g_alien=comparison(amat)["circuit_vs_base_gain_pct"]
    return {
        "native_gain_pct":g_native,
        "introduced_gain_pct":g_alien,
        "introduced_minus_native_gain_points":g_alien-g_native,
        "native_descriptive_5fold_interval":interval(native_draw),
        "introduced_descriptive_5fold_interval":interval(alien_draw),
        "paired_descriptive_5fold_interval":interval(difference),
        "bootstrap_draws":BOOT_N,
        "seed":BOOT_SEED,
        "statistical_significance_confirmed":False,
        "fold_independence_guaranteed":False,
    }

def main():
    parser=argparse.ArgumentParser()
    for name in ("plants","geometry","uniform-circuit","circuit-receipt","island-area","out-dir"):
        parser.add_argument("--"+name,type=Path,required=True)
    a=parser.parse_args()
    receipt=json.loads(a.circuit_receipt.read_text(encoding="utf-8"))
    if (receipt.get("status")!="SOURCE_VERIFIED_ORIGINAL_UNIFORM_MARINE_CIRCUIT_104_REGIONS_NOT_OCEAN_CURRENTS" or
        receipt.get("source_original_104_unordered_pair_resistances")!=5356 or
        receipt.get("source_component_marine_nodes")!=189550):
        raise RuntimeError("Original full-source true uniform circuit identity changed")
    entries,undefined=src.input_pair(SimpleNamespace(plants=a.plants,geometry=a.geometry))
    if len(entries)!=N_PAIRS or len(undefined)!=1:
        raise RuntimeError("Original Jaccard sample changed")
    areas={r["code"]:float(r["area_km2"]) for r in read_csv(a.island_area)}
    if len(areas)!=126:
        raise RuntimeError("Original region area source changed")
    original={
        "native":strict.evaluate(entries,"native",areas),
        "alien":strict.evaluate(entries,"alien",areas)
    }
    expected={"native":1.73510722916521,"alien":1.8576534605954167}
    for kind in TARGETS:
        obj=original[kind]
        if (obj["n_double_unseen_test_pairs"]!=N_STRICT or
            abs(obj["macro_equal_fold_relative_mse_improvement_percent"]["G_vs_B_pct"]-expected[kind])>1e-8):
            raise RuntimeError("Frozen shortest marine strict source did not reproduce")
    pair_data={}
    for row in read_csv(a.plants):
        key=tuple(sorted((row["source_region"],row["destination_region"])))
        if key in pair_data:raise RuntimeError("Duplicate WCVP pair in original source")
        pair_data[key]=row
    electr={}
    for row in read_csv(a.uniform_circuit):
        key=tuple(sorted((row["source_region"],row["destination_region"])))
        if key in electr:raise RuntimeError("Duplicate original electrical resistance pair")
        value=float(row["log1p_normalized_uniform_circuit"])
        if not (math.isfinite(value) and value>0):
            raise RuntimeError("Nonnumeric authentic circuit resistance")
        electr[key]=value
    if len(electr)!=5356:raise RuntimeError("Expected 5356 authentic circuit pairs")
    extra=[]
    for e in entries:
        key=tuple(sorted((e["island1"],e["island2"])))
        row=pair_data[key]
        for kind in TARGETS:
            if abs(float(row[kind+"_Jaccard"])-e[kind])>1e-10:
                raise RuntimeError("Original same-pair Jaccard source mismatch")
        if key not in electr:raise RuntimeError("Missing original pair circuit resistance")
        e["sea_detour"]=electr[key]
        extra.append(row)
    jaccard={}
    for kind in TARGETS:
        jaccard[kind]=evaluate(entries,kind,np.array([e[kind] for e in entries]),areas)
        full=original[kind]
        if jaccard[kind]["n_test_pairs"]!=N_STRICT:raise RuntimeError("Jaccard test population changed")
        for fold_a, fold_b in zip(jaccard[kind]["folds"],full["folds"]):
            if (fold_a["fold"]!=fold_b["fold"] or
                fold_a["n_train_pairs"]!=fold_b["n_train_pairs"] or
                fold_a["n_test_pairs"]!=fold_b["n_test_double_unseen_pairs"]):
                raise RuntimeError("Frozen Jaccard strict source fold membership changed")
            for model in ("B","P"):
                if abs(fold_a["MSE_"+model]-fold_b["MSE_"+model])>1e-12:
                    raise RuntimeError("Frozen source B/P baseline MSE mismatch")
        point=jaccard[kind]["equal_fold_relative_gain_pct"]["circuit_vs_base_gain_pct"]
        ref={"native":-0.04430509181116742,"alien":2.735040680490108}[kind]
        if abs(point-ref)>1e-8:raise RuntimeError("Original electrical circuit Jaccard numeric replay failed")
    # Freeze the originally defined 1,089 strict test pairs for zero-rate counts.
    strict_mask=np.array([strict.label(e["island1"])==strict.label(e["island2"]) for e in entries])
    if int(strict_mask.sum())!=N_STRICT:raise RuntimeError("Strict test pair mask drift")
    source_counts={}
    binary={};sorensen={};pos_only={}
    for kind in TARGETS:
        y=np.array([e[kind] for e in entries],dtype=float)
        if np.any(~np.isfinite(y)):raise RuntimeError("Source Jaccard absent")
        source_counts[kind]={
            "all_source_pairs":positive_quantiles(y),
            "strict_unseen_test_pairs":positive_quantiles(y[strict_mask])
        }
        by=(y>0).astype(float)
        binary[kind]=evaluate(entries,kind,by,areas,binary=True)
        if binary[kind]["n_test_pairs"]!=N_STRICT:raise RuntimeError("Binary original heldout sample differs")
        alt=np.array([numeric(row[kind+"_Sorensen"]) for row in extra],dtype=float)
        if np.any(~np.isfinite(alt)):raise RuntimeError("Sorensen undefined on the original Jaccard sample")
        sorensen[kind]=evaluate(entries,kind,alt,areas)
        pos_only[kind]=evaluate(entries,kind,y,areas,train_positive_only=True)
    paired_positive=np.array([e["native"]>0 and e["alien"]>0 for e in entries])
    matched_positive={}
    for kind in TARGETS:
        matched_positive[kind]=evaluate(entries,kind,np.array([e[kind] for e in entries]),areas,score_mask=paired_positive)
    # Simpson overlap requires a nonzero smaller flora. Missing is NOT zero.
    simpson_mask=np.array([all(numeric(row[k+"_Simpson_overlap"]) is not None for k in TARGETS) for row in extra])
    simpson_entries=[e for e,ok in zip(entries,simpson_mask) if ok]
    simpson_raw=[z for z,ok in zip(extra,simpson_mask) if ok]
    if len(simpson_entries)!=5151:raise RuntimeError("Original source-defined Simpson common population drift")
    simpson={}
    for kind in TARGETS:
        y=np.array([float(row[kind+"_Simpson_overlap"]) for row in simpson_raw],dtype=float)
        simpson[kind]=evaluate(simpson_entries,kind,y,areas)
    doc={
        "schema":"chocho_WCVP104_uniform_electrical_circuit_zero_floor_and_similarity_transform_20261118_v1",
        "status":"SOURCE_VERIFIED_EXPLORATORY_ZERO_FLOOR_AND_ALTERNATE_RESPONSE_DIAGNOSTICS",
        "source_original_CI_ids":{"flora":44010307551,"sea_geometry":44010673580,"authentic_circuit":44013658955,
                                  "reference_strict_electrical":44290887759},
        "protocol":"results/WCVP126_circuit_floor_diagnostic_20261118/PROTOCOL.md",
        "source_original_104_islands":N_ISLANDS,
        "source_both_origin_Jaccard_pairs":N_PAIRS,
        "strict_both_endpoints_unseen_test_pairs":N_STRICT,
        "source_original_alpha":src.RIDGE_ALPHA,
        "source_control_shortest_water_reproduced":{key:original[key]["macro_equal_fold_relative_mse_improvement_percent"] for key in TARGETS},
        "authentic_uniform_electrical_Jaccard_reproduced":jaccard,
        "zero_floor_source_descriptives":source_counts,
        "binary_any_taxon_shared_Brier":binary,
        "Sorensen_alternative_same_pair_source":sorensen,
        "positive_Jaccard_conditional_fit_only":pos_only,
        "same_positive_native_and_alien_test_pairs":matched_positive,
        "Simpson_common_defined_pair_n":len(simpson_entries),
        "Simpson_undefined_on_original_common_sample_n":int(N_PAIRS-len(simpson_entries)),
        "Simpson_overlap_common_defined_only":simpson,
        "paired_same_5fold_descriptive_intervals":{
            "Jaccard":paired_interval(jaccard["native"],jaccard["alien"]),
            "binary_any_taxon_shared_Brier":paired_interval(binary["native"],binary["alien"]),
            "Sorensen":paired_interval(sorensen["native"],sorensen["alien"]),
            "positive_Jaccard_fit":paired_interval(pos_only["native"],pos_only["alien"]),
            "both_originally_positive_test":paired_interval(matched_positive["native"],matched_positive["alien"]),
            "Simpson_common_defined_only":paired_interval(simpson["native"],simpson["alien"])
        },
        "limits":[
            "All alternate response definitions are exploratory extensions after positive source circuit results were known.",
            "Jaccard/Sorensen/Simpson are transforms or denominator choices on the SAME accepted WCVP catalogue, not three independent biological replications.",
            "Simpson has missing rather than zero where either source-accepted flora is empty; the common source-defined cohort is smaller.",
            "Binary sharing predicts catalogue co-occurrence, not plant exchange, propagule colonization or ecological equivalence of a rare shared taxon.",
            "Positive-only subset is defined with target response and must not be interpreted as forecasting whether sharing occurs.",
            "The five island bins share source observations and overlapping training data; any resampled intervals are descriptive, not independent geospatial inference.",
            "No directional currents, historical vessel OD or butterfly presence are outcomes of this test."
        ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    outfile=a.out_dir/"source_audited_zero_floor_and_origin_predictability.json"
    outfile.write_text(json.dumps(doc,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("WCVP104_NATIVE_ALIEN_JACCARD_FLOOR_DIAGNOSTIC",json.dumps({
        "status":doc["status"],"common_pairs":N_PAIRS,"strict_pairs":N_STRICT,
        "zero_native_pct":100*source_counts["native"]["strict_unseen_test_pairs"]["zero_fraction"],
        "zero_alien_pct":100*source_counts["alien"]["strict_unseen_test_pairs"]["zero_fraction"],
        "model_tests":{
            name:{
                "native_gain_pct":(result["native"]["equal_fold_relative_gain_pct"].get("circuit_vs_base_gain_pct")
                                   if result["native"]["status"]=="ESTIMATED_EXPLORATORY_SAME_SOURCE" else None),
                "alien_gain_pct":(result["alien"]["equal_fold_relative_gain_pct"].get("circuit_vs_base_gain_pct")
                                   if result["alien"]["status"]=="ESTIMATED_EXPLORATORY_SAME_SOURCE" else None)
            } for name,result in (
                ("binary",binary),("Sorensen",sorensen),
                ("positive_only",pos_only),("both_positive",matched_positive),
                ("Simpson",simpson))
        },
        "paired_difference_intervals":{
            k:v.get("paired_descriptive_5fold_interval") for k,v in doc["paired_same_5fold_descriptive_intervals"].items()
        },
        "simpson_common_n":len(simpson_entries)
    },ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
