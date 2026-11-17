#!/usr/bin/env python3
"""Falsifiable paired native/alien flora similarity prediction under island holdouts.

The 126-island strict WCVP source provides accepted species co-reported native
and introduced status. 104 of those have source-verified shoreline terminals on
the original GSHHG marine topology; 5356 unordered pairs are held as the
common sea-geometry/port-proxy eligibility cohort.

In each leave-one-ISLAND-out fold, fit on pairs involving neither the heldout
island, then predict its pair similarities. All four models share identical
train/test pairs and feature scaling learned ONLY from the training pairs.
At the test island we know the WCVP richness margins; this is a conditional
cross-region relationship diagnostic, NOT true prospective unseen flora.

Compare:
B = log great-circle distance + botanical Level1 equality + latitude
difference + WCVP richness and island-area margins
G = B plus GSHHG marine water-only shortest-path detour indicator
P = B plus product of contemporary WPI port counts (NOT vessel OD)
GP = B plus BOTH
The GSHHG scalar is NOT measured directional NOAA ocean currents, particle
arrival, uniform Circuitscape resistance, viable seed drift or butterfly use.
AIS2020 only has 7 geoverified source pairs; it is descriptive, never used
to construct a fabricated global shipping network predictor.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np

ORIGINAL_SHIPPING_CI=44010307551
ORIGINAL_SEA_CI=44010673580
SEED=20261117
N_BOOT=2499
RIDGE_ALPHA=20.
MODELS=("B","G","P","GP")
SCORE_GROUPS=("native","alien")

def rows(p):
    with p.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))

def positivefloat(s):
    v=float(s)
    if not math.isfinite(v) or v<=0:raise RuntimeError(f"Expected positive geographic value: {s}")
    return v

def input_pair(a):
    plants=rows(a.plants)
    geometry=rows(a.geometry)
    if len(plants)!=7875 or len(geometry)!=7875:
        raise RuntimeError("Original 126-island accepted WCVP and GSHHG pairs not preserved")
    geo={}
    for r in geometry:
        key=(r["source_region"],r["destination_region"])
        if key in geo:raise RuntimeError("Duplicate sea geometry pair")
        geo[key]=r
    join=[]
    unresolved=[]
    for r in plants:
        key=r["source_region"],r["destination_region"]
        if key not in geo:raise RuntimeError("Native/alien WCVP geographic island pair missing from GSHHG source")
        z=geo[key]
        if abs(float(r["centroid_geodesic_km"])-float(z["centroid_geodesic_km"]))>1e-5:
            raise RuntimeError("WCVP and GSHHG have different original representative-point distances")
        if z["source_sea_geometry_coverage_status"]!="BOTH_ANCHORED":
            continue
        if not r["alien_Jaccard"] or not r["native_Jaccard"]:
            unresolved.append({"region_A":key[0],"region_B":key[1],
                "native_Jaccard_defined":bool(r["native_Jaccard"]),
                "alien_Jaccard_defined":bool(r["alien_Jaccard"])})
            continue
        sea=positivefloat(z["GSHHG_original_0p0625_sea_shortest_graph_km"])
        d=positivefloat(r["centroid_geodesic_km"])
        # Original GSHHG approximate distance is not a lower bound on
        # geodesic because its block centers use discretization.
        geo_detour=math.log1p(sea)-math.log1p(d)
        port_proxy=float(r["modern_port_pair_proxy_log1p_product"])
        if not math.isfinite(port_proxy):raise RuntimeError("Missing modern port count proxy")
        n1,n2=float(r["source_native_WCVP_accepted"]),float(r["dest_native_WCVP_accepted"])
        a1,a2=float(r["source_alien_WCVP_accepted"]),float(r["dest_alien_WCVP_accepted"])
        # Source accepted richness is provided for all eligible regions;
        # physical area is not in original R pair CSV, so use regional
        # geography metadata separately. No missing values imputed.
        join.append({
          "island1":key[0],"island2":key[1],
          "n1":n1,"n2":n2,"a1":a1,"a2":a2,
          "native":float(r["native_Jaccard"]),
          "alien":float(r["alien_Jaccard"]),
          "log_geo":math.log1p(d),
          "sea_detour":geo_detour,
          "port_pair_proxy":port_proxy,
          "same_level1":float(r["same_plant_level1"]),
          "abs_lat":abs(float(r["source_representative_lat"])-float(r["dest_representative_lat"])),
          "has_2020_direct_OD":int(float(r["AIS2020_symmetric_direct_source_rows"])>0),
          "has_2020_OD_port_coverage":r["AIS2020_both_regions_have_verified_port_endpoint"]=="True",
          "share_native_species_count":int(r["native_shared_accepted_taxa"]),
          "share_alien_species_count":int(r["alien_shared_accepted_taxa"]),
        })
    if not join:raise RuntimeError("No common-source WCVP marine pair cohort")
    return join,unresolved

def make_X(entries,flora,area_by):
    # All feature slots are predeclared; no outcome-based feature selection.
    loggeo=np.array([e["log_geo"] for e in entries])
    equal=np.array([e["same_level1"] for e in entries])
    abs_lat=np.array([e["abs_lat"] for e in entries])
    minsp=np.array([min(e["n1"],e["n2"]) if flora=="native" else min(e["a1"],e["a2"]) for e in entries])
    maxsp=np.array([max(e["n1"],e["n2"]) if flora=="native" else max(e["a1"],e["a2"]) for e in entries])
    minarea=np.array([min(area_by[e["island1"]],area_by[e["island2"]]) for e in entries])
    maxarea=np.array([max(area_by[e["island1"]],area_by[e["island2"]]) for e in entries])
    base=np.column_stack([
      loggeo,loggeo**2,equal,abs_lat,
      np.log1p(minsp),np.log1p(maxsp),
      np.log1p(minarea),np.log1p(maxarea)
    ])
    geo=np.array([e["sea_detour"] for e in entries])[:,None]
    port=np.array([e["port_pair_proxy"] for e in entries])[:,None]
    # Model column counts are 8,9,9,10; no geography FE for novel island.
    return {"B":base,"G":np.column_stack([base,geo]),
            "P":np.column_stack([base,port]),"GP":np.column_stack([base,geo,port])}

def ridge(trainX,trainY,testX):
    # Mean/sd training-only. Ridge alpha is fixed BEFORE inspecting any
    # holdout score, not optimized to find a positive effect.
    means=trainX.mean(axis=0)
    sd=trainX.std(axis=0)
    sd=np.where(sd<1e-9,1.,sd)
    x=(trainX-means)/sd
    z=(testX-means)/sd
    xx=np.column_stack([np.ones(len(x)),x])
    zz=np.column_stack([np.ones(len(z)),z])
    penalty=np.eye(xx.shape[1]);penalty[0,0]=0.
    A=xx.T@xx+RIDGE_ALPHA*penalty
    beta=np.linalg.solve(A,xx.T@trainY)
    pred=np.clip(zz@beta,0.,1.)
    if not np.isfinite(pred).all():
        raise RuntimeError("Nonfinite held-island flora Jaccard predictions")
    return pred

def evaluate(entries,flora,area_by):
    islands=sorted({e["island1"] for e in entries}|{e["island2"] for e in entries})
    if len(islands)<85:raise RuntimeError(f"Common WCVP marine support too small: {len(islands)}")
    target=np.array([e[flora] for e in entries])
    if np.any((target<0)|(target>1)):raise RuntimeError("Original WCVP Jaccard outside [0,1]")
    X=make_X(entries,flora,area_by)
    foldrows={}
    all_test=[]
    for held in islands:
        testidx=np.array([i for i,e in enumerate(entries)
                          if held in (e["island1"],e["island2"])],dtype=int)
        trainidx=np.array([i for i,e in enumerate(entries)
                           if held not in (e["island1"],e["island2"])],dtype=int)
        if len(testidx)<70 or len(trainidx)<3000:
            raise RuntimeError(f"Held island has too little train/test source support {held}")
        y_train=target[trainidx]; y_test=target[testidx]
        er={}
        preds={}
        for model in MODELS:
            yh=ridge(X[model][trainidx],y_train,X[model][testidx])
            er[model]=float(np.mean((yh-y_test)**2))
            preds[model]=yh
        foldrows[held]={"n_test":len(testidx),**{f"MSE_{k}":v for k,v in er.items()}}
        all_test.append({
          "held":held,"n_test":len(testidx),
          "test_index":testidx,"target":y_test,"pred":preds
        })
    # Island equal-weighted: in leave-one-island-out each undirected pair
    # enters exactly two held-island test sets; do NOT call 2× pairs independent.
    macro={m:float(np.mean([f[f"MSE_{m}"] for f in foldrows.values()])) for m in MODELS}
    improvements={
      "G_vs_B_pct":(macro["B"]-macro["G"])/macro["B"]*100,
      "P_vs_B_pct":(macro["B"]-macro["P"])/macro["B"]*100,
      "GP_vs_B_pct":(macro["B"]-macro["GP"])/macro["B"]*100,
      "G_added_after_P_pct":(macro["P"]-macro["GP"])/macro["P"]*100,
      "P_added_after_G_pct":(macro["G"]-macro["GP"])/macro["G"]*100
    }
    rng=np.random.default_rng(SEED + (4 if flora=="alien" else 2))
    mat=np.array([[f[f"MSE_{m}"] for m in MODELS] for f in foldrows.values()])
    statkeys=["G_vs_B_pct","P_vs_B_pct","GP_vs_B_pct","G_added_after_P_pct","P_added_after_G_pct"]
    boot={k:[] for k in statkeys}
    for _ in range(N_BOOT):
        subset=mat[rng.integers(0,len(islands),len(islands))].mean(axis=0)
        base,sea,port,both=subset
        vals=[(base-sea)/base*100,(base-port)/base*100,
              (base-both)/base*100,(port-both)/port*100,(sea-both)/sea*100]
        for k,v in zip(statkeys,vals):boot[k].append(v)
    ci={k:[float(np.quantile(v,.025)),float(np.quantile(v,.975))] for k,v in boot.items()}
    byisland=[{"code":h,**v} for h,v in foldrows.items()]
    return {
      "response":flora+"_WCVP_accepted_species_Jaccard",
      "n_botanical_islands":len(islands),"n_unique_island_pairs":len(entries),
      "n_total_held_island_test_pair_predictions":sum(z["n_test"] for z in foldrows.values()),
      "held_island_exactly_once_folds":len(foldrows),
      "ridge_alpha_fixed":RIDGE_ALPHA,"no_tuning_to_holdout":True,
      "held_island_never_appears_in_training_pair":True,
      "conditional_predictive_marginal_richness_of_heldout_island_known":True,
      "macro_MSE_equal_held_island":macro,
      "relative_delta_MSE_percent_positive_is_improvement":improvements,
      "descriptive_island_bootstrap_95pct_MSE_delta":ci,
      "island_folds_improved_G_over_B":sum(r["MSE_G"]<r["MSE_B"] for r in foldrows.values()),
      "island_folds_improved_P_over_B":sum(r["MSE_P"]<r["MSE_B"] for r in foldrows.values()),
      "per_island_loss":byisland
    }

def main():
    ap=argparse.ArgumentParser()
    for k in ("plants","geometry","island-area","sea-audit","shipping-audit","out-dir"):
        ap.add_argument("--"+k,required=True,type=Path)
    args=ap.parse_args()
    one=json.loads(args.sea_audit.read_text())
    two=json.loads(args.shipping_audit.read_text())
    if one["source_original_island_pairs_total"]!=7875 or two["N_WCVP_native_alien_pairs"]!=7875:
        raise RuntimeError("Source geometry/plant crosswalk counts drift")
    if one["source_valid_original_island_pair_geometric_shortest_sea"]!=5356:
        raise RuntimeError("Original GSHHG source seawater pair source support changed")
    if two["strict126_unordered_plant_pairs_with_one_or_more_actual_source_2020_port_OD"]!=7:
        raise RuntimeError("AIS2020 independently verified real OD pair support changed")
    ar=rows(args.island_area)
    area_by={r["code"]:float(r["area_km2"]) for r in ar}
    if len(area_by)!=126 or len(ar)!=126 or min(area_by.values())<=0:
        raise RuntimeError("Original area and 126 WCVP island codes changed")
    entries,unresolved=input_pair(args)
    islands=set(e["island1"] for e in entries)|set(e["island2"] for e in entries)
    if (len(islands)!=104 or len(entries)+len(unresolved)!=5356 or
        len(unresolved)>2 or
        any(e["island1"] not in area_by or e["island2"] not in area_by for e in entries)):
        raise RuntimeError(f"Original source geometry eligible=5356 common={len(entries)} unresolved_WCVP={unresolved} islands={len(islands)}")
    native=evaluate(entries,"native",area_by)
    alien=evaluate(entries,"alien",area_by)
    if native["n_unique_island_pairs"]!=alien["n_unique_island_pairs"]:
        raise RuntimeError("Crossed native/alien not on exactly same pairs")
    crossed={
       "native_geometry_added_after_port_pct":native["relative_delta_MSE_percent_positive_is_improvement"]["G_added_after_P_pct"],
       "native_port_added_after_geometry_pct":native["relative_delta_MSE_percent_positive_is_improvement"]["P_added_after_G_pct"],
       "alien_geometry_added_after_port_pct":alien["relative_delta_MSE_percent_positive_is_improvement"]["G_added_after_P_pct"],
       "alien_port_added_after_geometry_pct":alien["relative_delta_MSE_percent_positive_is_improvement"]["P_added_after_G_pct"],
    }
    report={
       "schema":"chocho_WCVP126_104source_islands_native_vs_alien_source_geometry_port_proxy_island_holdout_v1",
       "status":"SOURCE_VERIFIED_NATIVE_ALIEN_WCVP_GEOMETRY_AND_MODERN_PORT_PROXY_CROSSED_PREDICTION_NOT_CURRENT_OR_SHIP_CAUSALITY",
       "input_original_WCVP_pair_source_CI":44010011288,
       "original_GSHHG_126island_geometry_source_CI":ORIGINAL_SEA_CI,
       "original_WPI_AIS2020_126island_coverage_source_CI":ORIGINAL_SHIPPING_CI,
       "n_original_WCVP_islands":126,"n_common_anchored_source_islands":len(islands),
       "n_common_valid_island_pairs":len(entries),
       "n_original_common_GSHHG_eligible_pairs_before_WCVP_response_missingness":5356,
       "original_WCVP_Jaccard_not_defined_in_one_or_both_origins":unresolved,
       "n_source_confirmed_real_ship_2020_OD_pairs_across_all126":7,
       "n_possible_2020_AIS_verified_endpoints_both_regions_pairs":45,
       "model":"same frozen source pairs, no island leakage, ridge alpha20, Level1 equality and source geography/botanical richness",
       "model_B":"geodesic distance+distance squared+Level1 equality+lat difference+native/alien richness min/max+area min/max",
       "model_G":"model_B+GSHHG shortest-water distance minus great-circle distance on log1p scale; symmetric sea geometry, NOT ocean current",
       "model_P":"model_B+WPI2023 logged port count product; NOT vessel voyages or route",
       "model_GP":"model_B+both marine shortest-path geometry and WPI port product",
       "native_flora":native,"introduced_flora":alien,"crossed_contrasts":crossed,
       "outcome_specification":"WCVP accepted taxa Jaccard contemporary co-recorded native or introduced state",
       "scientific_narrative_limits":[
         "The user-desired natural current versus human-mediated ship crossed-path comparison CANNOT be confirmed by this test: G is undirected shoreline geometry, and P is a modern port-number proxy.",
         "Correct real 2020 origin-destination port matched source has only 7 original 126-island cross-region pairs, too sparse for an independent global shipping model.",
         "The native-versus-alien Jaccard numerators and denominators have different richness and source/effort coverage; compare within-response relative MSE but not absolute MSE across responses.",
         "Islands are source units and dyads share endpoints; leave-one-island-out holds out all incidents, but original botanical geography is the same WCVP snapshot, not an independently replicated future period.",
         "Features include current flora richness of the held island, so these results condition on the observed marginal species pool and are not prospective new-region ecological forecasts.",
         "The main source graph topology is original GSHHG 0.0625° component-aware ocean graph, weighted coarse shortest-water distance, not the 51-island electrical resistance.",
         "Physical NOAA directional currents, source propagule viability, ocean drift times, plant arrival chronology, actual historical vessel cargo and butterfly colonisation remain separate unverified estimands."
       ]
    }
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"WCVP126_native_alien_crossed_held_island_prediction_audit.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("WCVP126_CROSSED_PREDICTION_SOURCE_GATES",json.dumps({
        "status":report["status"],"n_common_pairs":len(entries),
        "native_relative_MSE":native["relative_delta_MSE_percent_positive_is_improvement"],
        "alien_relative_MSE":alien["relative_delta_MSE_percent_positive_is_improvement"],
        "native_ci":native["descriptive_island_bootstrap_95pct_MSE_delta"],
        "alien_ci":alien["descriptive_island_bootstrap_95pct_MSE_delta"],
        "crossed":crossed
    }),flush=True)

if __name__=="__main__":main()
