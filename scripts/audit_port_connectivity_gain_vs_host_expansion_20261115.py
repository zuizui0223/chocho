#!/usr/bin/env python3
"""Does introduced-host resource expansion predict *species-specific* port-connectivity gain?

This is an observational SECONDARY analysis, *not* the unavailable species M3->M5
host-information effect. Instead it tests T1 port-connection added to the
already HOST-informed M5 model. Both per-species M5 and T1 squared errors
were originally committed in the DISTINCT port-transport forecast extension.

No actual host introduction date, shipborne larva, immigration, occupancy,
population establishment, or cross-validated geographic transfer is estimated.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

BASE="provenance/reviewer_defenses/results/transport_forecast_temporal_20261030"
SOURCES={
    "temporal155":f"{BASE}/transport_temporal_all155_species_errors.csv",
    "island95":f"{BASE}/transport_temporal_strict_islands_species_errors.csv",
    "island_nonGB95":f"{BASE}/transport_temporal_strict_islands_without_GB_Ireland_species_errors.csv",
    "comp":f"{BASE}/transport_temporal_comparison.csv",
    "traits":"data/frozen/figure_sources/anthropogenic_species_metrics.csv",
    "forecast_support":"provenance/reviewer_defenses/results/host_forecast_validation_20261029/forecast_all_225_final_coverage.csv",
    "island_support":"provenance/reviewer_defenses/results/island_forecast_overlap_20261030/island_forecast_overlap_all225.csv",
    "lock":f"{BASE}/transport_forecast_extension_lock.json",
}
SEED=20261115
R=4999

def read(root):
    docs={}
    shas={}
    for k,p in SOURCES.items():
        raw=(root/p).read_bytes()
        shas[k]=hashlib.sha256(raw).hexdigest()
        if p.endswith(".json"):
            docs[k]=json.loads(raw)
        else:
            with (root/p).open(newline="",encoding="utf-8-sig") as f:
                docs[k]=list(csv.DictReader(f))
    return docs,shas

def vec_rank(a):
    a=np.asarray(a,dtype=float)
    idx=np.argsort(a,kind="mergesort")
    z=np.empty(len(a),dtype=float)
    i=0
    while i<len(a):
        j=i+1
        while j<len(a) and a[idx[j]]==a[idx[i]]:
            j+=1
        z[idx[i:j]]=(i+j+1)*.5
        i=j
    return z

def corr(a,b):
    a=np.asarray(a,dtype=float)
    b=np.asarray(b,dtype=float)
    if np.std(a)<=0 or np.std(b)<=0:
        return float("nan")
    return float(np.corrcoef(a,b)[0,1])

def spearman(a,b):
    return corr(vec_rank(a),vec_rank(b))

def partial_rank(a,b,ctrl):
    a=vec_rank(a);b=vec_rank(b)
    c=np.column_stack((np.ones(len(a)),ctrl))
    ar=a-c@np.linalg.lstsq(c,a,rcond=None)[0]
    br=b-c@np.linalg.lstsq(c,b,rcond=None)[0]
    return corr(ar,br)

def controls(rows):
    # Observational covariate adjustment. Host-family, global native source
    # range, original known host count, GBIF source recording and M5 loss.
    family=["2_families","3_to_5_families","6plus_families"]
    return np.column_stack([
        vec_rank([np.log1p(z["native_global"]) for z in rows]),
        vec_rank([np.log1p(z["known_hosts"]) for z in rows]),
        vec_rank([np.log1p(z["source_records"]) for z in rows]),
        vec_rank([np.log1p(z["m5"]) for z in rows]),
        *[np.array([int(z["family"]==f) for z in rows]) for f in family],
    ])

def est(rows,seed=SEED,rep=R):
    n=len(rows)
    if n<20:
        raise ValueError("Too few source taxa for useful correlation")
    xx=np.array([z["added_global"] for z in rows],float)
    lx=np.log1p(xx)
    ratio=np.array([z["log_gain"] for z in rows],float)
    y=np.array([z["m5"]-z["t1"] for z in rows],float)
    net=np.array([z["t0"]-z["t1"] for z in rows],float)
    rel=np.array([100*(z["m5"]-z["t1"])/z["m5"] for z in rows])
    ctrl=controls(rows)
    rho=spearman(xx,y);partial=partial_rank(lx,y,ctrl)
    rng=np.random.default_rng(seed)
    boot=[]
    for _ in range(rep):
        sel=rng.integers(0,n,size=n)
        rho_b=spearman(xx[sel],y[sel])
        if np.isfinite(rho_b):
            boot.append(rho_b)
    boot=np.array(boot)
    return {
        "n_species":n,
        "baseline_M5_mse_species_equal":float(np.mean([z["m5"] for z in rows])),
        "T0_port_location_mse":float(np.mean([z["t0"] for z in rows])),
        "T1_port_connectivity_mse":float(np.mean([z["t1"] for z in rows])),
        "T1_vs_M5_relative_MSE_reduction_pct":float(100*np.mean(y)/np.mean([z["m5"] for z in rows])),
        "T1_vs_T0_relative_MSE_reduction_pct":float(100*np.mean(net)/np.mean([z["t0"] for z in rows])),
        "species_T1_improves_over_M5":int(np.sum(y>0)),
        "species_T1_worsens_from_M5":int(np.sum(y<0)),
        "primary_rho_added_units_vs_absolute_M5_minus_T1_mse":rho,
        "descriptive_4999_species_bootstrap_rho_CI":np.quantile(boot,[.025,.975]).tolist(),
        "partial_rank_rho_added_vs_M5minusT1_given_7_source_covariates":partial,
        "additional_rhos":{
            "global_log_resource_gain_vs_absolute_delta_MSE":spearman(ratio,y),
            "added_host_units_vs_species_relative_gain":spearman(xx,rel),
            "added_host_units_vs_T1_beyond_T0_port_location":spearman(xx,net),
        },
        "native_and_catalogue_proxy_covariates":"rank(log1p(global native resource), log1p(known host species), log1p(GBIF source records), log1p(M5 loss)) + 3 frozen diet-family indicator variables",
    }

def stratified_perm(rows,observed,reps=R,seed=SEED+17):
    n=len(rows)
    native=np.array([r["native_global"] for r in rows])
    qs=np.quantile(native,[.25,.5,.75])
    classkey=lambda z: (z["family"],sum(z["native_global"]>q for q in qs))
    strata=defaultdict(list)
    for i,z in enumerate(rows):strata[classkey(z)].append(i)
    y=np.array([z["m5"]-z["t1"] for z in rows])
    x=np.array([z["added_global"] for z in rows])
    rng=np.random.default_rng(seed)
    null=[]
    for _ in range(reps):
        out=y.copy()
        for ids in strata.values():
            if len(ids)>1:
                out[ids]=rng.permutation(y[ids])
        null.append(spearman(x,out))
    null=np.array(null)
    center=float(np.mean(null))
    p=(1+np.sum(np.abs(null-center)>=abs(observed-center)-1e-12))/(reps+1)
    return {
        "strata_definition":"frozen dietary-family class x quartile of global native host resource-region count",
        "strata":len(strata),
        "singleton_strata":sum(len(ids)==1 for ids in strata.values()),
        "sampled_species":n,
        "permutations":reps,
        "observed_rho":observed,
        "stratified_null_mean_rho":center,
        "stratified_null_q025_q975":np.quantile(null,[.025,.975]).tolist(),
        "exploratory_two_sided_p_centered_on_null":float(p),
        "not_causal_or_phylogenetically_independent":True,
    }

def quartiles(rows):
    arr=sorted(rows,key=lambda z:(z["added_global"],z["species"]))
    result=[]
    for i in range(4):
        sub=arr[int(i*len(arr)/4):int((i+1)*len(arr)/4)]
        metric=est(sub,seed=SEED+100*i,rep=499)
        result.append({
            "nominal_quartile":i+1, "source_species":len(sub),
            "min_added":min(z["added_global"] for z in sub),
            "max_added":max(z["added_global"] for z in sub),
            "mean_added":float(np.mean([z["added_global"] for z in sub])),
            "T1_over_M5_reduction_pct":metric["T1_vs_M5_relative_MSE_reduction_pct"],
            "n_T1_improved":metric["species_T1_improves_over_M5"],
        })
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument("--out",type=Path,default=Path("results/port_resource_predgain_species_20261115"))
    a=p.parse_args()
    data,digests=read(a.root)
    if data["lock"]["status"]!="EXPLORATORY_EXTENSION_FIXED_BEFORE_FEATURE_SUPPORT_AND_MODEL_RESULTS":
        raise ValueError("Original port extension stage lock altered")
    traits={z["species"]:z for z in data["traits"]}
    selections={z["species"]:z for z in data["forecast_support"]}
    islands={z["species"]:z for z in data["island_support"]}
    if (len(traits),len(selections),len(islands))!=(239,225,225):
        raise ValueError("Original 239 resource / 225 candidate source identity changed")
    source_tables={"all155":"temporal155","strict_islands":"island95",
                   "strict_islands_excluding_GB_Ireland":"island_nonGB95"}
    results={}
    joined=[]
    for name,key in source_tables.items():
        rows=[]
        for r in data[key]:
            species=r["species"]
            if species not in traits or species not in selections or species not in islands:
                raise ValueError("Unknown source species "+species)
            t=traits[species];sel=selections[species]
            if sel["final_disposition"]!="MODELED_ALL_PLANNED_VALIDATIONS":
                raise ValueError("Species in transport source never in originally fitted 155")
            z={
                "species":species,
                "group":name,
                "m5":float(r["M5_host_composition_sq"]),
                "t0":float(r["T0_port_location_sq"]),
                "t1":float(r["T1_port_connectivity_sq"]),
                "native_global":int(t["native_resource_units"]),
                "added_global":int(t["introduced_added_units"]),
                "log_gain":float(t["log_resource_expansion"]),
                "known_hosts":int(t["resolved_host_species"]),
                "family":t["host_breadth_stratum"],
                "source_records":int(float(sel["source_records"])),
                "temporal_test_pairs":int(sel["temporal_test_pairs"]),
                "strict_source_island_test_rows":int(islands[species]["strict_island_rows"]),
            }
            if z["m5"]<=0 or z["t0"]<=0 or z["t1"]<=0:
                raise ValueError("Nonpositive species MSE impossible")
            rows.append(z)
        if len({z["species"] for z in rows})!=len(rows):
            raise ValueError("Duplicate original species MSE")
        n_expect=155 if name=="all155" else 95
        if len(rows)!=n_expect:
            raise ValueError("Original port-forecast species count unexpected")
        orig={}
        for model,column in [
            ("M5_host_composition","m5"),("T0_port_location","t0"),
            ("T1_port_connectivity","t1")
        ]:
            source=[z for z in data["comp"] if z["cohort"]==name
                    and z["model"]==model]
            if len(source)!=1:
                raise ValueError("Original port comparison scope missing "+name+" "+model)
            mse=float(source[0]["species_equal_mse"])
            if abs(mse-np.mean([z[column] for z in rows]))>1e-12:
                raise ValueError("Original prediction model published MSE mismatch")
            orig[model]=mse
        result=est(rows,seed=SEED+len(name))
        result["source_mean_MSE"]=orig
        result["rho_stratified_permutation"]=stratified_perm(
            rows,result["primary_rho_added_units_vs_absolute_M5_minus_T1_mse"],
            seed=SEED+len(name)+100)
        result["quartiles_introduced_added_units"]=quartiles(rows)
        results[name]=result
        for z in rows:
            out=z.copy()
            out["M5_minus_T1_squared_error"]=z["m5"]-z["t1"]
            out["T0_minus_T1_squared_error"]=z["t0"]-z["t1"]
            out["T1_vs_M5_relative_MSE_reduction_pct"]=100*(z["m5"]-z["t1"])/z["m5"]
            joined.append(out)
    # The island records are a subset of the already-fit known butterfly model
    # but have uneven individual source sample sizes.
    allnames={z["species"] for z in data["temporal155"]}
    assert {z["species"] for z in data["island95"]} <= allnames
    assert {z["species"] for z in data["island_nonGB95"]} <= allnames
    if abs(results["all155"]["T1_vs_M5_relative_MSE_reduction_pct"]-.22476045477479722)>1e-9:
        raise ValueError("Original all155 port comparison not faithfully reconstructed")
    if abs(results["strict_islands"]["T1_vs_M5_relative_MSE_reduction_pct"]+.5953721285598235)>1e-9:
        raise ValueError("Original island subset port comparison not faithfully reconstructed")
    if abs(results["strict_islands_excluding_GB_Ireland"]["T1_vs_M5_relative_MSE_reduction_pct"]+.5056564926063656)>1e-9:
        raise ValueError("Original nonGBI insular port comparison not faithfully reconstructed")
    island_sensitivity=[]
    for threshold in (1,10,30,100):
        keep=[z for z in joined if z["group"]=="strict_islands"
              and z["strict_source_island_test_rows"]>=threshold]
        if len(keep)>=20:
            result=est(keep,seed=SEED+threshold,rep=1999)
            island_sensitivity.append({
                "min_strict_target51_rows_per_species":threshold,
                "important_scope_caveat":"Row threshold measured on 51 target islands, not all 25 transport strict islands",
                "n_species":len(keep),
                "rho_added_vs_delta_MSE":result["primary_rho_added_units_vs_absolute_M5_minus_T1_mse"],
                "T1_vs_M5_reduction_pct":result["T1_vs_M5_relative_MSE_reduction_pct"],
            })
    report={
        "study":"species level introduced-host opportunity vs additional maritime port feature prediction value",
        "input_sha256":digests,
        "original_analysis":"retrospective time held port T1 vs already host-informed M5 for same 155 fitted butterflies",
        "can_test_host_only_predictive_gain_M3_to_M5_from_these_files":False,
        "can_test_added_port_connectivity_T1_beyond_M5_from_these_files":True,
        "comparison_window":"2015-2018 predictors / 2019-2022 GBIF recording response; historical port coordinates and host maps not frozen in those periods",
        "results_by_scope":results,
        "island_sample_size_sensitivity":island_sensitivity,
        "critical_limitations":[
            "The test is *T1 port connectivity added to host model M5*, NOT the originally requested per species host information M3 to M5 effect",
            "Species MSE is measured in one selected 155-species original training cohort, no external colonization evidence",
            "Absolute MSE changes are small, may have outliers; rho and bootstrap intervals are descriptive, not phylogenetically adjusted",
            "Species are diet- and ancestry-dependent; the conditional stratified permutation does not address geography, observer effort, or species-specific climate",
            "Island portion includes 95 species but is dominated by Britain and Ireland; remaining direct indexed port support is 25 rows over 3 cells",
            "The port registry and geocoding are 2019/2026 snapshots; quarter indices do not guarantee historical trajectories",
            "Original M5 uses the current host catalogue and does not establish alien host arrival before butterfly response",
            "Missing transport predictors are retained under an unknown flag, so island-only changes may reflect training elsewhere",
            "Original T0 includes port location, T1 adds historical port connectivity; lack of positive island MSE gain does not disprove any biological dispersal",
            "No independent establishment data, plant-local host consumption, or actual ship transport exists in these result tables",
        ],
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"summary.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    with (a.out/"original_species_port_forecast_and_resources.csv").open("w",newline="",encoding="utf-8") as f:
        fields=list(joined[0]);w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(joined)
    print(json.dumps({
        "main":{"n":results["all155"]["n_species"],
                "M5_vs_T1_reduction_pct":results["all155"]["T1_vs_M5_relative_MSE_reduction_pct"],
                "rho_added_vs_connectivity_gain":results["all155"]["primary_rho_added_units_vs_absolute_M5_minus_T1_mse"],
                "rho_bootstrap_CI":results["all155"]["descriptive_4999_species_bootstrap_rho_CI"],
                "partial_rank_rho":results["all155"]["partial_rank_rho_added_vs_M5minusT1_given_7_source_covariates"],
                "stratified_permutation_p":results["all155"]["rho_stratified_permutation"]["exploratory_two_sided_p_centered_on_null"]},
        "islands":{"n":results["strict_islands"]["n_species"],
                   "gain_pct":results["strict_islands"]["T1_vs_M5_relative_MSE_reduction_pct"],
                   "rho_added":results["strict_islands"]["primary_rho_added_units_vs_absolute_M5_minus_T1_mse"]},
        "noGBI":{"n":results["strict_islands_excluding_GB_Ireland"]["n_species"],
                 "gain_pct":results["strict_islands_excluding_GB_Ireland"]["T1_vs_M5_relative_MSE_reduction_pct"],
                 "rho_added":results["strict_islands_excluding_GB_Ireland"]["primary_rho_added_units_vs_absolute_M5_minus_T1_mse"]},
        "output":str(a.out)},ensure_ascii=False,flush=True)

if __name__=="__main__":
    main()
