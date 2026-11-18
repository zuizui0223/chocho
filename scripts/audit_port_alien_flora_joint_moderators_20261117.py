#!/usr/bin/env python3
"""Identify whether port × island area and port × native flora moderators remain
separately distinguishable when entered in the SAME WCVP island regression.

This is POST-HOC diagnostic after observing both moderators separately.
No route inference, butterfly establishment, or host-specific port introduction.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import audit_island_port_alien_flora_native_area_moderators_20261117 as previous

N_PERM=1499
N_BOOT=999
SEED=20261117+250

def design(items,*,specific=False,effort=False):
    y,base,_,names,ma=previous.build_design(items,"island_area",specific=specific,effort=effort)
    _,base_n,_,names_n,mn=previous.build_design(items,"native_flora",specific=specific,effort=effort)
    if not (np.array_equal(base,base_n) and names==names_n):
        raise RuntimeError("Joint moderators do not share original nuisance source and port model")
    ports=np.array([r["ln2ports"] for r in items],dtype=float)
    Z=np.column_stack([base,ports*ma,ports*mn])
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        raise RuntimeError("Joint port×area/native fixed-effect source regression rank deficient")
    return y,base,Z,names

def bootstrap(items,specific,effort):
    rng=np.random.default_rng(SEED+len(items)+(1 if specific else 0)+(50 if effort else 0))
    bylevel=defaultdict(list)
    for i,r in enumerate(items):bylevel[r["level1"]].append(i)
    groups=[np.array(v,dtype=int) for v in bylevel.values()]
    vals=[]
    for _ in range(N_BOOT):
        idx=np.concatenate([rng.choice(g,len(g),replace=True) for g in groups])
        try:
            y,x0,x,names=design([items[int(i)] for i in idx],specific=specific,effort=effort)
        except RuntimeError:
            continue
        b,_=previous.fit(y,x)
        vals.append([float(b[-2]),float(b[-1])])
    if len(vals)<400:raise RuntimeError("Too few rank-complete joint regional bootstrap draws")
    arr=np.array(vals)
    return arr,len(vals)

def interaction_conditional_perm(items,*,specific,effort,which):
    # Freedman-Lane null keeps the other interaction PLUS all main effects.
    y,base,full,names=design(items,specific=specific,effort=effort)
    selected=full.shape[1]-2+which
    X0=np.delete(full,selected,axis=1)
    if np.linalg.matrix_rank(X0)!=X0.shape[1]:
        raise RuntimeError("Joint conditional null rank deficient")
    nullbeta,nullres=previous.fit(y,X0)
    expected=X0@nullbeta
    obs=float(previous.fit(y,full)[0][selected])
    groups=previous.blocks(items,within_area_half=True)
    rng=np.random.default_rng(SEED+which+len(items)*7+(13 if specific else 0)+(59 if effort else 0))
    ps=[]
    for _ in range(N_PERM):
        yp=expected+nullres
        yp=yp.copy()
        for g in groups:
            yp[g]=expected[g]+nullres[rng.permutation(g)]
        b,_=previous.fit(yp,full)
        ps.append(float(b[selected]))
    arr=np.array(ps)
    return {"p_two_sided_descriptive":float((1+np.sum(abs(arr)>=abs(obs)-1e-12))/(N_PERM+1)),
            "null_95pct":[float(np.quantile(arr,.025)),float(np.quantile(arr,.975))],
            "observed_beta":obs}

def fit_joint(items,*,cohort,specific=False,effort=False):
    y,base,full,names=design(items,specific=specific,effort=effort)
    b,res=previous.fit(y,full)
    n=len(items);pn=names.index("log2_1plus_WPI_ports")
    boot,nb=bootstrap(items,specific,effort)
    tests=[interaction_conditional_perm(items,specific=specific,effort=effort,which=i) for i in range(2)]
    b_area,b_native=map(float,b[-2:])
    # Cross-moderator collinearity is a quantifiable source selection limit.
    mx=np.array([r["lnarea"] for r in items]);mn=np.array([r["lnnative"] for r in items])
    partial=float(1-res@res/(previous.fit(y,base)[1]@previous.fit(y,base)[1]))
    loolevel={}
    for level in sorted(set(r["level1"] for r in items),key=int):
        sub=[r for r in items if r["level1"]!=level]
        try:
            yy,xx0,xx,nn=design(sub,specific=specific,effort=effort)
            coeff=previous.fit(yy,xx)[0]
            loolevel[level]={"n":len(sub),"beta_area":float(coeff[-2]),"beta_native":float(coeff[-1])}
        except RuntimeError as exc:
            loolevel[level]={"status":"rank_deficient","reason":str(exc)}
    return {
      "cohort":cohort,"response":"HOSTS introduced larval-host taxa conditional on all alien flora" if specific else "all accepted WCVP introduced flora",
      "n":n,"number_TDWG_Level1_groups":len(set(r["level1"] for r in items)),
      "botanical_GBIF_effort_adjusted":effort,
      "port_main_effect_at_median_area_and_native_richness":float(b[pn]),
      "moderator_coefficients_standardized":{"port_x_log_area":b_area,
            "port_x_log_native_flora":b_native},
      "interaction_conditional_Freedman_Lane_p_with_other_interaction_in_null":{
         "port_x_log_area":tests[0]["p_two_sided_descriptive"],
         "port_x_log_native_flora":tests[1]["p_two_sided_descriptive"]},
      "interaction_conditional_null_95pct":{
          "port_x_log_area":tests[0]["null_95pct"],
          "port_x_log_native_flora":tests[1]["null_95pct"]},
      "island_within_Level1_bootstrap_95pct":{
          "port_x_log_area":[float(np.quantile(boot[:,0],.025)),float(np.quantile(boot[:,0],.975))],
          "port_x_log_native_flora":[float(np.quantile(boot[:,1],.025)),float(np.quantile(boot[:,1],.975))]},
      "bootstrap_valid_draws":nb,
      "log_area_native_richness_correlation":float(np.corrcoef(mx,mn)[0,1]),
      "incremental_joint_moderators_partial_R2":partial,
      "leave_one_TDWG_Level1_out":loolevel,
      "source_controlled_outcome_noncausal":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--flora",required=True,type=Path)
    ap.add_argument("--hosts",required=True,type=Path)
    ap.add_argument("--effort",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    a=ap.parse_args()
    for name,p in (("flora",a.flora),("hosts",a.hosts),("effort",a.effort)):
        if previous.gitblob(p)!=previous.SOURCE_SHA[name]:
            raise RuntimeError("Original source git blob drift before joint model "+name)
    full=previous.vec_data(previous.rows(a.flora),previous.rows(a.hosts),previous.rows(a.effort))
    previous.REFERENCE={
      "area_median":float(np.median([x["lnarea"] for x in full])),
      "area_sd":float(np.std([x["lnarea"] for x in full],ddof=1)),
      "native_median":float(np.median([x["lnnative"] for x in full])),
      "native_sd":float(np.std([x["lnnative"] for x in full],ddof=1))}
    active=[r for r in full if r["has_host"]]
    efforts=[r for r in active if r["has_effort"]]
    cases=[("all_126_WCVP",full,False,False),
           ("same_123_WCVP_HOSTS",active,False,False),
           ("same_123_HOSTS_given_all_alien_flora",active,True,False),
           ("same_98_WCVP_with_GBIF_effort",efforts,False,True),
           ("same_98_HOSTS_given_alien_flora_and_effort",efforts,True,True)]
    reports=[]
    for label,items,specific,effort in cases:
        z=fit_joint(items,cohort=label,specific=specific,effort=effort)
        reports.append(z)
        print("PORT_FLORA_JOINT_MODERATOR",json.dumps({
          "cohort":label,"n":z["n"],
          "beta":z["moderator_coefficients_standardized"],
          "p_conditional":z["interaction_conditional_Freedman_Lane_p_with_other_interaction_in_null"],
          "boot_95":z["island_within_Level1_bootstrap_95pct"],
          "correlation":z["log_area_native_richness_correlation"]
        }),flush=True)
    if len(reports)!=5 or [r["n"] for r in reports]!=[126,123,123,98,98]:
        raise RuntimeError("Original three source datasets or five joint cohorts changed")
    report={
      "schema":"chocho_posthoc_independent_two_WPI_moderators_island_area_native_flora_v1",
      "status":"JOINT_AREA_NATIVE_FLORA_PORT_MODERATION_IS_OBSERVATIONAL_POSTHOC",
      "source_git_blobs":previous.SOURCE_SHA,
      "source_island_sample":[126,123,98],
      "original_source_main_model_beta":.23347121981392502,
      "Freedman_Lane_permutations_each_interaction_conditional_on_other":N_PERM,
      "Level1_island_bootstraps":N_BOOT,
      "posthoc_after_two_single_mods":"These joint tests were added after both single moderator analyses returned negative signs; not new preregistered tests",
      "results":reports,
      "limitations":[
        "Area and native WCVP flora are correlated proxies for island species pools and source/observer opportunity.",
        "No observational WPI port-plant coefficient is causal shipping traffic or independent of population, habitat and collection intensity.",
        "Original plant richness is present-day and not measured before indexed port construction.",
        "A significant interaction conditional on other moderator is still post-hoc and conditional on botanical regional observations.",
        "Original HOSTS literature is a frozen catalogue, not independent local larval host feeding or butterfly colonization.",
      ]
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"source_verified_joint_port_island_moderator_sensitivity.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("JOINT_PORT_MODERATION_SOURCE_GATE",json.dumps({"n_models":len(reports),"results":[{
       "cohort":r["cohort"],"area":r["moderator_coefficients_standardized"]["port_x_log_area"],
       "native":r["moderator_coefficients_standardized"]["port_x_log_native_flora"]}
       for r in reports]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
