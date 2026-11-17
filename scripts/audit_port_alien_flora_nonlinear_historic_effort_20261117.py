#!/usr/bin/env python3
"""Post-hoc falsification of port–alien-flora moderation in fixed island sources.

Reuses Git-SHA verified original 126/123/98 WGSRPD3 island sources from
PR #173. Tests flexible NONLINEAR port/area/native main functions versus
port×area/native interactions; checks 98-island contemporaneous / pre-1900
GBIF plant recording confounding and port-present support. Last test is a
HISTORICAL DOCUMENTATION DIAGNOSTIC, not a valid causal placebo.

The explicit goal is to find when earlier source-positive interactions become
model-dependent, not to pick a favourable specification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

import audit_island_port_alien_flora_native_area_moderators_20261117 as previous

SEED=20261117 + 900
N_PERM=1499
N_BOOT=699

def fit(y,x):
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    return b,y-x@b

def array(rows,k):
    return np.array([float(r[k]) for r in rows],dtype=float)

def matrix(rows,outcome,*,flexible,effort_controls=(),host_conditioned=False,
           include_interactions=True):
    n=len(rows)
    levs=sorted({x["level1"] for x in rows},key=int)
    if n<65 or len(levs)<6:
        raise RuntimeError(f"Small island geographic support {n}/{len(levs)}")
    p=array(rows,"ln2ports")
    area=array(rows,"lnarea")
    native=array(rows,"lnnative")
    aa=(area-previous.REFERENCE["area_median"])/previous.REFERENCE["area_sd"]
    nn=(native-previous.REFERENCE["native_median"])/previous.REFERENCE["native_sd"]
    y=array(rows,outcome)
    vals=[np.ones(n),area,native,p]
    labels=["intercept","log_area","log1p_native_flora","log2_1plus_ports"]
    if host_conditioned:
        vals.append(array(rows,"lnalien"))
        labels.append("all_WCVP_introduced_flora")
    for e in effort_controls:
        vals.append(array(rows,e))
        labels.append(e)
    if flexible:
        # Orthogonalize these only by regression fit (not outcome fitting).
        # The robust claim cannot rest on a mandatory linear port or area
        # response while other island plant opportunity scales nonlinearly.
        vals.extend([p*p,aa*aa,nn*nn,aa*nn])
        labels.extend(["port_log2_squared","log_area_standard_squared",
                       "log_native_standard_squared","area_x_native_main"])
    for lev in levs[1:]:
        vals.append(np.array([1. if r["level1"]==lev else 0. for r in rows]))
        labels.append("level1_"+lev)
    red=np.column_stack(vals)
    if np.linalg.matrix_rank(red,tol=1e-10)<red.shape[1]:
        raise RuntimeError(f"Main terms rank deficient n={n} X={red.shape}")
    if include_interactions:
        X=np.column_stack([red,p*aa,p*nn])
        labels+=["port_x_area","port_x_native"]
    else:
        X=red
    if np.linalg.matrix_rank(X,tol=1e-10)<X.shape[1]:
        raise RuntimeError(f"Interactions rank deficient n={n} X={X.shape}")
    return y,red,X,labels

def groups(items, *,by_area):
    by=defaultdict(list)
    if by_area:
        within=defaultdict(list)
        for i,r in enumerate(items):within[r["level1"]].append(i)
        median_by={k:float(np.median([items[i]["lnarea"] for i in idx]))
                   for k,idx in within.items()}
        for i,r in enumerate(items):
            by[r["level1"]+"|"+str(int(r["lnarea"]>median_by[r["level1"]]))].append(i)
    else:
        for i,r in enumerate(items):by[r["level1"]].append(i)
    return [np.array(z,dtype=int) for z in by.values() if len(z)>1]

def bootstrap(items,outcome,flexible,effort_controls,host_conditioned,rng):
    # Resample ISLAND botanical units within original Level1 groups, preserving
    # even singleton groups. Do not treat HOSTS species×island cells as trials.
    grouped=defaultdict(list)
    for i,r in enumerate(items):grouped[r["level1"]].append(i)
    draws=[]
    groups=[np.asarray(v,dtype=int) for v in grouped.values()]
    for k in range(N_BOOT):
        sel=np.concatenate([rng.choice(gr,len(gr),replace=True) for gr in groups])
        tmp=[items[int(i)] for i in sel]
        try:
            yy,_,xx,_=matrix(tmp,outcome,flexible=flexible,
                             effort_controls=effort_controls,host_conditioned=host_conditioned)
            b,_=fit(yy,xx)
            draws.append(b[-2:].tolist())
        except RuntimeError:continue
    if len(draws)<300:return None,len(draws)
    draws=np.asarray(draws)
    return {
      "area":[float(np.quantile(draws[:,0],.025)),float(np.quantile(draws[:,0],.975))],
      "native":[float(np.quantile(draws[:,1],.025)),float(np.quantile(draws[:,1],.975))]
    },len(draws)

def permute(items,y,X,which,rng):
    # Freedman-Lane: null has the other interaction, all controlled main
    # effects (including nonlinear terms), and the identical source regions.
    key=X.shape[1]-2+which
    XR=np.delete(X,key,axis=1)
    bnull,residual=fit(y,XR)
    mu=XR@bnull
    permgroups=groups(items,by_area=True)
    vals=[]
    for _ in range(N_PERM):
        yp=mu+residual
        yp=yp.copy()
        for g in permgroups:
            yp[g]=mu[g]+residual[rng.permutation(g)]
        b,_=fit(yp,X)
        vals.append(float(b[key]))
    vals=np.asarray(vals)
    return vals

def run(items,label,outcome,*,flexible=False,effort_controls=(),host_conditioned=False):
    rng=np.random.default_rng(SEED+sum(ord(ch) for ch in label)+len(items))
    try:
        Y,XR,X,names=matrix(items,outcome,flexible=flexible,
                            effort_controls=effort_controls,host_conditioned=host_conditioned)
    except RuntimeError as e:
        return {"case":label,"n":len(items),"status":"MODEL_UNIDENTIFIABLE",
                "message":str(e)}
    b,r=fit(Y,X)
    br,rr=fit(Y,XR)
    z={}
    conf,draws=bootstrap(items,outcome,flexible,effort_controls,host_conditioned,rng)
    for i,axis in enumerate(("area","native")):
        null=permute(items,Y,X,i,rng)
        beta=float(b[-2+i])
        z[axis]={
            "beta_port_by_standardized_moderator":beta,
            "conditional_permutation_descriptive_p":float((1+sum(abs(null)>=abs(beta)-1e-12))/(N_PERM+1)),
            "conditional_null95":[float(np.quantile(null,.025)),float(np.quantile(null,.975))],
            "island_stratified_bootstrap95":conf[axis] if conf else None,
            "frozen_same_parameter_inclusion":True
        }
    ps=np.array([z[k]["conditional_permutation_descriptive_p"] for k in ("area","native")])
    small=int(np.argmin(ps))
    large=1-small
    holm=np.zeros(2)
    holm[small]=min(1.,2*ps[small])
    holm[large]=max(holm[small],ps[large])
    z["area"]["exploratory_Holm2_p"]=float(holm[0])
    z["native"]["exploratory_Holm2_p"]=float(holm[1])
    l1={}
    for lev in sorted({x["level1"] for x in items},key=int):
        sub=[r for r in items if r["level1"]!=lev]
        try:
            _,_,x_sub,_=matrix(sub,outcome,flexible=flexible,effort_controls=effort_controls,
                               host_conditioned=host_conditioned)
            l1[lev]=list(map(float,fit(array(sub,outcome),x_sub)[0][-2:]))
        except RuntimeError:
            l1[lev]=None
    return {
       "status":"FIT_DESCRIPTIVE_CROSS_SECTIONAL_NOT_CAUSAL",
       "case":label,"n":len(items),"botanical_regions":len({x["level1"] for x in items}),
       "WPI_port_positive_islands":sum(x["ports"]>0 for x in items),
       "source_response":outcome,"flexible_nuisance_main_effects":bool(flexible),
       "effort_controls":list(effort_controls),"host_conditioned_on_total_alien_flora":host_conditioned,
       "model_df":X.shape[1],"main_port_beta_at_reference":float(b[names.index("log2_1plus_ports")]),
       "moderation":z,
       "joint_moderator_incremental_partial_R2":float(1-r@r/(rr@rr)),
       "ln_area_ln_native_flora_correlation":float(np.corrcoef(array(items,"lnarea"),array(items,"lnnative"))[0,1]),
       "valid_original_level1_bootstraps":draws,"leave_one_level1_beta_area_native":l1,
       "limits":[
          "This sensitivity was designed after initial port×area/native interactions were found; do not promote p to independent confirmation.",
          "With p^2, a^2, native^2 and area×native main terms, interactions are conditional on a flexible source-association curve.",
          "Ports-only changes the estimand from existence plus count of registered ports to within-ported-region port intensity.",
          "Plant GBIF recording effort is a count of documented occurrences, not a probability of observing a particular nonnative plant.",
          "Pre1900 GBIF plant indexing is a historical documentation diagnostic, not exogenous port assignment or an observed plant introduction date.",
          "The same WCVP species occurrences underlie both outcome and native/alien flora metrics and include historical botanical bias.",
          "No actual historical vessel OD, seed transport, butterfly larval food consumption or colonisation is observed."
       ],
    }

def main():
    a=argparse.ArgumentParser()
    for name in ("flora","hosts","effort","out-dir"):
        a.add_argument("--"+name,required=True,type=Path)
    args=a.parse_args()
    for name,path in (("flora",args.flora),("hosts",args.hosts),("effort",args.effort)):
        if previous.gitblob(path)!=previous.SOURCE_SHA[name]:
            raise RuntimeError(f"Wrong locked source original Git blob for {name}")
    z=previous.vec_data(previous.rows(args.flora),previous.rows(args.hosts),previous.rows(args.effort))
    previous.REFERENCE={
      "area_median":float(np.median([r["lnarea"] for r in z])),
      "area_sd":float(np.std([r["lnarea"] for r in z],ddof=1)),
      "native_median":float(np.median([r["lnnative"] for r in z])),
      "native_sd":float(np.std([r["lnnative"] for r in z],ddof=1))
    }
    h=[r for r in z if r["has_host"]]
    e=[r for r in h if r["has_effort"]]
    # Mandatory backward-compatibility controls with the previous exact
    # joint-results CI #44003404503 before a nonlinear refit.
    yy,red,xx,names=matrix(h,"lnalien",flexible=False)
    source_b=fit(yy,xx)[0][-2:]
    if not (abs(source_b[0]-(-.04254168022063994))<1e-9 and
            abs(source_b[1]-(-.13899715182792102))<1e-9):
        raise RuntimeError("Frozen 123-island joint moderations not reproduced")
    specs=[
      ("all126_nonlinear_main",z,"lnalien",True,(),False),
      ("host123_linear_reproduction",h,"lnalien",False,(),False),
      ("host123_nonlinear_main",h,"lnalien",True,(),False),
      ("host123_nonlinear_ports_only",[r for r in h if r["ports"]>0],"lnalien",True,(),False),
      ("effort98_nonlinear_modern",e,"lnalien",True,("lneffort",),False),
      ("effort98_nonlinear_pre1900",e,"lnalien",True,("lnpre1900effort",),False),
      ("effort98_nonlinear_both_efforts",e,"lnalien",True,("lneffort","lnpre1900effort"),False),
      ("effort98_nonlinear_ports_only",[r for r in e if r["ports"]>0],"lnalien",True,("lneffort",),False),
      ("effort98_historical_botanical_indexing_negative_timing_diagnostic",
       e,"lnpre1900effort",True,(),False),
      ("host123_nonlinear_hosts_given_all_flora",h,"lnhost",True,(),True),
      ("effort98_nonlinear_hosts_given_all_flora",e,"lnhost",True,("lneffort",),True)
    ]
    outcomes=[]
    for label,pop,response,nonlinear,effort,host in specs:
        z0=run(pop,label,response,flexible=nonlinear,effort_controls=effort,host_conditioned=host)
        outcomes.append(z0)
        print("PORT_NATIVE_FLORA_NONLINEAR_SENSITIVITY",json.dumps({
            "name":label,"n":z0["n"],"status":z0["status"],
            "area":z0.get("moderation",{}).get("area"),
            "native":z0.get("moderation",{}).get("native")
        }),flush=True)
    report={
      "schema":"chocho_original_island_port_alien_flora_moderator_nonlinear_and_historical_index_gate_v1",
      "status":"EXPLORATORY_CURVATURE_AND_HISTORIC_BOTANICAL_EFFORT_GATE_NOT_SEED_TRANSPORT",
      "source_Git_blobs":previous.SOURCE_SHA,
      "reproduced_linear_joint123_area_native":[float(x) for x in source_b],
      "source_baseline_single_moderator_CI":44003165119,
      "source_joint_moderator_CI":44003404503,
      "n_original_island_cohorts":[126,123,98],
      "permutation_draws":N_PERM,"bootstrap_draws":N_BOOT,
      "source_observational_GBIF_pre1900_not_true_pre_port_control":True,
      "cases":outcomes,
      "limits":[
        "All fits and falsifications are post-hoc sensitivity checks on the original 126/123/98 island WCVP source; no new ecological replication or preregistration.",
        "The 123-island observed port–alien flora baseline positive main effect and the original single moderator results remain archived independent of these functional-form diagnostics.",
        "A model of pre1900 botanical index density uses different response from WCVP accepted plant flora and is not a historical seed-arrival outcome.",
        "An association in an old plant-index observation is a warning for geographic sampling and long-term human access, not proof there was no later plant introduction.",
        "There is no direct link between contemporary WPI port counts and actual vessel origin-destination or butterfly colonisation.",
      ],
    }
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"nonlinear_port_native_flora_historic_indexing_audit.json").write_text(
      json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print("PORT_NONLINEAR_HISTORIC_EFFORT_SOURCE_GATE",json.dumps({
      "status":report["status"],"cases":[{
        "case":r["case"],"n":r["n"],"status":r["status"],
        "area":r.get("moderation",{}).get("area",{}).get("beta_port_by_standardized_moderator"),
        "native":r.get("moderation",{}).get("native",{}).get("beta_port_by_standardized_moderator")
      } for r in outcomes]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
