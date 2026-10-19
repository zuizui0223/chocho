#!/usr/bin/env python3
"""Follow-up macroecology of butterfly host-resource climate-analog gains.

Consumes only the derived, source-verified output of
analyze_butterfly_climate_analog_macro_20261019.py plus the
pinned original 239/339 species descriptors.
This is a post-hoc, negative-control analysis; no new individual-paper review.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from analyze_global_butterfly_habitat_filter_20261019 import (
    SEED, _columns, annotate, cluster_se, read_inputs, wild_genus_signflip
)
REGIONS={
    "1":"Europe","2":"Africa","3":"Temperate Asia","4":"Tropical Asia",
    "5":"Australasia","6":"Pacific","7":"Northern America",
    "8":"Southern America","9":"Antarctic"
}

def joined(der,macro,*,minK=5,quality=False,nonzero=False):
    mapped={r["species"]:r for r in der}
    out=[]
    for base in annotate(macro):
        source=mapped.get(base["species"])
        if not source or source["eligible_climate_analog_non_native_regions"]<minK:
            continue
        if quality and not base["adequate"]:
            continue
        if nonzero and source["added_regions"]==0:
            continue
        r={**base}
        r["climate_fill_logit"]=float(source["analog_fill_logit"])
        r["climate_enrichment"]=source["analog_excess_share"]
        r["dominant_continent"]=source["dominant_native_resource_level1"]
        r["dominant_share"]=float(source["dominant_native_level1_share"])
        r["center_bio1"]=float(source["native_host_region_mean_bio1"])
        r["center_bio12"]=float(source["native_host_region_mean_bio12"])
        r["analog_K"]=source["eligible_climate_analog_non_native_regions"]
        if r["center_bio12"]<0:
            raise ValueError("invalid negative annual precipitation")
        out.append(r)
    return out


def base_gear(rows,all_families,*,extra_geo=True,keep_k=False):
    assert rows
    cols=np.array([_columns(r,all_families) for r in rows],dtype=float)
    added=[]
    if extra_geo:
        present=set(r["dominant_continent"] for r in rows)
        l1=sorted(present)
        if len(l1)>1:
            for g in l1[1:]:
                added.append(np.array([r["dominant_continent"]==g for r in rows],dtype=float))
        # Additional continuous geography context (scaling fixed before fitting).
        added.extend([
            np.array([r["center_bio1"]/10. for r in rows]),
            np.array([math.log1p(r["center_bio12"])/6. for r in rows]),
            np.array([r["dominant_share"] for r in rows]),
        ])
    if keep_k:
        added.append(np.array([math.log1p(r["analog_K"]) for r in rows]))
    if added:
        cols=np.column_stack([cols,*added])
    # Drop linearly dependent nuisance columns, do not select by response.
    keep=[]
    for j in range(cols.shape[1]):
        subset=keep+[j]
        if np.linalg.matrix_rank(cols[:,subset])==len(subset):
            keep=subset
    return cols[:,keep]


def coefficient(rows,families,*,y_key="climate_fill_logit",geo=True,include_K=False,
                with_bootstrap=True,n_boot=1499,seed=SEED):
    y=np.array([r[y_key] for r in rows],dtype=float)
    if not np.all(np.isfinite(y)):raise ValueError("missing outcome")
    z=base_gear(rows,families,extra_geo=geo,keep_k=include_K)
    x=np.column_stack([z,np.array([r["exposure"] for r in rows],dtype=float)])
    rank=np.linalg.matrix_rank(x)
    if rank!=x.shape[1]:
        return {"status":"NOT_IDENTIFIABLE_EXPOSURE_RANK_LOST",
                "n":len(rows),"design_rank":int(rank),"columns":x.shape[1]}
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    resid=y-x@b
    se,ngen=cluster_se(x,resid,[r["genus"] for r in rows])
    estimate=float(b[-1])
    label="genus_cluster_approx"
    report={
        "status":"EXPLORATORY_ASSOCIATION",
        "n":len(rows),
        "genera":ngen,
        "associated_n":sum(r["exposure"]==1 for r in rows),
        "avoidant_n":sum(r["exposure"]==0 for r in rows),
        "design_rank":int(rank),
        "covariate_columns":x.shape[1],
        "with_geo_level1_and_native_bio1_bio12":geo,
        "with_climate_analogue_capacity_K_covariate":include_K,
        "response":y_key,
        "estimate_associated_minus_avoidant":estimate,
        "genus_cluster_se":se,
        "genus_cluster_approx_95ci":[estimate-1.96*se,estimate+1.96*se],
        "wild_genus_residual_signflip":wild_genus_signflip(
            z,y,x,[r["genus"] for r in rows],seed=seed,n=4999),
    }
    if with_bootstrap:
        by_gen=defaultdict(list)
        for r in rows:by_gen[r["genus"]].append(r)
        genus=sorted(by_gen)
        rng=np.random.default_rng(seed+73)
        vals=[]
        for _ in range(n_boot):
            draws=rng.integers(0,len(genus),size=len(genus))
            sample=[r for i in draws for r in by_gen[genus[int(i)]]]
            if len({r["exposure"] for r in sample})<2:continue
            controls=base_gear(sample,families,extra_geo=geo,keep_k=include_K)
            Xs=np.column_stack([controls,np.array([r["exposure"] for r in sample])])
            if np.linalg.matrix_rank(Xs)<Xs.shape[1]:
                continue
            Y=np.array([r[y_key] for r in sample],dtype=float)
            vals.append(float(np.linalg.lstsq(Xs,Y,rcond=None)[0][-1]))
        if len(vals)<n_boot*.90:
            report["bootstrap"]={"status":"INSUFFICIENT_IDENTIFIABLE_DRAW_SUPPORT",
                                 "valid_draws":len(vals),"requested":n_boot}
        else:
            report["bootstrap"]={
                "status":"EXPLORATORY",
                "valid_draws":len(vals),"requested":n_boot,
                "percentile_95ci":[float(np.quantile(vals,.025)),float(np.quantile(vals,.975))]
            }
    return report


def genus_cv(rows,families,*,response="climate_fill_logit",add_geo=False,
             include_K=False,seed=SEED):
    by_gen=defaultdict(list)
    for i,r in enumerate(rows):
        by_gen[r["genus"]].append(i)
    genera=sorted(by_gen)
    assert len(genera)>40
    scores=[]
    for rep in range(10):
        rng=np.random.default_rng(seed+rep*541)
        randomized=[genera[int(i)] for i in rng.permutation(len(genera))]
        folds=[[] for _ in range(5)];sizes=[0]*5
        for g in randomized:
            j=int(np.argmin(sizes))
            folds[j].append(g);sizes[j]+=len(by_gen[g])
        pred0=np.zeros(len(rows));pred1=np.zeros(len(rows))
        obs=np.array([r[response] for r in rows],dtype=float)
        # Use SAME full original set of predictors before reduced nuisance
        # rank-selection; unique train-only groups then remain fixed zeros.
        levels=sorted({r["dominant_continent"] for r in rows})
        def controls(v,geo):
            out=np.array([_columns(r,families) for r in v],dtype=float)
            if geo:
                extra=[]
                for g in levels[1:]:
                    extra.append(np.array([r["dominant_continent"]==g for r in v],dtype=float))
                extra.extend([
                    np.array([r["center_bio1"]/10 for r in v]),
                    np.array([math.log1p(r["center_bio12"])/6 for r in v]),
                    np.array([r["dominant_share"] for r in v]),
                ])
                out=np.column_stack([out,*extra])
            if include_K:
                out=np.column_stack([out,np.array([math.log1p(r["analog_K"]) for r in v])])
            return out
        for fg in folds:
            test_g=set(fg)
            ids_test=[i for i,r in enumerate(rows) if r["genus"] in test_g]
            ids_train=[i for i,r in enumerate(rows) if r["genus"] not in test_g]
            tr=[rows[i] for i in ids_train];te=[rows[i] for i in ids_test]
            Z0=controls(tr,add_geo)
            Zt=controls(te,add_geo)
            X0=np.column_stack([Z0,np.array([r["exposure"] for r in tr])])
            Xt=np.column_stack([Zt,np.array([r["exposure"] for r in te])])
            y=np.array([r[response] for r in tr])
            b0=np.linalg.lstsq(Z0,y,rcond=None)[0]
            b1=np.linalg.lstsq(X0,y,rcond=None)[0]
            pred0[ids_test]=Zt@b0
            pred1[ids_test]=Xt@b1
        e0=float(np.mean((obs-pred0)**2));e1=float(np.mean((obs-pred1)**2))
        scores.append({"baseline_MSE":e0,"plus_habitat_MSE":e1,
                       "relative_improvement":(e0-e1)/e0})
    return {
        "repeats":10,"folds":5,"heldout_unit":"butterfly_genus",
        "baseline_MSE":float(np.mean([s["baseline_MSE"] for s in scores])),
        "plus_disturbance_MSE":float(np.mean([s["plus_habitat_MSE"] for s in scores])),
        "mean_relative_MSE_improvement":float(np.mean([s["relative_improvement"] for s in scores])),
        "repeats_improved":sum(s["relative_improvement"]>0 for s in scores),
        "repeats_details":scores
    }


def native_region_checks(rows):
    out={}
    for g in sorted({r["dominant_continent"] for r in rows}):
        x=[r for r in rows if r["dominant_continent"]==g]
        A=[r for r in x if r["exposure"]==1]
        B=[r for r in x if r["exposure"]==0]
        out[REGIONS.get(g,g)]={
            "n":len(x),"assoc":len(A),"avoid":len(B),
            "unadjusted_logit_diff": (
                float(np.mean([r["climate_fill_logit"] for r in A])-
                np.mean([r["climate_fill_logit"] for r in B]))
                if len(A)>=1 and len(B)>=1 else None),
            "strict_min3_each":len(A)>=3 and len(B)>=3,
        }
    return out


def analyse(derived,macro):
    assert derived["scope"].startswith("239 focal butterflies")
    families=sorted({x["family"] for x in macro})
    results={}
    for tag,part in (
        ("main_5_95",derived["main_threshold_5_to_95"]),
        ("sensitivity_10_90",derived["sensitivity_threshold_10_to_90"]),
    ):
        v=part["species_derived"]
        selected=joined(v,macro,minK=5)
        quality=joined(v,macro,minK=5,quality=True)
        positive=joined(v,macro,minK=5,nonzero=True)
        results[tag]={
            "observed_versus_continent_benchmark":part["overview"]["by_level1_introduced_additions"],
            "primary_eligible_trait_count":len(selected),
            "geographic_continent_composition":native_region_checks(selected),
            "old_covariates_only":coefficient(selected,families,geo=False,n_boot=899),
            "plus_native_continent_and_bio1_bio12":coefficient(selected,families,geo=True,n_boot=1499),
            "plus_native_continent_climate_and_candidate_K":
                coefficient(selected,families,geo=True,include_K=True,n_boot=899),
            "plus_native_continent_climate_adequate_HOSTS_only":
                coefficient(quality,families,geo=True,n_boot=899),
            "climate_enrichment_nonzero_only_geo_adjusted":
                coefficient(positive,families,geo=True,y_key="climate_enrichment",n_boot=899),
            "out_of_genus_CV_baseline_vs_disturbance_including_geo":
                genus_cv(selected,families,add_geo=True),
            "note":"Posthoc geographic covariates and conditional climate enrichment. These tests cannot establish causal plant globalization."
        }
    return {
        "schema":"chocho_20261019_climate_analog_macro_geographic_control_v01",
        "status":"POSTHOC_239_BUTTERFLY_GEOGRAPHY_AND_CLIMATE_SENSITIVITY",
        "source":"Source-verified original 239 ×365 model, existing 339 habitat descriptors; derived JSON only",
        "baseline_output_schema":derived["schema"],
        "by_cutoff":results,
        "claims_not_identified":[
            "Actual butterfly population expansion",
            "Local viability or larval recruitment",
            "Introduction timing or anthropogenic causes",
            "Plant-network degree-preserving climatic null",
            "Trait representativeness of all Lepidoptera",
            "Area-weighted global habitat coverage"
        ]
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--derived-json",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    derived=json.loads(args.derived_json.read_text(encoding="utf-8"))
    macro,traits,digests=read_inputs()
    result=analyse(derived,macro)
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    for key,x in result["by_cutoff"].items():
        out={}
        for test in [
            "old_covariates_only",
            "plus_native_continent_and_bio1_bio12",
            "plus_native_continent_climate_and_candidate_K",
            "plus_native_continent_climate_adequate_HOSTS_only",
            "climate_enrichment_nonzero_only_geo_adjusted"
        ]:
            v=x[test]
            out[test]={
                "status":v["status"],
                "n":v["n"],
                "estimate":v.get("estimate_associated_minus_avoidant"),
                "CI":v.get("genus_cluster_approx_95ci"),
                "bootstrap_CI":v.get("bootstrap",{}).get("percentile_95ci"),
            }
        out["out_of_genus_cv"]={
            k:v for k,v in x["out_of_genus_CV_baseline_vs_disturbance_including_geo"].items()
            if k!="repeats_details"
        }
        out["continent_breakdown"]=x["observed_versus_continent_benchmark"]
        out["dominant_native_region_strata"]=x["geographic_continent_composition"]
        print(json.dumps({key:out},sort_keys=True))
