#!/usr/bin/env python3
"""Source-linked original 239 butterfly potential hosts vs same-site-year larval field use.

A matched-site/year exploratory source model, not biological proof of host
fitness, same-species invasion, competition, or larval colonization.
Field data are all Lepidoptera (moth/butterfly) caterpillars; host source is
the original 239 butterflies, so positive HOSTS exact linkage is a proxy.
"""
from __future__ import annotations
import argparse,json,hashlib,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

from audit_field_raw_site_year_caterpillar_20261027 import (
    raw_filter, read_inferred, read_origin,site_year_matched, aggregate_matched
)
from analyze_field_caterpillar_host_realization_filter_20261027 import (
    accepted_global_links,solve_logit,expit,
)
from analyze_hosts_location_wcvp_origin_20261027 import read_rows,sha256

SCHEMA="chocho_original_239_field_siteyear_HOSTS_portfolio_incremental_prediction_20261027_v01"
EXTERNAL_COMMIT="2c5690ac7c6dba7eca506da886c49781c37617bd"
SOURCE_SHA="408cb206d7ce4056cee1fbf614762e1a6efb3db4650ff7c2663a7d6d6003e323"


def aggregate_groups(surveys,sites_years):
    accepted={(x["site"],x["year"]) for x in sites_years}
    grouped=defaultdict(lambda:{"n":0,"pos":0})
    for r in surveys:
        if r["host_239"] is None or (r["site"],r["year"]) not in accepted:
            continue
        key=(r["site"],r["year"],r["plant"],r["origin"],r["host_239"])
        grouped[key]["n"]+=1
        grouped[key]["pos"]+=r["caterpillar"]
    rows=[
        {"site":key[0],"year":key[1],"plant":key[2],"origin":key[3],
         "source_239_butterfly_consumers":key[4],
         "n_field_surveys":v["n"],"positive_caterpillar_events":v["pos"]}
        for key,v in sorted(grouped.items())
    ]
    return rows


def assess_source_field_link_counts(grouped,*,require_source_gate=True):
    sites={r["site"] for r in grouped}
    taxa={r["plant"] for r in grouped}
    genera={r["plant"].split()[0] for r in grouped}
    nat={r["plant"] for r in grouped if r["origin"]=="native"}
    alien={r["plant"] for r in grouped if r["origin"]=="alien"}
    strata={(r["site"],r["year"]) for r in grouped}
    successes=sum(r["positive_caterpillar_events"] for r in grouped)
    total=sum(r["n_field_surveys"] for r in grouped)
    hist=Counter(r["source_239_butterfly_consumers"] for r in grouped)
    gate=(len(taxa)>=15 and len(alien)>=3 and len(nat)>=3 and
          len(genera)>=10 and len(sites)>=10 and len(strata)>=10 and
          successes>=20 and len(hist)>=3)
    return {
        "grouped_plant_site_year_records":len(grouped),
        "original_239_HOSTS_exact_linked_field_plants":len(taxa),
        "original_239_linked_field_native_plant_species":len(nat),
        "original_239_linked_field_alien_plant_species":len(alien),
        "source_plant_genera":len(genera),
        "source_field_sites":len(sites),
        "source_site_year_strata":len(strata),
        "source_plant_field_survey_events":total,
        "source_caterpillar_detection_events":successes,
        "original_239_source_consumer_count_histogram_by_plant_site_year":{
            str(k):v for k,v in sorted(hist.items())
        },
        "model_information_gate":(
            "SOURCE_LINKED_PLANT_AND_SITE_YEAR_SUPPORT" if gate
            else "INSUFFICIENT_MATCHED_BUTTERFLY_HOST_SOURCE_FIELD_SUPPORT"
        )
    }


def build_X(grouped,*,with_site_year=False):
    name_count={r["plant"]:r["source_239_butterfly_consumers"] for r in grouped}
    orig=np.array([math.log1p(n) for n in name_count.values()])
    mu=float(orig.mean());sd=float(orig.std(ddof=1))
    if sd<=1e-12:
        raise ValueError("Original 239 documented consumers lack variance")
    y=np.array([r["positive_caterpillar_events"] for r in grouped],float)
    n=np.array([r["n_field_surveys"] for r in grouped],float)
    nat=np.array([r["origin"]=="native" for r in grouped],float)
    h=np.array([(math.log1p(r["source_239_butterfly_consumers"])-mu)/sd
                for r in grouped])
    x=[np.ones(len(grouped)),nat,h]
    names=[]
    if with_site_year:
        names=sorted({(r["site"],r["year"]) for r in grouped})
        for label in names[1:]:
            x.append(np.array([(r["site"],r["year"])==label for r in grouped],float))
    return np.column_stack(x),y,n,names,(mu,sd)


def blocked_leave_site_out(grouped):
    X,y,n,_,normalized=build_X(grouped)
    names=sorted({r["site"] for r in grouped})
    result=Counter()
    for site in names:
        test=np.array([r["site"]==site for r in grouped])
        train=~test
        if np.linalg.matrix_rank(X[train])<3 or \
           np.linalg.matrix_rank(X[train,:2])<2:
            result["sites_skipped_rank"]+=1
            continue
        try:
            a=solve_logit(X[train],y[train],n[train])
            b=solve_logit(X[train,:2],y[train],n[train])
        except (ValueError,np.linalg.LinAlgError):
            result["sites_skipped_fit"]+=1
            continue
        full=expit(X[test]@a)
        base=expit(X[test,:2]@b)
        def score(p):
            return float(np.sum(y[test]*np.log(np.maximum(p,1e-12))+
                               (n[test]-y[test])*np.log(np.maximum(1-p,1e-12))))
        result["sites_validated"]+=1
        result["surveys_held_out"]+=int(n[test].sum())
        result["full_log_score"]+=score(full)
        result["origin_only_log_score"]+=score(base)
    return {
        **dict(result),
        "incremental_log_score_per_1000_out_of_site_surveys":(
            1000*(result["full_log_score"]-result["origin_only_log_score"])/
            result["surveys_held_out"] if result["surveys_held_out"] else None),
        "source_positive_gain_means_original239_portfolio_predicts_outside_training_sites":True,
    }


def fit_site_year_model(grouped):
    qc=assess_source_field_link_counts(grouped)
    if qc["model_information_gate"]!="SOURCE_LINKED_PLANT_AND_SITE_YEAR_SUPPORT":
        return {"source_gate":qc,"status":"NO_MODEL_INSUFFICIENT_MATCHED_SOURCE_SUPPORT"}
    # Group/site-year intercepts and fixed original plant origin/source count;
    # artificial ridge stabilisation is not a causal intervention.
    X,y,n,names,norm=build_X(grouped,with_site_year=True)
    out={
        "source_gate":qc,
        "original_source_site_year_strata":len(names),
        "fixed_effects_plant_site_year_group_rows":len(grouped),
        "source_full_design_rank":int(np.linalg.matrix_rank(X)),
        "source_full_design_n_predictors":X.shape[1],
        "condition_number":float(np.linalg.cond(X)),
    }
    if len(grouped)<X.shape[1]+5 or np.linalg.matrix_rank(X)<X.shape[1] or \
       out["condition_number"]>1e5:
        return {**out,"status":"NO_MODEL_RANK_OR_SMALL_STRATA"}
    try:
        betas=solve_logit(X,y,n)
        origin_only=solve_logit(np.delete(X,2,axis=1),y,n)
    except (ValueError,np.linalg.LinAlgError) as e:
        return {**out,"status":"NO_MODEL_NUMERICAL_IDENTIFIABILITY",
                "source_error_type":type(e).__name__}
    out.update({
        "status":"DESCRIPTIVE_SOURCE_SITE_YEAR_FIXED_EFFECT_MODEL_ONLY",
        "original_239_butterfly_host_consumer_log1p_SD":norm[1],
        "original_239_HOSTS_consumers_perSD_adjusted_log_odds_beta":float(betas[2]),
        "original_239_HOSTS_consumers_perSD_adjusted_odds_ratio":float(math.exp(betas[2])),
        "source_native_vs_alien_plants_adjusted_odds_ratio":float(math.exp(betas[1])),
        "source_origin_only_native_vs_alien_odds_ratio":float(math.exp(origin_only[1])),
        "leave_one_entire_survey_site_out":blocked_leave_site_out(grouped),
    })
    return out


def run(args):
    author=json.loads((args.author_source/"SOURCE_MANIFEST.json").read_text())
    if author.get("inputs_count")!=11:
        raise ValueError("Wrong original author-only source manifest")
    for path,expected in author["sha256"].items():
        if sha256(args.author_source/path)!=expected:
            raise ValueError("Original accepted HOSTS/WCVP source identity changed")
    rows=read_rows(args.author_source/"inputs/anthropogenic_species_metrics.csv")
    names={r["species"] for r in rows}
    if len(rows)!=239 or len(names)!=239:
        raise ValueError("Original 239 butterfly source roster changed")
    bname,ids,ambiguous,source_meta=accepted_global_links(
        args.author_source/"inputs/original/insect_host_accepted.csv",names)
    accepted={
        key:len(ids[next(iter(values))])
        for key,values in bname.items() if len(values)==1 and key not in ambiguous
    }
    plant_origin=read_origin(args.origin)
    inferred=read_inferred(args.inferred)
    surveys,qc,source,_=raw_filter(args.original_raw,inferred,plant_origin,accepted)
    source_paired=site_year_matched(surveys,min_survey=5,source_links=True)
    # Validate source original stage2 inventory source margins: prevents
    # losing same-site-year matched source positives through a bad join.
    count=aggregate_matched(source_paired)
    if count["matched_site_year_or_week_strata"]!=38 or \
       count["independent_field_sites_with_matching_both_plant_origins"]!=22 or \
       count["source_native_survey_events_in_matched_contexts"]!=1209 or \
       count["source_alien_survey_events_in_matched_contexts"]!=824:
        raise ValueError("Original source 239 site-year matched cohort changed")
    grouped=aggregate_groups(surveys,source_paired)
    model=fit_site_year_model(grouped)
    result={
        "schema":SCHEMA,
        "status":model["status"],
        "source":{
            "external_field_commit":EXTERNAL_COMMIT,
            "external_raw_csv_SHA256":sha256(args.original_raw),
            "original_author_WCVP_host_source_SHA256":SOURCE_SHA,
            "original_239_host_plant_links":source_meta["original_239_accepted_species_plant_pairs"],
        },
        "stage2_original_source_matched_site_year_population":count,
        "source_field_taxonomy_and_site_year_support":model["source_gate"],
        "stage3_model":model,
        "red_lines":[
            "The field caterpillar response is ALL Lepidoptera species, not the 239 original butterflies and not documented local invader feeding.",
            "Country-level WCVP native/introduced plant distributions cannot replace author's USDA site/plant native status.",
            "Original matched field site-year samples are observational and native/alien plant identities are not randomly assigned.",
            "Same site-year does not ensure same plant, canopy, morphology, phenology or survey day; host identity and branch availability are confounded.",
            "Original source HOSTS consumer count is a publication index, not a demonstrated local potential-host portfolio.",
            "Ridge logistic with site-year fixed effects and held-out site log score is exploratory, with survey IDs repeated by site and plant; do not compute naive organismal p values.",
            "No original field raw survey coordinates, site identifiers, exact plant species×site combinations, butterfly-host accepted IDs or private authored matrices exported."
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "status":result["status"],
        "siteyear_info":model["source_gate"],
        "siteyear_OR_per_SD":model.get("original_239_HOSTS_consumers_perSD_adjusted_odds_ratio"),
        "siteyear_origin_OR":model.get("source_native_vs_alien_plants_adjusted_odds_ratio"),
        "heldout_sites":model.get("leave_one_entire_survey_site_out"),
        "design_condition_number":model.get("condition_number")
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",required=True,type=Path)
    p.add_argument("--original-raw",required=True,type=Path)
    p.add_argument("--origin",required=True,type=Path)
    p.add_argument("--inferred",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    run(p.parse_args())
