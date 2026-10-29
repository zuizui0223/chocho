#!/usr/bin/env python3
"""Out-of-time UK butterfly detection: actual 1981–95 weather, 1970 wild hosts.

All original author-only HOSTS source, 3497 BSBI species×10km period records,
original UKBMS site-level visit/count data and NASA POWER daily weather are
reconstructed only in CI process memory. No site/grid/host data exported.
Source period flags 0 = not recorded, never source-verified plant absence.
"""
from __future__ import annotations

import argparse,collections,hashlib,json,math,random,time
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from pyproj import Transformer

from audit_bsbi_UK_host_hectad_coverage_20261029 import (
    site_historical_host_coverage,SCHEMA as HOST_SCHEMA)
from audit_uk_ltmn_historical_event_join_20261028 import (
    FIELDS_COUNTS,FIELDS_VISITS,get_source,excel_rows)
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk
from audit_uk_ltmn_butterfly_detection_panel_20261028 import (
    CROSSWALK,count_at_visits)
from audit_uk_ltmn_NASA_POWER_historical_baseline_20261028 import (
    nasa_site_effort,original_nasa_power_site_weather,site_weather_covariate_panel,
    score_site_species,train_eval as original_train_eval,SEED as OLD_SEED)

SCHEMA="chocho_uk_ltmn_NASA_POWER_bsbi_pre1970_host_source_hindcast_20261029_v01"
C=.1
SEED=20261029
N_PLACEBO=9
MIN_TEST_POSITIVE_SITE_FOLDS=4

def exposure(inputs,site,species,source_key="own"):
    historical=inputs["pre1970_observed_bsbi_plant_names_by_site"][site]
    if source_key=="own":
        hosts=inputs["exact_source_host_plant_names_by_butterfly"][species]
    else:
        hosts=source_key[species]
    return len(historical & hosts)

def unrelated_plant_portfolios(inputs,draw):
    """Loose, not prevalence- or family-matched, source-only alternative hosts."""
    pool=sorted(inputs["all_source_plant_name_pool"])
    out={}
    for species,accepted in sorted(inputs["exact_source_host_plant_names_by_butterfly"].items()):
        candidates=[p for p in pool if p not in accepted]
        if len(accepted)>len(candidates):
            raise ValueError("Insufficient distinct known nonhost source pool for breadth-compatible placebo")
        seed=int(hashlib.sha256(f"{SEED}|{draw}|{species}".encode()).hexdigest()[:16],16)
        out[species]=frozenset(random.Random(seed).sample(candidates,len(accepted)))
    return out

def covariates(raw,ids,inputs,host_source):
    species_early=inputs["exact_source_host_plant_names_by_butterfly"]
    X=[]
    for i in ids:
        site,sp,era,flag,features,weather=raw[i]
        if sp not in species_early:raise ValueError("Original focal UKBMS taxon changed")
        effort=inputs["original_pre1970_site_botanical_effort"][site]
        if effort<=0:raise ValueError("Original BSBI plant effort source has zero mapped taxa")
        row=list(features)+list(weather)+[math.log1p(effort)]
        if host_source is not None:
            row.append(math.log1p(exposure(inputs,site,sp,source_key=host_source)))
        X.append(row)
    out=np.asarray(X,float)
    if not np.isfinite(out).all():raise ValueError("Missing original climate or host predictor")
    return out

def site_species_weights(raw,indices):
    c=collections.Counter((raw[i][0],raw[i][1]) for i in indices)
    return np.asarray([len(indices)/(len(c)*c[(raw[i][0],raw[i][1])]) for i in indices],float)

def train_predict_score(raw,tr,te,inputs,host_source):
    xx=covariates(raw,tr,inputs,host_source)
    yy=np.asarray([raw[i][3] for i in tr],int)
    if set(yy)!={0,1}:raise ValueError("Original early training lacks detected and missed species visits")
    scaler=StandardScaler()
    model=LogisticRegression(C=C,solver="lbfgs",max_iter=1000,random_state=OLD_SEED)
    model.fit(scaler.fit_transform(xx),yy,sample_weight=site_species_weights(raw,tr))
    prediction=model.predict_proba(scaler.transform(covariates(raw,te,inputs,host_source)))[:,1]
    out=score_site_species(raw,te,prediction)
    return out,prediction

def butterfly_equal_scores(raw,ids,probs0,probs1):
    """All source taxa reported, without selecting impressive outcomes."""
    bucket=collections.defaultdict(lambda:collections.defaultdict(list))
    for idx,p0,p1 in zip(ids,probs0,probs1):
        site,sp,era,y,*unused=raw[idx]
        for key,prob in (("M0",p0),("M1",p1)):
            p=float(np.clip(prob,1e-7,1-1e-7))
            err=-(y*math.log(p)+(1-y)*math.log1p(-p))
            bucket[sp][key].append((site,err))
    result={}
    for sp in CROSSWALK:
        stats={}
        for key in ("M0","M1"):
            groups=collections.defaultdict(list)
            for site,loss in bucket[sp][key]:groups[site].append(loss)
            stats[key]=sum(sum(v)/len(v) for v in groups.values())/len(groups) if groups else None
        result[sp]={"source_late_site_equal_logloss_M0":stats["M0"],
                    "source_late_site_equal_logloss_M1":stats["M1"],
                    "M0_minus_M1_logloss_improvement":stats["M0"]-stats["M1"] if None not in stats.values() else None}
    return result

def test_contest(raw,sites,inputs,original):
    sites=set(sites)
    tr=[i for i,r in enumerate(raw) if r[2]=="early"]
    te=[i for i,r in enumerate(raw) if r[2]=="late"]
    if len(tr)<10000 or len(te)<10000:raise ValueError("Insufficient original source visit×species history")
    if len(raw)!=original["source"]["source_paired_visit_species_labels"]:
        raise ValueError("Original strict 27-species UKBMS source photo/visit labels drifted")
    if len(sites)!=original["source"]["n_NASA_POWER_source_transects"]:
        raise ValueError("Original NASA weather and source BSBI site cohort differ")
    # Compare old original model within real-data source SHA, not just code.
    original_ref=original_train_eval(raw,tr,te,True)
    old=original["validation"]["NASA_POWER_past30_temperature_precip"]
    if not math.isclose(original_ref["logloss"],old["logloss"],abs_tol=1e-8):
        raise ValueError("Frozen weather-only historical UKBMS result cannot be reproduced")
    base,p0=train_predict_score(raw,tr,te,inputs,None)
    host,p1=train_predict_score(raw,tr,te,inputs,"own")
    delta=base["logloss"]-host["logloss"]
    # Fixed 5 entire station fold groups (training early, testing late AND unseen).
    site_order=sorted(sites,key=lambda s:hashlib.sha256(f"{SEED}|{s}".encode()).hexdigest())
    partitions=[site_order[i::5] for i in range(5)]
    folds=[]
    for i,part in enumerate(partitions,1):
        if len(part)<3:raise ValueError("Few independent heldout UK transects")
        train=[j for j in tr if raw[j][0] not in part]
        test=[j for j in te if raw[j][0] in part]
        b,_=train_predict_score(raw,train,test,inputs,None)
        h,_=train_predict_score(raw,train,test,inputs,"own")
        folds.append({"fold":i,"source_heldout_transects":len(part),
                      "source_late_visit_species_labels":len(test),
                      "M0_climate_botanical_effort_site_equal_logloss":b["logloss"],
                      "M1_plus_pre1970_wild_hosts_site_equal_logloss":h["logloss"],
                      "M0_minus_M1_logloss_improvement":b["logloss"]-h["logloss"]})
    # Original host portfolio species specificity controls, all random draws
    # fixed from source taxon names independently of original butterfly status.
    placebo=[]
    for repeat in range(1,N_PLACEBO+1):
        fake=unrelated_plant_portfolios(inputs,repeat)
        sc,_=train_predict_score(raw,tr,te,inputs,fake)
        placebo.append({"source_only_predeclared_draw":repeat,
                        "M0_minus_loose_nonhost_control_logloss_improvement":base["logloss"]-sc["logloss"]})
    better_folds=sum(x["M0_minus_M1_logloss_improvement"]>0 for x in folds)
    wins_vs_all_placebo=all(delta>p["M0_minus_loose_nonhost_control_logloss_improvement"] for p in placebo)
    supported=(delta>0 and better_folds>=MIN_TEST_POSITIVE_SITE_FOLDS and wins_vs_all_placebo)
    return {
       "historical_ERA_NASA_POWER_original_reproduction":original_ref,
       "M0_climate_and_pre1970_botanical_effort":base,
       "M1_M0_plus_documented_source_native_or_introduced_wild_hosts":host,
       "overall_M0_minus_M1_site_species_equal_logloss":delta,
       "five_original_heldout_site_late_era_source_scores":folds,
       "five_folds_with_early_wild_host_predictive_gain":better_folds,
       "nine_frozen_nonhost_plant_portfolio_negative_controls":placebo,
       "true_source_host_improvement_larger_than_all_nine_loose_controls":wins_vs_all_placebo,
       "by_all_original_27_source_butterfly_taxa":butterfly_equal_scores(raw,te,p0,p1),
       "host_model_information_gate":
           "PASS_HISTORIC_WILD_HOST_OUTOFTIME_PREDICTIVE_INCREMENT_EXPLORATORY" if supported
           else "WILD_HOST_PREDICTIVE_INCREMENT_NOT_ESTABLISHED",
    }

def run(args):
    source_report,inputs=site_historical_host_coverage(
        args.original_source,args.bsbi_zip,return_internal=True)
    if source_report["pre_model_information_gate"]!="PASS_WILD_HOST_HISTORICAL_SOURCE_CONTRAST_CANDIDATE":
        raise ValueError("Original source historical UK host contrast insufficient for inferential model")
    original=json.loads(args.original_nasa_result.read_text())
    if original.get("schema")!="chocho_UK_LTMN_NASA_POWER_1981_2021_outoftime_detection_v01" or (
         original["validation"]["extrapolation_weather_prediction_gate"]!="WEATHER_SPATIOTEMPORAL_HINDCAST_NOT_RELIABLY_BETTER"):
        raise ValueError("Missing original frozen negative weather baseline; cannot alter original control")
    source_cross,_=read_crosswalk()
    visits,sites,by_site=nasa_site_effort(
        excel_rows(get_source("visits"),FIELDS_VISITS),source_cross)
    if set(sites)!=set(inputs["site_to_hectad"]):
        raise ValueError("Original UKBMS effort-only sites changed when loading BSBI")
    counted,qc=count_at_visits(excel_rows(get_source("counts"),FIELDS_COUNTS),visits,sites)
    src_weather={};replies=[];source_invalid=0
    transform=Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True)
    for i,site in enumerate(sorted(sites),1):
        lon,lat=transform.transform(*source_cross[site][:2])
        data,digest,invalid=original_nasa_power_site_weather(lat,lon)
        src_weather[site]=data;replies.append(digest);source_invalid+=invalid
        if i%5==0:print(json.dumps({"original_NASA_POWER_UK_source_sites_loaded":i}),flush=True)
        time.sleep(.35)
    previous=original["source"]["original_NASA_POWER_source_response_sha256_sorted_unlinked_to_locations"]
    original_sha_overlap=sum((collections.Counter(replies)&collections.Counter(previous)).values())
    if len(replies)!=19 or len(previous)!=19:
        raise ValueError("NASA original 19-site weather source cohort size changed")
    # NASA POWER original response body SHA may drift without a documented cause.
    # DO NOT automatically assert this is only metadata. Instead the next step
    # MUST reproduce original site/species/time heldout weather-model loss 1e-8
    # before *any* new host predictive performance result may be accepted.
    rows,source_qc=site_weather_covariate_panel(by_site,sites,counted,source_cross,src_weather)
    val=test_contest(rows,sites,inputs,original)
    result={
       "schema":SCHEMA,
       "status":"UK_SOURCE_VERIFIED_HISTORICAL_WILD_HOST_VS_ACTUAL_WEATHER_OUTOFTIME_CONTEST",
       "original_source":{
          "BSBI_3497_taxon_original_zip_MD5":source_report["source"]["bsbi_original_MD5"],
          "author_source_WCVP_HOSTS_SHA256":source_report["source"]["authorized_original_WCVP_host_links_sha256"],
          "original_NASA_historical_source_result_sha256":hashlib.sha256(args.original_nasa_result.read_bytes()).hexdigest(),
          "original_UKBMS_visits_SHA256":original["source"]["original_UKBMS_visits_sha256"],
          "original_UKBMS_source_counts_SHA256":original["source"]["original_UKBMS_counts_sha256"],
          "original_BNG_LTMN_site_crosswalk_SHA256":original["source"]["original_UKBMS_site_metadata_sha256"],
          "source_butterfly_species":len(CROSSWALK),
          "source_ukbms_repeated_transects":len(sites),
          "BSBI_source_10km_hectads":source_report["source_10km_hectad_cells"],
          "BSBI_source_exact_accepted_plant_host_taxa":source_report["focal_exact_accepted_BSBI_binomial_host_plant_taxa"],
          "source_labeled_visit_species_records":len(rows),
          "source_late_test_label_counts":sum(r[2]=="late" for r in rows),
          "source_historical_ERA_training_label_counts":sum(r[2]=="early" for r in rows),
          "NASA_original_daily_source_weather_1981to2021_all19_sites_hash_verified":original_sha_overlap==19,
          "original_NASA_prior_raw_body_SHA_match_count":original_sha_overlap,
          "original_NASA_prior_raw_body_SHA_count":len(previous),
          "original_NASA_climate_only_heldout_model_numerically_reproduced_to_1e_minus8":True,
          "original_weather_raw_body_drift_not_assumed_metadata_only":original_sha_overlap<19,
          "raw_source_drift_posthoc_amendment":"NASA_POWER_SOURCE_BODY_SHA_DRIFT_DIAGNOSTIC_20261029.md",
       },
       "model":{
          "primary_plant_period":"published BSBI wild plants recorded at 10km grid to 1969, no future class",
          "actual_plant_is_maintained_by_people":False,
          "BSBI_0_is_true_plant_absence":False,
          "host_network_source":"original approved WCVP accepted 27 butterfly×plant host pairs",
          "botanical_survey_effort_proxy_in_both_M0_M1":True,
          "frozen_baseline":"original 1981-95 UKBMS M0 climate T/P and source effort",
          "test_window":"2005-2021 original UKBMS visited survey detections",
          "no_original_UKBMS_site_or_BSBI_hectad_names_or_persons_exported":True,
          "no_human_maintained_larval_hosts_or_climate_mediation_claim":True,
       },
       "validation":val,
       "ecological_interpretation":
          "Potential out-of-time predictive value of previously documented WILD larval host distribution only; not a human maintained-host or warming-tracking causal estimate",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
       "status":result["status"],
       "sites":len(sites),"source_butterflies":len(CROSSWALK),
       "M0_logloss":val["M0_climate_and_pre1970_botanical_effort"]["logloss"],
       "M1_logloss":val["M1_M0_plus_documented_source_native_or_introduced_wild_hosts"]["logloss"],
       "M0_minus_M1":val["overall_M0_minus_M1_site_species_equal_logloss"],
       "site_fold_positive":val["five_folds_with_early_wild_host_predictive_gain"],
       "M1_better_than_all_placebo":val["true_source_host_improvement_larger_than_all_nine_loose_controls"],
       "original_NASA_source_response_body_hash_match_count":original_sha_overlap,
       "gate":val["host_model_information_gate"]},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-source",type=Path,required=True)
    p.add_argument("--bsbi-zip",type=Path,required=True)
    p.add_argument("--original-nasa-result",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
