#!/usr/bin/env python3
"""Actually test historical plant-HOSTS value in future-year UKBMS walk detection.

Freeze: 1987–95 training after BSBI prior1987 sources; test 2005–21, 19
same UKBMS transects, 16 botanical hectads. Compare static *documented* host
taxa vs climate+survey/plant effort and deterministic matched NONHOST controls.

This is an exploratory retrospective prediction, NOT human cultivation/host
colonisation/temperature causality. Raw original source IDs/site coordinates
and author-only host pairs are transient in-memory and NEVER serialized.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,random,statistics,time,zipfile
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from pyproj import Transformer

from audit_bsbi_historical_hectad_source_20261029 import obtain
from audit_bsbi_original_hosts_name_match_20261029 import (
    get_host_edges,canonical,AUTHOR_SHA)
from audit_bsbi_original_host_hectad_opportunity_20261029 import (
    to_osgb_hectad,recorded_flag,plant_lookup,PERIODS_PRIOR)
from audit_bsbi_nested_hectad_layout_20261029 import NESTED
from audit_uk_ltmn_historical_event_join_20261028 import (
    FIELDS_VISITS,FIELDS_COUNTS,get_source,excel_rows)
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk
from audit_uk_ltmn_butterfly_detection_panel_20261028 import (
    CROSSWALK,count_at_visits,original_taxonomy_ok)
import audit_uk_ltmn_NASA_POWER_historical_baseline_20261028 as nasa

SCHEMA="chocho_UK_BSBI_pre1987_host_late_2005to2021_survey_hindcast_v01"
TRAIN=(1987,1995)
TEST=(2005,2021)
SEED=20261029
PRIOR_HOST_SOURCE_SPATIAL_MIN=10
MIN_HOST_TAXA=2
MIN_HISTORIC_PLANT_GRID_VARIATION=10
N_PLACEBOS=19
KNOWN_HOST_GAP={"Aglais urticae","Lycaena dispar"}
SOURCE_19_SITE_SHA="051c82f9af21941bee3ea2fa9d8fa739990a36ed6e25e43441c25821e257b806"

def original_bsbi_features(original_source,sites,cross,author_artifact):
    """In-memory site×plant prior1987 presence, not declared absence."""
    host_edges,host_source=get_host_edges(author_artifact)
    name_id,raw_nested=plant_lookup(original_source)
    hectads={s:to_osgb_hectad(*cross[s][:2]) for s in sites}
    squares=set(hectads.values())
    if len(squares)<MIN_HISTORIC_PLANT_GRID_VARIATION:
        raise ValueError("Too few independent original pre1987 botanical hectads")
    documented={h:set() for h in squares}
    from audit_bsbi_original_host_hectad_opportunity_20261029 import PERIODS_LATER
    with zipfile.ZipFile(io.BytesIO(raw_nested)) as z:
        source_files=0
        for info in z.infolist():
            if not info.filename.endswith(".csv"):continue
            source_files+=1
            ident=Path(info.filename).stem
            with z.open(info) as stream:
                reader=csv.DictReader(io.TextIOWrapper(stream,encoding="utf-8-sig",newline=""))
                if not {"hectad",*PERIODS_PRIOR,*PERIODS_LATER}<=set(reader.fieldnames or []):
                    raise ValueError("Original BSBI 1970–86 source date class fields changed")
                for row in reader:
                    square=str(row["hectad"]).strip().upper()
                    if square not in squares:continue
                    if any(recorded_flag(row[x]) for x in PERIODS_PRIOR):
                        documented[square].add(ident)
    if source_files!=3497:
        raise ValueError("Source full 3497 atlas data files changed")
    host_ids={}
    for sp in CROSSWALK:
        ids=set()
        for name in host_edges.get(sp,[]):
            matched=name_id.get(canonical(name),set())
            if len(matched)==1:ids.update(matched)
        host_ids[sp]=ids
    early_counts={(s,sp):len(documented[hectads[s]] & host_ids[sp])
                  for s in sites for sp in CROSSWALK}
    eligible={
       sp for sp in CROSSWALK
       if sp not in KNOWN_HOST_GAP and len(host_ids[sp])>=MIN_HOST_TAXA
       and max(early_counts[s,sp] for s in sites)>0
       and len({early_counts[s,sp] for s in sites})>1
    }
    if len(eligible)<10:
        raise ValueError("Source-only documented historical host support below 10 eligible butterflies")
    all_prior={s:len(documented[hectads[s]]) for s in sites}
    pool=set().union(*documented.values()) | set().union(*host_ids.values())
    # Include names with no historic records as source absence-of-record counterfactual
    # but NEVER interpret source zero as actual botanical absence.
    universe={str(row["ddbidLookup"]).strip() for row in
              []}  # resolved below using source's entire 3497 distribution file list
    with zipfile.ZipFile(io.BytesIO(raw_nested)) as z:
        universe={Path(info.filename).stem for info in z.infolist() if info.filename.endswith(".csv")}
    prevalence={taxon:sum(taxon in documented[h] for h in squares) for taxon in universe}
    metadata={"n_original_sites":len(sites),"n_original_10km_plant_hectads":len(squares),
              "pre1987_available_plant_taxa_per_site_min":min(all_prior.values()),
              "pre1987_available_plant_taxa_per_site_max":max(all_prior.values()),
              "pre1987_available_plant_taxa_per_site_median":statistics.median(all_prior.values()),
              "all_source_UK_butterflies":len(CROSSWALK),
              "source_only_eligible_host_var_species":len(eligible),
              "excluded_documented_major_UK_host_caveat":sorted(KNOWN_HOST_GAP),
              "historical_date_classes_used":list(PERIODS_PRIOR),
              "later_botanical_date_classes_used_for_fit":False,
              "source_host_names_original":sum(len(v) for v in host_edges.values()),
              "host_source_SHA256":AUTHOR_SHA,
              "BSBI_archive_SHA256":hashlib.sha256(original_source).hexdigest(),
              "no_original_UKBMS_site_or_BSBI_hectad_names_exported":True}
    return hectads,documented,host_ids,early_counts,all_prior,eligible,prevalence,metadata

def placebo_ids(host_taxa,prevalence,draw,sp):
    """Equal host taxon count, exact pre1987 10km presence prevalence strata."""
    import random
    seed=int(hashlib.sha256(f"{SEED}:{sp}:{draw}".encode()).hexdigest()[:15],16)
    rng=random.Random(seed)
    available=defaultdict(list)
    for plant,p in prevalence.items():
        if plant not in host_taxa:available[p].append(plant)
    for key in available:available[key].sort()
    chosen=set()
    for host in sorted(host_taxa):
        target=prevalence[host]
        candidates=[p for p in available[target] if p not in chosen]
        if not candidates:
            raise ValueError("Not enough prevalence-matched original nonhost botanical controls")
        chosen.add(candidates[rng.randrange(len(candidates))])
    if len(chosen)!=len(host_taxa):
        raise ValueError("Placebo selected duplicate plant names")
    return chosen

def raw_design_with_era(rows,sites,selected,doc,hectads,host_ids,all_prior,placebo=None):
    qualified=[r for r in rows if r[1] in selected]
    out={"calendar":[],"weather":[],"host":[],"host_temperature":[]}
    for row in qualified:
        site,sp,era,status,base,weather=row
        old=base[:]
        # Existing source baseline site/year seasonal variables and historical
        # botanical richness are common to every comparator.
        common=old+[math.log1p(all_prior[site])]
        active=host_ids[sp] if placebo is None else placebo[sp]
        n=len(active & doc[hectads[site]])
        hostmath=math.log1p(n)
        out["calendar"].append((site,sp,era,status,common,weather))
        out["weather"].append((site,sp,era,status,common,weather))
        out["host"].append((site,sp,era,status,common+[hostmath],weather))
        out["host_temperature"].append((site,sp,era,status,common+[hostmath,hostmath*weather[0]],weather))
    return out

def run_models(raw_panel,source_sites,hectads,with_placebo=True):
    source=raw_panel["host"]
    tr=[i for i,r in enumerate(source) if r[2]=="early"]
    te=[i for i,r in enumerate(source) if r[2]=="late"]
    if len(tr)<3000 or len(te)<3000 or len(source_sites)<10:
        raise ValueError("Insufficient original UKBMS prior1987 host training/source response")
    outcomes={}
    for label,weather in (("calendar",False),("weather",True),("host",True),("host_temperature",True)):
        outcomes[label]=nasa.train_eval(raw_panel[label],tr,te,weather)
    def delta(a,b):return outcomes[a]["logloss"]-outcomes[b]["logloss"]
    outcomes["host_increment_vs_weather_logloss"]=delta("weather","host")
    outcomes["host_increment_vs_calendar_logloss"]=delta("calendar","host")
    outcomes["host_temperature_increment_vs_host_logloss"]=delta("host","host_temperature")
    # Independence is at the plant map 10km hectad, not 19 possibly shared transects.
    groups=[[] for _ in range(5)]
    uniq=sorted({hectads[s] for s in source_sites},key=lambda h:hashlib.sha256((str(SEED)+h).encode()).hexdigest())
    for i,h in enumerate(uniq):groups[i%5].append(h)
    if any(len(g)<2 for g in groups):
        raise ValueError("Historical botanical 10km held-out fold lacks independent hectads")
    folds=[]
    for rank,g in enumerate(groups,1):
        siteout={s for s in source_sites if hectads[s] in g}
        train_idx=[i for i in tr if source[i][0] not in siteout]
        test_idx=[i for i in te if source[i][0] in siteout]
        if len(train_idx)<1000 or len(test_idx)<300 or len({source[i][3] for i in test_idx})<2:
            raise ValueError("Historical source held-out hectad fold insufficient")
        res={name:nasa.train_eval(raw_panel[name],train_idx,test_idx,weather) for name,weather in
             (("calendar",False),("weather",True),("host",True),("host_temperature",True))}
        folds.append({"fold":rank,"heldout_independent_plant_hectads":len(g),
             "heldout_transects":len(siteout),
             "train_pre1987_source_plant_cutoff_visits_x_species":len(train_idx),
             "test_2005_2021_survey_visit_x_species":len(test_idx),
             "logloss":{k:res[k]["logloss"] for k in res},
             "host_vs_weather_logloss":res["weather"]["logloss"]-res["host"]["logloss"],
             "host_vs_calendar_logloss":res["calendar"]["logloss"]-res["host"]["logloss"],
             "host_Tinteraction_vs_host_logloss":res["host"]["logloss"]-res["host_temperature"]["logloss"]})
    outcomes["geographically_unseen_hectad_fold_scores"]=folds
    outcomes["geographically_unseen_hectad_folds_host_better_than_both_M0"]=sum(
        x["host_vs_weather_logloss"]>0 and x["host_vs_calendar_logloss"]>0 for x in folds)
    return outcomes

def run(args):
    taxon=original_taxonomy_ok(args.original_manifest)
    cross,crossqc=read_crosswalk()
    nasa.PERIODS={"early":TRAIN,"late":TEST}
    original,sites,events=nasa.nasa_site_effort(
        excel_rows(get_source("visits"),FIELDS_VISITS),cross)
    if len(sites)!=19:raise ValueError("New 1987–95 source visit-only qualified site panel differs from 19")
    counted,counts_qc=count_at_visits(excel_rows(get_source("counts"),FIELDS_COUNTS),original,sites)
    source,_,_=obtain()
    hectads,docs,hosts,hostcounts,botanical,eligible,prevalence,hostqc=original_bsbi_features(
         source,sites,cross,args.author_artifact_zip)
    # Original SOURCE water / this outcome untouched till source mapping fixed.
    coords=Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True)
    climate={};weather_shas=[]
    for i,site in enumerate(sorted(sites),1):
        east,north,habitat=cross[site];lon,lat=coords.transform(east,north)
        if not (49<lat<60 and -9<lon<3):raise ValueError("Original geocoded UKBMS route off Britain")
        weather,sha,invalid=nasa.original_nasa_power_site_weather(lat,lon)
        climate[site]=weather;weather_shas.append(sha)
        if i%5==0:print(json.dumps({"original_NASA_POWER_19_ERA_sites_retrieved":i}),flush=True)
        time.sleep(.25)
    aggregate_sha=hashlib.sha256("\n".join(sorted(weather_shas)).encode()).hexdigest()
    if aggregate_sha!=SOURCE_19_SITE_SHA:
        raise ValueError("Original 19-site 1981–2021 daily NASA POWER source changed")
    raw,qc=nasa.site_weather_covariate_panel(events,sites,counted,cross,climate)
    if set(raw[i][2] for i in range(len(raw)))!={"early","late"}:raise ValueError("Original two-era visit outcomes missing")
    selected=set(eligible)
    derived=raw_design_with_era(raw,sites,selected,docs,hectads,hosts,botanical)
    models=run_models(derived,sites,hectads)
    tr=sum(r[2]=="early" for r in derived["host"])
    te=sum(r[2]=="late" for r in derived["host"])
    # Original 19 prevalence-matched nonhost trial, temporal test only.
    placebo_gains=[]
    for draw in range(N_PLACEBOS):
        control={sp:placebo_ids(hosts[sp],prevalence,draw,sp) for sp in selected}
        pseudo=raw_design_with_era(raw,sites,selected,docs,hectads,hosts,botanical,control)
        tr_idx=[i for i,r in enumerate(pseudo["host"]) if r[2]=="early"]
        te_idx=[i for i,r in enumerate(pseudo["host"]) if r[2]=="late"]
        score=nasa.train_eval(pseudo["host"],tr_idx,te_idx,True)
        placebo_gains.append(models["weather"]["logloss"]-score["logloss"])
    n_better=sum(g>=models["host_increment_vs_weather_logloss"] for g in placebo_gains)
    models["source_pre1987_prevalence_matched_nonhost_placebo"]={
        "n_source_outcome_blind_plant_placebos":len(placebo_gains),
        "n_placebo_logloss_improvements_ge_host":n_better,
        "median_placebo_logloss_improvement":statistics.median(placebo_gains),
        "min_placebo_logloss_improvement":min(placebo_gains),
        "max_placebo_logloss_improvement":max(placebo_gains),
        "source_original_plant_presence_prevalence_matched":True,
        "not_a_causal_permutation_p_value":True}
    good=(models["host_increment_vs_weather_logloss"]>0 and
          models["host_increment_vs_calendar_logloss"]>0 and
          models["geographically_unseen_hectad_folds_host_better_than_both_M0"]>=4 and
          n_better<=1)
    result={
        "schema":SCHEMA,"status":"ORIGINAL_UK_1987TO1995_TRAIN_2005TO2021_TEST_HOST_GEOGRAPHY_PREDICTION_COMPLETE",
        "source":{
            "original_butterfly_manifest_sha256":hashlib.sha256(args.original_manifest.read_bytes()).hexdigest(),
            "original_source_27_butterflies":len(CROSSWALK),
            "historical_host_support_eligible_butterflies":len(selected),
            "source_only_excluded_primary_species":sorted(set(CROSSWALK)-selected),
            "original_repeat_monitoring_transects":len(sites),
            "source_botanical_hectads":len(set(hectads.values())),
            "original_training_1987to1995_qualified_visit_species_rows":tr,
            "original_testing_2005to2021_qualified_visit_species_rows":te,
            "original_NASA_POWER_source_SHA256_19_sites":aggregate_sha,
            "original_source_weather_19_site_daily_days_each":nasa.DAYS_EXPECTED,
            "source_original_botanical":hostqc,
            "taxon_key_qc":taxon,
            "butterfly_visit_source_qc":qc,
        },
        "comparisons":models,
        "original_atlas_observed_host_names_complete":False,
        "original_accepted_larval_host_use_documented_before_1987":False,
        "human_maintained_host_plants_identified":False,
        "future_butterfly_range_or_population_predicted":False,
        "host_plant_nonrecord_is_not_true_absence":True,
        "source_exact_survey_date_and_location_or_original_plant_records_exported":False,
        "information_gate":"PASS_HISTORICAL_HOST_PREDICTIVE_INCREMENT_NOT_CAUSAL"
                  if good else "FAIL_HISTORICAL_HOST_PREDICTIVE_INCREMENT_ROBUSTNESS"}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":result["status"],"information_gate":result["information_gate"],
        "host_eligible_butterfly_species":len(selected),"tr":tr,"te":te,
        "host_vs_weather":models["host_increment_vs_weather_logloss"],
        "host_vs_calendar":models["host_increment_vs_calendar_logloss"],
        "heldout_geography_folds_positive":models["geographically_unseen_hectad_folds_host_better_than_both_M0"],
        "placebos_ge_host":n_better},sort_keys=True))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-manifest",required=True,type=Path)
    p.add_argument("--author-artifact-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    run(p.parse_args())
