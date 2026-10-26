#!/usr/bin/env python3
"""Adult butterfly annual cofluctuation after coarse 14-bin Daymet forcing.

Weather is one historical Daymet modeled grid pixel at each of 14 PUBLIC
rounded source coordinate classes, NOT measured at actual Illinois transects
and NOT measured host plant leaf flush. 640 original fixed adult species pairs
are held constant under all 199 bipartite host-network nulls.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import frozen_focal,load_source,read_events
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import (
    field_panel,pair_correlations,meancontrast,null_test,
)
from analyze_sbus_illinois_host_degree_null_20261025 import (
    accepted_edges,butterfly_families,
)
from audit_sbus_daymet_14region_weather_source_20261026 import (
    original_coarse_coords,fetch_daymet_coarse_point
)

ROOT=Path(__file__).resolve().parents[1]
GATE=ROOT/"provenance/reviewer_defenses/results/sbus_daymet_14region_daily_source_gate_20261026.json"
QC=ROOT/"provenance/reviewer_defenses/results/sbus_illinois_10h_source_cofluctuation_20261026.json"
SOURCE_SAMPLE_CORRECTED=.14332588463806942


def fixed_source_panel():
    focal=frozen_focal()
    with tempfile.TemporaryDirectory() as tmp:
        source=Path(tmp)/"sbus_original.zip"
        load_source(source)
        events,_=read_events(source,focal)
        too_long=[key for key,e in events.items()
                  if key[0]=="Illinois Butterfly Monitoring Network"
                  and not e["bad_duration"] and e["duration"]>600]
        if len(too_long)!=5:
            raise ValueError("Original Henry et al source does not match verified five over-10h events")
        for key in too_long:
            events.pop(key)
        good,keys,sp,presence,field=field_panel(events)
        sitecoord,classes,info=original_coarse_coords(source,keys)
    if field["qualified_sites"]!=121 or field["qualified_site_years"]!=851:
        raise ValueError("Published effort filter changed 121-by-851 source cohort")
    pairdat=pair_correlations(good,keys,sp,presence)
    if pairdat["status"]!="FIXED_SPECIES_PAIRS_READY":
        raise ValueError("Original fixed temporal butterfly panel unavailable")
    temporal=pairdat["within_site_temporal_candidate_pairs"]
    if len(temporal)!=640:
        raise ValueError("Unexpected 640 source adult temporal pair count")
    return good,keys,field,temporal,sitecoord,classes


def verify_weather_replay(classes,receipt):
    if (receipt["status"]!="COARSE_REGION_YEAR_DAYMET_SOURCE_COMPLETE_NONLOCAL_PROXY"
        or len(classes)!=14 or receipt["public_rounded_coordinate_classes"]!=14
        or len(receipt["source_14_region_response_hash_receipts"])!=14):
        raise ValueError("Weather source gate failed")
    values={}
    for i,rounded in enumerate(classes):
        x,meta=fetch_daymet_coarse_point(rounded)
        record=receipt["source_14_region_response_hash_receipts"][i]
        if (record["nonidentifying_region_index"]!=i
            or record["response_sha256"]!=meta["sha256"]
            or record["weather_model_days"]!=meta["rows"]):
            raise ValueError("Weather source hash changed since independently verified acquisition")
        values[rounded]=x
    return values


def annual_partial_pairs(good,keys,temporal,sitecoord,weather):
    bysite=defaultdict(list)
    for st,yr in keys:
        bysite[st].append(yr)
    # cache site-specific observed source annual log adult encounter rate;
    # same records with original explicit-zero/implied non-detection convention.
    taxon=sorted({sp for a,b,r in temporal for sp in (a,b)})
    posat=defaultdict(set)
    for (site,yr),value in good.items():
        for sp,n in value[1].items():
            if n>0 and sp in taxon:
                posat[site].add(sp)
    weather_site={}
    for st,years in bysite.items():
        years=sorted(years)
        features=np.array([
            [weather[sitecoord[st]][y]["spring_GDD5"],
             weather[sitecoord[st]][y]["summer_precip_mm"]]
            for y in years],dtype=float)
        if len(years)<3 or not np.isfinite(features).all():
            raise ValueError("Incomplete gridded regional-season value for an eligible Illinois site-year")
        features-=features.mean(axis=0)
        weather_site[st]=(years,features)
    result=[]
    baseline_raw=[]
    weather_r2=[]
    lowrank=[]
    partial_differences=[]
    for a,b,ref in temporal:
        Xa=[];Ya=[];Weather=[]
        for st,ys in bysite.items():
            if a not in posat[st] or b not in posat[st]:
                continue
            years,W=weather_site[st]
            if sorted(ys)!=years:
                raise ValueError("Source site-year date mismatch")
            A=np.asarray([math.log1p(good[(st,y)][1].get(a,0)/good[(st,y)][0])
                          for y in years],dtype=float)
            B=np.asarray([math.log1p(good[(st,y)][1].get(b,0)/good[(st,y)][0])
                          for y in years],dtype=float)
            A-=np.mean(A)
            B-=np.mean(B)
            Xa.append(A)
            Ya.append(B)
            Weather.append(W)
        if len(Xa)<5 or sum(len(x) for x in Xa)<30:
            raise ValueError("Weather model changed original pair/site sample")
        aa=np.concatenate(Xa)
        bb=np.concatenate(Ya)
        ww=np.concatenate(Weather,axis=0)
        base=float((aa@bb)/math.sqrt((aa@aa)*(bb@bb)))
        if abs(base-ref)>1e-9:
            raise ValueError("Original 600-min source cofluctuation not reconstructed")
        beta_a,_,rank,_=np.linalg.lstsq(ww,aa,rcond=None)
        beta_b,_,rank2,_=np.linalg.lstsq(ww,bb,rcond=None)
        if rank!=2 or rank2!=2:
            lowrank.append((a,b))
            continue
        ra=aa-ww@beta_a
        rb=bb-ww@beta_b
        if (ra@ra)<=1e-10 or (rb@rb)<=1e-10:
            lowrank.append((a,b))
            continue
        pr=float((ra@rb)/math.sqrt((ra@ra)*(rb@rb)))
        if not math.isfinite(pr) or abs(pr)>1+1e-9:
            raise ValueError("Weather-adjusted source pair correlation impossible")
        pr=float(np.clip(pr,-1,1))
        result.append((a,b,pr))
        baseline_raw.append((a,b,base))
        partial_differences.append(pr-base)
        weather_r2 += [
            max(0.0, 1-float(ra@ra)/float(aa@aa)),
            max(0.0, 1-float(rb@rb)/float(bb@bb)),
        ]
    if lowrank or len(result)!=640:
        raise ValueError(f"Weather proxy correction cannot retain original all 640 temporal pairs: {len(lowrank)}")
    return result,baseline_raw,{
        "weather_coarse_source_position_classes":len(weather),
        "fixed_temporal_species_pairs":len(result),
        "exact_same_sites_years_species_pairs":True,
        "weather_variables":["spring_GDD5_MarJun","summer_precip_mm_JulSep"],
        "per_pair_two_independent_weather_terms":True,
        "source_weather_responses_per_pair":len(weather_r2),
        "median_fraction_variation_explained_per_partner_by_weather":float(np.median(weather_r2)),
        "fraction_partner_weather_r2_gt_0_1":float(np.mean(np.array(weather_r2)>.1)),
        "median_weather_partial_minus_unadjusted_pair_correlation":float(np.median(partial_differences)),
        "year_grouping_effective_spatial_weather_samples":len(weather),
        "observed_local_host_phenology":False,
    }


def run(accepted_host_csv):
    gate=json.loads(GATE.read_text(encoding="utf-8"))
    original=json.loads(QC.read_text(encoding="utf-8"))
    if original["status"]!="SOURCE_10H_CORRECTION_FULL_REPLAY_QUALIFIED":
        raise ValueError("Published source survey effort baseline correction missing")
    if not math.isclose(original["corrected_observed"]["difference"],SOURCE_SAMPLE_CORRECTED,abs_tol=1e-10):
        raise ValueError("Original site adult temporal effect changed")
    good,keys,field,temporal,sitecoord,classes=fixed_source_panel()
    weather=verify_weather_replay(classes,gate)
    adjusted,raw,diagnostics=annual_partial_pairs(good,keys,temporal,sitecoord,weather)
    focal=frozen_focal()
    fam=butterfly_families(focal)
    hosts,edges,groups,viable,plantfam=accepted_edges(accepted_host_csv,focal,fam)
    raw_result=meancontrast(raw,hosts,2)
    if (raw_result["shared"]!=62 or raw_result["nonshared"]!=578
        or abs(raw_result["difference"]-SOURCE_SAMPLE_CORRECTED)>1e-9):
        raise ValueError("Original source host sharing or baseline corrected contrast changed")
    observed=meancontrast(adjusted,hosts,2)
    if observed["shared"]!=62 or observed["nonshared"]!=578:
        raise ValueError("Weather correction changed fixed ecological host categories")
    null=null_test(adjusted,observed,hosts,edges,groups,viable,plantfam,fam,2)
    return {
        "schema":"chocho_sbus_illinois_coarse_weather_adjusted_host_synchrony_v01_20261026",
        "status":(
            "COARSE_14REGION_WEATHER_PROXY_SENSITIVITY_QUALIFIED_EXPLORATORY_ONLY"
            if null["status"]=="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
            else "NOT_EVALUABLE_UNDER_NETWORK_NULL"
        ),
        "original_US_butterfly_data_DOI":"10.6084/m9.figshare.27934602.v1",
        "independent_Daymet_daily_DOI":"10.3334/ORNLDAAC/2129",
        "original_source_repeated_adult_field_panel":field,
        "site_weather_accuracy":"14 public rounded regional proxy classes for 121 sites: NOT actual 1km transect temperature or plant leaf phenology",
        "verified_climate_source_gate_file":str(GATE.relative_to(ROOT)),
        "source_weather_quality":diagnostics,
        "original_600min_source_observed_pair_contrast":raw_result,
        "coarse_seasonally_weather_adjusted_observed_pair_contrast":observed,
        "observed_adjustment_change":observed["difference"]-raw_result["difference"],
        "weather_adjusted_bipartite_host_null":null,
        "scientific_limitations":[
            "Public monitored site coordinates are rounded to 14 classes; source actual site climate unresolved",
            "Modelled Daymet pixels are coarse regional proxies, not direct microhabitat weather",
            "GDD5 is climate-derived proxy, not observed larval host plant phenology",
            "No actual host plant availability, larval performance, competition, or historical host arrival date",
            "Prior outcome family contained spatial phi, Jaccard and time synchrony; all tests posthoc exploratory",
            "Species-pair correlations share sites and butterflies, and the null does not randomize butterfly census sampling",
        ],
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--accepted-host-csv",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    result=run(a.accepted_host_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    null=result["weather_adjusted_bipartite_host_null"]
    print(json.dumps({
        "status":result["status"],
        "old_600min_contrast":result["original_600min_source_observed_pair_contrast"]["difference"],
        "coarse_GDD5_plus_rain_adjusted_contrast":result["coarse_seasonally_weather_adjusted_observed_pair_contrast"]["difference"],
        "change":result["observed_adjustment_change"],
        "seasonal_weather_qc":result["source_weather_quality"],
        "host_network_null":{k:v for k,v in null.items() if k!="null_difference_history"}
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
