#!/usr/bin/env python3
"""Site-leave-out regional calendar-year adjustment to Illinois adult synchrony.

Not a weather measurement, not butterfly population interaction, and not causal.
Exactly original 640 adult temporal pairs, source-fixed 199 butterfly-host
family-constrained double-edge swaps. No site/individual data archived.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze_sbus_illinois_niche_counts_20261025 import frozen_focal,load_source,read_events
from analyze_sbus_illinois_host_degree_null_20261025 import (
    accepted_edges,butterfly_families
)
from analyze_sbus_illinois_adult_codetection_hostnull_20261026 import (
    field_panel,pair_correlations,meancontrast,null_test
)

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"provenance/reviewer_defenses/results/sbus_illinois_adult_codetection_hostnull_20261026.json"
EXPECTED_DELTA=.14323064389607087


def make_year_adjusted_pair_correlations(good,sitekeys,all_pair_records,
                                         fixed_temporal_records):
    # All original prevalence-eligible source butterflies, not outcome-chosen species.
    taxa=sorted({sp for p in all_pair_records for sp in p[:2]})
    lookup={sp:i for i,sp in enumerate(taxa)}
    bysite=defaultdict(list)
    for site,year in sitekeys:
        bysite[site].append(year)
    species_year_sum=defaultdict(float)
    species_year_n=Counter()
    baseline={}
    has_by_site={}
    annual_total=0
    for site,years in bysite.items():
        years=sorted(years)
        source=np.zeros((len(years),len(taxa)),dtype=float)
        seen=np.zeros(len(taxa),dtype=bool)
        for ti,y in enumerate(years):
            effort,counts=good[(site,y)]
            for sp,value in counts.items():
                i=lookup.get(sp)
                if i is not None and value>0:
                    seen[i]=True
                    source[ti,i]=math.log1p(value/effort)
        # The original correlation source centers each butterfly within a site
        # and only treats 'ever observed in site' taxa as a site-species panel.
        source[:,seen]-=source[:,seen].mean(axis=0)
        baseline[site]=(years,source)
        has_by_site[site]=seen
        for ti,y in enumerate(years):
            for i in np.flatnonzero(seen):
                key=(i,y)
                species_year_sum[key]+=float(source[ti,i])
                species_year_n[key]+=1
                annual_total+=1

    corrected={}
    adjusted=0
    nonadjusted=0
    sparse_by_species=Counter()
    for site,(years,src) in baseline.items():
        arr=src.copy()
        present=has_by_site[site]
        for ti,y in enumerate(years):
            for i in np.flatnonzero(present):
                key=(i,y)
                peers=species_year_n[key]
                if peers>=5:
                    # Site-excluded regional mean avoids directly reusing
                    # target site's own count when estimating annual shock.
                    leave_site_out=(species_year_sum[key]-float(src[ti,i]))/(peers-1)
                    arr[ti,i]-=leave_site_out
                    adjusted+=1
                else:
                    nonadjusted+=1
                    sparse_by_species[taxa[i]]+=1
        arr[:,present]-=arr[:,present].mean(axis=0)
        corrected[site]=arr
    assert adjusted+nonadjusted==annual_total
    base_idx={(p[0],p[1]):p[2] for p in fixed_temporal_records}
    fixed=[]
    raw_diffs=[]
    unidentifiable=[]
    change=[]
    for a,b,expected in fixed_temporal_records:
        ia=lookup[a];ib=lookup[b]
        xy=x2=y2=0.0
        ori_xy=ori_x2=ori_y2=0.0
        co_sites=0
        co_years=0
        for site,(years,src) in baseline.items():
            present=has_by_site[site]
            if not (present[ia] and present[ib]):
                continue
            co_sites+=1
            co_years+=len(years)
            xs=corrected[site][:,ia]
            ys=corrected[site][:,ib]
            ori_x=src[:,ia]
            ori_y=src[:,ib]
            xy+=float(xs@ys)
            x2+=float(xs@xs)
            y2+=float(ys@ys)
            ori_xy+=float(ori_x@ori_y)
            ori_x2+=float(ori_x@ori_x)
            ori_y2+=float(ori_y@ori_y)
        if co_sites<5 or co_years<30:
            raise ValueError("Fixed temporal species pair lost original site-year support")
        recon=ori_xy/math.sqrt(ori_x2*ori_y2)
        if not math.isclose(recon,expected,abs_tol=1e-9,rel_tol=0):
            raise ValueError(f"Original pair synchronous adult score source drift: {a},{b}")
        if x2<=1e-10 or y2<=1e-10:
            unidentifiable.append((a,b))
            continue
        val=xy/math.sqrt(x2*y2)
        if not math.isfinite(val) or abs(val)>1+1e-9:
            raise ValueError("Annual calendar-year-corrected pair correlation out of bounds")
        fixed.append((a,b,float(np.clip(val,-1,1))))
        change.append(val-expected)
    if unidentifiable or len(fixed)!=640:
        raise ValueError(f"Cannot retain all 640 original pair scores: {len(unidentifiable)} invalid")
    summary={
        "original_taxa_prevalence_eligible":len(taxa),
        "original_fixed_temporal_pairs":len(fixed),
        "eligible_species_site_year_anomaly_values":annual_total,
        "calendar_year_effect_identifiable_n_values":adjusted,
        "fraction_calendar_year_leave_siteout_adjusted":adjusted/annual_total,
        "calendar_year_n_other_sites_less_than_4":nonadjusted,
        "eligible_species_year_n_classes":len(species_year_n),
        "eligible_species_year_n_ge5":sum(v>=5 for v in species_year_n.values()),
        "source_original_minimum_pair_sample_held_exactly":True,
        "site_location_raw_values_not_persisted":True,
        "pairwise_correlation_source_difference_median":float(np.median(change)),
        "pairwise_correlation_source_difference_range":[min(change),max(change)],
    }
    return fixed,summary


def run(hostcsv:Path):
    orig=json.loads(RECEIPT.read_text(encoding="utf-8"))
    if (orig["adult_field_panel"]["qualified_sites"]!=121
        or orig["adult_field_panel"]["qualified_site_years"]!=851
        or orig["sample"]["within_site_temporal_pair_count"]!=640
        or abs(orig["within_site_temporal"]["observed"]["difference"]-EXPECTED_DELTA)>1e-10):
        raise ValueError("Original site/pair cohort or effect not source-fixed")
    f=frozen_focal()
    fam=butterfly_families(f)
    hosts,edgepairs,groups,viable,plantfam=accepted_edges(hostcsv,f,fam)
    with tempfile.TemporaryDirectory() as tmp:
        archive=Path(tmp)/"sbus.zip"
        load_source(archive)
        events,_=read_events(archive,f)
    good,sitekeys,species,presence,field=field_panel(events)
    original=pair_correlations(good,sitekeys,species,presence)
    if original["status"]!="FIXED_SPECIES_PAIRS_READY":
        raise ValueError("Original prevalence-eligible species-pair sample source gate failed")
    old=original["within_site_temporal_candidate_pairs"]
    baseline=meancontrast(old,hosts,2)
    if not (len(old)==640 and baseline["shared"]==62 and baseline["nonshared"]==578
            and math.isclose(baseline["difference"],EXPECTED_DELTA,abs_tol=1e-10)):
        raise ValueError("Previous adult cofluctuation source changed")
    corrected,qc=make_year_adjusted_pair_correlations(good,sitekeys,original["pairs"],old)
    new=meancontrast(corrected,hosts,2)
    if new["status"]!="OK" or new["shared"]!=62 or new["nonshared"]!=578:
        raise ValueError("Consumer host class changed after year anomaly adjustment")
    null=null_test(corrected,new,hosts,edgepairs,groups,viable,plantfam,fam,2)
    return {
        "schema":"chocho_illinois_species_year_common_shock_reduction_v01_20261026",
        "status":(
            "EXPLORATORY_YEAR_SHOCK_SENSITIVITY_SOURCE_QUALIFIED"
            if null["status"]=="EXPLORATORY_SOURCE_CONDITIONED_BIPARTITE_NULL_QUALIFIED"
            else "NOT_EVALUABLE_QUALIFICATION_OR_SOURCE_GATE"
        ),
        "source_doi":"10.6084/m9.figshare.27934602.v1",
        "daymet_site_match":"NOT_EVALUABLE_ROUNDED_14_UNIQUE_POSITIONS_FOR_121_SITES",
        "not_independent_weather_data":True,
        "field_panel":field,
        "calendar_year_correction":qc,
        "original_observed":baseline,
        "corrected_observed":new,
        "corrected_minus_original_delta":new["difference"]-baseline["difference"],
        "source_degree_graph_null":null,
        "inference_ceiling":"Observed adult count-rate correlations after species-specific leave-site-out shared calendar-year adjustment; cannot identify meteorology, botanical phenology, host use, competition, colonization, or true population growth.",
        "multiple_testing_note":"Exploratory follow-up to multiple already evaluated spatial/temporal original outcomes. Nominal host-network tail is not confirmatory p value.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--accepted-host-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=run(args.accepted_host_csv)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    null=result["source_degree_graph_null"]
    print(json.dumps({
        "status":result["status"],
        "original_delta":result["original_observed"]["difference"],
        "year_adjusted_delta":result["corrected_observed"]["difference"],
        "change":result["corrected_minus_original_delta"],
        "year_effect_qc":result["calendar_year_correction"],
        "network_null":{k:v for k,v in null.items() if k not in ("null_difference_history",)}
    },indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
