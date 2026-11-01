#!/usr/bin/env python3
"""Independent 224-species GBIF butterfly beta versus host and alien floras.

Separate response data: 2001–2024 GBIF original 14,637,152 indexed records
stream-aggregated by species×WGSRPD3×year and provider/date.

Potential host predictors: contemporaneous HOSTS and GloNAF/WCVP. Source
flora is NOT dated and may postdate GBIF records; do NOT claim that this is a
prospective historical hindcast, butterfly colonization or actual host use.

Primary target is 2018–2024 GBIF-recorded community similarity (not
true butterfly occupancy). Pre-2011 indexed butterfly similarity and the
2001–2010/2018–2024 provider-day observation opportunities are controls.
Every model gets the identical selected region pairs and source-geography
covariates. Spatial transfer is evaluated at complete held-out regions.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_colonial_legacy_resource_beta import (
    load as load_published_regional_controls,
    orthogonal_controls,residual,correlation,
)
from analyze_matched_glonaf_consumer_history import produce_matched_regions
from analyze_documented_host_specificity import (
    GLO_NAF,TDWG,HOSTS,glonaf_index,pair_jaccard,
)
from analyze_gbif_31_observed_host_beta import (
    canonical_exact_hosts,
    fit_ridge_predict,
)

ROOT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031")
ANNUAL=ROOT/"gbif_global224_annual_v1"
SOURCE_AUDIT=Path(str(ANNUAL)+".audit.json")
SOURCE_CELLS=Path(str(ANNUAL)+".species_region_year.csv.gz")
SOURCE_EFFORT=Path(str(ANNUAL)+".region_year_effort.csv")
LOCK=ROOT/"gbif_global224_observed_beta_lock_v1.json"
OUTPUT=ROOT/"gbif_global224_observed_beta_v1.json"
SCORES=ROOT/"gbif_global224_observed_beta_region_holdout_v1.csv"
PANEL=Path("data/frozen/figure_sources/anthropogenic_species_metrics.csv")
PERIODS={"prior":range(2001,2011),"middle":range(2011,2018),
         "later":range(2018,2025)}


def source_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path):
    if str(path).endswith(".gz"):
        file=gzip.open(path,"rt",encoding="utf8",newline="")
    else:
        file=Path(path).open("r",encoding="utf8",newline="")
    with file as f:
        yield from csv.DictReader(f)


def load_annual_data(audit=SOURCE_AUDIT,cells=SOURCE_CELLS,effort=SOURCE_EFFORT):
    source=json.loads(audit.read_text(encoding="utf8"))
    if source["schema"]!="chocho_global_GBIF_224_butterfly_369_region_year_effort_v1":
        raise RuntimeError("Wrong source annual data schema")
    if source["sampled_mode"] or source["records"]["rows"]!=14637152:
        raise RuntimeError("The full GBIF source has not been validated")
    if source["coverage"]["original_candidate_species"]!=225:
        raise RuntimeError("Different original butterfly cohort")
    for p in (cells,effort):
        if source["files"][p.name]["sha256"]!=source_sha(p):
            raise RuntimeError("Source-verified independent GBIF derived table checksum drift: "+p.name)
    snapshots={
        "prior":defaultdict(set),"middle":defaultdict(set),
        "later":defaultdict(set),"all":defaultdict(set)
    }
    year_counts=Counter()
    mapped_species=set()
    for r in read_rows(cells):
        name=r["butterfly_species"].strip()
        area=r["wgsrpd3_code"].strip()
        year=int(r["year"])
        n=int(r["indexed_gbif_records"])
        if not (2001<=year<=2024 and n>0 and area):
            raise RuntimeError("GBIF cell contains negative or outside-period source")
        mapped_species.add(name)
        snapshots["all"][area].add(name)
        for p,years in PERIODS.items():
            if year in years:
                snapshots[p][area].add(name)
        year_counts[year]+=n
    if sum(year_counts.values())!=source["records"]["candidate_records_mapped_to_369"]:
        raise RuntimeError("GBIF source group counts fail full row audit")
    base=defaultdict(lambda:defaultdict(lambda:{"records":0,"days":0,"providers":0}))
    for r in read_rows(effort):
        area=r["wgsrpd3_code"].strip()
        year=int(r["year"])
        for p,years in PERIODS.items():
            if year in years:
                b=base[area][p]
                b["records"]+=int(r["indexed_candidate_butterfly_records"])
                b["days"]+=int(r["distinct_gbif_datasetkey_calendar_days"])
                # This is an upper bound sum of dataset/year uniques, not
                # distinct dataset keys across the whole seven-year period.
                b["providers"]+=int(r["distinct_gbif_dataset_keys"])
    if len(mapped_species)<100:
        raise RuntimeError("Source lacks adequate original butterfly species")
    return snapshots,base,source,year_counts


def species_masks(by_region,selected):
    taxa=sorted({name for r in selected for name in by_region.get(r,set())})
    ix={name:j for j,name in enumerate(taxa)}
    masks=[]
    for region in selected:
        bits=0
        for name in by_region.get(region,set()):
            bits|=1<<ix[name]
        masks.append(bits)
    return masks,taxa


def jaccard(masks,i,j):
    a,n_missing=pair_jaccard(masks,i,j,undefined_zero=True)
    return np.asarray(a,float),int(n_missing)


def score_held_out(y,models,ia,ib,region_codes,world,lam=10.0):
    if not models:
        raise ValueError("No models")
    n=len(region_codes)
    if n<35:
        raise RuntimeError("Too few source-supported regions for leave-region-out tests")
    model_matrices={name:np.column_stack(values) for name,values in models.items()}
    scores={name:[] for name in models}
    # Region-held-out unit is the original WGSRPD3 region, not dependent pair.
    for heldout in range(n):
        testing=(ia==heldout)|(ib==heldout)
        training=~testing
        if testing.sum()<20 or training.sum()<100:
            raise RuntimeError("Poorly supported heldout-region test")
        for name,X in model_matrices.items():
            pred=fit_ridge_predict(training,testing,y,X,lam)
            err=pred-y[testing]
            scores[name].append({
                "wgsrpd3_code":region_codes[heldout],
                "world_division":str(world[heldout]),
                "test_pairs":int(testing.sum()),
                "training_pairs":int(training.sum()),
                "mse":float(np.mean(err*err)),
                "mae":float(np.mean(np.abs(err))),
            })
    return scores


def consolidate_scores(folds,main="M0_all_flora_effort_prior"):
    scores={}
    for key,rows in folds.items():
        scores[key]={
            "region_equal_MSE":float(np.mean([x["mse"] for x in rows])),
            "region_equal_MAE":float(np.mean([x["mae"] for x in rows])),
            "region_holdout_count":len(rows),
        }
    names=("M1_known_host","M2_other_aliens","M3_both")
    pair_comp={}
    for name in names:
        differences=np.asarray([b["mse"]-x["mse"] for b,x in zip(folds[main],folds[name])])
        pair_comp[name]={
            "MSE_reduction_vs_M0_percent":float(100*np.mean(differences)/scores[main]["region_equal_MSE"]),
            "regions_improved":int(sum(differences>0)),
            "regions_worsened":int(sum(differences<0)),
            "mse_diff_mean":float(np.mean(differences)),
            "mse_diff_median":float(np.median(differences)),
        }
    for numerator,denominator,key in [
        ("M3_both","M2_other_aliens","known_host_after_other_aliens"),
        ("M3_both","M1_known_host","other_aliens_after_known_host")
    ]:
        pair_comp[key]={
            "MSE_reduction_percent":100*(
                scores[denominator]["region_equal_MSE"]-
                scores[numerator]["region_equal_MSE"]
            )/scores[denominator]["region_equal_MSE"]
        }
    return scores,pair_comp


def world_holdout(y,models,ia,ib,codes,level1,lam):
    matrices={k:np.column_stack(v) for k,v in models.items()}
    results={}
    for world in sorted(set(level1)):
        is_held=np.asarray([w==world for w in level1],bool)
        testing=is_held[ia]|is_held[ib]
        training=~testing
        if testing.sum()<30 or training.sum()<100:
            continue
        result={"regions_heldout":int(is_held.sum()),
                "training_pairs":int(training.sum()),"test_pairs":int(testing.sum())}
        for name,X in matrices.items():
            pred=fit_ridge_predict(training,testing,y,X,lam)
            result[name]=float(np.mean((pred-y[testing])**2))
        result["known_host_vs_base_pct"]=100*(
            result["M0_all_flora_effort_prior"]-result["M1_known_host"]
        )/result["M0_all_flora_effort_prior"]
        result["other_aliens_vs_base_pct"]=100*(
            result["M0_all_flora_effort_prior"]-result["M2_other_aliens"]
        )/result["M0_all_flora_effort_prior"]
        result["known_host_given_other_pct"]=100*(
            result["M2_other_aliens"]-result["M3_both"]
        )/result["M2_other_aliens"]
        results[world]=result
    return {
        "world_divisions_evaluable":len(results),
        "world_divisions_known_host_improved":sum(v["known_host_vs_base_pct"]>0 for v in results.values()),
        "world_divisions_other_aliens_improved":sum(v["other_aliens_vs_base_pct"]>0 for v in results.values()),
        "world_divisions_known_host_given_other_improved":sum(v["known_host_given_other_pct"]>0 for v in results.values()),
        "by_world_division":results,
    }


REGION_GATES={
    "primary":(5,10,20),
    "relaxed":(5,5,10),
    "stringent":(10,20,50),
}


def analyse(eligibility="primary"):
    if eligibility not in REGION_GATES:
        raise ValueError("Unknown predeclared observational-effort eligibility gate")
    species_needed,prior_days_needed,later_days_needed=REGION_GATES[eligibility]
    lock=json.loads(LOCK.read_text(encoding="utf8"))
    if lock["schema"]!="chocho_global224_2001_2024_GBIF_observed_beta_lock_v1":
        raise RuntimeError("Study design lock changed")
    communities,effort,source,source_years=load_annual_data()
    beta=load_published_regional_controls()
    region_codes=beta["ids"]
    if len(region_codes)!=284:
        raise RuntimeError("Colonial matched source changed from 284 regions")
    def support(code,p):
        return effort[code][p]["days"]
    regions=[
        code for code in region_codes
        if len(communities["all"].get(code,set()))>=species_needed
        and support(code,"prior")>=prior_days_needed
        and support(code,"later")>=later_days_needed
    ]
    if len(regions)<35:
        raise RuntimeError("Too few prior-and-late provider-day-supported 284 regions: "+str(len(regions)))
    index={code:j for j,code in enumerate(regions)}
    f=(np.fromiter((region_codes[int(a)] in index and region_codes[int(b)] in index
                    for a,b in zip(beta["ia"],beta["ib"])),bool,count=len(beta["ia"])))
    original_i=beta["ia"][f]
    original_j=beta["ib"][f]
    pair_i=np.array([index[region_codes[int(i)]] for i in original_i],int)
    pair_j=np.array([index[region_codes[int(j)]] for j in original_j],int)
    if len(pair_i)<400:
        raise RuntimeError("Insufficient geographically matched independent butterfly pairs")
    y_by_period={}
    observed_taxa={}
    misses={}
    for period in ("prior","middle","later","all"):
        masks,sp=species_masks(communities[period],regions)
        y_by_period[period],misses[period]=jaccard(masks,pair_i,pair_j)
        observed_taxa[period]=len(sp)
    # The selected period windows should not produce both-empty pairs due to
    # the fixed pre/post source recording opportunity requirement. Report and
    # exclude any such pair on ALL models rather than code missing as 0.
    keep=np.fromiter((bool(communities["prior"].get(regions[int(a)])
                            or communities["prior"].get(regions[int(b)]))
                      and bool(communities["later"].get(regions[int(a)])
                            or communities["later"].get(regions[int(b)]))
                      for a,b in zip(pair_i,pair_j)),bool,count=len(pair_i))
    if int((~keep).sum())>0:
        pair_i=pair_i[keep];pair_j=pair_j[keep]
        f_indices=np.flatnonzero(f)[keep]
        for p in y_by_period:
            y_by_period[p]=y_by_period[p][keep]
    else:
        f_indices=np.flatnonzero(f)
    outcome=y_by_period["later"]
    prior=y_by_period["prior"]
    plant_masks,_,glo_audit=produce_matched_regions(GLO_NAF,TDWG,(.05,))
    fullflora=[plant_masks[.05][code] for code in regions]
    # Focal hosts must be restricted to the butterflies with independently
    # GBIF-indexed records, not the 239 GEB species nor 30 pilot alone.
    gbif_indexed_species=set().union(*(communities["all"].get(c,set()) for c in regions))
    butterfly_host_names=canonical_exact_hosts(
        {n.casefold() for n in gbif_indexed_species},HOSTS)
    aliases,glo_names=glonaf_index(GLO_NAF)
    known_host_names=set().union(*butterfly_host_names.values())
    known_ids={j for name in known_host_names&set(aliases) for j in aliases[name]}
    if len(known_ids)<80:
        raise RuntimeError("Insufficient GloNAF source-matched documented larval hosts for global panel")
    mask=sum(1<<j for j in known_ids)
    hosts=[m&mask for m in fullflora]
    other=[m&~mask for m in fullflora]
    all_alien_j,miss_all=jaccard(fullflora,pair_i,pair_j)
    host_j,miss_hosts=jaccard(hosts,pair_i,pair_j)
    other_j,miss_other=jaccard(other,pair_i,pair_j)
    covariates=beta["control_sets"]["colonial_vs_resource_beta_after_flora"]
    controls=[np.asarray(v)[f_indices] for v in covariates]
    if len(controls)!=9:
        raise RuntimeError("Original 284 matched geography controls changed")
    history=np.asarray(beta["exposure"])[f_indices]
    flora_rich=np.asarray([m.bit_count() for m in fullflora],float)
    def minmax(v):
        return (np.log1p(np.minimum(v[pair_i],v[pair_j])),
                np.log1p(np.maximum(v[pair_i],v[pair_j])))
    old_days=np.asarray([support(code,"prior") for code in regions],float)
    late_days=np.asarray([support(code,"later") for code in regions],float)
    prior_records=np.asarray([effort[c]["prior"]["records"] for c in regions],float)
    later_records=np.asarray([effort[c]["later"]["records"] for c in regions],float)
    # Records are a separate outcome-side nuisance, not an ecological
    # predictor available at the historical start of the later period.
    base=controls+[history,all_alien_j,*minmax(flora_rich),prior,
                   *minmax(old_days),*minmax(late_days),
                   *minmax(prior_records),*minmax(later_records)]
    if any(len(x)!=len(outcome) for x in base):
        raise RuntimeError("GBIF/GloNAF regional predictor length mismatch")
    models={
        "M0_all_flora_effort_prior":base,
        "M1_known_host":base+[host_j],
        "M2_other_aliens":base+[other_j],
        "M3_both":base+[host_j,other_j],
    }
    base_host=orthogonal_controls(base+[other_j])
    base_other=orthogonal_controls(base+[host_j])
    partials={
       "known_hosts_given_other_aliens":correlation(
           residual(host_j,base_host),residual(outcome,base_host)),
       "not_documented_aliens_given_known_hosts":correlation(
           residual(other_j,base_other),residual(outcome,base_other)),
       "not_iid_pair_significance":True,
    }
    lookup=json.loads(Path("data/frozen/wgsrpd3_level1_lookup_v0.1.json").read_text())["units"]
    world=[lookup[c] for c in regions]
    folds=score_held_out(outcome,models,pair_i,pair_j,regions,world,lam=10.)
    summary,comparisons=consolidate_scores(folds)
    worldscores=world_holdout(outcome,models,pair_i,pair_j,regions,world,lam=10.)
    source_rows={
        "region_count":len(regions),"evaluated_region_pairs":len(pair_i),
        "source_284_region_possible_pairs":len(beta["ia"]),
        "global_source_mapped_candidate_species":source["coverage"]["candidate_species_matched_to_a_369_region"],
        "source_225_candidate_names_all_subset_239":True,
        "positive_GBIF_species_in_matched_region_cohort":len(gbif_indexed_species),
        "unique_observed_taxa_by_period":observed_taxa,
        "valid_gbif_source_record_years":source_years,
        "World_division_region_counts":dict(Counter(world)),
        "region_recording_eligibility":eligibility,
        "region_recording_threshold":(
            f"all-time >={species_needed} indexed candidate butterfly species, "
            f"prior datasetKey-days>={prior_days_needed}, "
            f"late datasetKey-days>={later_days_needed}"
        ),
        "regions_with_prior_year_provider_days_at_threshold":sum(
            support(c,"prior")>=prior_days_needed for c in region_codes),
        "regions_with_late_year_provider_days_at_threshold":sum(
            support(c,"later")>=later_days_needed for c in region_codes),
        "region_pair_empty_butterfly_portfolios_by_period":misses,
        "pairs_both_empty_GloNAF_alien":miss_all,
        "pairs_both_empty_known_host":miss_hosts,
        "pairs_both_empty_not_documented_alien":miss_other,
        "focal_GBIF_butterfly_known_host_binomial_names":len(known_host_names),
        "GloNAF_documented_host_names":len(known_host_names&set(aliases)),
        "GloNAF_original_plant_raw_indices_known_hosts":len(known_ids),
        "GloNAF_source_taxa_count":glo_names["distinct_glonaf_names"],
    }
    return {
        "schema":"chocho_global224_GBIF_observed_WGSRPD3_2018_2024_beta_v1",
        "status":"RECONSTRUCTED_ORIGINAL_OBSERVED_BUTTERFLY_ASSOCIATION_NOT_ARRIVAL_OR_TRUE_OCCUPANCY",
        "sources":{
          "annual_GBIF_audit_sha256":source_sha(SOURCE_AUDIT),
          "annual_GBIF_species_region_year_sha256":source_sha(SOURCE_CELLS),
          "annual_GBIF_region_year_effort_sha256":source_sha(SOURCE_EFFORT),
          "method_lock_sha256":source_sha(LOCK),
          "GBIF_doi":"10.15468/dl.umqmsf",
          "GloNAF_zenodo_doi":"10.5281/zenodo.6640264",
          "original_TDWG_commit":"52da7828aba9d461dd133c27b3bd7a4407161f54",
          "source_WCVP_HOSTS_239_beta":"source-locked original GEB 239 resource panel, not 224-only native-resource control",
        },
        "cohort":source_rows,
        "outcome_Jaccard_summary":{
          "late_2018_2024_mean":float(np.mean(outcome)),
          "late_2018_2024_median":float(np.median(outcome)),
          "late_2018_2024_zero_pair_fraction":float(np.mean(outcome==0)),
          "prior_2001_2010_mean":float(np.mean(prior)),
          "prior_2001_2010_zero_pair_fraction":float(np.mean(prior==0)),
          "current_all_GloNAF_mean":float(np.mean(all_alien_j)),
          "documented_focal_HOSTS_GloNAF_mean":float(np.mean(host_j)),
          "unregistered_hosts_GloNAF_mean":float(np.mean(other_j)),
        },
        "partial_rank":partials,
        "heldout_region_model_scores":summary,
        "heldout_region_model_comparisons":comparisons,
        "heldout_world_division_results":worldscores,
        "limitations":[
          "The source 14,637,152 GBIF records are a preselected candidate butterfly cohort, not all Lepidoptera.",
          "GBIF absence is not established absence; provider-day counts only proxy recording opportunity.",
          "Data from after 2018 are used to select/adjust later sample: this is retrospective prediction only.",
          "Source-colonial variables measure area-weighted historical empire duration, not plant/voyage arrival dates.",
          "GloNAF and HOSTS are contemporary species records and may postdate GBIF early records.",
          "Known larval-host lists are catalogued relations, not independently demonstrated larval feeding in each region.",
          "Regional pairs are dependent; no independent pairwise p-values or causal mediated effects.",
          "The control native-resource geography is original frozen GEB 239, not source 224 alone.",
        ],
        "folds":folds,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=OUTPUT)
    ap.add_argument("--region-fold-csv",type=Path,default=SCORES)
    ap.add_argument("--eligibility",choices=tuple(REGION_GATES),default="primary")
    args=ap.parse_args()
    report=analyse(eligibility=args.eligibility)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    folds=report.pop("folds")
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    with args.region_fold_csv.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["model","wgsrpd3_code","world_division",
            "test_pairs","training_pairs","mse","mae"])
        w.writeheader()
        for name,rows in folds.items():
            for row in rows:
                w.writerow({"model":name,**row})
    print(json.dumps({
        "cohort":report["cohort"],"Jaccard":report["outcome_Jaccard_summary"],
        "partial":report["partial_rank"],
        "scores":report["heldout_region_model_scores"],
        "contrasts":report["heldout_region_model_comparisons"],
        "world":report["heldout_world_division_results"],
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
