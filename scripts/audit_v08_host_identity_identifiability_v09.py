#!/usr/bin/env python3
"""Independently audit frozen v0.8 known-host vs modern-climate-matched other-plant null.

No new plant choice or host-climate model fit. Source v0.8 is exploratory, after
prior host/climate findings, so uncertainty is descriptive and not confirmatory.

Unusable exact-four-way locality edges are NOT negative outcomes. Repeated plant
control groups are technical comparisons, never new independent insect species.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from scipy.stats import spearmanr

SOURCE_SHA={
  "species":"d5e5e374e89eeac701c449e05edd775e7b0f0014c714b965a03378af0c7fa4b4",
  "groups":"2574cd7c6b346720733c232b2044e9c4d88407b0d597264b5ffcf115ef7cb261",
  "summary":"7fcd67aa367d138960aa7d0d94e6482139a5678a1c2a4be08a38dcbe4264bc55"
}
SCENARIOS=("analog","extrapolated")
RADII=(50,100,150)
EPS=2e-12

def file_sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda:fh.read(1<<20),b""):h.update(block)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def number(value):
    if value is None or value=="":
        return None
    x=float(value)
    if not np.isfinite(x):
        raise ValueError("unexpected NaN or infinite source comparison")
    return x

def equal_numeric(actual,expected,field):
    if (actual is None)!=(expected is None):
        raise ValueError(f"missingness changed: {field}")
    if actual is not None and abs(float(actual)-float(expected))>EPS:
        raise ValueError(f"source v08 numerical identity changed: {field}")

def audit_group_row(row):
    n=int(row["fourway_exact_common_original_edges_n"])
    all_original=int(row["all_original_edges_n"])
    if not 0<=n<=all_original:
        raise ValueError("source same-edge denominator inconsistent")
    numeric=("actual_0BP_same_fraction","actual_21ka_same_fraction",
             "alternate_0BP_same_fraction","alternate_21ka_same_fraction",
             "actual_paired_delta","alternate_paired_delta",
             "actual_minus_alternate_paired_delta","modern_0BP_baseline_same_patch_gap")
    if n==0:
        if any(number(row[field]) is not None for field in numeric):
            raise ValueError("zero-overlap source-group fabricated a climate effect")
        return None
    actual0=number(row["actual_0BP_same_fraction"])
    actual21=number(row["actual_21ka_same_fraction"])
    alt0=number(row["alternate_0BP_same_fraction"])
    alt21=number(row["alternate_21ka_same_fraction"])
    for x in (actual0,actual21,alt0,alt21):
        if x is None or not 0<=x<=1:
            raise ValueError("invalid source same-patch fraction")
    actual_delta=actual21-actual0
    alt_delta=alt21-alt0
    diff=actual_delta-alt_delta
    for field,expected in {
        "actual_paired_delta":actual_delta,
        "alternate_paired_delta":alt_delta,
        "actual_minus_alternate_paired_delta":diff,
        "modern_0BP_baseline_same_patch_gap":actual0-alt0
    }.items():
        equal_numeric(number(row[field]),expected,field)
    for which in ("actual","alternate"):
        split=int(row[which+"_same_to_different_n"])
        merged=int(row[which+"_different_to_same_n"])
        if split<0 or merged<0 or split+merged>n:
            raise ValueError("source split/merge counts exceed same original edges")
        equal_numeric(number(row[which+"_paired_delta"]),(merged-split)/n,
                      which+" transition identity")
    return {"n":n,"diff":diff,"baseline":actual0-alt0,
            "modern_alternative_moment_gap":float(row["original_0BP_climate_moment_gap"])}

def audit_species_row(row,groups):
    n_good=[x for x in groups if int(x["fourway_exact_common_original_edges_n"])>0]
    n_stable=[x for x in n_good if int(x["fourway_exact_common_original_edges_n"])>=5]
    if len(groups)!=int(row["source_frozen_other_plant_group_n"]):
        raise ValueError("original pre-outcome donor group count drift")
    if len(n_good)!=int(row["other_plant_groups_with_common_original_edges_n"]):
        raise ValueError("original compatible donor group count drift")
    if len(n_stable)!=int(row["other_plant_groups_with_at_least_5_common_edges_n"]):
        raise ValueError("source five original edges denominator drift")
    if not n_good:
        if row["status"]!="NOT_EVALUABLE_FOURWAY_COMMON_EDGE_DENOMINATOR":
            raise ValueError("unevaluable source insect silently called evaluable")
        if number(row["median_actual_minus_otherplant_same_edge_difference"]) is not None:
            raise ValueError("missing modern anchor comparisons silently estimated")
        return
    if row["status"]!="COMPUTABLE_PAIRWISE_OTHER_PLANT_NEGATIVE_CONTROL":
        raise ValueError("source insect result status incorrect")
    metrics={
      "median_actual_minus_otherplant_same_edge_difference":np.median(
          [number(x["actual_minus_alternate_paired_delta"]) for x in n_good]),
      "median_0BP_actual_vs_otherplant_baseline_gap":np.median(
          [number(x["modern_0BP_baseline_same_patch_gap"]) for x in n_good]),
      "median_fourway_same_edge_n_per_group":np.median(
          [int(x["fourway_exact_common_original_edges_n"]) for x in n_good]),
      "median_difference_on_groups_with_5plus_edges":np.median(
          [number(x["actual_minus_alternate_paired_delta"]) for x in n_stable])
          if n_stable else None,
    }
    for field,val in metrics.items():
        equal_numeric(number(row[field]),val,field)

def grouped_cluster_bootstrap(species,field,seed=20261020,draws=2000):
    if draws<10 or not species:
        raise ValueError("need multiple independent genus groups")
    by_genus=defaultdict(list)
    for row in species:
        genus=row["insect_species"].split()[0]
        by_genus[genus].append(float(row[field]))
    clusters=[np.asarray(vals,float) for _,vals in sorted(by_genus.items())]
    rng=np.random.default_rng(seed)
    means=np.empty(draws);medians=np.empty(draws)
    for d in range(draws):
        inds=rng.integers(0,len(clusters),size=len(clusters))
        score=np.concatenate([clusters[k] for k in inds])
        means[d]=np.mean(score);medians[d]=np.median(score)
    return {
      "genus_clusters":len(clusters),"bootstrap_draws":draws,"bootstrap_seed":seed,
      "median_95pct_descriptive_interval":np.quantile(medians,[.025,.975]).tolist(),
      "mean_95pct_descriptive_interval":np.quantile(means,[.025,.975]).tolist(),
      "captures_shared_plant_donors_across_genera":False,
      "confirmatory_inference":False}

def describe_species(rows):
    if not rows:return {"species":0}
    diff=np.array([float(x["median_actual_minus_otherplant_same_edge_difference"]) for x in rows])
    gap=np.array([float(x["median_0BP_actual_vs_otherplant_baseline_gap"]) for x in rows])
    eligible_hosts={x["_hostset"] for x in rows}
    genus={x["insect_species"].split()[0] for x in rows}
    rho=float(spearmanr(diff,gap).statistic) if len(np.unique(diff))>1 and len(np.unique(gap))>1 else None
    return {"species":len(rows),"genus_n":len(genus),
        "different_known_host_combinations":len(eligible_hosts),
        "actual_host_more_split_species":int((diff<0).sum()),
        "other_plant_more_split_species":int((diff>0).sum()),
        "ties_no_difference_species":int((diff==0).sum()),
        "median_actual_minus_otherplant_delta":float(np.median(diff)),
        "mean_actual_minus_otherplant_delta":float(np.mean(diff)),
        "median_abs_modern_baseline_same_patch_gap":float(np.median(np.abs(gap))),
        "baseline_gap_nonzero_species":int((gap!=0).sum()),
        "baseline_delta_spearman_rho":rho,
        "baseline_vs_difference_mathematically_coupled":True,
        "species_with_exact_zero_modern_baseline_gap":int((gap==0).sum()),
        "genus_cluster_descriptive_bootstrap":grouped_cluster_bootstrap(rows,"median_actual_minus_otherplant_same_edge_difference")}

def validate_and_audit(species_file,groups_file,summary_file):
    for path,key in ((species_file,"species"),(groups_file,"groups"),(summary_file,"summary")):
        if file_sha256(path)!=SOURCE_SHA[key]:
            raise ValueError(f"source v08 SHA changed: {key}")
    source=json.loads(Path(summary_file).read_text())
    if (source["status"]!="EXPLORATORY_HOST_SPECIFIC_PATCH_CLIMATE_NEGATIVE_CONTROL_NOT_GENETIC_OR_OCCUPANCY"
        or source["source_qualified_actual_hosts_plus_four_otherplant_groups"]!=292
        or source["original_0BP_only_altplant_control_group_count"]!=2108):
        raise ValueError("source eligibility changed")
    species=read_csv(species_file); groups=read_csv(groups_file)
    if len(species)!=292*2*3 or len(groups)!=2108*2*3:
        raise ValueError("v08 original 292/2108 group x scenario/radius table incomplete")
    by_group=defaultdict(list)
    hostsets=defaultdict(set);plant_controls=defaultdict(set)
    per_group=[]
    source_scenario=Counter()
    for row in groups:
        sp=row["insect_species"];scenario=row["scenario"];radius=int(row["radius_km"])
        if scenario not in SCENARIOS or radius not in RADII:
            raise ValueError("source scenario or predeclared radius changed")
        key=(sp,scenario,radius)
        original= audit_group_row(row)
        by_group[key].append(row)
        hostsets[sp].add(row["known_actual_host_ids"])
        for donor in row["other_plant_ids"].split(";"):
            if donor:plant_controls[sp].add(donor)
        source_scenario[(scenario,radius)]+=1
        per_group.append({"insect_species":sp,"scenario":scenario,"radius_km":radius,
            "donor_group":int(row["other_plant_group_index"]),
            "fourway_original_edges":int(row["fourway_exact_common_original_edges_n"]),
            "paired_actual_minus_otherplant_delta":original["diff"] if original else None,
            "modern_0BP_baseline_gap":original["baseline"] if original else None})
    if not all(len(hostsets[k])==1 for k in hostsets):
        raise ValueError("known host plant set changed across source scenarios/groups")
    sp_index={}
    for row in species:
        key=(row["insect_species"],row["scenario"],int(row["radius_km"]))
        if key in sp_index:raise ValueError("same species scenario/radius duplicated")
        if key not in by_group:raise ValueError("species missing source selected donor plants")
        audit_species_row(row,by_group[key])
        sp_index[key]=row
    if set(sp_index)!=set(by_group) or len(hostsets)!=292:
        raise ValueError("not all source 292 insects included")
    output={}
    for mode in SCENARIOS:
        for radius in RADII:
            key=f"{mode}_{radius}km"
            candidates=[dict(r,_hostset=next(iter(hostsets[r["insect_species"]])))
                        for r in species if r["scenario"]==mode and int(r["radius_km"])==radius]
            contactable=[r for r in candidates if r["status"].startswith("COMPUTABLE")]
            stable=[r for r in contactable
                if int(r["other_plant_groups_with_at_least_5_common_edges_n"])>=4]
            gs=[g for g in groups if g["scenario"]==mode and int(g["radius_km"])==radius]
            with_pairs=[g for g in gs if int(g["fourway_exact_common_original_edges_n"])>0]
            with_five=[g for g in gs if int(g["fourway_exact_common_original_edges_n"])>=5]
            vals=np.array([float(g["actual_minus_alternate_paired_delta"]) for g in with_pairs])
            source_table=source["by_scenario_radius"][key]
            if (source_table["species_with_at_least_one_comparable_otherplant_control"]!=len(contactable)
                or source_table["species_with_at_least_four_groups_of_five_common_pairs"]!=len(stable)
                or source_table["source_groups_contactable"]!=len(with_pairs)):
                raise ValueError("per-group or per-insect denominator cannot reproduce original v08")
            b=describe_species(contactable)
            if (b["actual_host_more_split_species"]!=source_table["negative_species_actual_host_loses_more"]
                or b["other_plant_more_split_species"]!=source_table["positive_species_otherplant_loses_more"]
                or b["ties_no_difference_species"]!=source_table["equal_species"]):
                raise ValueError("source sign counts not reproduced")
            output[key]={
              "original_pretreatment_eligible_insects":len(candidates),
              "frozen_distinct_donor_plant_taxa":len(set().union(*[plant_controls[s["insect_species"]] for s in candidates])),
              "source_groups":len(gs),
              "groups_no_fourway_common_edge":len(gs)-len(with_pairs),
              "groups_at_least_one_common_edge":len(with_pairs),
              "groups_at_least_five_common_edges":len(with_five),
              "groups_exact_zero_effect_among_contactable":int(np.count_nonzero(vals==0)),
              "group_effect_median":float(np.median(vals)) if len(vals) else None,
              "contactable_fourway_edge_count_median":float(np.median(
                   [int(g["fourway_exact_common_original_edges_n"]) for g in with_pairs])) if with_pairs else None,
              "species_at_least_one_donor":b,
              "species_at_least_four_donors_with_five_pairs":describe_species(stable),
              "shared_plant_donor_across_insect_taxa_caveat":True}
    report={"schema":"chocho_v09_v08_source_bound_host_identity_climate_patches_identifiability_audit",
       "status":"EXPLORATORY_NEGATIVE_HOST_IDENTITY_CONTROL_WITH_LIMITED_COMMON_EDGE_INFORMATION",
       "source_sha256":dict(SOURCE_SHA),
       "source_preselected_species":292,"source_frozen_donor_groups":2108,
       "source_species_scenario_radius_records":len(species),
       "source_group_scenario_radius_records":len(groups),
       "by_scenario_radius":output,
       "limits":[
          "The primary v08 host-identity patch split finding remains negative or no positive advantage, not an equivalent-effect statistical proof.",
          "Many controls have zero common ancient/modern source edges and cannot identify any comparison; they remain unevaluable, not negative.",
          "Many nonempty exact-paired comparisons have zero contrast due discrete 0.5deg patch occupancy and small denominators.",
          "The 0BP same-patch gap and difference-in-differences share an algebraic term: their negative correlation is not a causal regression-to-mean estimate.",
          "Posthoc ≥5 matched edges / ≥4 donor group subset is a stability description, never retune of source cohort.",
          "Genus-level bootstrap cannot represent dependence from common host plant models, reused other plants or ancestral phylogeny.",
          "These are MODERN sample anchor counterfactuals on 21ka climate and PMIP4 landscapes, not observed ancestral host occupancy or insect genotypes.",
          "Do not reinterpret the original TTF v0.5 true-host versus matched-other-real-plant negative result as supported."]}
    return per_group,report

def main():
    a=argparse.ArgumentParser()
    for key in ("species","groups","summary","out"):
        a.add_argument("--"+key,type=Path,required=True)
    args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    _,report=validate_and_audit(args.species,args.groups,args.summary)
    (args.out/"identifiability_summary_v09.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:report["by_scenario_radius"][k] for k in
         ("analog_100km","extrapolated_100km")},sort_keys=True))
if __name__=="__main__":main()
