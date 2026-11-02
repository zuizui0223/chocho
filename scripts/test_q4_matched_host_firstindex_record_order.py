#!/usr/bin/env python3
"""Q4 descriptive 13 matched-source DOI cases: is a dated suitable plant
documented before the recorded butterfly's first DOI year?

This is an index chronology comparison, NOT a colonization, causal port effect,
true butterfly absence, or ecological time-to-arrival estimate.
"""
from __future__ import annotations

import argparse,csv,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

REPS=9999
SEED=20261102
GOOD=("DATED_OCCURRENCE_FOUND","NO_GBIF_INDEXED_PLANT_RECORD_IN_COHORT")

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def main():
    p=argparse.ArgumentParser()
    for k in ("source-cells","plant-years","q2-edges","out-dir"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    frame=rows(a.source_cells)
    plants=rows(a.plant_years)
    edges=rows(a.q2_edges)
    if len(frame)!=52 or len(plants)!=100 or len(edges)!=11573:
        raise RuntimeError("Fixed Q2/Q4 matched source cohorts drift")
    bycell=defaultdict(set)
    for r in edges:bycell[(r["species"],r["island_code"])].add(r["plant_id"])
    plant={(r["plant_id"],r["island_code"]):r for r in plants}
    if len(plant)!=100:raise RuntimeError("Duplicate plant query")
    classified=[]
    byblock=defaultdict(list)
    for r in frame:
        sp,code=r["species"],r["island_code"]
        key=(sp,code)
        docs=[plant[(pid,code)] for pid in bycell[key]]
        cut=int(r["case_reference_first_indexed_butterfly_year"])
        prior=[int(x["first_gbif_indexed_plant_year"]) for x in docs if
               x["source_date_status"]=="DATED_OCCURRENCE_FOUND" and int(x["first_gbif_indexed_plant_year"])<cut]
        unidentified=[x for x in docs if x["source_date_status"] not in GOOD]
        has_prior=int(bool(prior))
        status="DOCUMENTED_PRIOR" if has_prior else ("QUERY_INCOMPLETE_UNKNOWN" if unidentified else "NO_PRIOR_INDEXED_IN_SOURCE")
        case=int(r["recorded"])
        row={
          "species":sp,"island_code":code,"index_case":r["index_case"],
          "butterfly_indexed_in_original_DOI":case,
          "matched_control_basis":r["matched_block"],
          "cutoff_first_butterfly_GBIF_indexed_year":cut,
          "n_hosts_from_corrected_Q2":len(docs),
          "n_GBIF_indexed_host_dates":sum(x["source_date_status"]=="DATED_OCCURRENCE_FOUND" for x in docs),
          "n_GBIF_query_unresolved":len(unidentified),
          "GBIF_prior_plant_documentation":has_prior,
          "prior_documentation_classification":status,
          "earliest_indexed_prior_plant_year":min(prior) if prior else "",
          "positive_record_order_gap_max_years":cut-min(prior) if prior else "",
        }
        classified.append(row)
        byblock[(sp,r["index_case"])].append(row)
    if len(byblock)!=13:raise RuntimeError("Wrong 13 matched index case blocks")
    for block,rr in byblock.items():
        if len(rr)!=4 or sum(x["butterfly_indexed_in_original_DOI"] for x in rr)!=1:
            raise RuntimeError("Not 1 index + 3 nonindex in "+str(block))
    def contrast(blocks,maskunknown=False):
        a=[]
        for rr in blocks:
            if maskunknown and any(x["prior_documentation_classification"]=="QUERY_INCOMPLETE_UNKNOWN" for x in rr):
                continue
            case=[r for r in rr if r["butterfly_indexed_in_original_DOI"]==1][0]
            others=[r for r in rr if r["butterfly_indexed_in_original_DOI"]==0]
            score=case["GBIF_prior_plant_documentation"]-statistics.mean(x["GBIF_prior_plant_documentation"] for x in others)
            a.append(score)
        return a
    listed=list(byblock.values())
    selected=[rr for rr in listed if not any(x["prior_documentation_classification"]=="QUERY_INCOMPLETE_UNKNOWN" for x in rr)]
    values=contrast(selected)
    observed=float(statistics.mean(values)) if values else None
    rng=random.Random(SEED)
    null=[]
    if len(selected)>=4:
        for _ in range(REPS):
            y=[]
            for block in selected:
                tmp=[r["GBIF_prior_plant_documentation"] for r in block]
                chosen=rng.randrange(4)
                y.append(tmp[chosen]-sum(tmp[i] for i in range(4) if i!=chosen)/3.)
            null.append(statistics.mean(y))
    if null:
        mu=statistics.mean(null)
        test={
          "status":"DESCRIPTIVE_WITHIN_MATCHED_SET_INDEX_LABEL_PERMUTATION",
          "observed_mean_case_minus_controls":observed,
          "null_mean":mu,
          "null_95pct":[sorted(null)[int(.025*len(null))],sorted(null)[int(.975*len(null))]],
          "permutations":len(null),
          "p_two_sided_centered":
            (1+sum(abs(x-mu)>=abs(observed-mu)-1e-12 for x in null))/(len(null)+1),
          "not_independent_island_or_species_replicates":True
        }
    else:test={"status":"TOO_FEW_COMPLETE_MATCHED_SETS","n_complete":len(selected)}
    groups={role:[r for r in classified if r["butterfly_indexed_in_original_DOI"]==role] for role in (1,0)}
    summary={
      "schema":"chocho_Q4_31_butterfly_DOI_matched_source_first_host_doc_v0.1",
      "status":"DESCRIPTIVE_DOCUMENTATION_ASSOCIATION_NOT_COLONIZATION_OR_HALF_LAG",
      "n_observed_butterfly_index_case_cells":13,
      "n_nonindexed_DOI_control_cells_NOT_absent":39,
      "butterfly_species_in_cases":len({r["species"] for r in groups[1]}),
      "matched_groups":len(byblock),
      "source_GBIF_plant_firstyear_pairs":len(plants),
      "group_outcomes":{
        role: {
         "cells":len(group),
         "documented_prior_plant_index":sum(x["GBIF_prior_plant_documentation"] for x in group),
         "query_incomplete_with_no_prior_date":sum(x["prior_documentation_classification"]=="QUERY_INCOMPLETE_UNKNOWN" for x in group),
         "no_prior_date_with_all_queries_resolved":sum(x["prior_documentation_classification"]=="NO_PRIOR_INDEXED_IN_SOURCE" for x in group),
         "raw_documented_prior_fraction":sum(x["GBIF_prior_plant_documentation"] for x in group)/len(group),
        } for role,group in groups.items()
      },
      "blocked":{
        "n_complete_source_comparison_groups":len(selected),
        "n_incomplete_source_comparison_groups":len(listed)-len(selected),
        "test_on_complete_blocks":test,
        "raw_all_13_case_minus_39_control_documented_fraction":
          statistics.mean(x["GBIF_prior_plant_documentation"] for x in groups[1])-
          statistics.mean(x["GBIF_prior_plant_documentation"] for x in groups[0]),
      },
      "method_limits":[
        "Positive outcome is any original GBIF DOI observation of a butterfly and not introduced-butterfly arrival, population establishment or actual feeding on imported plants.",
        "The 39 no-index cells are NOT species absences; both groups come from nonrandom 31-species validation source.",
        "Plant first indexed dates are not actual introduction or naturalization years; local botanical collection and digitization can induce apparent temporal lags.",
        "The cutoff for every 1:3 block is its indexed butterfly case's earliest GBIF date; matched controls do not have butterfly first-arrival dates.",
        "Missing GBIF plant index/query errors are unknown and are not biological host absence; complete-block analysis is itself selected.",
        "Repeated butterfly species, ports and archipelago codes make the 13 blocks dependent; index-label randomization is descriptive under an unverified exchangeability assumption.",
        "The conditional p, when computed, is not a causal or independent establishment likelihood.",
        "Modern WPI and wing size were used to match source-index contexts, not to estimate a wing x port effect.",
      ]
    }
    with (a.out_dir/"q4_matched_host_documentation_audit.json").open("w",encoding="utf-8") as f:
        json.dump(summary,f,indent=2,ensure_ascii=False);f.write("\n")
    with (a.out_dir/"matched_species_island_doc_flags.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(classified[0]))
        w.writeheader();w.writerows(classified)
    print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=="__main__":main()
