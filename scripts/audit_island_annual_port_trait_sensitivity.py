#!/usr/bin/env python3
"""Offline sensitivity of source-locked WCVP×GIFT global island trait filtering.

Intentionally no external API/GIFT calls; source artifact from run 41946522894.
Do not promote nominal p from six exploratory tests to a verified effect.
"""
from __future__ import annotations
import argparse
import collections
import json
import statistics
from pathlib import Path
import numpy as np
import pandas as pd
from exploratory_island_annual_trait_selection import port_perm

def b_h(p):
    vals=np.asarray(p,float)
    n=len(vals)
    order=np.argsort(vals)
    out=np.zeros(n)
    m=1.0
    for rank in range(n,0,-1):
        i=int(order[rank-1])
        m=min(m,vals[i]*n/rank)
        out[i]=m
    return list(map(float,out))

def eval_subset(frame,*,min_cov,min_species,exclude_level1=None,seed=901):
    d=frame.copy()
    if exclude_level1 is not None:
        d=d[d.level1!=str(exclude_level1)].copy()
    i=d[
      (d.introduced_lifeform_both_known>=min_species)&
      (d.native_lifeform_both_known>=min_species)&
      (d.introduced_both_trait_coverage>=min_cov)&
      (d.native_both_trait_coverage>=min_cov)&
      d.GBIF_allplant_indexed.notna()
    ].copy()
    i["intro_minus_native"]=i.introduced_annual_herb_fraction_both_known-i.native_annual_herb_fraction_both_known
    controls=["introduced_both_trait_coverage","native_both_trait_coverage"]
    return {
       "cohort_n":len(i),"biogeographic_regions":int(i.level1.nunique()),
       "introduced_annual_herb":port_perm(i,"introduced_annual_herb_fraction_both_known",controls,B=999,seed=seed),
       "native_annual_herb":port_perm(i,"native_annual_herb_fraction_both_known",controls,B=999,seed=seed+101),
       "intro_minus_native":port_perm(i,"intro_minus_native",controls,B=999,seed=seed+202),
       "region_codes":sorted(i.code.to_list())
    }

def run(original):
    if original["schema"]!="chocho_WCVP_v13_GIFT3p2_lifecycle_annual_herb_island_port_trait_selection_v1":
        raise ValueError("Incorrect frozen WCVP×GIFT source")
    if (len(original["WCVP_island_trait_status"]),original["model_results"]["n_regions"],
        original["model_results"]["source_regions_with_botanical_effort"])!=(126,126,98):
        raise ValueError("Fixed 126 islands / 98 botanical-effort islands source changed")
    base=original["model_results"]["models"]
    assert (base["lifecycle"]["cohort_n"],base["both"]["cohort_n"])==(90,87)
    rows=pd.DataFrame(original["WCVP_island_trait_status"])
    for name,cols in (("introduced",["introduced_N","introduced_lifecycle_known","introduced_lifeform_both_known"]),
                      ("native",["native_N","native_lifecycle_known","native_lifeform_both_known"])):
        for c in cols:
            if (rows[c]<0).any():raise ValueError("Negative source counts")
        if (rows[cols[1]]>rows[cols[0]]).any() or (rows[cols[2]]>rows[cols[1]]).any():
            raise ValueError("Trait coverage counted outside source WCVP plants")
    measures={}
    ps=[];fields=[]
    for x in ("lifecycle","both"):
        m=base[x]
        for key in ("introduced_trait_selection","native_geographic_comparator","introduced_minus_native_selector"):
            fields.append(x+"__"+key)
            ps.append(m[key]["two_sided_descriptive_permutation_p"])
    qs=b_h(ps)
    adj={field:{"nominal_p":p,"BH_q_exploratory_six_tests":q}
         for field,p,q in zip(fields,ps,qs)}
    d=rows[rows.GBIF_allplant_indexed.notna()].copy()
    for origin in ("introduced","native"):
        for item in ("lifecycle_coverage","both_trait_coverage"):
            col=origin+"_"+item
            df=d[d[col].notna()]
            measures[col]={
              "n_regions":len(df),
              "mean":float(df[col].mean()),
              "median":float(df[col].median()),
              "q10":float(df[col].quantile(.10)),
              "q90":float(df[col].quantile(.90)),
              "rho_vs_ports":float(df[col].rank().corr(df.log2ports.rank())),
              "rho_vs_botanical_effort":float(df[col].rank().corr(df.GBIF_allplant_indexed.rank()))
            }
    analyses={}
    for label,cov,n in (("base_replay_20pct_min20",.20,20),
                        ("coverage30pct_min20",.30,20),
                        ("coverage40pct_min20",.40,20),
                        ("coverage20pct_min40",.20,40),
                        ("coverage30pct_min40",.30,40),
                        ("coverage50pct_min30",.50,30)):
        analyses[label]=eval_subset(rows,min_cov=cov,min_species=n,seed=20261105+len(label))
        print("GIFT_GLOBAL_TRAIT_COVERAGE_SENSITIVITY",label,
              json.dumps({k:v for k,v in analyses[label].items() if k!="region_codes"},
                         ensure_ascii=False)[:7000],flush=True)
    level1=sorted(rows.level1.unique())
    leave={}
    for lev in level1:
        z=eval_subset(rows,min_cov=.20,min_species=20,exclude_level1=lev,seed=20261211+int(lev))
        leave[lev]={"n":z["cohort_n"],"region_count":z["biogeographic_regions"],
                    "intro_beta":z["introduced_annual_herb"].get("observed_beta_per_log2_indexed_port"),
                    "native_beta":z["native_annual_herb"].get("observed_beta_per_log2_indexed_port"),
                    "relative_beta":z["intro_minus_native"].get("observed_beta_per_log2_indexed_port"),
                    "relative_permutation_p":z["intro_minus_native"].get("two_sided_descriptive_permutation_p")}
    report={
       "schema":"chocho_GIFT3p2_annual_herb_port_selection_WCVP_global_robustness_v1",
       "status":"EXPLORATORY_COMPOSITION_DIAGNOSTICS_NOT_CAUSAL_PORT_N_SELECTION",
       "source_run":41946522894,
       "biological_interpretation":"Port count not positively associated with absolute annual or annual-herb fraction of introduced flora in primary 90/87 regions. Relative introduced-native difference, if present, needs examination of native comparator and multiple comparisons.",
       "raw_six_tests_BH_q":adj,
       "allplant_trait_coverage_vs_port_reporting":measures,
       "coverage_minN_sensitivity":analyses,
       "leave_one_Level1_out":leave,
       "limitations":[
         "Each island is WGSRPD3 botanical region, not individual port or physical island.",
         "GIFT life-cycle and growth-form coverage not random by vegetation family, biogeography or popularity.",
         "Annual, herb and annual herb are NOT direct nitrophily, ruderal habitat fidelity or nutrient use.",
         "Native trait share is a comparator, not true counterfactual introduced taxon mix.",
         "Selection on trait-known plants and on GBIF effort may create collider bias; all estimates observational.",
         "Multiple exploratory contrasts (six primary), further sensitivity diagnostics; BH q descriptive not confirmatory.",
         "Modern 2023 WPI count cannot identify historic introduction or causal nutrient/maritime mechanisms."
       ]
    }
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    x=run(json.loads(a.source.read_text(encoding="utf-8")))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("GIFT_GLOBAL_PORT_ANNUAL_TRAIT_FINAL_AUDIT",
          json.dumps({"BH_q":x["raw_six_tests_BH_q"],
                      "coverage":x["allplant_trait_coverage_vs_port_reporting"],
                      "leave_level1_out":x["leave_one_Level1_out"]},ensure_ascii=False)[:19000],flush=True)

if __name__=="__main__":main()
