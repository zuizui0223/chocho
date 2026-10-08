#!/usr/bin/env python3
"""Bounds for adult emergence and parasitoid outcomes among collected Phoebis larvae.

This is a sensitivity analysis of rearing counts already reported by Koptur
et al. 2024 (Insects 15:123), not new field work. The input is the frozen
source-audited chocho empirical-quality benchmark.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def calculate(data: dict) -> dict:
    if data.get("schema") != "chocho_pieridae_fabaceae_empirical_quality_decision_v0.1":
        raise ValueError("unexpected source benchmark schema")
    paper=data["evidence"]["published_Phoebis_Senna_rearing"]
    assert (paper["native_site_sum_parasitoids"],
            paper["native_published_all_site_parasitoids"]) == (16, 17)
    assert (paper["native_counted_fates"],
            paper["introduced_site_sum_and_published_parasitoids"],
            paper["introduced_counted_fates"]) == (235, 49, 283)
    source={
        "native":{"counted":paper["native_counted_fates"],
                  "parasite":paper["native_site_sum_parasitoids"],
                  "unresolved":paper["unresolved_fates_native"]},
        "introduced":{"counted":paper["introduced_counted_fates"],
                      "parasite":paper["introduced_site_sum_and_published_parasitoids"],
                      "unresolved":paper["unresolved_fates_introduced"]}
    }
    out={}
    for label,row in source.items():
        n=row["counted"]+row["unresolved"]
        adult=row["counted"]-row["parasite"]
        assert n>0 and 0<=adult<=row["counted"] and row["unresolved"]>=0
        out[label]={
            "found":n,"resolved":row["counted"],"unresolved":row["unresolved"],
            "confirmed_adults":adult,"confirmed_parasitoid_emergences":row["parasite"],
            "conditional_parasitism_resolved_only":row["parasite"]/row["counted"],
            "adult_emergence_full_cohort_logical_bounds":[adult/n,(adult+row["unresolved"])/n],
            "parasitoid_emergence_full_cohort_logical_bounds":[row["parasite"]/n,(row["parasite"]+row["unresolved"])/n]
        }
    def contrast(key):
        n=out["native"][key];a=out["introduced"][key]
        return [a[0]-n[1],a[1]-n[0]]
    adult=contrast("adult_emergence_full_cohort_logical_bounds")
    parasitoid=contrast("parasitoid_emergence_full_cohort_logical_bounds")
    assert adult[0]<0<adult[1] and parasitoid[0]<0<parasitoid[1]
    assert sum(x["unresolved"] for x in out.values())==176
    return {
        "schema":"chocho_senna_rearing_fate_bounds_v0.1",
        "source":"benchmarks/exploratory/pieridae_fabaceae_empirical_quality_decision_v0.1.json",
        "source_doi":"10.3390/insects15020123",
        "table_inconsistency":{"native_site_sum":16,"native_printed_overall":17},
        "groups":out,
        "alien_minus_native_adult_emergence_risk_difference_bounds":adult,
        "alien_minus_native_parasitoid_risk_difference_bounds":parasitoid,
        "total_unresolved_fates":176,
        "sign_identified":False,
        "limitations":[
           "Upper/lower bounds assume no information about the 176 unresolved fates and are NOT estimated survival probabilities.",
           "Reported parasite contrast among resolved fates cannot be converted into a treatment causal effect.",
           "The same researchers cultivated native and nonnative Senna; this is NOT a cultivated-versus-wild comparison.",
           "Species identity, sites, host phenology and repeated observations per plant remain confounded.",
           "Study is previously published; this audit does not establish a new ecological mechanism."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-json",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    args=p.parse_args()
    x=calculate(json.loads(args.source_json.read_text(encoding="utf-8")))
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"unresolved":x["total_unresolved_fates"],
       "adult_difference_bounds":x["alien_minus_native_adult_emergence_risk_difference_bounds"],
       "parasitoid_difference_bounds":x["alien_minus_native_parasitoid_risk_difference_bounds"],
       "sign_identified":x["sign_identified"]},indent=2))
if __name__=="__main__":
    main()
