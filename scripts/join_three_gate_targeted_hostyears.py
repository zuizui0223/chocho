#!/usr/bin/env python3
"""Join newly measured historical plant *index* dates to pre-frozen 47-case three-gate source matrix.

The old 500-pair host pilot and the 47 literature cases are distinguishable.
Full proxy overlap still cannot identify the last causal gate because host
introduction year, local larvae, climate suitability and actual vessel arrivals
are all unobserved or incomparable across physical islands.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
from collections import Counter

def rows(p):
    with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))
def dump(p,items):
    with p.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-three-gate-v02",type=Path,required=True)
    ap.add_argument("--targeted-host-years",type=Path,required=True)
    ap.add_argument("--targeted-host-audit",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    cases=rows(a.original_three_gate_v02)
    target=rows(a.targeted_host_years)
    audit=json.loads(a.targeted_host_audit.read_text(encoding="utf8"))
    if len(cases)!=47 or len(target)!=47 or audit["source"]["candidate_literature_events"]!=47:
        raise RuntimeError("Expected exact 47-case source panels")
    tby={(r["species"],r["physical_island"],r["wgsrpd3_code"]):r for r in target}
    if len(tby)!=47:raise RuntimeError("Targeted botanical case identity conflict")
    result=[]
    for r in cases:
        key=r["species"],r["physical_island"],r["wgsrpd3_code"]
        if key not in tby:raise RuntimeError("No linked targeted host evidence for "+str(key))
        t=tby[key]
        source_host=int(r["modern_introduced_host_candidate_count"])
        if source_host!=int(t["modern_introduced_host_candidates"]):
            raise RuntimeError("Q2 host case changed "+str(key))
        nr=int(t["host_plant_region_pairs_with_GBIF_firstyear"])
        prior=int(t["host_plant_region_pairs_documented_before_butterfly"])
        if nr<prior or nr>source_host:
            raise RuntimeError("Targeted host first-year count impossible")
        year=str(t["earliest_preliterature_candidate_host_year"])
        first=int(r["literature_first_record_min_year"])
        if prior>0 and (not year or int(year)>=first):
            raise RuntimeError("Prior host year does not precede literature "+str(key))
        climate=bool(r["nasa_decade_prior_mean_T2M_C"])
        access=bool(r["WPI_indexed_port_count_2023"])
        wings=bool(r["leptraits_wingspan_mm"])
        allproxy=bool(source_host and prior>0 and climate and access and wings)
        result.append({
            **r,
            "targeted_host_pairs_documented":nr,
            "targeted_host_pairs_with_prebutterfly_plant_record":prior,
            "targeted_earliest_indexed_compatible_host_year":
                t["earliest_documented_candidate_host_year"],
            "targeted_earliest_prebutterfly_compatible_host_year":year,
            "targeted_record_order_lead_years":first-int(year) if year else "",
            "targeted_full_documentation_proxy_support":int(allproxy),
            "same_physical_island_date_precedence_verified":0,
            "biological_last_gate_identifiable":0,
        })
    full=[r for r in result if r["targeted_full_documentation_proxy_support"]]
    prior=[r for r in result if r["targeted_host_pairs_with_prebutterfly_plant_record"]]
    groups={(r["species"],r["wgsrpd3_code"]) for r in full}
    audit_out={
      "schema":"chocho_three_gates_enriched_47_case_documented_host_years_v0.3",
      "status":"OBSERVATIONAL_DOCUMENTATION_TIMELINE_NOT_LAST_CAUSAL_GATE",
      "source":{
         "literature_case_events":47,
         "original_old_500_sample":500,
         "new_targeted_plant_island_pairs":audit["targets"]["distinct_plant_island_pairs"],
         "new_targeted_plant_pairs_with_GBIF_index_year":audit["plant_pairs_with_first_index_year"],
      },
      "coverage":{
        "candidate_literature_events":len(result),
        "with_modern_additional_resource":sum(int(r["modern_introduced_host_new_butterfly_resource_cell"]) for r in result),
        "dated_compatible_plant_first_index":sum(int(r["targeted_host_pairs_documented"])>0 for r in result),
        "prebutterfly_compatible_plant_first_index":len(prior),
        "distinct_butterfly_region_cells_with_prebutterfly_plant_index":
            len({(r["species"],r["wgsrpd3_code"]) for r in prior if r["wgsrpd3_code"]}),
        "nasa_10yr_prefirst_available":sum(bool(r["nasa_decade_prior_mean_T2M_C"]) for r in result),
        "modern_WPI_port_count_available":sum(bool(r["WPI_indexed_port_count_2023"]) for r in result),
        "LepTraits_wingspan_available":sum(bool(r["leptraits_wingspan_mm"]) for r in result),
        "all_proxy_fields_and_prebutterfly_plant_index":len(full),
        "all_proxy_distinct_species":len({r["species"] for r in full}),
        "all_proxy_distinct_butterfly_WGSRPD3_cells":len(groups),
        "biological_last_gate_identified":0,
      },
      "cases_with_all_proxy_fields":[
       {"species":r["species"],"physical_island":r["physical_island"],
        "region":r["wgsrpd3_code"],
        "first_butterfly_documentation":r["literature_first_record_min_year"],
        "earliest_indexed_preliterature_plant_year":r["targeted_earliest_prebutterfly_compatible_host_year"],
        "modern_WPI_port_count":r["WPI_indexed_port_count_2023"],
        "nasa_prefirst_T2M":r["nasa_decade_prior_mean_T2M_C"],
        "wing_mm":r["leptraits_wingspan_mm"]}
       for r in full],
      "limits":[
        "Case-selected additional GBIF plant searches are exploratory, not a preregistered causal risk-set comparison.",
        "True introduction and butterfly arrival times may precede GBIF first indexed dates by unknown amounts.",
        "The plant and butterfly are indexed in the same broad WGSRPD3 region, which can span several physical islands.",
        "No independent local botanical naturalization, larval use, butterfly establishment year or year-specific vessel arrival is verified.",
        "NASA POWER region mean anomaly does not establish species-specific climate threshold crossing.",
        "All rows start from recorded butterfly cases; there are no never-colonized sampling units and thus a wing×port colonization rate is not identifiable.",
      ],
    }
    dump(a.out_dir/"three_gate_case_matrix_targeted_hostyears.csv",result)
    (a.out_dir/"three_gate_enrichment_audit.json").write_text(json.dumps(audit_out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("THREE_GATE_TARGETED_SUMMARY",json.dumps(audit_out,ensure_ascii=False))
if __name__=="__main__":
    main()
