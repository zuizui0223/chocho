#!/usr/bin/env python3
"""Q4 source-locked 31-butterfly island observation vs dated introduced hosts.

Answerable *documentation* question: among 920 butterfly x WGSRPD3 island
resource opportunities, how many have dated GBIF-indexed introduced hosts,
and do any indexed butterfly records occur after a plant record?

Source-respecting inference:
- original 500 host x island cases selected before response, plus three
  prespecified Canary cases; NOT a complete botanical cohort;
- separate butterfly-case-selected 84 pairs as explicitly NONCONFIRMATORY;
- 31 butterflies from a previously sampled DOI, not full 239/224 source;
- absence of GBIF indexing is NOT species absence or unexploited opportunity;
- temporal index order is NOT biological host-arrival / butterfly colonization.
"""
from __future__ import annotations

import argparse, csv, json, statistics
from collections import Counter, defaultdict
from pathlib import Path

def read_csv(p):
    with p.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def write_csv(p, values):
    if not values: raise RuntimeError("No data for "+str(p))
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer=csv.DictWriter(fh, fieldnames=list(values[0]))
        writer.writeheader();writer.writerows(values)

def idxyear(s):
    if not str(s).strip():return None
    y=int(str(s))
    if not (1700<=y<=2026):raise ValueError("Invalid GBIF year "+str(s))
    return y

def main():
    ap=argparse.ArgumentParser()
    for name in ("q2-edges","q2-audit","q3-risk-cells","q3-recorded",
                 "pilot-plant-years","canary-plant-years","targeted-plant-years",
                 "out-dir"):
        ap.add_argument("--"+name, type=Path, required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q2=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    edges=read_csv(a.q2_edges)
    risk=read_csv(a.q3_risk_cells)
    positive=read_csv(a.q3_recorded)
    pilot=read_csv(a.pilot_plant_years)
    canary=read_csv(a.canary_plant_years)
    targeted=read_csv(a.targeted_plant_years)
    if (q2["scope"]["added_strict_island_cells"],
        q2["scope"]["introduced_butterfly_plant_island_edges"],
        len(edges),len(risk),len(positive),len(pilot),len(canary),len(targeted)
        )!=(5741,11573,11573,920,13,500,3,84):
        raise RuntimeError("Source panel drift; never impute unknown Q2/Q3/Q4 data")
    plant_by_cell=defaultdict(set)
    for edge in edges:
        plant_by_cell[(edge["species"],edge["island_code"])].add(edge["plant_id"])
    if len(plant_by_cell)!=5741:raise RuntimeError("Q2 corrected 5741 cells drift")
    observed={(r["species"],r["area_code_l3"]):r for r in positive}
    if len(observed)!=13:raise RuntimeError("Duplicate or missing recorded cases")
    panel={(r["species"],r["wgsrpd3_code"]) for r in risk}
    if len(panel)!=920 or not observed.keys()<=panel:
        raise RuntimeError("Source risk cells not unique or Q3 positives outside population")

    source=defaultdict(dict)
    for r in pilot:
        key=(r["plant_id"],r["island_code"])
        if key in source["RESP_BLIND_PILOT"]:
            raise RuntimeError("Duplicate blind pilot plant-island key")
        source["RESP_BLIND_PILOT"][key]={
          "status":r["q4_pilot_status"],
          "year":idxyear(r["earliest_gbif_indexed_year"])
          if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND" else None}
    for r in canary:
        key=(r["plant_id"],r["island_code"])
        if key in source["RESP_BLIND_CANARY"]:
            raise RuntimeError("Duplicate prespecified Canary plant")
        source["RESP_BLIND_CANARY"][key]={
           "status":r["GBIF_query_status"],
           "year":idxyear(r["earliest_gbif_plant_indexed_year"])
           if r["GBIF_query_status"]=="DATED_OCCURRENCE_FOUND" else None}
    source["RESP_BLIND_PILOT_PLUS_CANARY"]={**source["RESP_BLIND_PILOT"],**source["RESP_BLIND_CANARY"]}
    if len(source["RESP_BLIND_PILOT_PLUS_CANARY"])!=503:
        raise RuntimeError("Unexpected duplicate in two prespecified cohorts")
    for r in targeted:
        key=(r["plant_id"],r["wgsrpd3_code"])
        source["CASE_TARGETED_ENRICHMENT"][key]={
            "status":r["status"],
            "year":idxyear(r["first_gbif_indexed_plant_year"])
            if r["status"]=="DATED_OCCURRENCE_FOUND" else None}
    if len(source["CASE_TARGETED_ENRICHMENT"])!=84:
        raise RuntimeError("Targeted host chronology source drift")
    source["EXPLORATORY_POOLED"]={**source["RESP_BLIND_PILOT_PLUS_CANARY"],
                                    **source["CASE_TARGETED_ENRICHMENT"]}

    outcome_rows=[];cohorts={}
    for cohort_name in ("RESP_BLIND_PILOT","RESP_BLIND_PILOT_PLUS_CANARY",
                        "CASE_TARGETED_ENRICHMENT","EXPLORATORY_POOLED"):
        seen=source[cohort_name]
        n_plants=0;dated=0;pre2001=0;pre2010=0;pre2018=0;obs_covered=0
        obs_year_order=0;obs_preyear=0;obs_nocov=0
        years_order=[]
        status=Counter()
        regional=defaultdict(lambda:{"risk":0,"obs":0,"pre2001":0})
        for r in risk:
            cell=(r["species"],r["wgsrpd3_code"])
            if cell not in plant_by_cell:raise RuntimeError("Q3 cell outside audited Q2")
            candidates=plant_by_cell[cell]
            measured=[seen[(pid,cell[1])] for pid in candidates if (pid,cell[1]) in seen]
            dated_years=[m["year"] for m in measured if m["status"]=="DATED_OCCURRENCE_FOUND"]
            if any(y is None for y in dated_years):raise RuntimeError("Incomplete dated plant year")
            earliest=min(dated_years) if dated_years else None
            earliest_butterfly=idxyear(observed[cell]["first_indexed_year_in_DOI_subset"]) if cell in observed else None
            earlier=[y for y in dated_years if earliest_butterfly is not None and y<earliest_butterfly]
            n_plants+=int(bool(measured))
            dated+=int(bool(earliest is not None))
            pre2001+=int(any(y<=2000 for y in dated_years))
            pre2010+=int(any(y<=2009 for y in dated_years))
            pre2018+=int(any(y<=2017 for y in dated_years))
            obs_covered+=int(bool(measured) and cell in observed)
            obs_year_order+=int(bool(earlier))
            obs_preyear+=int(cell in observed and bool(earliest is not None))
            obs_nocov+=int(cell in observed and not measured)
            if earlier:years_order.append(earliest_butterfly-min(earlier))
            for item in measured:status[item["status"]]+=1
            regional[r["level1"]]["risk"]+=1
            regional[r["level1"]]["obs"]+=int(cell in observed)
            regional[r["level1"]]["pre2001"]+=int(any(y<=2000 for y in dated_years))
            if cohort_name=="RESP_BLIND_PILOT_PLUS_CANARY":
                outcome_rows.append({
                  "species":cell[0],
                  "island_code":cell[1],
                  "level1":r["level1"],
                  "butterfly_doi_recorded":int(cell in observed),
                  "q2_introduced_hosts_in_cell":len(candidates),
                  "indexed_pilot_hosts_in_cell":len(measured),
                  "indexed_pilot_hosts_with_gbif_date":len(dated_years),
                  "any_plant_indexed_pre2001":int(any(y<=2000 for y in dated_years)),
                  "any_plant_indexed_pre2010":int(any(y<=2009 for y in dated_years)),
                  "earliest_pilot_plant_indexed_year":earliest if earliest is not None else "",
                  "first_butterfly_doi_indexed_year":earliest_butterfly if earliest_butterfly else "",
                  "indexed_plant_earlier_than_butterfly_index":
                     int(bool(earlier)) if cell in observed else "",
                  "indexed_record_order_lag_years_max":earliest_butterfly-min(earlier) if earlier else "",
                  "WPI_ports_2023":r["indexed_WPI_port_count_2023"],
                  "wing_proxy":r["wing_size_proxy_original"],
                  "gbif_other_butterfly_indexed_records":r["other_30_butterfly_record_count"],
                  "IMPORTANT_not_observed_is_not_absent":1,
                })
        cohorts[cohort_name]={
           "source_plant_island_pairs":len(seen),
           "butterfly_island_opportunity_cells":len(panel),
           "cells_with_one_or_more_queried_host":n_plants,
           "cells_with_one_or_more_dated_host":dated,
           "cells_with_host_index_pre2001":pre2001,
           "cells_with_host_index_pre2010":pre2010,
           "cells_with_host_index_pre2018":pre2018,
           "butterfly_indexed_positive_cells":13,
           "positive_cells_with_a_queried_host":obs_covered,
           "positive_cells_with_a_dated_host":obs_preyear,
           "positive_cells_with_dated_host_index_earlier_than_butterfly_index":obs_year_order,
           "positive_cells_with_no_queried_host":obs_nocov,
           "record_order_difference_min_years":min(years_order) if years_order else None,
           "record_order_difference_median_years":statistics.median(years_order) if years_order else None,
           "record_order_difference_max_years":max(years_order) if years_order else None,
           "year_cohort_by_level1":dict(regional),
           "sample": "selected before butterfly case outcome" if "BLIND" in cohort_name else
                "selected based on literature butterfly cases; NOT a confirmatory risk cohort",
        }
    report={
       "schema":"chocho_q4_dated_host_q3_independent_butterfly_doi_screen_v0.1",
       "status":"DOCUMENTATION_ORDER_NOT_BIOLOGICAL_HOST_ARRIVAL_OR_ESTABLISHMENT",
       "source":{
         "q2_corrected_host_edges":11573,"q2_corrected_butterfly_island_cells":5741,
         "q3_31species_butterfly_DOI":"10.15468/dl.pp5nc9",
         "q3_31_species_added_resource_opportunities":920,
         "q3_species_island_recorded_cells":13,
         "q4_response_blind_pilot_pairs":500,
         "q4_prespecified_canary_pairs":3,
         "q4_case_targeted_plant_pairs":84,
       },
       "cohorts":cohorts,
       "limits":[
          "Even among 920 opportunity cells, 13 observed butterfly cases are not biological colonization successes and the other 907 cells are not established absences.",
          "A randomization test is not interpretable if the dated-plant risk cohort has too few butterfly indexed positives; record support and missingness before model fitting.",
          "First GBIF indexed host dates are not plant import or naturalization years, and selected host is not verified as actually used by larvae on the same island.",
          "A WGSRPD3 archipelago may contain separate plant and butterfly occurrences on different physical islands.",
          "All 31 butterflies are a previously selected, incomplete DOI sample, not 224 or 239 independent taxa.",
          "The case-targeted extra plant-year search is conditioned on documented butterfly detections and must not be pooled into a confirmatory effect estimate.",
          "Record order is not a butterfly host-lag or a demographic time-to-half-occupancy estimator.",
       ],
    }
    if len(outcome_rows)!=920:
        raise RuntimeError("Incorrect 920 fixed opportunity cell rows")
    write_csv(a.out_dir/"dated_host_x_archived_butterfly_31species.csv",outcome_rows)
    (a.out_dir/"q4_hostyear_vs_archived_butterfly_audit.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
