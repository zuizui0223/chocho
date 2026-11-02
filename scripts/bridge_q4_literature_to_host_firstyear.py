#!/usr/bin/env python3
"""Bridge 47 literature butterfly-island leads to Q2 novel host opportunities and Q4 GBIF first indexed plant years.

Do NOT interpret as plant arrival, island-specific larval feeding, butterfly arrival,
colonization, pathway of transport, or opportunity occupation rate. WGSRPD3 is a
region often containing multiple islands; same-region botanical/animal index overlap
is not guaranteed on the same physical island.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from collections import defaultdict,Counter

def records(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def output_csv(p, data):
    if not data:raise RuntimeError(f"empty output {p}")
    with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)

def main():
    ap=argparse.ArgumentParser()
    for k in ("seed","q2-edges","q2-audit","q4-plant-years","q4-audit","out-dir"):
        ap.add_argument("--"+k,type=Path,required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    events=records(a.seed)
    if len(events)!=47:raise RuntimeError("47 candidate literature entries changed")
    q2a=json.loads(a.q2_audit.read_text())
    edges=records(a.q2_edges)
    if (len(edges),q2a["scope"]["added_strict_island_cells"])!=(11455,5674):
        raise RuntimeError("Q2 fixed WCVP HOSTS panel changed")
    q4a=json.loads(a.q4_audit.read_text())
    indexed=records(a.q4_plant_years)
    if len(indexed)!=500 or q4a["pilot_dated_pairs"]!=219:
        raise RuntimeError("Q4 fixed historical plant-index pilot changed")
    by_cell=defaultdict(set)
    for r in edges:
        by_cell[(r["species"],r["island_code"])].add(r["plant_id"])
    dated={(r["plant_id"],r["island_code"]):r for r in indexed}
    outcomes=[]
    for event in events:
        sp=event["species"]; region=event["wgsrpd3_code"]
        first=int(event["first_evidence_earliest_year"])
        last=int(event["first_evidence_latest_year"])
        if first>last:raise RuntimeError("inverted first-record interval")
        candidates=by_cell.get((sp,region),set()) if region else set()
        reviewed=[dated[(p,region)] for p in candidates if (p,region) in dated]
        observed_years=[int(r["earliest_gbif_indexed_year"]) for r in reviewed
            if r["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND"]
        before=[y for y in observed_years if y<first]
        after=[y for y in observed_years if y>last]
        within=[y for y in observed_years if first<=y<=last]
        outcomes.append({
          **event,
          "q2_source_has_new_host_resource_cell":int(bool(candidates)),
          "q2_introduced_hosts_in_region":len(candidates),
          "q4_pilot_hosts_in_cell":len(reviewed),
          "q4_pilot_hosts_with_dated_plant_record":len(observed_years),
          "q4_pilot_hosts_with_indexed_plant_before_first_literature_record":len(before),
          "q4_pilot_hosts_with_indexed_plant_after_last_literature_record":len(after),
          "q4_pilot_hosts_with_indexed_plant_in_literature_interval":len(within),
          "earliest_indexed_plant_year_among_pilot_hosts":min(observed_years) if observed_years else "",
          "earliest_pre_butterfly_indexed_plant_year":min(before) if before else "",
          "observed_documentation_lead_years_at_most":first-min(before) if before else "",
          "first_butterfly_documentation_pre_POWER":int(first<1981),
          "same_physical_island_botanical_match_verified":0,
          "interpretation":"REGIONAL_COINDEX_NOT_ARRIVAL_OR_COLONIZATION"
        })
    if len(outcomes)!=47:raise RuntimeError("literature seed did not fully map")
    groups=defaultdict(list)
    for row in outcomes:
        if row["wgsrpd3_code"]:groups[(row["species"],row["wgsrpd3_code"])].append(row)
    cellrows=[]
    for (sp,reg),rr in sorted(groups.items()):
        cellrows.append({
          "species":sp,"wgsrpd3_code":reg,
          "physical_literature_cases_collapsed":len(rr),
          "first_documented_year_lower_bound":min(int(q["first_evidence_earliest_year"]) for q in rr),
          "first_documented_year_upper_bound":min(int(q["first_evidence_latest_year"]) for q in rr),
          "q2_added_resource_cell":max(q["q2_source_has_new_host_resource_cell"] for q in rr),
          "q2_introduced_hosts_in_region":max(q["q2_introduced_hosts_in_region"] for q in rr),
          "q4_pilot_dated_introduced_hosts":max(q["q4_pilot_hosts_with_dated_plant_record"] for q in rr),
          "q4_pilot_hosts_documented_before_butterfly":
            max(q["q4_pilot_hosts_with_indexed_plant_before_first_literature_record"] for q in rr),
        })
    stats={
      "all_literature_physical_island_or_locality_events":len(outcomes),
      "unique_species":len({r["species"] for r in outcomes}),
      "unique_mapped_species_x_WGSRPD3_cells":len(cellrows),
      "events_with_Q2_new_resource_cell":sum(r["q2_source_has_new_host_resource_cell"] for r in outcomes),
      "distinct_cells_with_Q2_new_resource_cell":sum(r["q2_added_resource_cell"] for r in cellrows),
      "events_with_any_host_in_fixed_500_pair_first_year_pilot":sum(bool(r["q4_pilot_hosts_in_cell"]) for r in outcomes),
      "events_with_dated_indexed_host_in_pilot":sum(bool(r["q4_pilot_hosts_with_dated_plant_record"]) for r in outcomes),
      "events_with_dated_pilot_host_observed_before_literature_first_record":sum(bool(r["q4_pilot_hosts_with_indexed_plant_before_first_literature_record"]) for r in outcomes),
      "cells_with_indexed_host_earlier_than_literature":sum(bool(r["q4_pilot_hosts_documented_before_butterfly"]) for r in cellrows),
      "cases_first_evidence_1981plus":sum(int(r["first_evidence_earliest_year"])>=1981 for r in outcomes),
      "cases_with_subsequent_establishment_or_breeding_year":sum(bool(r["established_by_year"]) for r in outcomes),
    }
    summary={
      "schema":"chocho_q4_literature_host_observation_bridge_v0.1",
      "status":"EXPLORATORY_COVERAGE_AND_RECORD_ORDER_ONLY",
      "source": {
        "q2_edge_count":len(edges),"q2_added_island_cells":5674,
        "q4_fixed_sample_plant_island_pairs":len(indexed),
        "q4_original_dated_pairs":q4a["pilot_dated_pairs"],
        "literature_seed":"47 curated leads with separate first record and establishment statuses"
      },
      "coverage":stats,
      "limitations":[
        "Q4 plants were originally chosen in a response-blind fixed 500-pair pilot, not chosen to enrich literature species-island events.",
        "A Q2 'added opportunity' merely indicates introduced larval host presence in the modern broad WGSRPD3 region. It does not confirm that the butterfly is invasive or absent natively.",
        "The earliest recorded plant year is biased by historical botanical sampling and is not introduction or establishment year.",
        "The first published butterfly record is not arrival and may be a preexisting native or overlooked population.",
        "WGSRPD3 regional overlap can reflect different physical islands within archipelagos.",
        "No population survival, breeding, ship arrival, butterfly detection-effort model or climate eligibility threshold is estimated.",
        "Same biological entity can contribute multiple documented physical islands, so physical events are not independent WGSRPD3 observations."
      ]
    }
    output_csv(a.out_dir/"literature_firstrecord_x_host_source_pilot.csv",outcomes)
    output_csv(a.out_dir/"literature_source_cells_collapsed.csv",cellrows)
    (a.out_dir/"literature_host_firstyear_coverage_audit.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
