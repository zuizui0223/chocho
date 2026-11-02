#!/usr/bin/env python3
"""Assemble a source-audited three-gate island butterfly case/trait evidence matrix.

A first indexed record is NOT introduction/establishment year; a present-day port
count is NOT historical ship exposure; an annual temperature anomaly is NOT a
species climate suitability crossing. This is a feasibility and record-order
integration, not a three-way attribution or a wing×port colonization model.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from collections import defaultdict,Counter
from pathlib import Path

def rows(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def csvwrite(path,values):
    if not values:raise RuntimeError("No case rows")
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(values[0]))
        writer.writeheader();writer.writerows(values)
def numeric(x):
    try:
        if x is None or not str(x).strip():return None
        val=float(x)
        if not math.isfinite(val):return None
        return val
    except (ValueError,TypeError):return None
def best_wingspan(rows):
    cols=("WS_L","WS_U","WS_L_Fem","WS_U_Fem","WS_L_Mal","WS_U_Mal")
    out=[]
    for rr in rows:
        low=numeric(rr.get("WS_L"))
        high=numeric(rr.get("WS_U"))
        if low is not None and high is not None and low>0 and high>=low:
            out.append((low+high)/2)
        elif low is not None and low>0:out.append(low)
        elif high is not None and high>0:out.append(high)
    if out:return float(statistics.median(out)),"LEPTR_WS_L_U"
    for col in cols[2:]:
        numbers=[numeric(rr.get(col)) for rr in rows]
        valid=[n for n in numbers if n is not None and n>0]
        if valid:return float(statistics.median(valid)),"LEPTR_SEX_SPECIFIC"
    return None,"MISSING"
def main():
    p=argparse.ArgumentParser()
    for k in ("literature-cases","q2-edges-v02","q2-audit-v02","q4-firstyear-pilot",
              "nasa-power-case-climate","port-counts-v02","leptraits","geb-species-metrics",
              "out-dir"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    literature=rows(a.literature_cases)
    edges=rows(a.q2_edges_v02)
    q2=json.loads(a.q2_audit_v02.read_text(encoding="utf8"))
    q4=rows(a.q4_firstyear_pilot)
    nasa=rows(a.nasa_power_case_climate)
    ports=rows(a.port_counts_v02)
    leptraits=rows(a.leptraits)
    geb=rows(a.geb_species_metrics)
    if len(literature)!=47 or len(ports)!=369 or len(q4)!=500 or len(geb)!=239:
        raise RuntimeError("Q1–Q4 and 47-case input population drift")
    if len(edges)!=q2["scope"]["introduced_butterfly_plant_island_edges"]:
        raise RuntimeError("Q2 corrected edge source drift")
    if q2["scope"]["added_strict_island_cells"]!=5741:
        raise RuntimeError("Requires corrected 126-island CNY host geography")
    source=defaultdict(set)
    for r in edges:source[(r["species"],r["island_code"])].add(r["plant_id"])
    if len(source)!=5741 or len(set(r["island_code"] for r in edges))!=123:
        raise RuntimeError("Q2 corrected unique butterfly–island resource membership drift")
    byplant={(r["plant_id"],r["island_code"]):r for r in q4}
    bynasa={(r["species"],r["physical_island"],r["wgsrpd3_code"]):r for r in nasa}
    if len(bynasa)!=47 or len(nasa)!=47:
        raise RuntimeError("POWER case study record keys changed")
    byport={r["area_code_l3"]:r for r in ports}
    if sum(x["geographic_group"]=="island" for x in ports)!=126 or byport["CNY"]["geographic_group"]!="island":
        raise RuntimeError("WPI 126-island geography not yet corrected")
    traitby=defaultdict(list)
    for r in leptraits:traitby[r.get("Species","").strip()].append(r)
    metrics={r["species"]:r for r in geb}
    details=[]
    for entry in literature:
        sp=entry["species"];code=entry["wgsrpd3_code"]
        linked=source.get((sp,code),set())
        plantdata=[byplant[(h,code)] for h in linked if (h,code) in byplant]
        documented=[p for p in plantdata if p["q4_pilot_status"]=="DATED_OCCURRENCE_FOUND"]
        first=int(entry["first_evidence_earliest_year"])
        firstyears=[int(p["earliest_gbif_indexed_year"]) for p in documented]
        before=[v for v in firstyears if v<first]
        nasa_key=(sp,entry["physical_island"],code)
        cc=bynasa[nasa_key]
        pp=byport.get(code)
        proxy,origin=best_wingspan(traitby.get(sp,[]))
        # Never combine a prescaled GEB proxy with raw LepTraits wingspans.
        geb_wing_proxy=numeric(metrics[sp].get("wing_size_proxy")) if sp in metrics else None
        diapause=[r.get("DiapauseStage") for r in traitby.get(sp,[]) if r.get("DiapauseStage")]
        ovip=[r.get("OvipositionStyle") for r in traitby.get(sp,[]) if r.get("OvipositionStyle")]
        port_n=int(pp["WPI_port_count"]) if pp else None
        climate_10yr=numeric(cc.get("T2M_mean_decade_before_first_record_C"))
        details.append({
          "species":sp,
          "physical_island":entry["physical_island"],
          "wgsrpd3_code":code,
          "literature_first_record_min_year":first,
          "literature_first_record_max_year":int(entry["first_evidence_latest_year"]),
          "documented_establishment_by_year":entry["established_by_year"],
          "literature_first_record_grade":entry["evidence_grade"],
          "species_in_fixed_GEB239":int(sp in metrics),
          "modern_introduced_host_new_butterfly_resource_cell":int(bool(linked)),
          "modern_introduced_host_candidate_count":len(linked),
          "hosts_in_response_blind_500_firstyear_pilot":len(plantdata),
          "hosts_with_any_earliest_GBIF_plant_index":len(documented),
          "hosts_with_GBIF_index_before_butterfly_literature":len(before),
          "earliest_indexed_candidate_plant_year":min(firstyears) if firstyears else "",
          "earliest_preliterature_indexed_candidate_host_year":min(before) if before else "",
          "record_order_indexed_host_years_prior_to_butterfly":first-min(before) if before else "",
          "nasa_1981_2025_status":cc["climate_coverage"],
          "nasa_decade_prior_mean_T2M_C":climate_10yr if climate_10yr is not None else "",
          "nasa_decade_minus_1981_2000_T2M_C":cc["T2M_decade_before_minus_1981_2000_C"],
          "nasa_decadal_T2M_trend_C_per_decade":cc["T2M_decadal_trend_at_first_record_C_per_decade"],
          "WPI_indexed_port_count_2023":port_n if port_n is not None else "",
          "leptraits_wingspan_source_units":proxy if proxy is not None else "",
          "leptraits_wingspan_source_column":origin,
          "geb_prescaled_wing_size_proxy":geb_wing_proxy if geb_wing_proxy is not None else "",
          "diapause_stage_LepTraits_raw":diapause[0] if diapause else "",
          "oviposition_style_LepTraits_raw":ovip[0] if ovip else "",
          "full_temporal_proxy_context":
             int(bool(linked) and bool(before) and climate_10yr is not None
                 and port_n is not None and proxy is not None),
          "last_satisfied_gate_identified":0,
          "ship_or_flight_pathway_identified":0,
        })
    if len(details)!=47:raise RuntimeError("Case evidence alignment did not preserve source")
    joint={
        "schema":"chocho_3gate_literature_host_climate_port_wing_integrated_v0.1",
        "status":"SAME_WGSRPD3_MACROREGION_PROXY_SUPPORT_NOT_CAUSAL_GATE_IDENTIFICATION",
        "source_cohorts":{
            "literature_leads":47,
            "GEB_species":239,
            "Q2_geo_corrected_resource_cells":len(source),
            "Q2_source_host_edges":len(edges),
            "first_year_pilot_pairs":len(q4),
            "NASA_case_climate_rows":len(nasa),
            "WPI_regions":len(ports)
        },
        "observed":{
            "literature_rows_in_GEB":sum(r["species_in_fixed_GEB239"] for r in details),
            "literature_rows_at_Q2_new_resource_cell":sum(r["modern_introduced_host_new_butterfly_resource_cell"] for r in details),
            "literature_rows_with_dated_host_in_old_fixed_pilot":sum(bool(r["hosts_with_any_earliest_GBIF_plant_index"]) for r in details),
            "literature_rows_with_host_index_prior_to_butterfly_record":sum(bool(r["hosts_with_GBIF_index_before_butterfly_literature"]) for r in details),
            "literature_rows_with_10yr_prefirst_T2M":sum(bool(r["nasa_decade_prior_mean_T2M_C"]!="") for r in details),
            "literature_rows_with_WPI_port_count":sum(r["WPI_indexed_port_count_2023"]!="" for r in details),
            "literature_rows_with_wing_measure":sum(r["leptraits_wingspan_source_units"]!="" for r in details),
            "literature_rows_with_all_recorded_proxies":sum(r["full_temporal_proxy_context"] for r in details),
            "distinct_species_with_all_recorded_proxies":len({r["species"] for r in details if r["full_temporal_proxy_context"]}),
        },
        "critical_limits":[
           "The 47 literature rows are preselected known detections and cannot form a risk-set logistic regression of arrival/establishment with meaningful non-detection controls.",
           "Most physical island cases reside within WGSRPD3 archipelagos; a host record in one member and butterfly first report in another is not temporal evidence for a single island.",
           "A first record is not introduction/colonization time; WCVP current alien-status is not a dated plant arrival; NASA warming is not a climate suitability threshold.",
           "2023 WPI modern indexed port count is not shipping frequency at historical year; big wings do not unambiguously imply long-distance flight.",
           "Only four empire histories were published in Lenzner; no empire-origin mechanism fitted here.",
           "This table explicitly leaves all three bottleneck gate labels unidentified. It provides source-provenance coverage for future causal comparison only.",
        ],
    }
    csvwrite(a.out_dir/"island_species_three_gate_case_matrix.csv",details)
    (a.out_dir/"three_gate_case_integration_audit.json").write_text(
         json.dumps(joint,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(joint,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
