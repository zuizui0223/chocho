#!/usr/bin/env python3
"""Time-consistent GBIF index-risk test for historical port hosts x ocean access.

Unlike the 2001–2010-reference pilot, this analysis indexes historical
butterfly records in 1700–2000 for EVERY original 225 x 51 risk cell, and
excludes previously indexed island/species from candidate subsequent first
indices. This is still only indexing, NOT verified colonization or breeding.

A compatible host unrecorded around a port = unknown presence, NOT absence.
Original positive global resource and marine community results remain intact.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_port_host_ocean_interaction_20261114 import (
    FILES, source_join, model_sets,
)


def read_rows(path):
    with path.open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument("--history-csv",type=Path,required=True)
    p.add_argument("--history-summary",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()

    oldrisk, source_hashes, coeff, original_population = source_join(a.source_root)
    history=read_rows(a.history_csv)
    summary=json.loads(a.history_summary.read_text(encoding="utf8"))
    if len(history)!=11475 or summary["status"]!="SOURCE_COMPLETE":
        raise ValueError("Historical full-panel risk provenance is incomplete")
    if summary["original_panel_sha256"] != source_hashes["butterflies"]:
        raise ValueError("The original 225 x 51 source changed")
    histmap={(x["speciesKey"],x["landmass_id"]):x for x in history}
    if len(histmap)!=11475 or any(x["pre2001_source_complete"]!="1" for x in history):
        raise ValueError("Historical source has incomplete/duplicate 225 x 51 cells")

    butterflies=read_rows(a.source_root/FILES["butterflies"])
    host=read_rows(a.source_root/FILES["hosts"])
    plant=read_rows(a.source_root/FILES["plant_effort"])
    marine=read_rows(a.source_root/FILES["circuits"])
    hostmap={(x["key"],x["island"]):x for x in host}
    effort=defaultdict(int)
    for x in plant:
        effort[x["island"]]+=int(x["records"])
    marine_map={}
    for x in marine:
        for key in ((x["island_a"],x["island_b"]),(x["island_b"],x["island_a"])):
            marine_map[key]={
                "distance_km":float(x["port_distance_km"]),
                "uniform_resistance":float(x["uniform_distributed_port_resistance"]),
            }

    by_species=defaultdict(list)
    for x in butterflies:
        by_species[x["speciesKey"]].append(x)
        if (x["speciesKey"],x["landmass_id"]) not in histmap:
            raise ValueError("Missing historical species-island key")

    modeled={}
    audits={}
    def prior_index(other,after_2017):
        h=histmap[other["speciesKey"],other["landmass_id"]]
        return (int(h["pre2001_indexed_records"])>0 or
                int(other["pre2001_2010_indexed_records"])>0 or
                (after_2017 and int(other["future2011_2017_indexed_records"])>0))

    for period,after_2017 in (("2011_2017",False),("2018_2022",True)):
        group=[]
        status=defaultdict(int)
        expected_new=sum(int(x[f"indexed_first_{period}_longer_history"]) for x in history)
        for x in butterflies:
            key=(x["speciesKey"],x["landmass_id"])
            historic=histmap[key]
            if int(historic[f"eligible_{period}_longer_history"])!=1:
                continue
            status["source_eligible"]+=1
            h=hostmap.get(key)
            if h is None or effort[x["landmass_id"]]<=0:
                status["not_historical_port_plant_covered"]+=1
                continue
            status["botanical_covered"]+=1
            prior=[]
            for other in by_species[x["speciesKey"]]:
                if other["landmass_id"]==x["landmass_id"] or not prior_index(other,after_2017):
                    continue
                edge=marine_map.get((x["landmass_id"],other["landmass_id"]))
                if edge is not None:
                    prior.append((edge["distance_km"],other["landmass_id"],edge["uniform_resistance"]))
            if not prior:
                status["no_indexed_other_original_island"]+=1
                continue
            dist,source,resistance=min(prior)
            logd=np.log1p(dist)
            fit_resistance=float(coeff[0] + coeff[1]*logd + coeff[2]*logd**2)
            access=float(-(np.log1p(resistance)-fit_resistance))
            target=int(historic[f"indexed_first_{period}_longer_history"])
            group.append({
                "species":x["speciesKey"],
                "island":x["landmass_id"],
                "earlier_island":source,
                "event":target,
                "native":int(int(h["native"])>0),
                "introduced":int(int(h["introduced"])>0),
                "effort":effort[x["landmass_id"]],
                "distance_km":dist,
                "access_residual":access,
            })
        status["source_total_extended_first_event"]=expected_new
        status["current_model_risk_rows"]=len(group)
        status["current_model_first_index_rows"]=sum(x["event"] for x in group)
        status["current_model_introduced_host_indexed_event_rows"]=sum(
            x["event"] and x["introduced"] for x in group
        )
        status["current_model_distinct_islands"]=len({x["island"] for x in group})
        status["current_model_distinct_species"]=len({x["species"] for x in group})
        modeled[period]=group
        audits[period]=dict(status)

    results=model_sets(modeled)
    manuscript_notice={
        "interpretation":"Extended-time GBIF observation-only first indexing, not documented new arrival",
        "do_not_use_as_colonization_effect":True,
        "model_basis":"Species+island FE LPM, native and candidate introduced host indices, distance, uniform-ocean-residual access, and their host × access interaction",
        "source_only_comparison":"All species-island cells screened for 1700-2000 index; unlike pre-2001 positive-only exclusion",
        "source_historical_biases":[
            "A zero GBIF index (even 1700–2010) does not imply genuine butterfly absence",
            "Nearest historical conspecific source is restricted to original 51 selected islands, not a proven exporting population",
            "2001–2010 plant records at 5 km ports represent candidate regional known hosts, not historic island-local introduction or larval use",
            "Uneven historical butterfly recording effort and species/port selection remain",
            "Observational first indexing is not verified colonization or establishment",
        ],
        "do_not_change":[
            "239 butterfly resource geography +54.9%",
            "126 island vs 237 mainland resource reconstruction (+118.9% vs +40.9%)",
            "51-island retrospective Circuitscape community held-island prediction +5.082%",
        ],
    }

    a.out_dir.mkdir(parents=True,exist_ok=True)
    report={
        "analysis":"source-uniform-pre2001-full11475-port-host-ocean-model",
        "source_panel_sha256":source_hashes,
        "historical_full_panel_csv_sha256":digest(a.history_csv),
        "historical_full_panel_summary_sha256":digest(a.history_summary),
        "historical_risk_stats":summary["stats"],
        "original_index_pilot_population":{
            k:{"n":v["plant_records_ge_10"]["n"],
               "events":v["plant_records_ge_10"]["events"],
               "introduced_event_cells":v["plant_records_ge_10"]["introduced_event_cells"],
               "interaction_pp_perSD":v["plant_records_ge_10"]["introduce_x_sea_access_effect_pp_per_SD"]}
            for k,v in model_sets(oldrisk).items()
        },
        "extended_time_risk_support":audits,
        "extended_time_model_results":results,
        "inference":manuscript_notice,
    }
    out=a.out_dir/"long_history_port_host_ocean_model.json"
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("HISTORICALLY_COHERENT_SOURCE_INDEX_MODEL",json.dumps({
        "history_stats":summary["stats"],
        "status":audits,
        "model":{key:dict(
          early_index_risk=val["plant_records_ge_10"]["n"],
          observed_first_index=val["plant_records_ge_10"]["events"],
          host_positive_events=val["plant_records_ge_10"]["introduced_event_cells"],
          host_x_sea_access_beta_pp=val["plant_records_ge_10"].get("introduce_x_sea_access_effect_pp_per_SD"),
          approx_95pct_interval_pp=val["plant_records_ge_10"].get("approximate_95pct_interval_pp"),
        ) for key,val in results.items()},
        "file":str(out)
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
