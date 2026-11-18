#!/usr/bin/env python3
"""Match independently identified endemic butterflies with source-observed
named host-plant species on individual named islands, without projecting a
coarse TDWG3 botanical range to local caterpillars.

22 independently checklisted butterfly×island combinations for 8 source-known
butterflies form the denominator. Only the observed plant-with-caterpillar
records at source-locality grain are accepted as direct named-host use.

Do not call a garden caterpillar self-sustaining, or a host from a nonlocal
NHM database locally used. Oviposition-only evidence is not completed larval
development; unknown host species is not zero native host resources.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

ROOT=Path("provenance/reviewer_defenses/results")
INBO=ROOT/"independent_island_butterfly_host_coverage_20261117"/"species_island_coverage.csv"
FOCAL=ROOT/"island_endemic_synonym_origin_20261118"/"independent_butterfly_host_plant_species_events_v1.csv"
NEW=ROOT/"island_endemic_synonym_origin_20261118"/"E_grancanariensis_H_tamadabae_source_evidence_20261118.csv"

def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def run(source_edges,output):
    edges=read(source_edges)
    original=read(INBO)
    focus=read(FOCAL)
    new=read(NEW)
    taxa={x["butterfly"] for x in edges}
    if (len(edges),len(taxa),len(original),len(focus),len(new))!=(13,8,227,13,4):
        raise ValueError("Frozen independent checklist/focal evidence changed")
    focal_events=[]
    allowed={
        "wild_larva_photo":"wild_named_plant",
        "garden_larva_photo":"garden_named_plant",
        "wild_eggs_larvae_report":"wild_named_plant"
    }
    for r in focus:
        if r["butterfly_name"] not in taxa:continue
        grade=allowed.get(r["source_evidence_type"])
        if grade is None:continue
        island=r["named_island"]
        if not island or not r["reported_host_binomial"].strip():
            raise ValueError("Direct named-plant source field lacks site or species")
        focal_events.append({
            "butterfly":r["butterfly_name"],
            "island_name":"Canary Islands, "+island,
            "plant_species":r["reported_host_binomial"],
            "plant_source_origin":r["reported_plant_origin"],
            "evidence_class":grade,
            "event_id":r["event_id"],
            "source_url":r["butterfly_host_source_url"]
        })
    for r in new:
        if r["butterfly"] not in taxa or r["evidence_class"]!="field_larvae_on_named_plant_and_laboratory_pupation":
            continue
        if r["plant_name_resolution"]!="species":
            raise ValueError("New source larval host not species-identified")
        focal_events.append({
            "butterfly":r["butterfly"],
            "island_name":"Canary Islands, "+r["island"],
            "plant_species":r["host_plant"],
            "plant_source_origin":"native", # independently validated with WCVP in #44213255146
            "evidence_class":"wild_named_plant",
            "event_id":r["event_id"],
            "source_url":r["source_url"]
        })
    if len(focal_events)!=4 or len({(r["butterfly"],r["island_name"]) for r in focal_events})!=4:
        raise ValueError("Original four distinct island-local named-larva reports not rebuilt")
    inbo={(r["scientific_name"],r["island"]) for r in original if r["scientific_name"] in taxa}
    if len(inbo)!=22:
        raise ValueError("Original independent source-known eight butterfly×island denominator changed")
    direct={(e["butterfly"],e["island_name"]):e for e in focal_events}
    if not set(direct)<=inbo:raise ValueError("Direct field events have named island not in independently checklisted butterfly")
    if set(direct)!={
       ("Pieris cheiranthi","Canary Islands, Tenerife"),
       ("Pieris cheiranthi","Canary Islands, La Palma"),
       ("Gonepteryx cleobule","Canary Islands, La Gomera"),
       ("Euchloe grancanariensis","Canary Islands, Gran Canaria")
    }:
        raise ValueError("Direct butterfly+island names are different")
    rows=[]
    for name,island in sorted(inbo):
        e=direct.get((name,island))
        rows.append({
            "butterfly":name,
            "independent_INBO_named_island":island,
            "source_level_direct_named_plant_larva":int(e is not None),
            "source_level_wild_named_plant_larva":int(bool(e) and e["evidence_class"]=="wild_named_plant"),
            "source_level_garden_named_plant_larva":int(bool(e) and e["evidence_class"]=="garden_named_plant"),
            "plant_species_direct_named":e["plant_species"] if e else "",
            "source_plant_origin":e["plant_source_origin"] if e else "",
            "evidence_class":e["evidence_class"] if e else "NOT_OBSERVED_IN_CHECKED_SOURCES",
            "original_source_ID":e["event_id"] if e else "",
            "original_url":e["source_url"] if e else "",
            "caterpillar_survived_to_adult_confirmed":0,
            "new_butterfly_island_establishment_confirmed":0
        })
    if sum(x["source_level_direct_named_plant_larva"] for x in rows)!=4 or sum(x["source_level_wild_named_plant_larva"] for x in rows)!=3:
        raise ValueError("Wild versus garden caterpillar evidence accidentally combined")
    if sum(x["source_level_garden_named_plant_larva"] for x in rows)!=1:
        raise ValueError("Garden-only source not explicit")
    source_distinct_taxa={r["butterfly"] for r in rows if r["source_level_direct_named_plant_larva"]==1}
    if len(source_distinct_taxa)!=3:
        raise ValueError("Named field larval evidence expected for three selected endemic butterfly species")
    origin_counts={
      "introduced":sum(x["source_plant_origin"]=="introduced" for x in rows),
      "native":sum(x["source_plant_origin"]=="native" for x in rows)
    }
    if origin_counts!={"introduced":2,"native":2}:
        raise ValueError("Direct field alien/native host plant split drift")
    res={
      "schema":"chocho_independent8_endemic_source_grain_local_larval_identifiability_v1",
      "independent_endemic_butterflies_global_list_n":18,
      "source_selected_exact_named_host_species_butterflies_n":8,
      "independent_named_butterfly_island_pairs_in_selected8_n":22,
      "named_island_direct_larval_host_species_observation_n":4,
      "named_island_strict_wild_larval_host_species_observation_n":3,
      "named_island_garden_larval_host_species_observation_n":1,
      "butterfly_species_with_named_island_direct_plant_larval_observation_n":3,
      "named_island_direct_larval_evidence_by_origin":origin_counts,
      "sources_confirmed_local_plant_butterfly_events":focal_events,
      "limits":[
        "Four measured species×island cells out of 22 from an evidence-selected eight butterfly subset; this is NOT a survey-derived occupancy or host-use prevalence.",
        "INBO records of butterfly on an island do not themselves confirm reproduction.",
        "Other 18 butterfly×island pairs lack verified source-local plant use in this review; not evidence of no larval hosts.",
        "P. cheiranthi garden larvae are physically on plant but self-sustaining wild population success is unknown.",
        "The WCVP v13 source-origin codes CNY/AZO/MDR are archipelago-level, not species-level native/alien status on each island.",
        "A globally expanded WCVP host plant does not demonstrate a species is feeding or dispersing in the remote added region.",
        "Larva on plant/observed egg is not completion to adult; no demographic or colonial dispersal mechanism measured."
      ]
    }
    output.mkdir(parents=True,exist_ok=True)
    with (output/"independent8_butterfly_named_island_host_presence22.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/"independent8_named_island_direct_larvae_summary.json").write_text(
       json.dumps(res,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("INDEPENDENT_8_BUTTERFLY_ISLAND_LOCAL_HOST_SOURCE_GRAIN_AUDITED="+json.dumps({
       "eight_species":8,"source_checklisted_island_pairs":22,"direct_species_named_local_larvae":4,
       "wild_local_larvae":3,"garden_local_larvae":1,
       "actual_named_local_plant_butterfly_taxa":sorted(source_distinct_taxa),
       "direct_local_origin":origin_counts,
       "no_new_occupancy_or_colonization_assumed":True
    },ensure_ascii=False),flush=True)
    return res
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-edges",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    run(a.source_edges,a.out_dir)
