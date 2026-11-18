#!/usr/bin/env python3
"""Distinguish added host species from added island-scale host-presence units.

Only two source-verified Pieris cheiranthi island records qualify. The native
host references and alien host larvae are from different dates/sites, so
co-presence is an island-wide union of historical reports, NOT simultaneous
local resource accessibility, larval performance, or butterfly range expansion.
"""
import csv
import json
from pathlib import Path

ROOT=Path("results/focal_island_larval_plant_species_20261118")
EVIDENCE=ROOT/"focal_host_plant_species_evidence_checked.csv"
WCVP=ROOT/"focal_island_named_host_plant_species_WCVP_CNY_v13.csv"
SITES={"Tenerife":{"native":"Crambe strigosa","alien":"Tropaeolum majus"},
       "La Palma":{"native":"Descurainia millefolia","alien":"Tropaeolum majus"}}

def load(path):
    with path.open(encoding="utf-8",newline="") as h:
        return list(csv.DictReader(h))

def audit():
    e=load(EVIDENCE);p=load(WCVP)
    if len(e)!=13 or len(p)!=13:raise ValueError("Original source selected evidence cohort changed")
    em={r["event_id"]:r for r in e}
    pm={r["event_id"]:r for r in p}
    if set(em)!=set(pm):raise ValueError("Species level plant status must be measured on original evidence identities")
    site_records=[]
    for island,plants in SITES.items():
        support_native=[x for x in e if
           x["named_island"]==island and x["butterfly_name"]=="Pieris cheiranthi" and
           x["reported_host_binomial"]==plants["native"] and
           pm[x["event_id"]]["admissible_species_origin"]=="native_CNY_exact_species"]
        support_alien=[x for x in e if
           x["named_island"]==island and x["butterfly_name"]=="Pieris cheiranthi" and
           x["reported_host_binomial"]==plants["alien"] and
           x["strict_direct_named_island_larva_use"]=="1" and
           pm[x["event_id"]]["trusted_island_specific_direct_introduced_use"]=="1"]
        if len(support_native)!=1 or len(support_alien)!=1:
            raise ValueError("Focal native native-host vs alien caterpillar comparison unsupported on "+island)
        source_native=support_native[0]
        source_alien=support_alien[0]
        if source_native["event_id"]==source_alien["event_id"]:
            raise ValueError("Source independence accidentally collapsed")
        site_records.append({
            "butterfly":"Pieris cheiranthi",
            "island":island,
            "documented_native_larval_host_species":plants["native"],
            "documented_exotic_larval_host_species":plants["alien"],
            "native_evidence_type":source_native["source_evidence_type"],
            "exotic_evidence_type":source_alien["source_evidence_type"],
            "exotic_photo_date":source_alien["date_text"],
            "native_host_reference_event":source_native["event_id"],
            "exotic_host_photo_event":source_alien["event_id"],
            "native_host_doc_n_min":1,
            "introduced_host_doc_n_min":1,
            "documented_host_species_options_increment_min":1,
            "native_larval_host_documented_island_unit":1,
            "after_introduction_documented_larval_host_island_unit":1,
            "documented_new_island_host_presence_unit":0,
            "time_coincident_resources_verified":0,
            "population_benefit_verified":0
        })
    if {x["island"] for x in site_records}!={"Tenerife","La Palma"}:
        raise ValueError("Unexpected island comparison")
    if sum(x["documented_new_island_host_presence_unit"] for x in site_records)!=0:
        raise ValueError("Host-species options wrongly converted to expanded geographic units")
    if sum(x["documented_host_species_options_increment_min"] for x in site_records)!=2:
        raise ValueError("Source documented added host species unavailable")
    result={
      "schema":"chocho_focal_endemic_island_host_portfolio_vs_geo_unit_v1",
      "underlying_butterfly":"Pieris cheiranthi",
      "butterfly_island_units_source_documented":len(site_records),
      "minimum_new_documented_introduced_host_species_butterfly_island_links":2,
      "newly_gained_island_units_under_native_host_already_documented":0,
      "documented_native_host_island_units":2,
      "documented_larval_host_presence_contemporary_units":2,
      "islands":site_records,
      "logic":"Named alien host added to a butterfly's documented plant species repertoire in each source island where a named native larval host was already independently documented. Does not newly create the island as a geographic host-presence unit.",
      "strict_boundary":[
        "An island-level union across source years is not simultaneous co-occurrence in a habitat or during a butterfly's generation.",
        "One native plus one alien documented host species does NOT mean equal abundance/viability, increased usable biomass, or reproductive success.",
        "One butterfly on two islands is neither a global sample nor a causal estimate of GEB +118.9% resource geography.",
        "An island-level named flora record is a different measurement scale than the original butterfly x WGSRPD3 global range units.",
        "Original fossil/historical extinct records never used; plant name conflicts excluded."
      ]
    }
    with (ROOT/"P_cheiranthi_native_alien_host_portfolio_two_islands.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(site_records[0]));w.writeheader();w.writerows(site_records)
    (ROOT/"P_cheiranthi_native_alien_host_portfolio_two_islands.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("P_CHEIRANTHI_PORTFOLIO_VS_NEW_GEOGRAPHY_REVIEWED="+json.dumps({
        "islands":len(site_records),"new_documented_host_species_island_links":2,
        "new_documented_island_presence_units":0,
        "no_claim_of_reproduction_or_range_expansion":True
    },ensure_ascii=False))
    return result

if __name__=="__main__":
    audit()
