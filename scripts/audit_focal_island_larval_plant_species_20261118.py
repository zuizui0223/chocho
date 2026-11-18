#!/usr/bin/env python3
"""Validate independently sourced butterfly -> plant SPECIES larval-host reports.

A named plant congener is NOT assumed to be the actual larval host. Only
directly named plant binomials are retained; genus-only source reports cannot
support a species-level island resource expansion estimate. Photo/locality
observations, regional literature statements, and historical-extinct species
are strictly separated.
"""
import csv
import json
from collections import defaultdict, Counter
from pathlib import Path

ROOT=Path("provenance/reviewer_defenses/results")
SOURCE=ROOT/"island_endemic_synonym_origin_20261118"/"independent_butterfly_host_plant_species_events_v1.csv"
INBO=ROOT/"independent_island_butterfly_host_coverage_20261117"/"species_island_coverage.csv"
OUT=Path("results/focal_island_larval_plant_species_20261118")

def read(path):
    with path.open(encoding="utf-8-sig",newline="") as h:
        return list(csv.DictReader(h))

def run():
    events=read(SOURCE)
    checklist=read(INBO)
    assert len(events)==13 and len(checklist)==227
    checklist_names={r["scientific_name"] for r in checklist}
    checklist_pairs={(r["scientific_name"],r["island"]) for r in checklist}
    madeira_extinct={"Pieris wollastoni"}
    seen=set()
    checked=[]
    for e in events:
        event=e["event_id"]
        if not event or event in seen:
            raise ValueError("Missing or repeated event ID: "+str(event))
        seen.add(event)
        taxon=e["butterfly_name"]
        host=e["reported_host_binomial"]
        island=e["named_island"]
        if taxon not in checklist_names:
            raise ValueError("Taxon absent from independent butterfly denominator: "+taxon)
        if len(host.split(" "))!=2:
            raise ValueError("Plant is not a named species binomial: "+host)
        if island and island!="Madeira":
            if (taxon,"Canary Islands, "+island) not in checklist_pairs:
                raise ValueError("Literature locality is outside source-checklisted butterfly islands: "+taxon+", "+island)
        if island=="Madeira" and (taxon,"Madeira, Madeira") not in checklist_pairs:
            raise ValueError("Madeira original butterfly checklist disagreement: "+taxon)
        if not e["butterfly_host_source_url"].startswith("https://"):
            raise ValueError("Missing transparent original host report source: "+event)
        if not e["plant_origin_source_url"].startswith("https://"):
            raise ValueError("Missing independent origin source: "+event)
        if e["reported_plant_origin"] not in {"native","introduced"}:
            raise ValueError("Unverified host plant origin: "+event)
        if e["name_resolution_grade"] not in {
            "high","medium","taxon_unresolved","historical_extinct"
        }:
            raise ValueError("Unknown taxonomic reliability: "+event)
        if not e["scope_note"].strip():
            raise ValueError("Evidence interpretation notes required: "+event)
        historical=taxon in madeira_extinct
        if historical and e["name_resolution_grade"]!="historical_extinct":
            raise ValueError("Historical extinct classification not separated")
        source_type=e["source_evidence_type"]
        if source_type not in {
            "wild_larva_photo","garden_larva_photo","wild_eggs_larvae_report",
            "wild_larva_report","field_host_report","literature_host_report",
            "curated_host_report","historical_host_report","historical_breeding"
        }:
            raise ValueError("Unknown larval evidence type: "+source_type)
        direct=source_type in {
            "wild_larva_photo","garden_larva_photo",
            "wild_eggs_larvae_report","wild_larva_report"
        }
        if direct and not island:
            raise ValueError("Direct larval event must have identified island: "+event)
        if direct and not e["date_text"].strip():
            raise ValueError("Direct larval event must have source date: "+event)
        output={
          **e,
          "named_island_matches_independent_butterfly_checklist":int(bool(island)),
          "evidence_direct_egg_or_larva_with_host":int(bool(direct)),
          "historically_extinct_butterfly":int(historical),
          "usable_current_species_host_identity":int(
             not historical and e["name_resolution_grade"] in {"high","medium"}
          ),
          "island_specific_extant_host_association":int(
             not historical and bool(island) and
             e["name_resolution_grade"] in {"high","medium"}
          ),
          "strict_direct_named_island_larva_use":int(
             not historical and bool(island) and direct and
             e["name_resolution_grade"] in {"high","medium"}
          ),
          "confirmed_larval_performance":0,
          "independently_verified_actual_butterfly_colonization":0
        }
        checked.append(output)
    alien_direct=[
        x for x in checked if x["strict_direct_named_island_larva_use"]
        and x["reported_plant_origin"]=="introduced"]
    if {(x["butterfly_name"],x["named_island"],x["reported_host_binomial"]) for x in alien_direct} != {
        ("Pieris cheiranthi","Tenerife","Tropaeolum majus"),
        ("Pieris cheiranthi","La Palma","Tropaeolum majus")
    }:
        raise ValueError("Alien larval use field localities must be independently documented")
    native_direct=[
        x for x in checked if x["strict_direct_named_island_larva_use"]
        and x["reported_plant_origin"]=="native"]
    if len(native_direct)!=1 or native_direct[0]["butterfly_name"]!="Gonepteryx cleobule":
        raise ValueError("Native direct event stratum changed")
    if any(x["island_specific_extant_host_association"] for x in checked if x["historically_extinct_butterfly"]):
        raise ValueError("Extinct species leaked into current resource estimate")
    pchei=[x for x in checked if x["butterfly_name"]=="Pieris cheiranthi"]
    if not {x["reported_plant_origin"] for x in pchei}=={"native","introduced"}:
        raise ValueError("Native and alien specific-host literature contrast lost")
    if len(checked)!=13: raise ValueError("Evidence ledger size changed")
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/"focal_host_plant_species_evidence_checked.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(checked[0]))
        w.writeheader(); w.writerows(checked)
    summary={
      "schema":"chocho_island_butterfly_larval_host_species_locality_evidence_gate_v1",
      "denominator":"13 source-selected evidence records, NOT all island butterflies, NOT a prevalence sample",
      "source_events":len(checked),
      "unique_butterflies_with_named_host_binomial":len({x["butterfly_name"] for x in checked}),
      "nonextinct_named_host_species_records":sum(not x["historically_extinct_butterfly"] for x in checked),
      "taxonomic_concept_unresolved":sum(x["name_resolution_grade"]=="taxon_unresolved" for x in checked),
      "historically_extinct_records_excluded":sum(x["historically_extinct_butterfly"] for x in checked),
      "strict_direct_named_island_host_use_count":sum(x["strict_direct_named_island_larva_use"] for x in checked),
      "strict_direct_introduced_host_use":{
          "events":len(alien_direct),
          "distinct_butterflies":sorted({x["butterfly_name"] for x in alien_direct}),
          "different_islands":sorted({x["named_island"] for x in alien_direct}),
          "host_plants":sorted({x["reported_host_binomial"] for x in alien_direct}),
          "source_event_ids":[x["event_id"] for x in alien_direct]
      },
      "strict_direct_native_host_use_count":len(native_direct),
      "focal_endemic_butterfly_native_and_alien_host_species":"Pieris cheiranthi",
      "cooccurrence_on_butterfly_named_island_vs_documented_larval_use":
        "Separate: source confirms introduced T. majus availability in Tenerife and La Palma, but local larva records are only two events, not population establishment",
      "taxon_reconciliation_status":{
        "Brachypodium arbusculum_to_arbuscula":"quarantined",
        "Crambe santosii_to_gigantea":"quarantined"
      },
      "important_negative_control":"Tenerife 1984 collected P. cheiranthi larvae fed on T. majus but ~12 died before pupation due to parasitoid Hymenoptera (Acosta-Fernandez 2009); not evidence that T. majus caused mortality or that all larvae fail to develop",
      "unknown":"Butterfly population establishment, self-sustaining growth, larval quality, species-specific alien host availability on every island, and endemics' mean resource expansion ratio cannot be computed from this selected evidence ledger",
      "original_GEB_known_host_opportunity_preserved":"+54.9% global / +118.9% original selected island potential resources, not recomputed"
    }
    (OUT/"focal_host_plant_species_evidence_summary.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("DIRECT_ALIEN_HOST_SPECIES_ISLANDS_SOURCE_VERIFIED="+json.dumps({
      "alien_field_islands":summary["strict_direct_introduced_host_use"]["different_islands"],
      "alien_direct_events":len(alien_direct),
      "native_direct_events":len(native_direct),
      "nonextinct":summary["nonextinct_named_host_species_records"],
      "taxon_unresolved":summary["taxonomic_concept_unresolved"]
    },ensure_ascii=False))

if __name__=="__main__":
    run()
