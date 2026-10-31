#!/usr/bin/env python3
"""Reconcile focal butterfly foodplant 5-km historical source vs 1-km coastal samples.

Four *host × island* combinations support six nonindependent butterfly × island
endpoint mentions. A GBIF nonindex is UNKNOWN plant availability, not absence.
No ship introduction, island-native status, larval use, breeding or flights
can be inferred from these observational source indexes.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
FILES={
 "focal_years":("provenance/reviewer_defenses/results/dated_focal_host_taxonomy_20261031/focal_host_prior_to_butterfly_20261031.json",
                "94a82adf74eff89ee073fb2016613dda2aee7a26"),
 "five_km":("provenance/reviewer_defenses/results/dated_focal_host_taxonomy_20261031/focal_historical_host_taxon_year_audit.json",
            "3048c769d250f1b2f8df97c45aa20c8b91bc573b"),
 "one_km":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_plant_site_facets.csv",
           "eb09a3f3ec053c00b43e8bc0678d6643ce02697b"),
 "land":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/land_fraction_sites.csv",
         "b2f7a7597550b61abbc182b6c5c6f2f22ae6a898"),
 "prior_5links":("provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_links.csv",
                "88db773c74dea9cd2d58bde0b67a394d93226869"),
 "migrant_proxy":("provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
                  "439dd1034ca16f96ff94222695d98dc39747c2e4")
}

def source(name):
    path,expect=FILES[name]
    b=(ROOT/path).read_bytes()
    digest=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    if digest!=expect:
        raise ValueError("Pinned Git source changed: "+path+" "+digest)
    if path.endswith(".json"):
        return json.loads(b)
    return list(csv.DictReader(io.StringIO(b.decode("utf8"))))

def run():
    years=source("focal_years")
    five=source("five_km")
    one=source("one_km")
    land=source("land")
    links=source("prior_5links")
    migration=source("migrant_proxy")
    migration_by_species={x["species"]:x for x in migration}
    assert len(migration)==225
    assert len(one)==114 and len(land)==114 and len(links)==5
    assert years["butterfly_links_checked"]==3
    assert years["butterfly_species_island_endpoints"]==6
    assert years["host_species"]==2
    assert len(five["island_taxon_observations"])==4
    look={(x["gbif_species_key"],x["island_id"]):x for x in five["island_taxon_observations"]}
    landby={x["site_id"]:x for x in land}
    assert len(landby)==114
    # All shared endpoint records must belong to previously independently
    # published *butterfly × foodplant × island* pairs, not a selected new host.
    edge={(x["species"],x["host_name"]) for x in links}
    pairs={}
    butterfly=[]
    for q in years["focal_links"]:
        b=q["butterfly"]
        plant=q["source_catalogue_compatible_host"]
        key=q["gbif_accepted_plant_speciesKey"]
        assert (b,plant) in edge
        for ep in q["island_record_comparisons"]:
            island=str(ep["island"])
            before=int(ep["host_first_indexed_plant_year"])
            after=int(ep["butterfly_first_indexed_year"])
            assert before<after
            original=look[(key,island)]
            assert original["host_name"]==plant
            assert int(original["first_indexed_year"])==before
            assert original["number_ports_with_matching_species_key"]==len(original["port_ids"])
            assert len(original["port_ids"])>0
            pkey=(key,island)
            focal_rows=[x for x in one if x["landmass_id"]==island]
            if not focal_rows:raise ValueError("Missing all 1-km island floral samples")
            p=[x for x in focal_rows if x["site_type"]=="port"]
            off=[x for x in focal_rows if x["site_type"]=="offport"]
            assert len(p)==len(original["port_ids"])
            assert len(off)==2*len(p)
            assert {x["port_id"] for x in p}==set(original["port_ids"])
            assert all(x["query_status"]=="COMPLETE" for x in focal_rows)
            def f(r):
                return key in set(filter(None,r["species_keys"].split(";")))
            if pkey not in pairs:
                positives_port=[x for x in p if f(x)]
                positives_off=[x for x in off if f(x)]
                pairs[pkey]={
                    "plant":plant,"plant_accepted_speciesKey":key,"island":island,
                    "first_indexed_year_in_original_5km":before,
                    "five_km_source_ports_indexing_plant":len(original["port_ids"]),
                    "five_km_source_indexed_record_count_sum_may_overlap":original["sum_indexed_records_across_potentially_overlapping_port_buffers"],
                    "one_km_port_sites":len(p),
                    "one_km_port_plant_index_sites":len(positives_port),
                    "one_km_offport_sites":len(off),
                    "one_km_offport_plant_index_sites":len(positives_off),
                    "one_km_port_species_background_min":min(int(x["plant_indexed_species_keys"]) for x in p),
                    "one_km_offport_species_background_min":min(int(x["plant_indexed_species_keys"]) for x in off),
                    "one_km_all_sites_have_positive_botanical_index":all(int(x["plant_indexed_species_keys"])>0 for x in focal_rows),
                    "one_km_positive_offport_ids":[x["site_id"] for x in positives_off],
                    "one_km_positive_port_ids":[x["site_id"] for x in positives_port],
                    "individual_one_km_sites":[{
                        "site_id":x["site_id"],
                        "port_id":x["port_id"],
                        "type":x["site_type"],
                        "host_indexed":f(x),
                        "all_indexed_plant_species":int(x["plant_indexed_species_keys"]),
                        "total_plant_records":int(x["plant_occurrence_index_count"]),
                        "land_fraction_1km":float(landby[x["site_id"]]["land_area_fraction_1km"]),
                    } for x in sorted(focal_rows,key=lambda y:y["site_id"])]
                }
            if b not in migration_by_species:
                raise ValueError("Focal butterfly missing from source 225 migration candidate panel")
            migration_label=migration_by_species[b]["author_2025_GBIF_requested_taxon_key_proxy"]
            if migration_label!="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST":
                raise ValueError("Focal butterfly lacks independent author-review positive migration proxy: "+b)
            butterfly.append({
                "butterfly":b,"migrant_proxy":"REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST",
                "butterfly_migration_proxy_is_positive_only_not_island_direct_flight":True,
                "historical_foodplant":plant,
                "foodplant_speciesKey":key,"island":island,
                "foodplant_index_year_5km":before,"butterfly_first_index_year":after,
                "indexed_lead_years":after-before,
                "host_island_unit":plant+"|"+island,
                "one_km_port_host_index":pairs[pkey]["one_km_port_plant_index_sites"],
                "one_km_offport_host_index":pairs[pkey]["one_km_offport_plant_index_sites"],
            })
    values=list(pairs.values())
    assert len(values)==4 and len(butterfly)==6
    assert len({x["butterfly"] for x in butterfly})==3
    assert len({x["plant"] for x in values})==2
    assert sum(x["one_km_port_sites"] for x in values)==11
    assert sum(x["one_km_offport_sites"] for x in values)==22
    totalP=sum(x["one_km_port_plant_index_sites"] for x in values)
    totalO=sum(x["one_km_offport_plant_index_sites"] for x in values)
    # In new data, only Alnus incana at island 276 is indexed off harbor.
    # This fact is not a hypothesis test; original 5-km source recorded other
    # plant points from nonidentical buffers and earlier GBIF indexing runs.
    assert totalP==0 and totalO==2
    assert sum(x["one_km_all_sites_have_positive_botanical_index"] for x in values)==4
    return {
      "schema":"chocho_focal_butterfly_foodplant_harbor_5km_vs_1km_v1",
      "status":"HISTORICAL_FOODPLANT_BUFFER_SCALE_DISCREPANCY_WITH_NONPORT_ALNUS_RECORDS",
      "source_git_blobs":{k:{"path":v[0],"sha":v[1]} for k,v in FILES.items()},
      "plant_years_2001_2010_prior_to_butterfly_years":True,
      "source_original_port_distance_km":5,
      "new_port_and_nonport_distance_km":1,
      "independent_host_island_units":len(values),
      "butterfly_species_island_mentions_not_independent":len(butterfly),
      "butterfly_species":3,
      "butterfly_species_matched_to_published_author_selected_migrants":len({x["butterfly"] for x in butterfly}),
      "indirect_migrant_roster_not_local_displacement_evidence":True,
      "plant_species":2,
      "five_km_original_foodplant_island_units_indexed":len(values),
      "one_km_port_plant_indexes_per_host_island":totalP,
      "one_km_offport_plant_indexes_per_host_island":totalO,
      "total_one_km_port_times_host_island_combinations":sum(x["one_km_port_sites"] for x in values),
      "total_one_km_offport_times_host_island_combinations":sum(x["one_km_offport_sites"] for x in values),
      "host_island":values,
      "butterfly_endpoints":butterfly,
      "limitations":[
         "Original positive 5km plant records do not imply occurrence within 1km of the harbor, and cannot be localized by this aggregate alone.",
         "A GBIF plant nonindex is not evidence of biological absence, failed establishment or a fixed dispersal barrier.",
         "New 1km and earlier 5km facets were indexed in different GBIF acquisitions, so reindexing effects remain possible.",
         "2019 WPI port locations are retrospective sites for 2001-2010.",
         "The 1km harbor and nonport circles have different land fractions; use land-area QC where available.",
         "Colias species share Medicago host, two islands and identical plant opportunity, not independent replication.",
         "Historic host-native/introduced classification on specific Danish islands is not established.",
         "Plant records lack ship/seed pathway, larval diet verification or observed local breeding.",
         "Six butterfly island endpoints were identified from the earlier focal chain, not a predeclared unconditional risk set.",
         "No inference about positive/negative true butterfly colonization or island stepping-stone routes."
      ]
    }

def main():
    obj=run()
    target=OUT/"focal_harbor_host_1km_vs_5km_results.json"
    target.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:obj[k] for k in ("status","independent_host_island_units",
      "butterfly_species_island_mentions_not_independent","one_km_port_plant_indexes_per_host_island",
      "one_km_offport_plant_indexes_per_host_island")},indent=2))

if __name__=="__main__":
    main()
