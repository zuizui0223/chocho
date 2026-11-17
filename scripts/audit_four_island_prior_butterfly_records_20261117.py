#!/usr/bin/env python3
"""Four-island prior-publication audit for later butterfly GBIF index entries.

Historical source lists have unequal resolution:
Phu Quoc, Mauritius, Cape Breton: island-specific reports.
Sri Lanka: 2006 NATIONAL island checklist, lacking exact encounter coordinates.
Neither is an estimator of native origin, first colonization, endemicity or real
population turnover. Keep all 317 original first-index events/225x51 source fixed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "gbif": ("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv", "b6a3c98eb76864052d45d520a55b7d830816126c"),
    "geo": ("results/island51_host_circuit_atlas_20261106/island_geography.csv", "f0c5b1139e54b4872d28682f0991f89f2368ab43"),
    "phu": ("results/butterfly_GBIF_effort_origin_gate_20261117/PHUQUOC_PRE2011_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv", "c7dcc85a216ff5f76339f277bdeb6329d522a611"),
    "mauritius": ("results/butterfly_GBIF_effort_origin_gate_20261117/MAURITIUS_1993_INDEPENDENT_CHECKLIST_MATCHES_20261117.csv", "e06fb0032b8dd2a9bd6344e2e88db0f201f3a4ee"),
    "cape": ("results/butterfly_GBIF_effort_origin_gate_20261117/CAPE_BRETON_1948_1984_INDEPENDENT_FIELD_MATCHES_20261117.csv", "7bf77344e0be697ff4a2497a36b91a3fd577d78c"),
    "srilanka": ("results/butterfly_GBIF_effort_origin_gate_20261117/SRILANKA_2006_IUCN_CHECKLIST_15_CANDIDATE_CROSSWALK_20261117.csv", "ba952428113ccb56936c2cf62e3f808ff328c7fc"),
}


def blob(p:Path):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode("ascii")+b"\0"+b).hexdigest()


def source(k):
    rel,sha=PINS[k]
    f=ROOT/rel
    if blob(f)!=sha:
        raise ValueError("Original pinned Git blob drift: "+k)
    with f.open(newline="",encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def jaccard(ids,sets):
    sum_j,n=0.,0
    for x,i in enumerate(ids):
        for j in ids[x+1:]:
            a,b=sets[i],sets[j]
            u=len(a|b)
            if not u:raise ValueError("Undefined early Jaccard")
            sum_j+=len(a&b)/u
            n+=1
    return sum_j/n,n


def run(raw,geo,phu,mauritius,cape,srilanka):
    if len(raw)!=11475:raise ValueError("GBIF original panel changed")
    ids={int(z["landmass_id"]):z for z in geo}
    if len(ids)!=51:raise ValueError("Original landmass count changed")
    for i,country,port in [(471,"VN","DUONG DONG"),(218,"MU","PORT LOUIS"),
                           (79,"CA","SYDNEY"),(31,"LK","KANKESANTURAI")]:
        if ids[i]["example_port_country_code"]!=country or ids[i]["example_port_name"]!=port:
            raise ValueError("Original port/geography ID mismatch "+str(i))
    if json.loads(ids[31]["source_port_names"])!=["KANKESANTURAI","COLOMBO","GALLE HARBOR"]:
        raise ValueError("Original Sri Lanka main-island port sample changed")
    first={}
    early={i:set() for i in ids}
    late={i:set() for i in ids}
    species=set()
    for row in raw:
        name=row["species"]
        taxon=int(row["speciesKey"]);island=int(row["landmass_id"])
        a,b,c=map(int,(row["pre2001_2010_indexed_records"],
                    row["future2011_2017_indexed_records"],
                    row["future2018_2022_indexed_records"]))
        species.add(taxon)
        if island not in ids or min(a,b,c)<0:raise ValueError("Original source row invalid")
        if a:early[island].add(taxon)
        if a+b+c:late[island].add(taxon)
        if not a and (b+c):
            key=(taxon,island)
            if key in first:raise ValueError("Duplicate original species-island risk event")
            first[key]={
                "species":name,"speciesKey":taxon,"landmass_id":island,
                "GBIF_first_index_window":"2011-2017" if b else "2018-2022",
                "GBIF_early_count":a,"GBIF_mid_count":b,"GBIF_late_count":c,
                "historical_source_status":"not_reviewed",
                "source_year_upper_bound":"","source_geographic_resolution":"",
                "historical_evidence_url":"","native_status":"unknown",
                "reproductive_establishment_status":"unknown",
                "island_endemic_status":"unassessed"
            }
    if len(species)!=225 or len(first)!=317:
        raise ValueError("Changed source butterfly species or first-index count")
    seen=set()
    source_count={}
    summary=[]

    def assign(source_kind, rows, island, year, geo_resolution, url_field, validator):
        nonlocal seen
        found=set()
        for row in rows:
            key=(int(row["speciesKey"]),int(row["landmass_id"]))
            if key in seen or key in found or key not in first or key[1]!=island:
                raise ValueError("Duplicate/unknown/out-of-island historical record "+str(key))
            if first[key]["species"]!=row["species"]:
                raise ValueError("Unreconciled historical taxonomic name "+str(key))
            if int(year(row))>=2011:
                raise ValueError("Historical evidence publication after start of test period")
            validator(row,key)
            found.add(key)
            v=first[key]
            v.update(historical_source_status=source_kind,
                     source_year_upper_bound=int(year(row)),
                     source_geographic_resolution=geo_resolution,
                     historical_evidence_url=row[url_field])
            if v["native_status"]!="unknown" or v["reproductive_establishment_status"]!="unknown" or v["island_endemic_status"]!="unassessed":
                raise ValueError("Historical checklists are not native or endemicity evidence")
        seen|=found
        return found

    ph=assign("prior_publication_exact_named_island",phu,471,
              lambda x:x["earliest_independent_record_year_upper_bound"],
              "named_island_field_report","source_url",
              lambda row,key: (
                  None if (
                    row["butterfly_origin_status"]=="unknown" and
                    row["butterfly_establishment_status"]=="unknown" and
                    row["butterfly_endemicity_status"]=="unassessed" and
                    int(row["earliest_independent_record_year_upper_bound"]) in (2006,2007) and
                    first[key]["GBIF_first_index_window"]=="2018-2022"
                  ) else (_ for _ in ()).throw(ValueError("Phu Quoc prior-source validation failed"))))
    if len(ph)!=12:raise ValueError("Phu Quoc historical count changed")

    ma=assign("prior_publication_exact_named_island",mauritius,218,
              lambda x:x["publication_year_upper_bound"],
              "named_island_annotated_species_list","source_url",
              lambda row,key:(
                  None if (row["name_match_strength"]=="source_exact_species" and
                           int(row["publication_year_upper_bound"])==1993 and
                           row["origin_status"]=="unknown" and
                           row["establishment_status"]=="unknown" and
                           row["endemicity_status"]=="unassessed" and
                           first[key]["GBIF_first_index_window"]=="2011-2017")
                  else (_ for _ in ()).throw(ValueError("Mauritius prior-source validation failed"))))
    if len(ma)!=3:raise ValueError("Mauritius source count changed")

    ca=assign("prior_publication_exact_named_island",cape,79,
              lambda x:x["observation_year_upper_bound"],
              "named_island_original_field_localities","source_url",
              lambda row,key: (
                  None if (int(row["observation_year_upper_bound"]) in (1948,1983) and
                           int(row["publication_year_upper_bound"]) in (1950,1984) and
                           row["origin_status"]=="unknown" and
                           row["establishment_status"]=="unknown" and
                           row["endemicity_status"]=="unassessed" and
                           first[key]["GBIF_first_index_window"]=="2011-2017")
                  else (_ for _ in ()).throw(ValueError("Cape Breton historical source validation failed"))))
    if len(ca)!=4:raise ValueError("Cape Breton source count changed")

    verified_sri=[x for x in srilanka if x["checklist_exact_species_match"]=="yes"]
    unknown_sri=[x for x in srilanka if x["checklist_exact_species_match"]=="not_found_in_2006_checklist"]
    if len(srilanka)!=15 or len(verified_sri)!=14 or len(unknown_sri)!=1 or unknown_sri[0]["species"]!="Ypthima huebneri":
        raise ValueError("2006 IUCN source-match scope mismatch")
    si=assign("prior_publication_national_island_checklist",verified_sri,31,
              lambda x:x["checklist_publication_year"],
              "national_island_not_exact_geo_point","source_url",
              lambda row,key: (
                  None if (int(row["checklist_publication_year"])==2006 and
                           row["island_name"]=="Sri Lanka" and
                           int(row["printed_page"]) in (59,60,61,63,64) and
                           row["butterfly_origin_status"]=="unknown" and
                           row["butterfly_establishment_status"]=="unknown" and
                           row["butterfly_endemicity_status"]=="unassessed" and
                           row["source_url"]=="https://portals.iucn.org/library/efiles/documents/2006-030.pdf")
                  else (_ for _ in ()).throw(ValueError("Sri Lanka national checklist source validation failed"))))
    if len(si)!=14:raise ValueError("Sri Lanka 14 first-index matches changed")
    unresolved=unknown_sri[0]
    ukey=(int(unresolved["speciesKey"]),31)
    if ukey not in first or unresolved["printed_page"].strip():
        raise ValueError("Unmatched Sri Lanka binomial cannot be classified")
    first[ukey]["historical_source_status"]="searched_2006_checklist_no_exact_binomial_unresolved"

    old={i:set(s) for i,s in early.items()}
    J0,p0=jaccard(sorted(ids),old)
    accumulations=[]
    for island_source,entries in [("Phu Quoc",ph),("Mauritius",ma),("Cape Breton",ca),("Sri Lanka 2006",si)]:
        for taxon,island in entries:old[island].add(taxon)
        z,n=jaccard(sorted(ids),old)
        accumulations.append({"source":island_source,"added_cells":len(entries),
            "early_indexed_or_prior_source_cells":sum(len(s) for s in old.values()),
            "mean_early_Jaccard":z,"delta_to_2022_cumulative_GBIF":None,"pairs":n})
    Jlate,pl=jaccard(sorted(ids),late)
    if (p0,pl)!=(1275,1275):raise ValueError("Original source pair count changed")
    for row in accumulations:row["delta_to_2022_cumulative_GBIF"]=Jlate-row["mean_early_Jaccard"]
    if [row["early_indexed_or_prior_source_cells"] for row in accumulations] != [748,751,755,769]:
        raise ValueError("Previous three-island historical ledger or new 14-source count drift")
    if abs(accumulations[2]["mean_early_Jaccard"]-0.1270751832600804)>1e-9:
        raise ValueError("Previously verified three-island source mismatch")
    if abs(accumulations[3]["mean_early_Jaccard"]-0.12806690628478226)>1e-9:
        raise ValueError("Four-island historical expected baseline mismatch")
    if len(seen)!=33:raise ValueError("Historical four-island source overlap or missingness")
    first_index_unknown=len(first)-len(seen)
    if first_index_unknown!=284:raise ValueError("Historical first-index count changed")
    return {
        "schema":"chocho_original225x51_four_island_historical_butterfly_GBIF_first_index_gates_20261117_v1",
        "status":"PUBLISHED_SOURCE_RECORD_PRECEDES_GBIF_INDEX_NOT_FIRST_ARRIVAL_OR_NATIVE_RANGE",
        "source_sha_blobs":{k:sha for k,(path,sha) in PINS.items()},
        "cohort":{"species":225,"landmasses":51,"species_island_cells":11475,
            "later_first_index_cells":317,"independent_sources_before_GBIF":33,
            "as_yet_unverified_first_indices":284,"pairs":1275,
            "geographically_named_island_source_matches":19,
            "national_island_2006_checklist_matches":14,
            "unmatched_2006_exact_binomial_but_not_absent":["Ypthima huebneri"]},
        "source_coverage":{"Phu_Quoc":12,"Mauritius":3,"Cape_Breton":4,"Sri_Lanka_national":14},
        "recorded_community_Jaccard":{
            "original_early_GBIF_only_cells":sum(len(x) for x in early.values()),
            "original_early_mean_Jaccard":J0,
            "cumulative_through_2022_indexed_cells":sum(len(x) for x in late.values()),
            "cumulative_through_2022_mean_Jaccard":Jlate,
            "GBIF_only_delta":Jlate-J0,
            "historical_source_additions":accumulations,
            "all_33_source_augmented_delta":accumulations[-1]["delta_to_2022_cumulative_GBIF"],
            "difference_in_descriptive_delta":accumulations[-1]["delta_to_2022_cumulative_GBIF"]-(Jlate-J0)},
        "claims_do_not_make":[
            "The GBIF first-index date is an ecological first-arrival or breeding-colonization date",
            "The 33/317 source-selected fraction estimates the population-wide first-index error rate",
            "A national island checklist proves exact within-island sighting locality or native origin",
            "These checklists demonstrate island endemicity, native range, occupancy or reproductive establishment",
            "Historical source-augmented recorded Jaccard is an unbiased true assemblage change or complete effort correction",
            "Unknown source records imply ecological absence or new colonisation"
        ]
    }, list(first.values())


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir",type=Path,default=Path("_four_island_history_index"))
    a=p.parse_args()
    tables=[source(k) for k in PINS]
    result,records=run(*tables)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    (a.output_dir/"FOUR_ISLAND_HISTORY_VS_GBIF_317_SOURCE_AUDIT_20261117.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    with (a.output_dir/"FOUR_ISLAND_317_FIRST_INDICES_SOURCE_STATUS_20261117.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(sorted(records,key=lambda x:(x["landmass_id"],x["speciesKey"])))
    print(json.dumps({"cohort":result["cohort"],"Jaccard":result["recorded_community_Jaccard"]},indent=2))
if __name__=="__main__":
    main()
