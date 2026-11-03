#!/usr/bin/env python3
"""Q2-F. Whose reconstructed *native host-resource geography* gains on islands?

This is NOT butterfly native range, butterfly origin, local presence, or proof
of colonial import. The geographic source label is an outcome-blind, pinned
WCVP × HOSTS native-host footprint definition, fixed before island additions.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

REGION_NAMES={
 "1":"Europe","2":"Africa","3":"Temperate Asia","4":"Tropical Asia",
 "5":"Australasia","6":"Pacific","7":"Northern America","8":"Southern America",
 "9":"Antarctic"
}

def records(path):
    with path.open(newline="",encoding="utf-8") as handle:
        yield from csv.DictReader(handle)

def savecsv(path,values):
    if not values:
        raise RuntimeError(f"No data to save: {path}")
    with path.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)

def main():
    p=argparse.ArgumentParser()
    for name in ("panel-csv","host-csv","native-distribution","q2-edges","q2-audit",
                 "level1-json","q1-region-metrics","empire-islands","output-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    m=list(records(a.panel_csv))
    if len(m)!=239 or len({r["species"] for r in m})!=239:
        raise RuntimeError("Wrong frozen 239-species panel")
    q2=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    edges=list(records(a.q2_edges))
    if q2["scope"]["introduced_butterfly_plant_island_edges"] != len(edges):
        raise RuntimeError("Q2 plant-link artifact drift")
    native=defaultdict(set)
    for r in records(a.native_distribution):
        native[r["accepted_plant_name_id"]].add(r["area_code_l3"])
    hosts=defaultdict(set)
    for r in records(a.host_csv):
        hosts[r["insect_species"]].add(r["accepted_plant_name_id"])
    level1=json.loads(a.level1_json.read_text(encoding="utf-8"))["units"]
    resource_island_codes={
      r["area_code_l3"] for r in records(a.q1_region_metrics) if r["class"]=="island"
    }
    if len(resource_island_codes)!=125:raise RuntimeError("Q1 island identification drift")

    island_cells=defaultdict(set)
    island_species=defaultdict(set)
    for r in edges:
        sp, island=r["species"],r["island_code"]
        if island not in resource_island_codes:raise RuntimeError("Non-island cell in Q2 edges")
        island_cells[sp].add(island)
        island_species[island].add(sp)
    if sum(map(len,island_cells.values())) != 5674:
        raise RuntimeError("Q2 introduced host edges do not reconstruct 5674 unique butterfly–island cells")

    species_rows=[]
    n_native=0
    byregion=Counter()
    for entry in m:
        sp=entry["species"]
        footprint=set()
        for h in hosts[sp]:
            footprint.update(native.get(h,()))
        native_len=int(entry["native_resource_units"])
        total_added=int(entry["introduced_added_units"])
        if len(footprint)!=native_len or native_len<=0:
            raise RuntimeError(f"Unmatched frozen native footprint: {sp} n={len(footprint)} expected={native_len}")
        groups=Counter(level1[x] for x in footprint if x in level1)
        if sum(groups.values())!=native_len:
            raise RuntimeError(f"Unmapped native WGSRPD region for {sp}")
        most,count=sorted(groups.items(),key=lambda x:(-x[1],x[0]))[0]
        fraction=count/native_len
        island_added=len(island_cells[sp])
        mainland_or_mixed_added=total_added-island_added
        if mainland_or_mixed_added<0:
            raise RuntimeError(f"Negative non-island added units for {sp}")
        byregion[most]+=1
        species_rows.append({
          "species":sp,"native_resource_units":native_len,
          "dominant_native_host_resource_level1":most,
          "dominant_native_host_resource_name":REGION_NAMES[most],
          "dominant_native_resource_fraction":fraction,
          "strict_dominance_70pct":int(fraction>=.7),
          "europe_share_of_native_resource":groups.get("1",0)/native_len,
          "global_added_resource_regions":total_added,
          "island_added_resource_regions":island_added,
          "other_added_resource_regions":mainland_or_mixed_added,
        })
    if sum(int(r["global_added_resource_regions"]) for r in species_rows)!=14553:
        raise RuntimeError("Global added opportunity mismatch")
    region_summary=[]
    for code in sorted(REGION_NAMES):
        sample=[r for r in species_rows if r["dominant_native_host_resource_level1"]==code]
        if not sample:continue
        region_summary.append({
            "region":REGION_NAMES[code],"level1":code,
            "eligible_butterflies":len(sample),
            "strict_native_host_resource_70pct_butterflies":
                sum(int(r["strict_dominance_70pct"]) for r in sample),
            "added_island_butterfly_region_units":sum(r["island_added_resource_regions"] for r in sample),
            "added_global_butterfly_region_units":sum(r["global_added_resource_regions"] for r in sample),
            "expanded_butterflies_with_any_island_resource":
                sum(r["island_added_resource_regions"]>0 for r in sample),
        })
    all_island=sum(r["added_island_butterfly_region_units"] for r in region_summary)
    if all_island!=5674:raise RuntimeError("Region contributions do not sum to Q2 island additions")

    empire={r["area_code_l3"]:r for r in records(a.empire_islands)}
    if len(empire)!=125:raise RuntimeError("Historical spatial overlay drift")
    region_by_species={r["species"]:r for r in species_rows}
    island_summary=[]
    for island, species in sorted(island_species.items()):
        historic=empire[island]
        n=len(species)
        european=sum(region_by_species[sp]["dominant_native_host_resource_level1"]=="1" for sp in species)
        europe_weighted=sum(region_by_species[sp]["europe_share_of_native_resource"] for sp in species)
        strict_eu=sum(region_by_species[sp]["dominant_native_host_resource_level1"]=="1"
                      and region_by_species[sp]["strict_dominance_70pct"] for sp in species)
        years=[float(historic[e+"_weighted_exposure_years"]) for e in
               ("Great_Britain","Spain","Portugal","Netherlands")]
        island_summary.append({
          "island_code":island,"level1":historic["level1"],
          "added_resource_butterflies":n,
          "europe_dominant_resource_butterflies":european,
          "europe_dominant_resource_fraction":european/n,
          "europe_strict_70pct_resource_butterflies":strict_eu,
          "mean_share_native_resources_in_europe":europe_weighted/n,
          "lenzner_polygon_coverage":float(historic["lenzner_polygon_coverage"]),
          "max_four_european_empire_years":max(years),
        })
    savecsv(a.output_dir/"beneficiary_species_native_resource_geography.csv",species_rows)
    savecsv(a.output_dir/"beneficiary_native_resource_regions.csv",region_summary)
    savecsv(a.output_dir/"beneficiary_island_geography.csv",island_summary)
    total_europe=sum(r["island_added_resource_regions"] for r in species_rows
                     if r["dominant_native_host_resource_level1"]=="1")
    european_strict=sum(r["island_added_resource_regions"] for r in species_rows
                        if r["dominant_native_host_resource_level1"]=="1"
                        and r["strict_dominance_70pct"])
    diagnostic={
      "schema":"chocho_q2_resource_beneficiary_native_host_geography_v0.1",
      "status":"RECONSTRUCTED_AND_GLOBAL_COUNTS_AUDITED",
      "population":{"eligible_butterflies":len(m),"total_native_resource_cells":sum(r["native_resource_units"] for r in species_rows),
                    "global_added_resource_cells":14553,"island_added_resource_cells":5674,
                    "islands_with_added_opportunity":len(island_summary)},
      "dominant_native_host_resource_regions_by_species":dict(byregion),
      "region_contributions":region_summary,
      "europe_beneficiaries":{
          "europe_dominant_species":byregion.get("1",0),
          "europe_dominant_new_island_resources":total_europe,
          "europe_dominant_fraction_of_all_new_island_resources":total_europe/5674,
          "europe_strict_70pct_new_island_resources":european_strict,
          "europe_strict_70pct_fraction_of_all_new_island_resources":european_strict/5674,
      },
      "cautions":[
         "The dominant native host-resource region is NOT the butterfly's native geographic range or origin.",
         "Classification is obtained from plant-native WCVP ranges and HOSTS links, and does not prove host use on a destination island.",
         "This is a non-causal descriptive stratification, without a history hypothesis test.",
         "Different native-resource regions have unequal numbers of sampled butterfly species, so raw cell counts are not rate-standardized.",
         "Native plant-region prevalence and introduced geography are correlated; this does not establish colonial or shipping origins."
      ],
    }
    if diagnostic["population"]["total_native_resource_cells"]!=26530:
        raise RuntimeError("Q2 native resource footprint reconciliation failed")
    (a.output_dir/"beneficiary_geography_audit.json").write_text(
        json.dumps(diagnostic,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(diagnostic,indent=2,ensure_ascii=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
