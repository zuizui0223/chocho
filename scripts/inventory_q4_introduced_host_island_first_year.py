#!/usr/bin/env python3
"""Audited Q4 inventory of introduced host plant × island opportunity cells.

Input: pinned Q2 artifact with 11,455 butterfly–host–island edges and 5,674
butterfly × WGSRPD3 additions. Distinct plant-island observations are the
unit of the earliest-GBIF-botanical-record audit, not butterfly×island events.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,collections
from pathlib import Path

def read(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def write(path,records):
    if not records:raise RuntimeError("empty")
    with path.open("w",newline="",encoding="utf-8") as f:
        out=csv.DictWriter(f,fieldnames=list(records[0]))
        out.writeheader();out.writerows(records)
def digest(s):return hashlib.sha256(s.encode("utf-8")).hexdigest()
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--edges",type=Path,required=True)
    p.add_argument("--audit",type=Path,required=True)
    p.add_argument("--island-metrics",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--sample-size",type=int,default=500)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    edges=read(a.edges)
    audit=json.loads(a.audit.read_text())
    island={r["area_code_l3"]:r for r in read(a.island_metrics) if r["class"]=="island"}
    if len(island)!=125 or len(edges)!=11455 or audit["scope"]["added_strict_island_cells"]!=5674:
        raise RuntimeError("Q1/Q2 frozen source mismatch")
    perplant=collections.defaultdict(set)
    percell=collections.defaultdict(set)
    per_pair=collections.defaultdict(set)
    for r in edges:
        sp,isle,plant=r["species"],r["island_code"],r["plant_id"]
        if isle not in island:raise RuntimeError("Unknown island")
        perplant[plant].add(isle)
        percell[(sp,isle)].add(plant)
        per_pair[(plant,isle)].add(sp)
    if len(perplant)!=560 or len(percell)!=5674 or len({b for b,i in percell})>239:
        raise RuntimeError(f"Q2 population mismatch {len(perplant)=} {len(percell)=}")
    pairs=[]
    for (plant,isle),butterflies in sorted(per_pair.items()):
        pairs.append({"plant_id":plant,"island_code":isle,
                      "level1":island[isle]["level1"],
                      "native_baseline_butterfly_cells":island[isle]["native"],
                      "added_butterfly_resource_cells_island":island[isle]["added"],
                      "butterflies_benefiting_from_this_plant":len(butterflies),
                      "butterfly_resource_species_overlap_from_multi_hosts":int(len(butterflies)>1),
                      "host_island_sha256":digest(plant+"|"+isle)})
    counts=collections.Counter(r["island_code"] for r in pairs)
    priorities={}
    for x in pairs:
        key=x["island_code"]
        priorities.setdefault(key,[]).append(x)
    fixed=[]
    for isle, items in sorted(priorities.items()):
        fixed.extend(sorted(items,key=lambda x:x["host_island_sha256"])[:3])
    selected={(r["plant_id"],r["island_code"]) for r in fixed}
    for r in sorted(pairs,key=lambda x:x["host_island_sha256"]):
        if len(selected)>=max(a.sample_size,len(fixed)):break
        selected.add((r["plant_id"],r["island_code"]))
    sample=[r for r in pairs if (r["plant_id"],r["island_code"]) in selected]
    sample.sort(key=lambda r:(r["island_code"],r["host_island_sha256"]))
    write(a.output_dir/"all_introduced_host_island_pairs.csv",pairs)
    write(a.output_dir/"fixed_first_year_sample_pairs.csv",sample)
    size_by_island=collections.Counter(r["island_code"] for r in sample)
    distribution=sorted(counts.values())
    summary={
      "schema":"chocho_q4_host_island_first_record_inventory_v0.1",
      "status":"AUDITED_SOURCE_ONLY_NO_GBIF_DATES_YET",
      "frozen_input":{"butterfly_host_island_edges":len(edges),
                      "new_butterfly_island_resource_cells":len(percell),
                      "distinct_introduced_host_plant_ids":len(perplant),
                      "islands_with_added_resources":len(priorities)},
      "analysis_unit":{
        "distinct_introduced_host_x_island_pairs":len(per_pair),
        "host_island_pairs_by_biogeographic_level1":dict(collections.Counter(r["level1"] for r in pairs)),
        "unique_host_islands_from_one_host":sum(len(isles)==1 for isles in perplant.values()),
        "median_plant_islands":sorted(map(len,perplant.values()))[len(perplant)//2],
        "min_host_pairs_per_island":min(distribution),
        "median_host_pairs_per_island":distribution[len(distribution)//2],
        "max_host_pairs_per_island":max(distribution),
      },
      "fixed_pilot":{
         "sample_size":len(sample),
         "islands_in_pilot":len(size_by_island),
         "minimum_three_pairs_per_island":len(fixed),
         "selection":"3 deterministic SHA256 smallest host-island pairs per island plus SHA256 global until 500; do not select by GBIF first year, port data or butterfly detection",
         "no_overlapping_pairs":len(sample)==len(selected),
      },
      "unresolved": [
        "WCVP accepted plant ID → unambiguous GBIF accepted species key mapping.",
        "GBIF first indexed year per plant × entire WGSRPD3 island, potentially left-censored.",
        "Era-specific all-botanical-observation effort by the same island geometry.",
        "Species × island historic butterfly first records and observation effort over same years."
      ],
      "claims":["Source host-island pairs are WCVP contemporary introduced-host potential resources, not dated introductions.",
                "GBIF non-detections cannot be imputed as plant absence or unused opportunity.",
                "The pilot's 500 sampled pairs cannot establish coverage of the entire pair panel."]
    }
    (a.output_dir/"first_year_inventory_audit.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=="__main__":main()
