#!/usr/bin/env python3
"""Count *distinct introduced larval host plant taxa* on corrected islands and test modern ports.

Unlike the host-resource *fraction* screen, this counts unique accepted WCVP
introduced plant IDs within the pinned HOSTS butterfly panel. Source-preserved
Q2 v0.2 butterfly×host×island edges are deduplicated on plant×island, so hosts
serving many butterfly species are not counted multiple times.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from test_three_gate_port_resource_island_association import regression,group_moments

def read(p):
    with p.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open(encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    for field in ("q2-edges","q2-audit","q1-regions","q1-audit","wpi-regions","island-area","out-dir"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q1=json.loads(a.q1_audit.read_text(encoding="utf-8"))
    q2=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    if q1["groups"]["island"]["added_butterfly_region_units"]!=5741:
        raise RuntimeError("Q1 must be corrected 126-island geography")
    edges=read(a.q2_edges)
    if len(edges)!=11573 or q2["scope"]["introduced_butterfly_plant_island_edges"]!=11573:
        raise RuntimeError("Q2 source introduced host-edge count mismatch")
    pairs={(r["plant_id"],r["island_code"]) for r in edges}
    speciescells={(r["species"],r["island_code"]) for r in edges}
    if (len(pairs),len(speciescells),len({p for p,c in pairs}))!=(5611,5741,564):
        raise RuntimeError("WCVP×HOSTS introduction geography changed")
    plants=defaultdict(set)
    spcount=Counter()
    edgecount=Counter()
    for plant,code in pairs:plants[code].add(plant)
    for sp,code in speciescells:spcount[code]+=1
    for edge in edges:edgecount[edge["island_code"]]+=1
    regions={r["area_code_l3"]:r for r in read(a.q1_regions) if r["class"]=="island"}
    ports={r["area_code_l3"]:r for r in read(a.wpi_regions)}
    areas={r["code"]:r for r in read(a.island_area)}
    if len(regions)!=126 or len(plants)!=123:
        raise RuntimeError("Expected 126 strict islands, 123 active introduced-host ones")
    joined=[]
    for code,r in sorted(regions.items()):
        cont=int(r["contemporary"])
        if cont<=0:
            if code in plants:raise RuntimeError("Introduced plants despite no contemporary resource "+code)
            continue
        if code not in plants or code not in areas or code not in ports:
            raise RuntimeError("Missing source island "+code)
        source=n=len(plants[code])
        if spcount[code]!=int(r["added"]):
            raise RuntimeError("Butterfly resource count mismatch for "+code)
        w=int(ports[code]["WPI_port_count"])
        area=float(areas[code]["area_km2"])
        if area<=0:raise RuntimeError("Wrong island area")
        joined.append({
          "code":code,"level1":r["level1"],"region_name":r["region_name"],
          "native":int(r["native"]), "contemporary":cont,
          "added":int(r["added"]),"added_share":int(r["added"])/cont,
          "area_km2":area,"log_area":math.log(area),
          "ports":w,"has_port":int(w>0),"log2ports":math.log2(1+w),
          "introduced_larval_host_taxa":n,
          "log_introduced_host_taxa":math.log1p(n),
          "butterfly_host_island_edges":edgecount[code],
          "log_added":math.log1p(int(r["added"])),
        })
    if len(joined)!=123 or sum(x["introduced_larval_host_taxa"] for x in joined)!=5611:
        raise RuntimeError("Every plant-island pair must enter the island taxon count once")
    groups={
       "all_active_123":joined,
       "native_resource_positive":[r for r in joined if r["native"]>0],
       "excluding_pacific_antarctic":[r for r in joined if r["level1"] not in ("6","9")],
       "contemporary_at_least_75":[r for r in joined if r["contemporary"]>=75],
       "excluding_biggest_port_region":[r for r in joined
          if r["ports"]<max(v["ports"] for v in joined)],
    }
    outcome={
       "schema":"chocho_three_gate_port_vs_introduced_larval_host_taxa_v0.1",
       "status":"DESCRIPTIVE_MODERN_PORTS_VS_HOST_TAXA_NOT_INTRODUCTION_PATHWAY",
       "source":{"Q1":"fixed 126-region 239 butterfly WCVP HOSTS source",
          "Q2_edges":len(edges),
          "distinct_plant_island_pairs":len(pairs),
          "Q2_species_island_resources":len(speciescells),
          "global_plant_taxa":564,
          "public_port_inventory":"2023 WPI published mirror",
          "units":"WGSRPD3 strict island level3 regions, not individual physical islands"},
       "main_response":"ln(1 + distinct introduced WCVP larval host plant taxa in selected HOSTS butterfly panel)",
       "exposure":"log2(1+2023 indexed WPI count)",
       "seed":20261102,"permutations":4999,
       "baseline":{
           "active_islands":len(joined),
           "ports_positive":sum(r["ports"]>0 for r in joined),
           "mean_host_taxa_if_ports":group_moments(
               [r for r in joined if r["ports"]>0],"introduced_larval_host_taxa"),
           "mean_host_taxa_if_no_ports":group_moments(
               [r for r in joined if r["ports"]==0],"introduced_larval_host_taxa"),
           "naive_spearman_port_count_vs_introduced_host_taxa":{
              "rho":float(spearmanr([r["ports"] for r in joined],
                                     [r["introduced_larval_host_taxa"] for r in joined]).statistic),
              "p_nonspatial":float(spearmanr([r["ports"] for r in joined],
                                             [r["introduced_larval_host_taxa"] for r in joined]).pvalue)
           }
       },
       "cohorts":{},
       "limits":[
          "Modern WPI ports and modern naturalized larval host geography do not indicate historical shipping route or causal introduction.",
          "All plant taxa are restricted to the HOSTS network used by 239 butterflies. This is not total introduced flora diversity.",
          "Port absence in an index does not imply a plant-free or traffic-free island.",
          "Area, Level1 and native resource abundance cannot control climate, human population, colonial history, horticultural imports or sampling bias.",
          "Count responses within one WGSRPD3 region do not represent actual host suitability on each constituent island.",
          "An associated port count alone cannot identify direct maritime transport."
       ]
    }
    for i,(name,subset) in enumerate(groups.items()):
        source=regression(subset,"log_introduced_host_taxa",native_control=False,
                          exposure="log2ports",rngseed=20261102+i*109)
        adjusted=regression(subset,"log_introduced_host_taxa",native_control=True,
                            exposure="log2ports",rngseed=20261102+i*109+43)
        out2=regression(subset,"log_added",native_control=True,
                        exposure="log2ports",rngseed=20261102+i*109+71,permutations=False)
        outcome["cohorts"][name]={
           "n":len(subset),
           "indexed_port_regions":sum(r["has_port"] for r in subset),
           "host_taxon_count_primary":source,
           "host_taxon_native_baseline_adjusted":adjusted,
           "resource_added_count_secondary":out2,
        }
    write(a.out_dir/"distinct_introduced_host_taxa_and_port_by_island.csv",joined)
    (a.out_dir/"introduced_host_taxon_port_association_audit.json").write_text(
         json.dumps(outcome,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "baseline":outcome["baseline"],
        "cohorts":{name:{
          "n":x["n"],"primary_beta":x["host_taxon_count_primary"].get("coefficient"),
          "adjusted_beta":x["host_taxon_native_baseline_adjusted"].get("coefficient"),
          "secondary_added_beta":x["resource_added_count_secondary"].get("coefficient"),
          "permutation_p_within_level1_area":x["host_taxon_native_baseline_adjusted"].get("permutations",{}).get("within_level1_and_within_level_area_half",{}).get("two_sided_centered_p")
        } for name,x in outcome["cohorts"].items()}
    },indent=2,ensure_ascii=False))
if __name__=="__main__":main()
