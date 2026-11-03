#!/usr/bin/env python3
"""Independent botanical validation: all WCVP accepted introduced plants × 2023 WPI ports.

The response is all introduced plants, not selected HOSTS butterfly resource
taxa. Pinned WCVP v13, official 126 WGSRPD3 island polygons, fixed WPI snapshot.
Descriptive association; marine introductions, causal path and ship arrivals
are not inferred.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from analyze_wgsrpd3_island_area_adjusted import area_by_code
from test_three_gate_port_resource_island_association import regression,group_moments

def load(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))
def csvout(path, rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    for field in ("all-flora","all-flora-audit","q1-regions","q1-audit","wpi-regions",
                  "wgsrpd3-geojson","out-dir"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q1=json.loads(a.q1_audit.read_text(encoding="utf-8"))
    if len(q1["classification"]["strict_island_codes"])!=126 or "CNY" not in q1["classification"]["strict_island_codes"]:
        raise RuntimeError("Wrong Q1 island geometry census")
    flora={r["area_code_l3"]:r for r in load(a.all_flora)}
    ports={r["area_code_l3"]:r for r in load(a.wpi_regions)}
    q1reg={r["area_code_l3"]:r for r in load(a.q1_regions)}
    src={r["metric"]:int(r["value"]) for r in load(a.all_flora_audit)}
    if len(flora)!=126 or src["strict_island_regions"]!=126:
        raise RuntimeError("All-flora botanical census not 126 regions")
    if src["conflicting_native_alien_statuses"]!=0:
        raise RuntimeError("Ambiguous native and introduced status in WCVP")
    area=area_by_code(a.wgsrpd3_geojson)
    data=[]
    for code,r in sorted(flora.items()):
        if code not in q1reg or q1reg[code]["class"]!="island" or ports[code]["geographic_group"]!="island":
            raise RuntimeError("Geography class mismatch: "+code)
        km2=float(area[code])
        if km2<=0:raise RuntimeError("Bad geographic area: "+code)
        native=int(r["native_flora_species_count"])
        introduced=int(r["introduced_flora_species_count"])
        total=int(r["contemporary_flora_species_count"])
        port=int(ports[code]["WPI_port_count"])
        if native+introduced!=total:raise RuntimeError("Floristic census not additive for "+code)
        data.append({
          "code":code,"region_name":q1reg[code]["region_name"],
          "level1":q1reg[code]["level1"],
          "area_km2":km2,"log_area":math.log(km2),
          "native":native,"introduced":introduced,"total":total,
          "wpi_port_count":port,"ports":port,"has_port":int(port>0),
          "log2ports":math.log2(1+port),
          "log_introduced_flora":math.log1p(introduced),
          "intro_flora_fraction":(introduced/total if total>0 else None),
          "contemporary_butterfly_resource_units":int(q1reg[code]["contemporary"]),
        })
    if len(data)!=126:raise RuntimeError("Joined island sample not 126")
    if sum(x["introduced"] for x in data)!=src["introduced_species_island_records"]:
        raise RuntimeError("Plant count sum vs full WCVP audit disagrees")
    active=[x for x in data if x["contemporary_butterfly_resource_units"]>0]
    if len(active)!=123:raise RuntimeError("Source 123 butterfly-active island sample drift")
    groups={
      "all_126_islands":data,
      "same_123_resource_active_islands":active,
      "native_flora_positive":[x for x in data if x["native"]>0],
      "exclude_pacific_antarctic":[x for x in data if x["level1"] not in ("6","9")],
      "exclude_largest_indexed_port":[x for x in data if x["ports"]<max(y["ports"] for y in data)],
    }
    result={
      "schema":"chocho_all_WCVP_introduced_flora_vs_island_ports_v0.1",
      "status":"INDEPENDENT_ALL_FLORA_ASSOCIATION_NOT_CAUSAL_TRANSPORT",
      "input":{"wcvp_source_commit":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
               "wgsrpd_source_commit":"52da7828aba9d461dd133c27b3bd7a4407161f54",
               "all_flora_accepted_species":src["WCVP_accepted_species_global"],
               "regions":126,
               "native_flora_species_region_occurrences":src["native_species_island_records"],
               "introduced_flora_species_region_occurrences":src["introduced_species_island_records"],
               "nonconflicting_wcvp_native_introduced_assignment":True},
      "response":"ln(1+all WCVP introduced accepted plant species count per WGSRPD3 region)",
      "exposure":"log2(1+modern WPI port count)",
      "controls":"log official WGSRPD3 polygon area, botanical Level1 fixed effects; native flora richness sensitivity",
      "seed":20261102,
      "permutations":4999,
      "baseline":{
        "n":len(data),"indexed_port_regions":sum(x["ports"]>0 for x in data),
        "mean_intro_flora_with_port":group_moments([x for x in data if x["ports"]>0],"introduced"),
        "mean_intro_flora_without_port":group_moments([x for x in data if x["ports"]==0],"introduced"),
        "naive_spearman":float(spearmanr([x["ports"] for x in data],[x["introduced"] for x in data]).statistic),
        "naive_nonspatial_p":float(spearmanr([x["ports"] for x in data],[x["introduced"] for x in data]).pvalue)
      },
      "cohorts":{},
      "limits":[
        "Count all WCVP introduced plants, independently of HOSTS; no proof any plant is a larval host or suitable food.",
        "2023 WPI port counts do not reconstruct historical shipping volume, plant import routes or original introduction date.",
        "Current WCVP introduced status has unequal documentation quality among regions.",
        "Area and botanical biogeographic region do not fully control population, land use, economy, botanical effort or historical exposure.",
        "Spatial dependence remains and permutations are descriptive not causal randomization.",
      ],
    }
    for i,(name,subset) in enumerate(groups.items()):
        base=regression(subset,"log_introduced_flora",native_control=False,
                        exposure="log2ports",rngseed=20261102+i*117)
        n=regression(subset,"log_introduced_flora",native_control=True,
                     exposure="log2ports",rngseed=20261102+i*117+63)
        result["cohorts"][name]={"n":len(subset),"base":base,"native_flora_adjusted":n}
    csvout(a.out_dir/"all_flora_intro_native_wpi_islands.csv",data)
    (a.out_dir/"all_flora_port_association_audit.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"input":result["input"],"baseline":result["baseline"],
        "cohorts":{k:{"n":v["n"],"base_beta":v["base"].get("coefficient"),
                       "native_beta":v["native_flora_adjusted"].get("coefficient"),
                       "native_p_lev_area":v["native_flora_adjusted"].get("permutations",{}).get("within_level1_and_within_level_area_half",{}).get("two_sided_centered_p")}
                   for k,v in result["cohorts"].items()}},indent=2,ensure_ascii=False))
if __name__=="__main__":main()
