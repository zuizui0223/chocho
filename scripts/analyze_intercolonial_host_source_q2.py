#!/usr/bin/env python3
"""Exploratory Q2-D: Are new butterfly-resource hosts native to other former colonies?

This is a possible plant-source channel, not proof of actual intra-imperial
transport, port arrival, original butterfly geography, or host use.
Thresholds are source-geography-based; never select by outcome direction.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import defaultdict,Counter
from pathlib import Path

EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")
METROPOLE_CODES={"GRB","SPA","POR","NET"}


def rows(path):
    with path.open(newline="",encoding="utf-8") as h:
        yield from csv.DictReader(h)


def dumpcsv(path,vals):
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=list(vals[0]))
        w.writeheader();w.writerows(vals)


def main():
    p=argparse.ArgumentParser()
    for k in ("q2-edges","native-plant-regions","global-empire-crosswalk",
              "island-empire-crosswalk","q2-original-audit","output-dir"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    original=json.loads(a.q2_original_audit.read_text(encoding="utf-8"))
    edges=list(rows(a.q2_edges))
    native=defaultdict(set)
    for rec in rows(a.native_plant_regions):
        native[rec["plant_id"]].add(rec["area_code_l3"])
    global_hist={r["area_code_l3"]:r for r in rows(a.global_empire_crosswalk)}
    dest_hist={r["area_code_l3"]:r for r in rows(a.island_empire_crosswalk)}
    if len(global_hist)!=369 or len(dest_hist)!=125:
        raise RuntimeError("Not a complete pinned global/source island historical crosswalk")
    scope=original["scope"]
    if len(edges)!=scope["introduced_butterfly_plant_island_edges"] or len({r["plant_id"] for r in edges})!=560:
        raise RuntimeError("Host edge source artifact does not equal Q2-A audited cohort")
    cells=defaultdict(set)
    for r in edges:
        cells[(r["species"],r["island_code"])].add(r["plant_id"])
    if len(cells)!=5674:
        raise RuntimeError("Introduced resource cells not reconstructed exactly")
    variants={
       "non_european_other_colonies":lambda r:r["level1"]!="1",
       "all_other_empire_territories":lambda r:True,
       "other_colonial_islands_only":lambda r:r["geographic_group"]=="island",
    }
    audits={}
    for variant,check in variants.items():
        donors=defaultdict(lambda:defaultdict(set))
        donor_regions=defaultdict(set)
        # Strict rule for native source: majority mapped geometry and at least 50
        # area-weighted historical years in that empire; home metropoles excluded.
        for host,regset in native.items():
            for code in regset:
                region=global_hist.get(code)
                if not region or code in METROPOLE_CODES:continue
                if float(region["lenzner_polygon_coverage"])<.50 or not check(region):continue
                for empire in EMPIRES:
                    if float(region[empire+"_weighted_exposure_years"])>=50:
                        donors[host][empire].add(code)
                        donor_regions[empire].add(code)
        per_island=defaultdict(lambda:{"n":0,"matches":Counter(),"host_edges":0})
        for (butterfly,island),hostids in cells.items():
            target=per_island[island]
            target["n"]+=1
            target["host_edges"]+=len(hostids)
            for empire in EMPIRES:
                # A native source elsewhere in empire's past territories is
                # only *consistent* with such a corridor, not a route observation.
                if any(donors[h][empire] for h in hostids):
                    target["matches"][empire]+=1
        island_output=[]
        for island,obj in sorted(per_island.items()):
            island_output.append({
                "island_code":island,
                "level1":dest_hist[island]["level1"],
                "new_butterfly_resource_cells":obj["n"],
                "introduced_host_edges":obj["host_edges"],
                **{e+"_native_supported_cells":obj["matches"][e] for e in EMPIRES}
            })
        if len(island_output)!=122 or sum(r["new_butterfly_resource_cells"] for r in island_output)!=5674:
            raise RuntimeError("Island/source cohort drift")
        totals={e:sum(r[e+"_native_supported_cells"] for r in island_output) for e in EMPIRES}
        filename=f"colonial_network_{variant}_island_source_summary.csv"
        dumpcsv(a.output_dir/filename,island_output)
        audits[variant]={
            "source_membership":"WCVP native region + GloNAF territory with >=50 weighted historic empire years, >=50% polygon overlap; 4 European metropoles excluded",
            "islands_with_new_resources":len(island_output),
            "introduced_host_taxa":len(native),
            "source_territory_regions_by_empire":{e:len(donor_regions[e]) for e in EMPIRES},
            "native_source_hosts_by_empire":{e:sum(bool(donors[h][e]) for h in native) for e in EMPIRES},
            "new_opportunity_cells_supported_by_native_territory_sources":totals,
            "new_opportunity_cell_source_support_share":{e:totals[e]/5674 for e in EMPIRES},
            "island_csv":filename
        }
    payload={
       "schema":"chocho_q2_other_colony_native_host_source_v0.1",
       "status":"EXPLORATORY_ALTERNATIVE_AFTER_METROPOLE_NULL",
       "source":{"WCVP":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
                 "Lenzner_GIS":"7aae72aadd39ae3827ece6935c6740dc0aace030",
                 "WGSRPD3":"52da7828aba9d461dd133c27b3bd7a4407161f54"},
       "scope":scope,
       "variants":audits,
       "limits":[
         "This is an exploratory test developed after the simpler metropole-native source hypothesis did not have convincing positive support.",
         "Native plant distribution within an empire's territories does not imply human transport or intra-imperial colonization.",
         "Host records are literature/observation-derived and are not demonstrated local use.",
         "Donor identities and former empire memberships overlap; percentages for empires cannot be added.",
         "French and other empires are absent from this published four-empire source.",
         "WDG/GLONAF polygon crosswalk is incomplete or uncertain for many very small oceanic islands.",
       ]
    }
    (a.output_dir/"intercolonial_host_source_audit.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2,ensure_ascii=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
