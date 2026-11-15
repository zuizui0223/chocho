#!/usr/bin/env python3
"""Source-locked geographic sensitivity of potential island larval-host gains.

Filters are on the original botanical WGSRPD3 region's *native HOSTS resource*
count, broad Level-1 geographic group, or actual WGSRPD3 polygon area.
These are not field-verified butterfly occupancy or climatic habitat filters.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path

SOURCES={
    "regions":"results/island-mainland-resource-q1/region_metrics.csv",
    "areas":"results/island_resource_gap_filling/areas.csv",
    "original_pooled":"results/island_resource_gap_filling/pooled.csv",
}
def read(path):
    b=path.read_bytes()
    with path.open(newline="",encoding="utf8") as f:rows=list(csv.DictReader(f))
    return rows,hashlib.sha256(b).hexdigest()

def metrics(rows):
    n=len(rows)
    N=sum(x["native"] for x in rows)
    A=sum(x["added"] for x in rows)
    T=239*n
    return {"n_regions":n,"butterfly_region_cells":T,"native_covered":N,
        "introduced_added":A,"native_uncovered":T-N,
        "native_coverage_pct":100*N/T if T else None,
        "fill_of_native_gaps_pct":100*A/(T-N) if T>N else None,
        "relative_added_over_native_pct":100*A/N if N else None,
        "share_original_island_added_pct":None}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo-root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir",type=Path,default=Path("results/island_resource_extreme_sensitivity_20261115"))
    a=parser.parse_args()
    raw={};hashes={}
    for key,path in SOURCES.items():
        raw[key],hashes[key]=read(a.repo_root/path)
    area={x["area_code_l3"]:float(x["area_km2"]) for x in raw["areas"]}
    if len(area)!=369:raise ValueError("Original 369 WGSRPD3 area polygons changed")
    data=[]
    for x in raw["regions"]:
        if x["area_code_l3"] not in area:raise ValueError("Missing source WGSRPD3 polygon "+x["area_code_l3"])
        row=dict(x)
        row["native"]=int(x["native"])
        row["added"]=int(x["added"])
        row["contemporary"]=int(x["contemporary"])
        row["area_km2"]=area[x["area_code_l3"]]
        if row["contemporary"]!=row["native"]+row["added"]:
            raise ValueError("Regional resource non-additivity "+row["area_code_l3"])
        data.append(row)
    if len(data)!=369 or len({x["area_code_l3"] for x in data})!=369:
        raise ValueError("Original region cohort changed")
    island=[x for x in data if x["class"]=="island"]
    mainland=[x for x in data if x["class"]=="mainland"]
    mixed=[x for x in data if x["class"]=="mixed"]
    if not (len(island)==126 and len(mainland)==237 and len(mixed)==6):
        raise ValueError(f"Original 126/237/6 WGSRPD3 strict classification changed: "
                         f"{len(island)}/{len(mainland)}/{len(mixed)}")
    original={x["geography"]:x for x in raw["original_pooled"]}
    observed_island=metrics(island);observed_mainland=metrics(mainland)
    for group,actual in (("island",observed_island),("mainland",observed_mainland)):
        for k,l in (("cells","butterfly_region_cells"),("native","native_covered"),("added","introduced_added")):
            if int(float(original[group][k]))!=actual[l]:
                raise ValueError(f"Published pooled {group}/{k} changed")
    assert (observed_island["native_covered"],observed_island["introduced_added"],
            observed_mainland["native_covered"],observed_mainland["introduced_added"])==(4829,5741,21215,8686)

    filters={
        "all_126_strict_islands":lambda r:True,
        "exclude_all_zero_native_resource_regions":lambda r:r["native"]>=1,
        "native_resource_regions_ge_5":lambda r:r["native"]>=5,
        "native_resource_regions_ge_10":lambda r:r["native"]>=10,
        "native_resource_regions_ge_20":lambda r:r["native"]>=20,
        "native_resource_regions_ge_30":lambda r:r["native"]>=30,
        "exclude_antarctic_Level1_9":lambda r:r["level1"]!="9",
        "whole_WGSRPD3_region_area_ge_1000_km2":lambda r:r["area_km2"]>=1000,
        "whole_WGSRPD3_region_area_ge_10000_km2":lambda r:r["area_km2"]>=10000,
    }
    panels={k:metrics([x for x in island if condition(x)]) for k,condition in filters.items()}
    for v in panels.values():
        v["share_original_island_added_pct"]=100*v["introduced_added"]/observed_island["introduced_added"]
        v["original_mainland_gain_pct_for_context_not_matched_comparison"]=(observed_mainland["relative_added_over_native_pct"])
    zero=[x for x in island if x["native"]==0]
    if len(zero)!=13 or sum(x["added"] for x in zero)!=432:
        raise ValueError("Original no-native-host regions changed")
    report={
        "analysis":"original_126_strict_insular_wgsrpd3_native_host_scarcity_20261115",
        "source_sha256":hashes,
        "butterflies_per_region":239,
        "original_strict_islands":observed_island,
        "original_strict_mainland":observed_mainland,
        "geography_filters_are_nested_descriptive_sensitivities_not_independent_tests":panels,
        "original_zero_native_regions":[
            {"code":x["area_code_l3"],"name":x["region_name"],"added":x["added"],"level1":x["level1"]}
            for x in sorted(zero,key=lambda z:z["area_code_l3"])
        ],
        "zero_native_island_added_total":sum(x["added"] for x in zero),
        "interpretation":[
            "The +118.9% gain is not solely a zero-baseline or Antarctic Level1 region artefact",
            "Higher native-host thresholds select botanically richer regions and can change ecological and geographical support",
            "The 239 HOSTS butterflies are potential consumers of listed plants; zero native-host index does NOT establish actual butterfly absence",
            "WGSRPD3 polygon area >=1000 km2 is not equal to island land area, nor evidence of resident butterfly population",
            "The filtered island panels are not matched to a filtered mainland panel; comparisons to the old mainland ratio are context only",
            "The original insular gap-filling fraction is not generally higher than mainland, even if relative native-based gain is higher",
            "The available species-paired source independently tests host family breadth and multiple-candidate resource heterogeneity",
            "No climatic, historical establishment, breeding, host-quality, or ship-route causality is inferred"
        ],
    }
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"summary.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    with (a.out_dir/"island_regions.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=["area_code_l3","region_name","level1","area_km2","native","added","contemporary"])
        writer.writeheader()
        for x in island:writer.writerow({k:x[k] for k in writer.fieldnames})
    print(json.dumps({"island":observed_island,"mainland":observed_mainland,
                      "all_pooled_filters":{k:{"regions":v["n_regions"],
                         "gain_pct":v["relative_added_over_native_pct"],
                         "fill_pct":v["fill_of_native_gaps_pct"]}
                        for k,v in panels.items()},
                      "zero_native_regions":len(zero),"zero_native_added":sum(x["added"] for x in zero)
                     },ensure_ascii=False))
if __name__=="__main__":
    main()
