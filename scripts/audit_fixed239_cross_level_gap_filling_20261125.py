#!/usr/bin/env python3
"""Descriptive cross-level baseline/gap-filling diagnostic, not occupancy."""
import argparse,csv,hashlib,json
from pathlib import Path

BLOBS={
 "results/island_resource_gap_filling/by_species.csv":"a230f6d214c96f874a3ab184ee9823cc9023c1bb",
 "data/frozen/figure_sources/anthropogenic_species_metrics.csv":"9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
}
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(root,path):
    b=(root/path).read_bytes()
    h=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    require(h==BLOBS[path],"Frozen source blob mismatch "+path)
    with (root/path).open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))
def totals(rows):
    d={k:sum(int(r[k]) for r in rows) for k in ("cells","native","empty","added","unfilled")}
    require(d["native"]+d["empty"]==d["cells"] and d["added"]+d["unfilled"]==d["empty"],"Nested opportunity accounting")
    return {"n_butterflies":len(rows),**d,"current":d["native"]+d["added"],
        "gap_filling_fraction":d["added"]/d["empty"],"proportional_gain_relative_to_native":d["added"]/d["native"],
        "mean_native_regions":d["native"]/len(rows),"mean_current_regions":(d["native"]+d["added"])/len(rows)}
def gini(a):
    a=sorted(a);n=len(a)
    return sum((2*(i+1)-n-1)*x for i,x in enumerate(a))/(n*sum(a))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    by=read(a.repo,"results/island_resource_gap_filling/by_species.csv")
    metrics=read(a.repo,"data/frozen/figure_sources/anthropogenic_species_metrics.csv")
    island=[r for r in by if r["group"]=="island"]
    mainland=[r for r in by if r["group"]=="mainland"]
    require(len(island)==len(mainland)==len(metrics)==239,"Fixed panel drift")
    d={r["species"]:r for r in island}; m={r["species"]:r for r in metrics}
    require(len(d)==len(m)==239 and set(d)==set(m)=={r["species"] for r in mainland},"Taxon mismatch")
    require(all(int(r["cells"])==126 for r in island) and all(int(r["cells"])==237 for r in mainland),"Region counts")
    require(sum(int(r["native_resource_units"]) for r in metrics)==26530,"Global native resource total")
    q=sorted([{**r,"global_native":int(m[r["species"]]["native_resource_units"])} for r in island],
        key=lambda r:(r["global_native"],r["species"]))
    blocks=[q[:60],q[60:120],q[120:180],q[180:]]
    groups=[{"group":i+1,"minimum_global_native":g[0]["global_native"],"maximum_global_native":g[-1]["global_native"],**totals(g)} for i,g in enumerate(blocks)]
    require([g["native"] for g in groups]==[339,710,1238,2542] and [g["added"] for g in groups]==[644,1077,1749,2271],"Frozen quartiles")
    its,mts=totals(island),totals(mainland)
    require((its["native"],its["added"],mts["native"],mts["added"])==(4829,5741,21215,8686),"source totals")
    report={
      "schema":"chocho_cross_level_resource_baseline_and_gap_filling_v1","status":"DESCRIPTIVE_FROZEN_PANEL_ARITHMETIC_NOT_ESTABLISHMENT",
      "source_git_blobs":BLOBS,"islands":its,"mainlands":mts,"groups_by_global_native_resource_extent":groups,
      "upper_boundary_ties":[{"species":r["species"],"global_native":r["global_native"]} for r in q if r["global_native"]==169],
      "top_group_tie_sensitivity":{"both_tied_species_lower":totals([r for r in q if r["global_native"]>169]),"both_tied_species_upper":totals([r for r in q if r["global_native"]>=169])},
      "island_absolute_mean_top_minus_bottom":{"native":groups[-1]["mean_native_regions"]-groups[0]["mean_native_regions"],"current":groups[-1]["mean_current_regions"]-groups[0]["mean_current_regions"]},
      "island_relative_inequality_gini":{"native":gini([int(r["native"]) for r in island]),"current":gini([int(r["native"])+int(r["added"]) for r in island])},
      "interpretation_limits":["Fixed documented-host coverage, not feeding, realized range, abundance or establishment.",
        "Group definitions use global native resource extent (not butterfly dispersal ability or true ecological specialization); equal-count groups with deterministic alphabetical tie resolution.",
        "Gap-filling is pooled added/initially empty butterfly-region units, not the mean species ratio. No independence or test of causal port-mediated selection is assumed.",
        "The absolute top-bottom gap grows, but relative Gini declines: no claim of increased overall inequality. Post-result descriptive accounting, not a new independent biological test."]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("CROSS_LEVEL_GAP_ACCOUNTING_PASS",json.dumps(report),flush=True)
if __name__=="__main__":main()
