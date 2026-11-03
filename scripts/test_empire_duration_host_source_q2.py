#!/usr/bin/env python3
"""Q2-C: Does historic empire exposure covary with empire-home native-host opportunity?

Uses equal island-weighted 4-empire crossed panels and two-way-de-meaned
descriptive alignment. Island-wide and empire-wide baselines are removed; a
within-biogeographic permutation shuffles *entire* empire-history vectors,
preserving empire co-history and each donor's native-host availability scores.
An association does not demonstrate historical transport, occupied butterflies,
colonial causes, or harbor introduction.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")
SEED=20261102


def read_csv(path):
    with path.open(encoding="utf-8",newline="") as h:
        yield from csv.DictReader(h)


def doubly_center(mat):
    cols=len(mat[0])
    row_means=[sum(row)/cols for row in mat]
    col_means=[sum(row[j] for row in mat)/len(mat) for j in range(cols)]
    grand=sum(row_means)/len(row_means)
    return [[x-row_means[i]-col_means[j]+grand for j,x in enumerate(row)]
            for i,row in enumerate(mat)]


def slope(xs,ys):
    a=doubly_center(xs); b=doubly_center(ys)
    cross=sum(x*y for ax,by in zip(a,b) for x,y in zip(ax,by))
    norm=sum(x*x for row in a for x in row)
    return cross/norm if norm>1e-12 else None


def permute(xs, ys, blocks, draws=4999):
    obs=slope(xs,ys)
    groups=defaultdict(list)
    for i,b in enumerate(blocks):groups[b].append(i)
    movable=[ind for ind in groups.values()
             if len({tuple(round(v,6) for v in xs[i]) for i in ind})>1]
    if obs is None or not movable:
        return {"status":"UNIDENTIFIABLE_OR_NO_WITHIN_STRATUM_VARIATION",
                "observed_slope":obs,
                "effective_islands":sum(map(len,movable))}
    rng=random.Random(SEED)
    draws_value=[]
    for _ in range(draws):
        new=[r[:] for r in xs]
        for inds in movable:
            vals=[new[i][:] for i in inds]
            rng.shuffle(vals)
            for i,row in zip(inds,vals):new[i]=row
        fitted=slope(new,ys)
        if fitted is not None:draws_value.append(fitted)
    draws_value.sort()
    mean_null=sum(draws_value)/len(draws_value)
    return {
       "status":"DESCRIPTIVE_PERMUTATION",
       "observed_slope_fraction_per_century":obs,
       "null_mean_slope":mean_null,
       "observed_minus_null_slope":obs-mean_null,
       "null_interval_q025_q975":[draws_value[int(.025*(len(draws_value)-1))],
                                draws_value[int(.975*(len(draws_value)-1))]],
       "p_high":(1+sum(x>=obs-1e-12 for x in draws_value))/(len(draws_value)+1),
       "p_two":(1+sum(abs(x-mean_null)>=abs(obs-mean_null)-1e-12 for x in draws_value))/(len(draws_value)+1),
       "draws":len(draws_value),
       "effective_islands":len({i for group in movable for i in group}),
       "exchangeable_blocks":len(movable),
    }


def analyze(sample, min_years):
    xs=[];ys=[]
    for r in sample:
        history=[max(0,r["years"][e]) if r["years"][e]>=min_years else 0 for e in EMPIRES]
        xs.append([v/100 for v in history])
        ys.append([r["supports"][e] for e in EMPIRES])
    native_possible=sum(int(any(v>0 for v in row)) for row in ys)
    area_logs=sorted(math.log(r["area_km2"]) for r in sample)
    quartiles=[area_logs[int((len(area_logs)-1)*q)] for q in [.25,.5,.75]]
    def area_cat(v):return sum(math.log(v)>x for x in quartiles)
    groupings={
       "globally_unstratified":["all"]*len(sample),
       "within_biogeographic_level1":[r["level1"] for r in sample],
       "within_level1_and_global_area_quartile":[r["level1"]+"|"+str(area_cat(r["area_km2"])) for r in sample]
    }
    stats={k:permute(xs,ys,grp) for k,grp in groupings.items()}
    years_total={e:sum(r["years"][e] for r in sample) for e in EMPIRES}
    return {
       "islands":len(sample),
       "total_added_resource_cells":sum(r["new_cells"] for r in sample),
       "islands_with_any_metropole_native_source":native_possible,
       "islands_with_any_historic_empire":sum(any(r["years"][e]>=min_years for e in EMPIRES) for r in sample),
       "total_empire_years_geography_weighted":years_total,
       "area_quartiles_ln_km2":quartiles,
       "history_exposure_threshold_yrs":min_years,
       "test":stats,
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("host-source-islands","host-source-audit","empire-islands","empire-overlay-audit",
                 "area-metrics","output-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--source-type", choices=[
        "metropole", "non_european_other_colonies", "all_other_empire_territories",
        "other_colonial_islands_only"
    ], default="metropole")
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    hosts={r["island_code"]:r for r in read_csv(a.host_source_islands)}
    hist={r["area_code_l3"]:r for r in read_csv(a.empire_islands)}
    area={r["code"]:r for r in read_csv(a.area_metrics)}
    host_audit=json.loads(a.host_source_audit.read_text())
    overlay_audit=json.loads(a.empire_overlay_audit.read_text())
    if host_audit["scope"]["added_strict_island_cells"]!=5674 or overlay_audit["frozen_island_regions"]!=125:
        raise RuntimeError("Precommitted source panel changed")
    panel=[]
    for code,h in sorted(hosts.items()):
        if code not in hist or code not in area:raise RuntimeError(f"Unmapped donor island: {code}")
        n=int(h["new_butterfly_resource_cells"])
        if n<=0:continue
        source=hist[code]
        panel.append({
          "code":code,
          "level1":source["level1"],
          "cover":float(source["lenzner_polygon_coverage"]),
          "area_km2":float(area[code]["area_km2"]),
          "new_cells":n,
          "years":{e:float(source[e+"_weighted_exposure_years"]) for e in EMPIRES},
          "supports":{e:int(h[e+"_native_supported_cells"])/n for e in EMPIRES},
        })
    if len(panel)!=122:raise RuntimeError(f"Frozen 122 active island resources expected: {len(panel)}")
    def subset(cov, mincell=1, ex_europe=False):
        return [x for x in panel if x["cover"]>=cov and x["new_cells"]>=mincell
                and (not ex_europe or x["level1"]!="1")]
    variants={
      "primary_cover_ge80_min50yr":analyze(subset(.8),50),
      "full_history_cover_ge80_min1yr":analyze(subset(.8),1),
      "sensitivity_cover_ge50_min50yr":analyze(subset(.5),50),
      "sensitivity_cover_ge80_min50yr_ncells25":analyze(subset(.8,25),50),
      "sensitivity_non_european_cover80_50yr":analyze(subset(.8,1,True),50)
    }
    output={
      "schema":"chocho_empire_duration_metropole_host_alignment_v0.1" if a.source_type=="metropole"
          else "chocho_empire_duration_other_colonial_host_alignment_v0.1",
      "source_type":a.source_type,
      "status":"OBSERVATIONAL_COLONIAL_HISTORY_SOURCE_ASSOCIATION_TEST",
      "units":"Equal island-weighted array across four European historical empires, not butterfly abundances",
      "source_panel":host_audit["scope"],
      "all_active_islands":len(panel),
      "primary_coverage_ge80":len(subset(.8)),
      "source_refs":{"lenzner":overlay_audit["source"],"wcvp_hosts":host_audit["source_commits"]},
      "variants":variants,
      "controls":{
        "island_baseline":"Subtract per-island four-metropole source-opportunity mean",
        "metropole_baseline":"Subtract per-metropole mean across included islands",
        "permutation":"Permute full empire-year vectors among islands, holding the native-source opportunity matrix fixed, within strata",
        "duration":"Ignore a historical empire's fractional area-weighted occupancy shorter than 50 years for primary",
      },
      "caution":[
        "Descriptive cross-classified correlation and restricted exchangeability null; spatial and unobserved history remain uncontrolled.",
        "The published historical data covers four European empires, not France or every historical colonizer.",
        "Native plant source regions or native populations within other former colonies do not establish import origin, pathway, period or port involvement.",
        "Species gain is potential larval-host resource geography, not butterfly colonization or established population.",
        "The four source fractions overlap, since host species may have native distributions in more than one metropole.",
      ]
    }
    filename="empire_duration_source_alignment_audit.json" if a.source_type=="metropole" else (
        a.source_type+"_duration_source_alignment_audit.json")
    (a.output_dir/filename).write_text(
       json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(output,indent=2,ensure_ascii=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
