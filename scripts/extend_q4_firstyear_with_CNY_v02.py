#!/usr/bin/env python3
"""Extend rather than replace the response-blind Q4 GBIF botanical pilot for corrected Canary Island CNY scope.

Original fixed 500 plant×island pilot and 219 dated index cells remain unchanged.
Add the first 3 SHA256-ranked CNY host plant×island pairs for the newly included
126th strict island, without selecting based on butterfly or plant observations.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import defaultdict,Counter

def rr(path):
    with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def writecsv(path,rrr):
    if not rrr:raise RuntimeError("Cannot write empty output")
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rrr[0]))
        w.writeheader();w.writerows(rrr)
def sha(p,i):
    return hashlib.sha256((p+"|"+i).encode()).hexdigest()
def main():
    p=argparse.ArgumentParser()
    for arg in ("corrected-q1","q2-v02-audit","q2-v02-edges",
                "original-500-plant-years","original-500-audit","out-dir"):
        p.add_argument("--"+arg,required=True,type=Path)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q1=json.loads(a.corrected_q1.read_text(encoding="utf8"))
    q2=json.loads(a.q2_v02_audit.read_text(encoding="utf8"))
    olda=json.loads(a.original_500_audit.read_text(encoding="utf8"))
    edges=rr(a.q2_v02_edges)
    previous=rr(a.original_500_plant_years)
    if (q1["groups"]["island"]["all_regions"],q1["groups"]["island"]["added_butterfly_region_units"])!=(126,5741):
        raise RuntimeError("Q1 CNY correction absent")
    if len(edges)!=q2["scope"]["introduced_butterfly_plant_island_edges"]:
        raise RuntimeError("Q2 v02 source edge count mismatch")
    if q2["scope"]["added_strict_island_cells"]!=5741 or len(previous)!=500 or olda["pilot_dated_pairs"]!=219:
        raise RuntimeError("Q2 or historic pilot frozen source drift")
    allpairs=defaultdict(set)
    cells=defaultdict(set)
    for row in edges:
        p,i,s=row["plant_id"],row["island_code"],row["species"]
        allpairs[(p,i)].add(s)
        cells[(s,i)].add(p)
    if len(cells)!=5741 or "CNY" not in {i for _,i in allpairs}:
        raise RuntimeError("CNY added butterfly host cell absent")
    old={(r["plant_id"],r["island_code"]) for r in previous}
    if len(old)!=500 or not old.issubset(allpairs):
        missing=sorted(old-set(allpairs))[:10]
        raise RuntimeError(f"Old Q4 pilot pairs invalidated: {missing}")
    new=[(p,i) for p,i in sorted(allpairs) if i=="CNY" and (p,i) not in old]
    if len(new)<3:raise RuntimeError("Not enough newly added CNY plants for fixed 3")
    new=sorted(new,key=lambda p:sha(*p))[:3]
    selected=old|set(new)
    all_rows=[]
    for p,i in sorted(allpairs):
        all_rows.append({"plant_id":p,"island_code":i,
                         "butterflies_gaining_new_resource":len(allpairs[(p,i)]),
                         "plant_island_sha256":sha(p,i),
                         "q4_pilot_20261102_v01":int((p,i) in old),
                         "q4_geo_corrected_canary_new":int((p,i) in new)})
    newrows=[x for x in all_rows if x["q4_geo_corrected_canary_new"]]
    sampled=[x for x in all_rows if (x["plant_id"],x["island_code"]) in selected]
    if len(newrows)!=3 or len(sampled)!=503:
        raise RuntimeError("Response-blind CNY expansion failed")
    writecsv(a.out_dir/"all_corrected_host_island_pairs.csv",all_rows)
    writecsv(a.out_dir/"new_CNY_host_island_pairs.csv",newrows)
    writecsv(a.out_dir/"continued_500_plus_CNY_503.csv",sampled)
    summary={
       "schema":"chocho_corrected_Q4_126_island_host_pairs_preserve_original_500_v0.2",
       "status":"CNY_EXTENDED_PILOT_3_NEW_UNDATED_PLANT_ISLANDS",
       "source":{
         "Q1_strict_island_WGSRPD3":126,
         "Q1_new_butterfly_resource_cells":5741,
         "Q2_introduced_butterfly_host_island_edges":len(edges),
         "Q2_distinct_intro_plant_taxa":q2["scope"]["distinct_intro_host_taxa"],
       },
       "corrected_population":{
         "introduced_host_plant_x_island_pairs":len(allpairs),
         "distinct_introduced_host_taxa":len({p for p,i in allpairs}),
         "islands_with_new_butterfly_resources":len({i for p,i in allpairs}),
         "canary_new_butterfly_cells":sum(i=="CNY" for s,i in cells),
         "canary_introduced_host_plant_pairs":sum(i=="CNY" for p,i in allpairs),
       },
       "fixed_pilot":{
          "original_pilot_pairs_preserved":len(old),
          "original_pilot_with_dated_GBIF_plants":olda["pilot_dated_pairs"],
          "new_response_blind_CNY_host_pairs":len(new),
          "new_pilot_pairs":len(sampled),
          "canary_plant_ids":[p for p,i in new],
       },
       "limits":[
         "The three CNY pairs were selected solely by source-locked SHA256, not by GBIF or butterfly historical results.",
         "Original 219 dated plant indexes remain in the unrevised 500-pair pilot; 3 new CNY pairs have not yet been queried.",
         "This corrects geographic opportunity classification without re-estimating old literature or historic port introduction dates.",
         "GBIF botanical first-record dates are sampling-dependent observation dates, not true plant introduction or butterfly colonization times."
       ],
    }
    (a.out_dir/"Q4_CNY_extension_audit.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
