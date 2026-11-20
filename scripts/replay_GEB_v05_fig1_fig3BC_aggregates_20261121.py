#!/usr/bin/env python3
"""Independent standard-library GEB Fig1 and Fig3B/C AGGREGATE-table replay.

This can reproduce native/introduced known-host resource counts and
butterfly breadth rank patterns from three SHA-pinned derived tables.
Fig3A plant contribution is NOT reproduced without upstream original
species×region×host-plant fractional-credit edges; preserve a separate
original numerical receipt, but never call it an independent replay.

Not a substitute for source WCVP/HOSTS extraction or an anonymous data deposit.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path

PIN={
 "regions":"4aba89565354edc6c44d566523864cec66a30cbe",
 "butterflies":"a230f6d214c96f874a3ab184ee9823cc9023c1bb",
 "metrics":"9f83d2d9ac8e46f1cfd4ded572c7a023273503dd",
}
EXPECTED_GROUP={
 "island":(126,4829,10570,5741),
 "mainland":(237,21215,29901,8686),
 "mixed":(6,486,612,126)
}
def read(path,sha):
    raw=path.read_bytes()
    git=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if git!=sha:raise AssertionError("Source-pinned aggregate Git blob drift: "+path.name)
    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
def rank(x):
    order=sorted(range(len(x)),key=lambda k:x[k])
    out=[0.0]*len(x)
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and x[order[j]]==x[order[i]]:j+=1
        midpoint=(i+j+1)/2
        for idx in order[i:j]:out[idx]=midpoint
        i=j
    return out
def rho(a,b):
    x=rank(a);y=rank(b)
    mx=sum(x)/len(x);my=sum(y)/len(y)
    numer=sum((v-mx)*(w-my) for v,w in zip(x,y))
    denom=math.sqrt(sum((v-mx)**2 for v in x)*sum((w-my)**2 for w in y))
    if not denom:raise AssertionError("Rank correlation degenerate")
    return numer/denom
def run(regions,butterflies,metrics):
    q=read(regions,PIN["regions"])
    per=read(butterflies,PIN["butterflies"])
    m=read(metrics,PIN["metrics"])
    if len(q)!=369 or len({r["area_code_l3"] for r in q})!=369:
        raise AssertionError("Expected 369 unique WGSRPD3 region aggregate records")
    grouped={}
    for group,(count,native,current,added) in EXPECTED_GROUP.items():
        a=[x for x in q if x["class"]==group]
        got=(len(a),sum(int(x["native"]) for x in a),
             sum(int(x["contemporary"]) for x in a),
             sum(int(x["added"]) for x in a))
        if got!=(count,native,current,added):
            raise AssertionError(f"Fig1 {group} source counts mismatch: {got}")
        if any(int(x["native"])+int(x["added"])!=int(x["contemporary"]) for x in a):
            raise AssertionError(f"Fig1 {group} region resource counts not additive")
        obj={"regions":count,"native":native,"contemporary":current,"added":added}
        if group!="mixed":
            total=239*count
            obj.update({"native_coverage_pct":100*native/total,
              "relative_increase_pct":100*added/native,
              "missing_resource_gap_filled_pct":100*added/(total-native),
              "introduced_attributable_contemporary_pct":100*added/current})
        grouped[group]=obj
    assert sum(x["native"] for x in grouped.values())==26530
    assert sum(x["contemporary"] for x in grouped.values())==41083
    assert sum(x["added"] for x in grouped.values())==14553
    if len(m)!=239 or len(per)!=478:
        raise AssertionError("Expected fixed 239 butterfly metrics and two geographic responses")
    lookup={r["species"]:r for r in m}
    if len(lookup)!=239 or any(float(v["host_family_count"])<1 for v in m):
        raise AssertionError("Duplicated or invalid original butterfly resource metrics")
    isle=[r for r in per if r["group"]=="island"]
    mainland=[r for r in per if r["group"]=="mainland"]
    if len(isle)!=239 or len(mainland)!=239 or set(x["species"] for x in isle)!=set(lookup) or set(x["species"] for x in mainland)!=set(lookup):
        raise AssertionError("239 source butterfly IDs must occur exactly once in each geography")
    assert sum(int(x["added"]) for x in isle)==5741
    assert sum(int(x["added"]) for x in mainland)==8686
    assert sum(int(x["introduced_added_units"]) for x in m)==14553
    rows=[{"name":z["species"],"island_added":int(z["added"]),
           "native_global":int(lookup[z["species"]]["native_resource_units"]),
           "family":float(lookup[z["species"]]["host_family_count"]),
           "global_added":int(lookup[z["species"]]["introduced_added_units"])} for z in isle]
    r=rho([x["native_global"] for x in rows],[x["island_added"] for x in rows])
    if abs(r-0.5788964068803357)>1e-12:
        raise AssertionError(f"Fig3B native-host-resource breadth Spearman drift: {r}")
    # Tied cutoff issue: original published value uses native breadth descending,
    # then descending species lexical name. Ascending names yields 5 units more.
    source_rank=sorted(rows,key=lambda x:(x["native_global"],x["name"]),reverse=True)
    alternate=sorted(rows,key=lambda x:(-x["native_global"],x["name"]))
    top59=sum(x["island_added"] for x in source_rank[:59])
    bottom60=sum(x["island_added"] for x in source_rank[-60:])
    alternate59=sum(x["island_added"] for x in alternate[:59])
    if (top59,bottom60,alternate59)!=(2271,644,2276):
        raise AssertionError("Fig3B top/bottom quartile source tie-break drift")
    fam=[x for x in rows if x["family"]==1]
    c={"butterflies":len(rows),
       "globally_gaining_butterflies":sum(x["global_added"]>0 for x in rows),
       "island_gaining_butterflies":sum(x["island_added"]>0 for x in rows),
       "one_documented_host_plant_family":len(fam),
       "one_family_global_gainers":sum(x["global_added"]>0 for x in fam)}
    if (c["butterflies"],c["globally_gaining_butterflies"],c["island_gaining_butterflies"],
        c["one_documented_host_plant_family"],c["one_family_global_gainers"])!=(239,206,198,97,82):
        raise AssertionError("Fig3C source butterfly grouping drift")
    return {
      "schema":"GEB_v05_Fig1_and_Fig3BC_derived_aggregate_replay_v1",
      "status":"FIG1_AND_FIG3_B_C_SOURCE_AGGREGATE_REPLAY_PASS_FIG3_A_NOT_REPLAYED",
      "fig1":grouped,
      "fig3B":{"n_butterflies":len(rows),"Spearman_native_global_to_island_added":r,
               "top59_added":top59,"bottom60_added":bottom60,
               "top59_lexical_ascending_tie_alternative":alternate59,
               "boundary_native_global_resource_units":source_rank[58]["native_global"],
               "ambiguous_cutoff_tie_species":sorted(x["name"] for x in rows if x["native_global"]==source_rank[58]["native_global"])},
      "fig3C":c,
      "fig3A":"Only pinned numerical source receipt checked separately; plant-by-butterfly-by-region credit allocation NOT independently recomputed",
      "source_scope":"Previously derived WCVP v13 and NHM HOSTS geographic resource aggregates, not raw original accepted host edge records",
      "main_paper_ecological_estimates_modified":False}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--region-csv",required=True,type=Path)
    ap.add_argument("--by-butterfly-csv",required=True,type=Path)
    ap.add_argument("--species-metrics-csv",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    result=run(a.region_csv,a.by_butterfly_csv,a.species_metrics_csv)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("FIG1_REGION_NATIVE_ALIEN_KNOWN_HOST_239_AGGREGATE_REPLAY_PASS")
    print("FIG3_B_C_239_BUTTERFLY_SOURCE_RANK_AND_GROUPS_REPLAY_PASS")
    print("FIG3_A_PLANT_CREDIT_SOURCE_RECEIPT_ONLY_NOT_INDEPENDENT_REPLAY")
    print(json.dumps({"island":result["fig1"]["island"],
       "mainland":result["fig1"]["mainland"],
       "Fig3B":result["fig3B"],"Fig3C":result["fig3C"]},ensure_ascii=False))
if __name__=="__main__":
    main()
