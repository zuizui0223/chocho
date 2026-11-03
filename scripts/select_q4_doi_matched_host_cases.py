#!/usr/bin/env python3
"""Q4 dated host source verification for observed (13) vs non-indexed butterfly regions.

This is a *case-control retrieval design* to fill plant chronology gaps, not a
colony-establishment experiment. Cases are selected by DOI butterfly indexing
(which is subject to effort). Controls were non-indexed in a previously selected
31-species DOI, NOT verified absent butterflies. Each case picks up to 3
within-species (prefer same Level1) controls matched on log WGSRPD area, modern
WPI ports and records of other DOI butterflies. Matching never uses plant date.
"""
from __future__ import annotations
import argparse,csv,json,math,hashlib,statistics
from collections import defaultdict
from pathlib import Path

def read(p):
    with p.open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))
def write(p, rows):
    if not rows:raise RuntimeError("No rows")
    with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    for x in ("q3-risk","q3-records","q2-edges","q2-audit","out-dir"):
        ap.add_argument("--"+x,type=Path,required=True)
    a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q3=read(a.q3_risk)
    pos=read(a.q3_records)
    edge=read(a.q2_edges)
    audit=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    if (len(q3),len(pos),len(edge),audit["scope"]["added_strict_island_cells"])!=(920,13,11573,5741):
        raise RuntimeError("Q2/Q3 pinned cohort drift")
    observed={(r["species"],r["wgsrpd3_code"]) for r in q3 if int(r["archival_butterfly_record_at_added_resource"])}
    positive={(r["species"],r["area_code_l3"]) for r in pos}
    if len(observed)!=13 or observed!=positive:raise RuntimeError("Q3 positivity mismatch")
    bycell=defaultdict(set)
    for r in edge:
        bycell[(r["species"],r["island_code"])].add(r["plant_id"])
    scored=[]
    for r in q3:
        s=r["species"];reg=r["wgsrpd3_code"];cell=(s,reg)
        if cell not in bycell:raise RuntimeError("Q3 member not in corrected Q2: "+str(cell))
        scored.append({
            "species":s,"island_code":reg,"level1":r["level1"],
            "recorded":int(cell in positive),
            "ports":int(r["indexed_WPI_port_count_2023"]),
            "lnarea":math.log(float(r["WGSRPD3_area_km2"])),
            "lnother":math.log1p(int(r["other_30_butterfly_record_count"])),
            "n_hosts":len(bycell[cell]),
        })
    positives=[r for r in scored if r["recorded"]]
    nonindexed=[r for r in scored if not r["recorded"]]
    sample=[]
    selected=set()
    matching={}
    for p in sorted(positives,key=lambda r:(r["species"],r["island_code"])):
        key=(p["species"],p["island_code"])
        sample.append({**p,"role":"BUTTERFLY_DOI_INDEXED","index_case":p["island_code"],
                       "matched_block":"CASE","matched_rank":0})
        selected.add(key)
        near=[x for x in nonindexed if x["species"]==p["species"] and x["level1"]==p["level1"]]
        distant=[x for x in nonindexed if x["species"]==p["species"] and x["level1"]!=p["level1"]]
        def ranking(x):
            # Distances are source-fixed, do not use host year or outcome.
            dx=abs(x["lnarea"]-p["lnarea"])/3+abs(math.log1p(x["ports"])-math.log1p(p["ports"]))
            de=abs(x["lnother"]-p["lnother"])
            return (dx+de/2,x["island_code"])
        same=sorted(near,key=ranking)
        fallback=sorted(distant,key=ranking)
        chosen=[]
        for pool,reason in ((same,"WITHIN_SPECIES_LEVEL1"),(fallback,"WITHIN_SPECIES_OUTSIDE_LEVEL1")):
            for x in pool:
                if len(chosen)>=3:break
                cid=(x["species"],x["island_code"])
                if cid in selected:continue
                chosen.append((x,reason))
                selected.add(cid)
            if len(chosen)>=3:break
        for i,(c,why) in enumerate(chosen,1):
            sample.append({**c,"role":"NOT_INDEXED_IN_ARCHIVED_DOI_NOT_ABSENT",
                          "index_case":p["island_code"],"matched_block":why,"matched_rank":i})
        matching[str(key)]={"control_n":len(chosen),"same_Level1":sum(x[1]=="WITHIN_SPECIES_LEVEL1" for x in chosen)}
    if len(positives)!=13 or len(sample)<27:raise RuntimeError("Too few matched DOI comparison cells")
    plant_pairs={}
    for r in sample:
        for h in sorted(bycell[(r["species"],r["island_code"])]):
            key=(h,r["island_code"])
            if key not in plant_pairs:
                plant_pairs[key]={"plant_id":h,"island_code":r["island_code"],
                  "matched_case_source_species":r["species"],
                  "selected_case_or_control":r["role"],
                  "matched_on":"species+prefer_level1+logarea+logWPI+other_DOI_butterfly_records"}
            else:
                plant_pairs[key]["selected_case_or_control"]="SHARED_ACROSS_CASE_CONTROL_CELLS"
    write(a.out_dir/"matched_species_island_record_index_cells.csv",sample)
    write(a.out_dir/"all_matched_candidate_introduced_host_island_pairs.csv",list(plant_pairs.values()))
    summary={
      "schema":"chocho_q4_indexed_butterfly_hostcase_control_source_selection_v0.1",
      "status":"MATCHED_RETRIEVAL_DESIGN_NOT_BIOLOGICAL_ARRIVAL_OR_ESTABLISHMENT",
      "source":{"q3_31_species_sampled_opportunity_cells":920,"q3_index_positive_cells":13,
                "q2_126island_resource_cells":5741,"q2_introduced_host_edges":11573},
      "case_control_population":{
        "positive_DOI_butterfly_index_cells":13,
        "unique_nonindexed_control_cells":sum(r["role"]!="BUTTERFLY_DOI_INDEXED" for r in sample),
        "total_case_control_butterfly_region_cells":len(sample),
        "species_in_butterfly_records":len({r["species"] for r in positives}),
        "distinct_Q2_introduced_host_x_island_pairs":len(plant_pairs),
        "distinct_host_plants":len({h for h,code in plant_pairs}),
        "matched_control_blocks":dict({name:sum(r["matched_block"]==name for r in sample)
                                       for name in ("WITHIN_SPECIES_LEVEL1","WITHIN_SPECIES_OUTSIDE_LEVEL1")}),
        "per_case_control_support":matching,
      },
      "limitations":[
       "Controls have no indexed record in this particular earlier 31-species GBIF DOI, NOT verified butterfly absence.",
       "Case status is a modern observational index rather than independent island naturalization or arrival date.",
       "Selection uses butterfly DOI outcome to fill plant index dates; case-control comparisons are exploratory and require observation-effort correction.",
       "Even with matching, source GBIF butterfly sample is nonrandom and case-control samples are not a complete year-specific risk set.",
       "Matching on modern port count describes record access and might overmatch a true access effect; no port causal effect can be inferred.",
       "Plants' native-versus-introduced WGSRPD3 status is contemporary, and original year of introduction is unknown.",
       "WGSRPD3 areas can be archipelagos and botanical/entomological records can originate on different islands."
      ]
    }
    (a.out_dir/"matched_source_inventory.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
