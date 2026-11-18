#!/usr/bin/env python3
"""Immediate frozen-artifact check: indexed HOSTS locality alien share minus
accepted regional vascular alien species share, on *islands only*.
This does not condition on family, and is not evidence of field preference.
"""
from pathlib import Path
import csv
import json
import math
import random
import sys
from statistics import median

def read(p):
    with open(p,newline="",encoding="utf8") as f:return list(csv.DictReader(f))

def main(artifact,root,out):
    host=read(artifact)
    flora={r["code"]:r for r in read(root/"results/three-gate-port-all-flora/all_flora_intro_native_wpi_islands.csv")}
    loc={r["HOSTS_location"]:r for r in read(root/"data/frozen/HOSTS_locality_WGSRPD3_strict_crosswalk_20261116.csv")}
    outrows=[]
    for r in host:
        if r["threshold"]!="10":continue
        if not r["analysis_key"].startswith("HOSTS_Papilionoidea__strict_exact_or_composite__"):continue
        l=loc.get(r["location"])
        if not l or l["comparison_group"]!="island" or l["mapping_tier"]!="exact":continue
        code=l["WGSRPD3_codes"]
        if "|" in code or code not in flora:continue
        f=flora[code]
        hostshare=float(r["introduced_share"])
        florashare=float(f["intro_flora_fraction"])
        outrows.append({"scope":r["analysis_key"],"location":r["location"],
          "code":code,"n_HOSTS":int(r["n"]),"introduced_HOSTS":int(r["introduced"]),
          "reported_share":hostshare,"flora_share":florashare,
          "excess_pp":100*(hostshare-florashare)})
    allscopes={}
    for scope in sorted({r["scope"] for r in outrows}):
        rows=sorted((r for r in outrows if r["scope"]==scope),key=lambda x:x["location"])
        ds=[r["excess_pp"] for r in rows]
        rng=random.Random(20261116)
        boot=sorted(median([ds[rng.randrange(len(ds))] for _ in ds]) for _ in range(4999))
        allscopes[scope]={"islands":len(rows),"median_host_over_flora_pp":median(ds),
        "positive_islands":sum(x>0 for x in ds),
        "median_bootstrap_95_descriptive":[boot[124],boot[4874]],
        "source_islands":rows}
    if len(allscopes)!=4 or not allscopes:
        raise RuntimeError("Incomplete Papilionoidea fixed geographic source comparison")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps({"status":"WHOLE_FLORA_ISLAND_SOURCE_ONLY_NOT_FAMILY_MATCH",
       "input_source_run":43881198606,
       "all_floristic_reference":"Original independent 126 strict-island WCVP flora counts",
       "scopes":allscopes,
       "limits":"Source locality bootstrap descriptive, plant species shares are NOT actual larval host abundances; no mainland or family control"},indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("HOSTS_WHOLE_FLORA_ISLAND_REFERENCE",json.dumps(
        {k:{"islands":v["islands"],"median_over_flora_pp":v["median_host_over_flora_pp"],
        "positive_islands":v["positive_islands"],"interval":v["median_bootstrap_95_descriptive"],
        "source_islands":v["source_islands"]} for k,v in allscopes.items()},ensure_ascii=False),flush=True)

if __name__=="__main__":
    if len(sys.argv)!=4:raise SystemExit("Usage: script <source per_locality CSV> <repo root> <out.json>")
    main(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]))
