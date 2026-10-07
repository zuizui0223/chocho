#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resource-cells-csv",type=Path,required=True)
    ap.add_argument("--pseudoresource-matches-csv",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    a=ap.parse_args()

    single=set()
    with a.resource_cells_csv.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if int(float(r["actual_host_count"]))==1:
                single.add((str(r["species"]).strip(),str(r["wgsrpd3_code"]).strip()))
    if len(single)!=1143:
        raise RuntimeError(f"expected 1,143 single-host cells; got {len(single)}")

    by_cell=defaultdict(list)
    with a.pseudoresource_matches_csv.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key=(str(r["species"]).strip(),str(r["wgsrpd3_code"]).strip())
            if key in single:
                by_cell[key].append(r)

    out=[]; excluded=[]
    for key in sorted(single):
        vals=by_cell.get(key,[])
        actual=[r for r in vals if r["assignment_type"]=="actual"]
        pseudo=[r for r in vals if r["assignment_type"]=="pseudo" and int(float(r["match_rank"]))==1]
        if len(actual)!=1 or len(pseudo)!=1:
            excluded.append({
                "species":key[0],"wgsrpd3_code":key[1],
                "actual_rows":len(actual),"rank1_pseudo_rows":len(pseudo),
                "reason":"MISSING_MATCHABLE_ACTUAL_OR_RANK1_PSEUDO"
            })
            continue
        ar=dict(actual[0]); pr=dict(pseudo[0])
        ar["primary_pair_role"]="actual"
        pr["primary_pair_role"]="pseudo"
        out.extend([ar,pr])

    if not out:
        raise RuntimeError("no matched single-host primary cells")
    a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with a.output_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),extrasaction="ignore"); w.writeheader(); w.writerows(out)

    payload={
        "schema":"chocho_single_host_temporal_pair_panel_v0.2",
        "status":"FROZEN_PAIRED_ACTUAL_VS_NEAREST_PSEUDORESOURCE",
        "single_host_resource_cells_total":len(single),
        "matched_primary_cells":len(out)//2,
        "excluded_cells":len(excluded),
        "coverage_fraction":(len(out)//2)/len(single),
        "rows":len(out),
        "rule":"For every outcome-blind introduced-only cell with exactly one actual known host, retain that actual host and the rank-1 matched alien non-host pseudo-resource selected without local chronology.",
        "excluded":excluded,
        "claim_boundary":"Coverage exclusions arise only from resource/pseudo chronology matching availability, not from butterfly outcomes or local first-record dates."
    }
    a.output_json.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in payload.items() if k!="excluded"},indent=2))

if __name__=="__main__":
    main()
