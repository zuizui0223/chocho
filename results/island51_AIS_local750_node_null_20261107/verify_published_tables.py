#!/usr/bin/env python3
"""Compare committed geographic-source AIS placebo scores and original 18 capped swaps to fresh run."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent

def read(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f)
        return list(rd.fieldnames or []),list(rd)

def verify_table(name,n,results):
    pa,a=read(HERE/name);pb,b=read(Path(results)/name)
    if pa!=pb or len(a)!=n or len(b)!=n:
        raise ValueError("Published source placebo table columns or rows differ: "+name)
    for index,(x,y) in enumerate(zip(a,b)):
        for col in pa:
            if x[col]==y[col]:
                continue
            try:
                v,w=float(x[col]),float(y[col])
            except (ValueError,TypeError):
                raise ValueError(f"Source identity changed in {name} row {index} field {col}")
            if not all(map(math.isfinite,(v,w))) or not math.isclose(v,w,rel_tol=1e-10,abs_tol=1e-10):
                raise ValueError(f"Source numeric result differs in {name} row {index} field {col}")
    return {"rows":n,"columns":len(pa),"matches":True}

def verify(directory):
    result={
        "original_vs_local_placebo_rank":verify_table("observed_vs_local750_node_null.csv",12,directory),
        "geography_only_selected_swap_pairs":verify_table("local_swap_pairs.csv",18,directory),
    }
    print("AIS_LOCAL_NULL_PUBLISHED_TABLES_MATCH_FRESH_RESULTS="+str(result),flush=True)
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-dir",type=Path,required=True)
    verify(parser.parse_args().results_dir)
