#!/usr/bin/env python3
"""Verify published 120 predictive comparisons and 51 original port ranks against a fresh run."""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path

PUBLISHED=Path(__file__).resolve().parent

def rows(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        r=csv.DictReader(f);cols=list(r.fieldnames or []);items=list(r)
    if not cols or not items:
        raise ValueError("Published/fresh table empty: "+str(path))
    return cols,items

def match(name,count,output):
    fields1,a=rows(PUBLISHED/name)
    fields2,b=rows(Path(output)/name)
    if fields1!=fields2 or len(a)!=count or len(b)!=count:
        raise ValueError("Published/fresh table columns/row count mismatch: "+name)
    for i,(x,y) in enumerate(zip(a,b)):
        for field in fields1:
            if x[field]==y[field]:
                continue
            try:
                v=float(x[field]);w=float(y[field])
            except (ValueError,TypeError):
                raise ValueError(f"Source identity mismatch in {name} row {i} field {field}")
            if not math.isfinite(v) or not math.isfinite(w) or not math.isclose(v,w,rel_tol=1e-10,abs_tol=1e-11):
                raise ValueError(f"Changed numeric result in {name} row {i} field {field}")
    return {"rows":len(a),"columns":len(fields1),"matches":True}

def verify(output):
    result={
        "scope_comparisons":match("contrast_scope_summary.csv",120,output),
        "source_port_rank":match("original_port_roster_rank.csv",51,output),
    }
    print("TOP_PORT_PUBLISHED_SCORES_MATCH_FRESH_RESULTS="+str(result),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--results-dir",type=Path,required=True)
    verify(p.parse_args().results_dir)
