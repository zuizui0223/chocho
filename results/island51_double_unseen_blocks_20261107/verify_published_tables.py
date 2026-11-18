#!/usr/bin/env python3
"""Verify published double-unseen summary and geographical assignments from a fresh run."""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path

PUB=Path(__file__).resolve().parent

def read(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f);rows=list(reader)
    if not rows:
        raise ValueError("Published or regenerated table empty: "+str(path))
    return list(rows[0]),rows

def match(name,n,output):
    cols1,a=read(PUB/name);cols2,b=read(output/name)
    if cols1!=cols2 or len(a)!=n or len(b)!=n:
        raise ValueError("Published table shape/column mismatch: "+name)
    for i,(x,y) in enumerate(zip(a,b)):
        for k in cols1:
            if x[k]==y[k]:continue
            try:
                ax=float(x[k]);by=float(y[k])
            except (TypeError,ValueError):
                raise ValueError(f"Geographical identity or score mismatch: {name} row {i} {k}")
            if not math.isfinite(ax) or not math.isfinite(by) or not math.isclose(ax,by,rel_tol=1e-10,abs_tol=1e-11):
                raise ValueError(f"Published numeric mismatch: {name} row {i} {k}")
    return {"rows":n,"columns":len(cols1),"matches":True}

def verify(out):
    a=match("longitude_island_blocks.csv",51,out)
    b=match("contrast_scope_summary.csv",12,out)
    result={"longitude_groups":a,"scope_comparisons":b}
    print("DOUBLE_UNSEEN_PUBLISHED_TABLES_MATCH="+str(result),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--results-dir",type=Path,required=True)
    verify(p.parse_args().results_dir)
