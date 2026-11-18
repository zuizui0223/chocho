#!/usr/bin/env python3
"""Check 12 permanent AIS placebo comparison rows and 51 source latitude-band rows against a fresh run."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent

def load(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def check(filename,n,generated):
    a_fields,a=load(HERE/filename)
    b_fields,b=load(Path(generated)/filename)
    if a_fields!=b_fields or len(a)!=n or len(b)!=n:
        raise ValueError("Published vs fresh fields or row count changed: "+filename)
    for index,(saved,real) in enumerate(zip(a,b)):
        for field in a_fields:
            if saved[field]==real[field]:
                continue
            try:
                x,y=float(saved[field]),float(real[field])
            except (ValueError,TypeError):
                raise ValueError(f"Published source identity changed: {filename} / {index} / {field}")
            if not math.isfinite(x) or not math.isfinite(y) or not math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-11):
                raise ValueError(f"Published null result changed: {filename} / {index} / {field}")
    return {"rows":n,"columns":len(a_fields),"matches":True}

def verify(out):
    result={"AIS_source_placebo_comparisons":check("observed_vs_latitude_node_null.csv",12,out),
            "original_island_latitude_groups":check("latitude_bands.csv",51,out)}
    print("AIS_NODE_NULL_PUBLISHED_TABLES_MATCH_FRESH_RESULTS="+str(result),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,required=True)
    verify(p.parse_args().results_dir)
