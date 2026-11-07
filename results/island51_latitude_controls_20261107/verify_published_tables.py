#!/usr/bin/env python3
"""Compare all 72 published latitude-control summary rows to freshly verified results."""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

PUB=Path(__file__).resolve().parent/"contrast_scope_summary.csv"

def load(path):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or ()),list(reader)

def verify(results_dir):
    fields_a,expected=load(PUB)
    fields_b,actual=load(Path(results_dir)/"contrast_scope_summary.csv")
    if fields_a!=fields_b or len(expected)!=72 or len(actual)!=72:
        raise ValueError("Published/fresh source comparison shape mismatch")
    for index,(x,y) in enumerate(zip(expected,actual)):
        for key in fields_a:
            if x[key]==y[key]:
                continue
            try:
                a,b=float(x[key]),float(y[key])
            except (ValueError,TypeError):
                raise ValueError("Comparison source identity differs at "+str((index,key)))
            if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a,b,abs_tol=1e-11,rel_tol=1e-10):
                raise ValueError("Comparison result differs at "+str((index,key)))
    print("LATITUDE_CONTROL_PUBLISHED_72_COMPARISONS_MATCH_FRESH_RUN=True",flush=True)
    return True

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--results-dir",type=Path,required=True)
    verify(p.parse_args().results_dir)
