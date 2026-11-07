#!/usr/bin/env python3
"""Compare all thirty published AIS/wind summaries against fresh source-locked model outputs."""
from __future__ import annotations
import argparse,csv,math
from pathlib import Path
PUBLISHED=Path(__file__).resolve().parent/"contrast_scope_summary.csv"

def load(p):
    with Path(p).open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        return list(reader.fieldnames or []),list(reader)

def verify(results_dir):
    fields1,a=load(PUBLISHED);fields2,b=load(Path(results_dir)/"contrast_scope_summary.csv")
    if fields1!=fields2 or len(a)!=30 or len(b)!=30:
        raise ValueError("Archived source-factor summary shape/column mismatch")
    for idx,(x,y) in enumerate(zip(a,b)):
        for key in fields1:
            if x[key]==y[key]:
                continue
            try:
                u=float(x[key]);v=float(y[key])
            except (ValueError,TypeError):
                raise ValueError(f"Different model/target/scope ID at row {idx} field {key}")
            if not math.isfinite(u) or not math.isfinite(v) or not math.isclose(u,v,abs_tol=1e-11,rel_tol=1e-10):
                raise ValueError(f"Changed AIS/wind result at row {idx} field {key}")
    print("AIS_WIND_PUBLISHED_30_COMPARISONS_MATCH_FRESH_RUN=True",flush=True)
    return True

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--results-dir",required=True,type=Path)
    verify(parser.parse_args().results_dir)
