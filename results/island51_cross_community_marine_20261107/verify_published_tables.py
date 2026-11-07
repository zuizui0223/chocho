#!/usr/bin/env python3
"""Verify every committed model/pair/score CSV against fresh source-only offline refit.

Compares 7 tables including all 15,840 original source-pair forecast rows;
never requires GBIF network data or current latest sources.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import math
from pathlib import Path

PUBLISHED=Path(__file__).resolve().parent
ROWS={
    "source_pair_community_panel.csv":1275,
    "source_eligible_geo_folds.csv":15,
    "model_fold_losses.csv":240,
    "comparison_fold_losses.csv":300,
    "comparison_scope_summary.csv":80,
    "heldout_predictions.csv.gz":15840,
}

def open_table(path):
    path=Path(path)
    return (gzip.open(path,"rt",encoding="utf-8",newline="") if path.suffix==".gz"
            else path.open("rt",encoding="utf-8",newline=""))

def compare(name,rows,derived):
    original=PUBLISHED/name
    fresh=Path(derived)/name
    with open_table(original) as left,open_table(fresh) as right:
        a=csv.DictReader(left);b=csv.DictReader(right)
        if a.fieldnames!=b.fieldnames:
            raise ValueError("Published/derived source table columns changed: "+name)
        count=0
        while True:
            l=next(a,None);r=next(b,None)
            if l is None and r is None:
                break
            if l is None or r is None:
                raise ValueError("Published/derived source table length differs: "+name)
            for field in a.fieldnames or ():
                if l[field]==r[field]:
                    continue
                try:
                    va=float(l[field]);vb=float(r[field])
                except (TypeError,ValueError):
                    raise ValueError(f"Source ID/field differs: {name} row {count} field {field}")
                if not (math.isfinite(va) and math.isfinite(vb) and
                        math.isclose(va,vb,rel_tol=1e-10,abs_tol=1e-10)):
                    raise ValueError(f"Source numeric comparison differs: {name} row {count} field {field}")
            count+=1
        if count!=rows:
            raise ValueError("Published/derived source table rows incorrect for "+name)
        return {"rows":rows,"columns":len(a.fieldnames or []),"equal":True}

def verify(output):
    verified={name:compare(name,n,output) for name,n in ROWS.items()}
    print("CROSS_COMMUNITY_PUBLISHED_TABLES_MATCH_FRESH_SOURCE_ONLY_RESULTS="+str(verified),flush=True)
    return verified

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--results-dir",type=Path,required=True)
    verify(ap.parse_args().results_dir)
