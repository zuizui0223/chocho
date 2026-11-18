#!/usr/bin/env python3
"""Verify all published early general-flora source-covariates and scores match offline refits.

The immutable original194 source snapshot is loaded locally. No network required.
The independent audit of the refitted source pairs runs separately in CI.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import math
from pathlib import Path

PUBLISHED=Path(__file__).resolve().parent
FILES={
    "heldout_predictions.csv.gz":15300,
    "model_fold_losses.csv":180,
    "comparison_fold_losses.csv":180,
    "comparison_scope_summary.csv":48,
    "original51_early_general_flora_coverage.csv":51,
    "source_geographic_folds.csv":15,
    "source_pair_general_flora_covariates.csv":1275,
}

def open_rows(filename,base):
    path=Path(base)/filename
    if filename.endswith(".gz"):
        file=gzip.open(path,"rt",encoding="utf-8",newline="")
    else:
        file=path.open("rt",encoding="utf-8",newline="")
    return file

def verify(result_dir):
    verified={}
    for filename,expected in FILES.items():
        with open_rows(filename,PUBLISHED) as previous,open_rows(filename,result_dir) as current:
            old=csv.DictReader(previous);new=csv.DictReader(current)
            if old.fieldnames!=new.fieldnames:
                raise ValueError("Offline published source/result fieldnames differ: "+filename)
            count=0
            while True:
                a=next(old,None);b=next(new,None)
                if a is None and b is None:
                    break
                if a is None or b is None:
                    raise ValueError("Published/fresh result row counts differ: "+filename)
                for k in old.fieldnames or []:
                    if a[k]==b[k]:
                        continue
                    try:
                        x=float(a[k]);y=float(b[k])
                    except (TypeError,ValueError):
                        raise ValueError(f"Published identity changed at {filename} row {count} field {k}")
                    if not (math.isfinite(x) and math.isfinite(y) and
                            math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-10)):
                        raise ValueError(f"Published numeric value changed at {filename} row {count} field {k}")
                count+=1
            if count!=expected:
                raise ValueError("Published/fresh source required row count incorrect: "+filename)
            verified[filename]={"rows":count,"fields":len(old.fieldnames or []),"matches":True}
    print("EARLY_FLORA_PUBLISHED_OFFLINE_RECONSTRUCTION_JSON="+str(verified),flush=True)
    return verified

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-dir",type=Path,required=True)
    verify(parser.parse_args().results_dir)
