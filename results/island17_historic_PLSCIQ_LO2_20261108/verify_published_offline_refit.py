#!/usr/bin/env python3
"""Full numerical offline refit consistency against published old-PLSCIQ butterfly tables.

Runs on newly fit source-matched 17-island LO2 outcomes, not a network.
Compares all 2048 original pair predictions as decompressed data, all 20
comparisons and all source index/port/fold manifests to saved published rows.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path

ORIGINAL=Path(__file__).resolve().parent
FILES={
    "original24_geoverified_ports_historical_PLSCIQ_quality.csv":24,
    "original51_historical_PLSCIQ_island_eligibility.csv":51,
    "original17_136_source_pair_manifest.csv":136,
    "double_unseen_source_folds.csv":256,
    "comparison_scope_summary.csv":20,
    "LO2_heldout_predictions.csv.gz":2048
}

def text_source(path):
    p=Path(path)
    return gzip.open(p,"rt",encoding="utf-8",newline="") if str(p).endswith(".gz") else p.open("rt",encoding="utf-8",newline="")

def verify(fresh):
    compared={}
    for name,count in FILES.items():
        with text_source(ORIGINAL/name) as original,text_source(Path(fresh)/name) as newest:
            old=csv.DictReader(original);new=csv.DictReader(newest)
            if old.fieldnames!=new.fieldnames:
                raise ValueError("Published source old/historic index field list changed "+name)
            n=0
            while True:
                a=next(old,None);b=next(new,None)
                if a is None and b is None:break
                if a is None or b is None:
                    raise ValueError("Published old-port index source file row count differs "+name)
                for field in old.fieldnames or []:
                    if a[field]==b[field]:continue
                    try:
                        x=float(a[field]);y=float(b[field])
                    except (TypeError,ValueError):
                        raise ValueError(f"Original historical source identity changed in {name} at row {n} field {field}")
                    if not(math.isfinite(x) and math.isfinite(y) and
                           math.isclose(x,y,abs_tol=1e-10,rel_tol=1e-10)):
                        raise ValueError(f"Original historical index numeric calculation changed in {name} row {n} field {field}: {x} vs {y}")
                n+=1
        if n!=count:
            raise ValueError("Published original historical PLSCIQ source file has wrong number of rows: "+name)
        compared[name]={"rows":n,"fields":len(old.fieldnames or []),"matched":True}
    print("HISTORICAL_PLSCIQ_PUBLISHED_OFFLINE_REFIT_VERIFIED="+str(compared),flush=True)
    return compared

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-dir",required=True,type=Path)
    verify(parser.parse_args().results_dir)
