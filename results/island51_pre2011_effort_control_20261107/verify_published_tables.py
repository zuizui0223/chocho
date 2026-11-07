#!/usr/bin/env python3
"""Check publication summaries against independently verified fresh CI outputs.

Uses Python standard library only; no model runner imports or statistical refits.
"""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path

PUB=Path(__file__).resolve().parent

def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        rd=csv.DictReader(f)
        columns=tuple(rd.fieldnames or ())
        data=list(rd)
    if not columns or not data:
        raise ValueError("Empty published or fresh table: "+str(path))
    return columns,data

def same_tables(published,actual,expected_rows,kind):
    fields1,a=rows(published)
    fields2,b=rows(actual)
    if fields1!=fields2 or len(a)!=expected_rows or len(b)!=expected_rows:
        raise ValueError("Published fields or row counts disagree with fresh "+kind)
    for index,(old,new) in enumerate(zip(a,b)):
        for key in fields1:
            x,y=old[key],new[key]
            if x==y:
                continue
            try:
                f,g=float(x),float(y)
            except (ValueError,TypeError):
                raise ValueError(f"Changed publication identity {kind}/{index}/{key}")
            if not math.isfinite(f) or not math.isfinite(g) or not math.isclose(f,g,rel_tol=1e-10,abs_tol=1e-11):
                raise ValueError(f"Changed publication value {kind}/{index}/{key}: {f} != {g}")
    return {"rows":len(a),"fields":len(fields1),"matches":True}

def verify(results_dir):
    results=Path(results_dir)
    report={
        "comparison_summary":same_tables(PUB/"comparison_summary.csv",results/"comparison_summary.csv",14,"summary"),
        "historical_botanical_coverage":same_tables(
            PUB/"historical_botanical_coverage_by_island.csv",
            results/"historical_botanical_coverage_by_island.csv",51,"coverage"),
    }
    print("PUBLISHED_TABLES_MATCH_FRESH_RESULTS="+str(report),flush=True)
    return report

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-dir",required=True,type=Path)
    verify(parser.parse_args().results_dir)
