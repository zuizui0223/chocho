#!/usr/bin/env python3
"""Correct modern WPI exposure labels for official Canary Islands WGSRPD3 region CNY.

No new port coordinates are fetched. The original pinned 3,630-port source
assignment remains unchanged, and only the island-vs-mainland geographic
grouping is repaired for CNY. Do not interpret indexed ports as ship traffic.
"""
from __future__ import annotations
import argparse,csv,json,runpy
from pathlib import Path
def rows(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def csvwrite(p,rr):
    with p.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]))
        w.writeheader();w.writerows(rr)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-counts",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--corrected-q1-classifier",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    orig=json.loads(a.source_audit.read_text(encoding="utf-8"))
    if orig["frozen_strict_island_regions"]!=125 or orig["all_public_WPI_ports"]!=3630:
        raise RuntimeError("Original source audit version drift")
    source=rows(a.source_counts)
    if len(source)!=369 or sum(r["geographic_group"]=="island" for r in source)!=125:
        raise RuntimeError("Original public WPI region table drift")
    classifier=runpy.run_path(str(a.corrected_q1_classifier))
    islands=set(classifier["STRICT_ISLAND_CODES"])
    mixed=set(classifier["MIXED_CODES"])
    if len(islands)!=126 or "CNY" not in islands or "RHO" in islands:raise RuntimeError("New island classifier wrong")
    changes=[]
    for r in source:
        code=r["area_code_l3"]
        new="island" if code in islands else ("mixed" if code in mixed else "mainland")
        if r["geographic_group"]!=new:
            changes.append({"code":code,"before":r["geographic_group"],
                            "after":new,"n_WPI_ports":int(r["WPI_port_count"])})
        r["geographic_group"]=new
    if changes!=[{"code":"CNY","before":"mainland","after":"island","n_WPI_ports":4}]:
        raise RuntimeError(f"Unexpected island port source changes {changes}")
    islands=[r for r in source if r["geographic_group"]=="island"]
    report={
        "schema":"chocho_WPI_126_island_CNY_corrected_source_v0.2",
        "original_WPI_audit_schema":orig["schema"],
        "original_indexed_WPI_ports":orig["all_public_WPI_ports"],
        "original_matched_WPI_positions":orig["matched_port_positions"],
        "changes":changes,
        "corrected_island_regions":126,
        "corrected_islands_with_at_least_one_indexed_WPI_port":sum(int(r["WPI_port_count"])>0 for r in islands),
        "corrected_indexed_WPI_island_ports":sum(int(r["WPI_port_count"]) for r in islands),
        "not_a_shipping_volume_or_intro_date_proxy":True,
        "source_2019_vs_2023_mirror_mismatch":True,
    }
    csvwrite(a.out_dir/"wgsrpd3_port_counts_corrected_v02.csv",source)
    (a.out_dir/"wpi_port_exposure_audit_v02.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
