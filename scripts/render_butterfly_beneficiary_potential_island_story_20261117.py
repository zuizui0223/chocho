#!/usr/bin/env python3
"""Render a source-verified butterfly-first two-panel figure, without inventing
island-endemic status or claiming GBIF butterfly community homogenization."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--source-dir",type=Path,required=True)
 ap.add_argument("--out-prefix",type=Path,required=True)
 a=ap.parse_args()
 inp=a.source_dir/"BENEFICIARY_AND_POTENTIAL_HOMOGENIZATION_20261117.json"
 d=json.loads(inp.read_text())
 q=pd.read_csv(a.source_dir/"beneficiary_four_quartiles.csv")
 assert (d["cohort"]["butterflies"],d["cohort"]["strict_islands"],d["cohort"]["original_island_added_cells"])==(239,126,5741)
 assert len(q)==4 and q.new_island_resource_cells.sum()==5741
 assert [int(x) for x in q.butterflies]==[60,60,60,59]
 obs=d["potential_island_consumer_composition"]["observed"]
 assert obs["pairs_comparable_both_native_and_contemporary_nonempty"]==7797
 native,now=obs["mean_potential_butterfly_Jaccard_native"],obs["mean_potential_butterfly_Jaccard_contemporary"]
 assert abs(native-.1401370010070764)<1e-9 and abs(now-.3675355469925412)<1e-9
 within=d["potential_island_consumer_composition"]["within_vs_between_Level1"]
 assert within["between_Level1"]["pairs"]==6808 and within["within_Level1"]["pairs"]==989

 plt.rcParams.update({
  "font.family":"DejaVu Sans","font.size":10.5,
  "axes.linewidth":.7,"axes.spines.top":False,"axes.spines.right":False,
  "savefig.facecolor":"white","figure.facecolor":"white"
 })
 fig,(left,right)=plt.subplots(1,2,figsize=(11.4,5.15),gridspec_kw={"width_ratios":[1.05,.95]})
 fig.subplots_adjust(left=.075,right=.97,top=.74,bottom=.24,wspace=.35)
 green="#1d7467";blue="#618697";gray="#707983";orange="#b8793c";color1="#9bb7b2"
 vals=(100*q.fraction_all_island_resource_additions).to_numpy()
 bars=left.bar(np.arange(4),vals,color=[color1,"#82a8a0","#4d9587",green],width=.69)
 left.axhline(25,c=gray,linestyle=(0,(3,2)),linewidth=1,alpha=.8)
 left.set_xticks(range(4),["Q1\nNarrow","Q2","Q3","Q4\nBroad"])
 left.set_ylim(0,46)
 left.set_yticks([0,10,20,30,40])
 left.set_ylabel("Share of new island host opportunities (%)",fontsize=10)
 for i,v in enumerate(vals):
  left.text(i,v+1.0,f"{v:.1f}%",ha="center",va="bottom",fontsize=11,weight="bold")
 left.set_title("A   Which butterflies gain resources?",loc="left",fontsize=12,weight="bold",pad=13)
 left.text(.02,.96,"Quartiles of native HOSTS resource-region breadth",
          transform=left.transAxes,ha="left",va="top",fontsize=8.7,color=gray)
 left.text(-.32,27.0,"25% equal-share reference",ha="left",fontsize=8,color=gray)

 p=right.bar([0,1],[native,now],width=.52,color=[blue,orange])
 right.set_ylim(0,.43)
 right.set_ylabel("Mean potential butterfly Jaccard similarity",fontsize=10)
 right.set_xticks([0,1],["Native hosts\nonly","Native +\nintroduced hosts"])
 right.yaxis.set_major_formatter(FuncFormatter(lambda v,pos: f"{v:.2f}"))
 for i,v in enumerate([native,now]):
  right.text(i,v+.012,f"{v:.3f}",ha="center",fontsize=12,weight="bold")
 right.set_title("B   A shared potential consumer pool",loc="left",fontsize=12,weight="bold",pad=13)
 right.text(.5,.91,"Same 7,797 island pairs",ha="center",fontsize=9,color=gray,
            transform=right.transAxes)
 fig.text(.075,.94,"Human-assembled plants redistribute potential food for butterfly consumers",
          fontsize=16,weight="bold",ha="left")
 fig.text(.075,.86,"239 butterflies · 126 island botanical regions · fixed known larval-host repertoire",
          fontsize=10.1,color=gray,ha="left")
 fig.text(.075,.095,
   "Potential only: a butterfly is counted where ≥1 documented larval host is regionally available.",
   fontsize=9.2,color="#344148",ha="left")
 fig.text(.075,.065,
   "Native-host breadth is NOT the butterfly's actual range. Similarity gain largely follows resource-addition margins.",
   fontsize=8.6,color=gray,ha="left")
 a.out_prefix.parent.mkdir(parents=True,exist_ok=True)
 for ext in ["png","pdf","svg"]:
  fig.savefig(a.out_prefix.with_suffix("."+ext),dpi=280)
 plt.close(fig)
 outputs={x:{"path":str(a.out_prefix.with_suffix("."+x)),"sha256":sha(a.out_prefix.with_suffix("."+x)),
             "bytes":a.out_prefix.with_suffix("."+x).stat().st_size} for x in ["png","pdf","svg"]}
 receipt={"schema":"chocho_butterfly_resource_beneficiary_potential_Jaccard_main_figure_v1",
  "status":"SOURCE_VERIFIED_MAIN_FIGURE_ONLY_POTENTIAL_BUTTERFLY_COMMUNITIES",
  "source_file":str(inp),"source_sha256":sha(inp),
  "values":{"quartile_percents":vals.tolist(),"native_potential_Jaccard":native,
            "contemporary_potential_Jaccard":now,"pair_count":7797},
  "null_warning":"Two degree-preserving nulls in source JSON; Jaccard increase is mostly explained by potential resource expansion and marginals, not independently observed butterfly communities.",
  "output_files":outputs}
 (a.out_prefix.parent/(a.out_prefix.name+"_receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
 print("MAIN_FIGURE_SOURCE_VERIFIED="+json.dumps(receipt,sort_keys=True),flush=True)
if __name__=="__main__":main()
