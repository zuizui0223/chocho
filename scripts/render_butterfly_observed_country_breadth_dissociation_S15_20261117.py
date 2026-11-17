#!/usr/bin/env python3
"""Source-verified supplementary Figure S15: country-record breadth of butterfly
itself versus potential island host gains and later first GBIF indexing.

Important: x-axis are NOT native butterfly ranges; 222 species with >=30
2001–2010 georeferenced GBIF indexed occurrences, quartiles of exact expected
number of indexed countries among 30 records. Panels have different responses,
different geographic supports and different sampling processes. Not causal.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
QSOURCE="results/butterfly_global_GBIF_country_breadth_20261117/independent_GBIF_country_breadth_quartile_beneficiaries.csv"
REPORT="results/butterfly_global_GBIF_country_breadth_20261117/BUTTERFLY_GBIF_COUNTRY_BREADTH_VS_ISLAND_HOST_GAIN_20261117.json"
OOF="results/butterfly_global_GBIF_country_breadth_20261117/GLOBAL_COUNTRY_BREADTH_AND_HOST_RESOURCE_FIRSTINDEX_OOF_20261117.json"
PIN_Q="7cb264b03d9bff365b6d9be07c78656c952a5ac8"
PIN_R="c2136a2fbed02c5a6ffa213bfdd90651d90dbc1c"
PIN_O="bd615d68cfb1588ee65770379ecb9721a874c503"

def gitblob(path):
 b=Path(path).read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument("--out-prefix",type=Path,required=True);a=p.parse_args()
 for n,pin in ((QSOURCE,PIN_Q),(REPORT,PIN_R),(OOF,PIN_O)):
  if gitblob(ROOT/n)!=pin:raise ValueError("Frozen original source drift "+n)
 source=pd.read_csv(ROOT/QSOURCE)
 report=json.loads((ROOT/REPORT).read_text())
 oof=json.loads((ROOT/OOF).read_text())
 assert len(source)==4 and source.species.sum()==222
 assert int(source.new_island_host_resource_opportunities.sum())==5617
 assert int(source.new_later_51_indexed_species_island_cells.sum())==317
 assert oof["n_risk_cells"]==10586 and oof["n_future_first_index_events"]==317
 assert report["cohort"]["eligible_30_record_rarefaction_species"]==222
 a1=100*source.share_island_added_covered_by_rarefaction_eligible_species.to_numpy()
 a2=100*source.fraction_new_51_index_events_in_rarefied_sample.to_numpy()
 assert abs(sum(a1)-100)<1e-8 and abs(sum(a2)-100)<1e-8

 plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10.5,
     "axes.linewidth":.75,"axes.spines.top":False,"axes.spines.right":False,
     "figure.facecolor":"white","savefig.facecolor":"white"})
 fig,axes=plt.subplots(1,2,figsize=(12.0,5.75))
 fig.subplots_adjust(left=.065,right=.98,top=.71,bottom=.28,wspace=.31)
 green="#287a6a";amber="#c27a3a";gray="#73808a"
 x=np.arange(4)
 def panel(ax,vals,col,head,ylabel):
  ax.bar(x,vals,color=col,width=.64)
  ax.axhline(25,color="#8c969d",linewidth=1,linestyle=(0,(3,3)))
  ax.set_xticks(x,["Q1\nNarrow","Q2","Q3","Q4\nBroad"])
  ax.set_ylim(0,58);ax.set_yticks([0,10,20,30,40,50]);ax.set_ylabel(ylabel,fontsize=10.5)
  ax.set_title(head,loc="left",fontsize=12,weight="bold",pad=10)
  ax.text(.98,.97,"Dashed line: equal share",transform=ax.transAxes,
          ha="right",va="top",fontsize=8.6,color=gray)
  for i,v in enumerate(vals):
   ax.text(i,v+1.1,f"{v:.1f}%",ha="center",fontsize=12,weight="bold")
 panel(axes[0],a1,green,"A  Potential HOSTS resource additions","Share of 5,617 newly covered\nbutterfly × botanical-island pairs (%)")
 panel(axes[1],a2,amber,"B  Later recorded island butterflies","Share of 317 first GBIF-indexed\nbutterfly × sampled-landmass pairs (%)")
 fig.text(.065,.94,"New potential resources are not reserved for already widely recorded butterflies",
          ha="left",fontsize=16,weight="bold")
 fig.text(.065,.855,
          "222 butterflies · quartiles by 2001–2010 GBIF countries rarefied to exactly 30 indexed records",
          fontsize=10.5,ha="left",color=gray)
 fig.text(.065,.153,
          "Resource additions: 126 WCVP island regions (A). Butterfly first-indexing: 51 DIFFERENT WPI-linked landmasses (B).",
          ha="left",fontsize=9.4)
 fig.text(.065,.116,
          "Country breadth is observed GBIF distribution, NOT native geographic range or true endemism; first indexing is NOT colonisation.",
          ha="left",fontsize=9.1,color=gray)
 fig.text(.065,.077,
          "Country breadth vs new HOSTS resource gain: ρ≈0.053 (partial rank≈−0.002). OOF first-index model gains: negative.",
          ha="left",fontsize=9.0,color=gray)
 a.out_prefix.parent.mkdir(parents=True,exist_ok=True)
 for ext in ("png","pdf","svg"):
  fig.savefig(a.out_prefix.with_suffix("."+ext),dpi=280)
 plt.close(fig)
 files={ext:{"sha256":sha(a.out_prefix.with_suffix("."+ext)),"bytes":a.out_prefix.with_suffix("."+ext).stat().st_size} for ext in ("png","pdf","svg")}
 result={
  "schema":"chocho_independent_actual_indexed_country_breadth_resource_vs_first_island_index_FigS15_v1",
  "status":"SOURCE_VERIFIED_SUPPLEMENT_FIGURE_NOT_NATIVE_BUTTERFLY_RANGE_OR_ESTABLISHMENT",
  "inputs_gitblobs":{QSOURCE:PIN_Q,REPORT:PIN_R,OOF:PIN_O},
  "butterflies":222,"resource_panel_total_new_island_opportunities":5617,
  "record_index_panel_total_later_first_index":317,
  "quartile_num_butterflies":source.species.to_list(),
  "panel_A_resource_pct":a1.tolist(),"panel_B_index_pct":a2.tolist(),
  "raw_country_breadth_vs_host_gain_rho":report["primary_correlations"]["rarefied_30_indexed_records"]["Spearman_countrybreadth_vs_island_host_resource_gains"],
  "partial_country_breadth_vs_host_gain_rank":report["primary_correlations"]["rarefied_30_indexed_records"]["partial_rank_with_native_host_footprint_catalogued_HOSTS_count_early_51_index_and_country_observation_effort"],
  "geographic_OOF_country_only_brier_gain_pct":oof["heldout_scores"]["source_frozen_5_geographic_groups"]["incremental_over_base"]["B_plus_global_pre2011_country_breadth"]["Brier_error_reduction_pct"],
  "files":files,
  "scope":"A tests potential known HOSTS resources in botanical regions. B tests GBIF FIRST INDEX in selected landmasses; neither is population arrival or true butterfly community homogenization. A and B cannot be numerically pooled."}
 (a.out_prefix.parent/(a.out_prefix.name+"_receipt.json")).write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
 print("BUTTERFLY_OBSERVED_COUNTRY_BREADTH_FIGS15="+json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__":main()
