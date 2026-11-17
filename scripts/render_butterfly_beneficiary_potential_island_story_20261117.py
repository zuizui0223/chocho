#!/usr/bin/env python3
"""Source-verified GEB Figure 3:
A: 239 butterfly beneficiaries based on known HOSTS native plant-resource width
B: potential consumer-pool similarity across 126 WCVP botanical island regions
C: independent GBIF-record-index similarity for 225 butterflies × 51 distinct
   original WPI-linked landmasses, with four-chain fixed record-margin null.

Figure panels have distinct sources, populations and estimands. C is indexed
GBIF records, NOT measured biological butterfly homogenization or colonization.
"""
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
 p=a.source_dir/"BENEFICIARY_AND_POTENTIAL_HOMOGENIZATION_20261117.json"
 r=Path("results/butterfly_record_index_similarity_20261117/CUMULATIVE_GBIF_BUTTERFLY_ISLAND_INDEX_SIMILARITY_20261117.json")
 m=Path("results/butterfly_record_index_similarity_20261117/MULTICHAIN_225_51_RECORD_SIMILARITY_SOURCE_NULL_20261117.json")
 d=json.loads(p.read_text())
 h=json.loads(r.read_text())
 c=json.loads(m.read_text())
 q=pd.read_csv(a.source_dir/"beneficiary_four_quartiles.csv")
 assert (d["cohort"]["butterflies"],d["cohort"]["strict_islands"],d["cohort"]["original_island_added_cells"])==(239,126,5741)
 assert len(q)==4 and q.new_island_resource_cells.sum()==5741
 assert list(map(int,q.butterflies))==[60,60,60,59]
 pobs=d["potential_island_consumer_composition"]["observed"]
 assert pobs["pairs_comparable_both_native_and_contemporary_nonempty"]==7797
 potential0,potential1=pobs["mean_potential_butterfly_Jaccard_native"],pobs["mean_potential_butterfly_Jaccard_contemporary"]
 assert abs(potential0-.1401370010070764)<1e-9 and abs(potential1-.3675355469925412)<1e-9
 recorded=h["cumulative_recording_Jaccard"]["all_51_original_landmass_pairs"]
 obs0,obs1=recorded["mean_Jaccard_2001_2010"],recorded["mean_Jaccard_by_2022"]
 assert (h["cohort"]["butterflies"],h["cohort"]["fixed_landmass_groups"],recorded["pairs"])==(225,51,1275)
 assert (h["cohort"]["earlier_indexed_cells"],h["cohort"]["new_first_indexed_2011_2022"])==(736,317)
 assert (c["chains_per_null"],c["draws_per_chain"])==(4,199)
 z=c["nulls"]["geography_preserving"]["diagnostic"]["all_pairs"]
 assert z["chain_means_stable"] and abs(z["observed_increment"]-recorded["mean_increment"])<1e-10
 expected_end=obs0+z["null_mean"]
 null_lo, null_hi=[obs0+k for k in z["null_descriptive_2p5_97p5"]]
 assert expected_end<obs1 and .17<expected_end<.172

 plt.rcParams.update({
  "font.family":"DejaVu Sans","font.size":11,
  "axes.linewidth":.7,"axes.spines.top":False,"axes.spines.right":False,
  "savefig.facecolor":"white","figure.facecolor":"white"
 })
 fig=plt.figure(figsize=(12.8,7.0))
 gs=fig.add_gridspec(2,2,width_ratios=[1.02,1.0],height_ratios=[1,1],
                     left=.085,right=.975,top=.79,bottom=.175,
                     wspace=.34,hspace=.62)
 axa=fig.add_subplot(gs[:,0]);axb=fig.add_subplot(gs[0,1]);axc=fig.add_subplot(gs[1,1])
 green="#1d7467";blue="#618697";gray="#707983";orange="#b8793c";color1="#9bb7b2"
 vals=100*q.fraction_all_island_resource_additions.to_numpy()
 axa.bar(np.arange(4),vals,color=[color1,"#82a8a0","#4d9587",green],width=.66)
 axa.axhline(25,c=gray,linestyle=(0,(3,2)),linewidth=1,alpha=.75)
 axa.set_xticks(range(4),["Q1\nNarrow","Q2","Q3","Q4\nBroad"])
 axa.set_ylim(0,47);axa.set_yticks([0,10,20,30,40])
 axa.set_ylabel("Share of new island host opportunities (%)",fontsize=10)
 for i,v in enumerate(vals):
  axa.text(i,v+1.2,f"{v:.1f}%",ha="center",va="bottom",fontsize=12,weight="bold")
 axa.set_title("A  Who gains potential resources?",loc="left",fontsize=13,weight="bold",pad=14)
 axa.text(.02,.94,"Native host-resource geography quartiles\n239 butterflies × 126 island plant regions",va="top",
          ha="left",transform=axa.transAxes,fontsize=9,color=gray)
 axa.text(-.31,26.0,"Equal-share 25%",ha="left",va="bottom",fontsize=8,color=gray)

 axb.bar([0,1],[potential0,potential1],color=[blue,orange],width=.46)
 axb.set_ylim(0,.43)
 axb.set_xticks([0,1],["Native host\nranges","Native + introduced\nhost ranges"])
 axb.set_ylabel("Potential-pool Jaccard",fontsize=9.4)
 axb.yaxis.set_major_formatter(FuncFormatter(lambda y,pos:f"{y:.2f}"))
 for i,v in enumerate((potential0,potential1)):
  axb.text(i,v+.017,f"{v:.3f}",ha="center",va="bottom",fontsize=11.6,weight="bold")
 axb.set_title("B  Potential butterfly consumer sets",loc="left",fontsize=12,weight="bold",pad=13)
 axb.text(.98,1.035,"239 butterflies · 7,797 island pairs",transform=axb.transAxes,
          ha="right",va="bottom",fontsize=8.6,color=gray)

 axc.bar([0,1,2],[obs0,expected_end,obs1],
         color=[blue,"#aeb3b4",orange],width=.52)
 axc.errorbar([1],[expected_end],yerr=[[expected_end-null_lo],[null_hi-expected_end]],
              color="#424e53",fmt="none",capsize=3,elinewidth=1.25)
 axc.set_ylim(0,.215)
 axc.set_xticks([0,1,2],["Early index\n2001–2010","2022 null\nrecord margins","Observed\nby 2022"])
 axc.set_ylabel("Indexed-butterfly Jaccard",fontsize=9.4)
 axc.yaxis.set_major_formatter(FuncFormatter(lambda y,pos:f"{y:.2f}"))
 for i,v in enumerate((obs0,expected_end,obs1)):
  axc.text(i,v+.006,f"{v:.3f}",ha="center",va="bottom",fontsize=11,weight="bold")
 axc.set_title("C  Independently indexed butterfly records",loc="left",fontsize=12,weight="bold",pad=13)
 axc.text(.98,1.035,"225 butterflies · 51 landmasses · 1,275 pairs",transform=axc.transAxes,
          ha="right",va="bottom",fontsize=8.5,color=gray)

 fig.text(.085,.958,"Introduced plants assemble shared potential resources for butterflies",
          fontsize=17,weight="bold",ha="left")
 fig.text(.085,.895,
          "Two independently measured layers: potential resources (A–B) and cumulative GBIF species indexing (C)",
          fontsize=10.3,ha="left",color=gray)
 fig.text(.085,.098,
          "A: plant-derived host breadth, not butterfly native range. B: potential resources. C: GBIF indexing, not arrival or abundance.",
          fontsize=9.2,ha="left",color="#344148")
 fig.text(.085,.060,
          "C null holds each butterfly's and island's new index counts and butterfly×5-geographic-group counts (4 chains).",
          fontsize=8.8,ha="left",color=gray)
 a.out_prefix.parent.mkdir(parents=True,exist_ok=True)
 for ext in ("png","pdf","svg"):
  fig.savefig(a.out_prefix.with_suffix("."+ext),dpi=260)
 plt.close(fig)
 outputs={ext:{"path":str(a.out_prefix.with_suffix("."+ext)),"sha256":sha(a.out_prefix.with_suffix("."+ext)),
               "bytes":a.out_prefix.with_suffix("."+ext).stat().st_size} for ext in ("png","pdf","svg")}
 receipt={
  "schema":"chocho_fig3_butterfly_beneficiary_potential_and_independent_indexed_communities_v2",
  "status":"SOURCE_VERIFIED_THREE_PANEL_POTENTIAL_VS_GBIF_RECORD_INDEX_NOT_REAL_BUTTERFLY_HOMOGENIZATION",
  "source_file":str(p),"source_sha256":sha(p),
  "record_index_source_file":str(r),"record_index_source_sha256":sha(r),
  "fourchain_null_source_file":str(m),"fourchain_null_source_sha256":sha(m),
  "values":{"resource_quartile_shares_pct":vals.tolist(),
            "potential_consumer_native_Jaccard":potential0,
            "potential_consumer_contemporary_Jaccard":potential1,"potential_pairs":7797,
            "cumulative_GBIF_index_native_Jaccard":obs0,
            "cumulative_GBIF_index_2022_Jaccard":obs1,
            "cumulative_GBIF_index_expected_2022_under_group_null":expected_end,
            "cumulative_index_pairs":1275,
            "geographic_null_extra_Jaccard":z["observed_excess"],
            "fourchain_source_stability_pass":z["chain_means_stable"]},
  "null_warning":"Panels A and B are reconstructed potential host-supported consumer resources on 126 botanical regions. C is actually indexed GBIF record membership on a DIFFERENT 51-landmass sample; expected indexed similarity preserves per butterfly/island/year pooled new indexing counts and butterfly-by-group margins, not all observer bias. Neither is observed ecological butterfly-community homogenization.",
  "output_files":outputs
 }
 (a.out_prefix.parent/(a.out_prefix.name+"_receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
 print("BUTTERFLY_GEB_THREE_PANEL_SOURCE_FIG3_VERIFIED="+json.dumps(receipt,sort_keys=True),flush=True)
if __name__=="__main__":main()
