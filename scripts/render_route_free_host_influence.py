"""Plot all single deletions and the uncertainty of the strongest attenuations."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from analyze_route_free_host_influence import OUT,ROOT


def draw():
    fits=pd.read_csv(OUT/"panel_estimates.csv").set_index("panel")
    cuts=pd.read_csv(OUT/"deletions.csv")
    fig,axes=plt.subplots(2,2,figsize=(12,7.25),sharex="row",gridspec_kw={"height_ratios":[1.22,1]})
    fig.subplots_adjust(left=.17,right=.98,top=.79,bottom=.16,hspace=.5,wspace=.65)
    kinds=["Region","Island","Butterfly","Introduced host"]
    for col,panel in enumerate(fits.index):
        original=fits.loc[panel]; subset=cuts[cuts.panel.eq(panel)]
        top,low=axes[:,col]
        color=["#096851","#83552d"][col]
        top.set_title(f"{['A','B'][col]}   {['Full comparison','Additional islands'][col]}\n{int(original.islands)} islands · {int(original.species)} spp. · {int(original.events)} records",loc="left",fontsize=12,weight="bold")
        weakest=[subset[subset.kind.eq(k)].nsmallest(1,"diff_pp").iloc[0] for k in kinds]
        for y,row in enumerate([original]+weakest):
            top.errorbar(row.diff_pp,y,xerr=[[row.diff_pp-row.lo],[row.hi-row.diff_pp]],fmt="o",capsize=3,color=color if y==0 else "#566d75")
        top.set_yticks(range(5),["All retained"]+[r.label for r in weakest],fontsize=9)
        top.set_ylim(4.6,-.6)
        low.set_title("All single-removal estimates",loc="left",size=10,weight="bold")
        low.axvline(original.diff_pp,color=color,lw=1,ls=":")
        for y,k in enumerate(kinds):
            data=subset[subset.kind.eq(k)].diff_pp
            low.plot([data.min(),data.max()],[y,y],lw=1,color=color)
            low.scatter(data,np.full(len(data),y),s=14,color=color,alpha=.32)
        low.set_yticks(range(4),kinds,fontsize=9)
        low.set_ylim(3.6,-.6)
        for ax in [top,low]:
            ax.axvline(0,ls="--",lw=.8,color="#899591")
            ax.spines[["top","right","left"]].set_visible(False)
            ax.tick_params(axis="y",length=0)
    for ax in axes[0]:ax.set_xlim(-.75,1.1)
    for ax in axes[1]:ax.set_xlim(-.06,.90);ax.set_xlabel("Annual rate difference (percentage points)",fontsize=9)
    fig.text(.035,.96,"Which islands, butterflies and introduced hosts influence the rate difference?",weight="bold",fontsize=15)
    fig.text(.035,.895,"Top: full estimate and the single removal in each category giving the smallest estimate; bars show 95% intervals.",size=10)
    fig.text(.035,.069,"Bottom: each dot is a refitted estimate; the dotted line marks the full estimate. Overlapping dots can represent several removals.",size=9)
    fig.text(.035,.032,"Both minus native only, adjusted for distance, year, island and butterfly. Host removal changes recorded status only where no other introduced host remains.",size=8.8)
    stem=ROOT/"figures/integrated_island_port_butterfly_v03/FigS10_overlap_influence"
    for ext in ("png","pdf"):fig.savefig(stem.with_suffix("."+ext),dpi=195)
    plt.close(fig)


if __name__=="__main__":draw()
