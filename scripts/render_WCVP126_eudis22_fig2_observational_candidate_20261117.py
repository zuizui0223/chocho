#!/usr/bin/env python3
"""Evidence-accurate 2-panel *candidate*, NOT authoritative GEB Fig.2 replacement.

Panel A: whole WCVP alien plant richness vs *registered* WPI 2023 port
index, with original 126-island Level1+area+native-richness adjusted
partial residuals and exact original slope. This is observational, not
ships, cargo, dated introductions, or larval use.

Panel B: 22 eligible original same-region EuDiS v5 matched pairs for
morphological ocean-dispersal trait; accepted native vs alien WCVP taxa.
It is a European-ascertainment-biased TRAIT pattern, not real drift.

No connection between the two panels should be portrayed as an
identified causal path: numerical tests show geometry relates to both
floras and verified historical ship OD is insufficient.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics
from pathlib import Path
import numpy as np

SOURCE_WCVP126_SHA="53b618fa51a1a418ccc6e4363f4c866c389c90cc"
ORIGINAL_126_SOURCE_BETA=0.22140536684115839
ORIGINAL_EUDIS_22_MEDIAN=2.9515

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def read(path):
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def one(root,name):
    p=list(root.rglob(name))
    if len(p)!=1:raise RuntimeError(f"Exactly one original EuDiS source file {name} expected; found {len(p)}")
    return p[0]

def residualize(x,control):
    Q,R=np.linalg.qr(control,mode="reduced")
    if np.linalg.matrix_rank(R)!=R.shape[0]:
        raise RuntimeError("Original WCVP nuisance controls not identifiable")
    return x-Q@(Q.T@x)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--flora",type=Path,required=True)
    p.add_argument("--floristic-audit",type=Path,required=True)
    p.add_argument("--original-EuDiS-artifact",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    if gitblob(a.flora)!=SOURCE_WCVP126_SHA:
        raise RuntimeError("WCVP original fixed 126 accepted botanical region source SHA mismatch")
    flora=read(a.flora)
    if len(flora)!=126 or len({x["code"] for x in flora})!=126:
        raise RuntimeError("Whole accepted flora 126 islands not found")
    original=json.loads(a.floristic_audit.read_text(encoding="utf8"))
    source=original["cohorts"]["all_126_islands"]["native_flora_adjusted"]
    expected=float(source["coefficient"])
    if abs(expected-ORIGINAL_126_SOURCE_BETA)>1e-11:
        raise RuntimeError(f"Original pre-existing WPI/all flora regression changed: {expected}")
    levels=sorted({z["level1"] for z in flora},key=int)
    if len(levels)!=9:raise RuntimeError("Original nine broad botanical regions changed")
    x=np.array([math.log2(1+int(z["wpi_port_count"])) for z in flora],dtype=float)
    y=np.array([math.log1p(int(z["introduced"])) for z in flora],dtype=float)
    Z=np.column_stack([
        np.ones(126),
        np.array([math.log(float(z["area_km2"])) for z in flora]),
        np.array([math.log1p(int(z["native"])) for z in flora]),
        *[np.array([float(z["level1"]==level) for z in flora]) for level in levels[1:]]
    ])
    xr=residualize(x,Z);yr=residualize(y,Z)
    slope=float((xr@yr)/(xr@xr))
    if abs(slope-ORIGINAL_126_SOURCE_BETA)>1e-10:
        raise RuntimeError("Original full-flora port association not recovered from source partial residuals")
    matched=read(one(a.original_EuDiS_artifact,"within_region_eudis_native_introduced_pairs.csv"))
    candidates=[r for r in matched if r["class"]=="island"]
    if len(candidates)!=22 or len({x["code"] for x in candidates})!=22:
        raise RuntimeError("Original 22 EuDiS within-region paired island comparisons changed")
    diffs=[float(z["native_minus_introduced_thalasso_pp"]) for z in candidates]
    med=float(statistics.median(diffs))
    if abs(med-ORIGINAL_EUDIS_22_MEDIAN)>0.0001:
        raise RuntimeError(f"Original EuDiS paired composition median differs: {med}")
    plus=sum(q>0 for q in diffs)
    minus=sum(q<0 for q in diffs)
    ties=sum(q==0 for q in diffs)
    if (plus,minus,ties)!=(17,4,1):
        raise RuntimeError("Source original 22 EuDiS effect directions changed")
    paired=json.loads(one(a.original_EuDiS_artifact,"paired_source_trait_composition.json").read_text())
    original_report=paired["islands_all"]
    ci=original_report["descriptive_region_bootstrap_95pp"]
    if (paired["status"]!="PAIRED_EUDIS_TRAIT_COMPOSITION_NOT_OCEAN_TRANSPORT" or
        original_report["same_regions"]!=22 or
        abs(original_report["median_native_minus_introduced_trait_pp"]-med)>1e-10):
        raise RuntimeError("Original EuDiS source pair receipt not reproduced")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(13.2,8.25),layout="constrained")
    spec=fig.add_gridspec(nrows=1,ncols=2,width_ratios=[1.05,1.30],wspace=.17)
    ax=fig.add_subplot(spec[0,0])
    ax.scatter(xr,yr,s=32,alpha=.58,edgecolors="none")
    lim=[float(xr.min()),float(xr.max())]
    ax.plot(lim,np.array(lim)*slope,linewidth=2)
    ax.axhline(0,linewidth=.7,alpha=.35)
    ax.axvline(0,linewidth=.7,alpha=.35)
    ax.set_title("A  Ports and all introduced flora",loc="left",fontweight="bold",fontsize=13)
    ax.set_xlabel("Registered port index (residualized)")
    ax.set_ylabel("ln(1 + introduced plant species) (residualized)")
    ax.text(.03,.97,f"126 island regions\nAdjusted slope = +{slope:.3f}",
        ha="left",va="top",transform=ax.transAxes,fontsize=11,
        bbox={"facecolor":"white","alpha":.85,"edgecolor":"none","pad":5})
    ax.grid(alpha=.16)
    ax2=fig.add_subplot(spec[0,1])
    candidates.sort(key=lambda r:float(r["native_minus_introduced_thalasso_pp"]))
    ya=np.arange(len(candidates))
    native=np.array([float(q["native_fraction"])*100 for q in candidates])
    alien=np.array([float(q["introduced_fraction"])*100 for q in candidates])
    for yy,pa,pb in zip(ya,native,alien):
        ax2.plot([pa,pb],[yy,yy],linewidth=1.25,alpha=.4)
    ax2.scatter(native,ya,marker="o",s=30,label="WCVP native")
    ax2.scatter(alien,ya,marker="s",s=28,label="WCVP introduced")
    ax2.set_yticks(ya,[r["code"] for r in candidates])
    ax2.set_xlim(min(0,float(np.min([*native,*alien]))-1),max(25,float(np.max([*native,*alien]))+3))
    ax2.set_xlabel("EuDiS thalassochory trait among matched plants (%)")
    ax2.set_title("B  Ocean-dispersal morphology",loc="left",fontweight="bold",fontsize=13)
    ax2.legend(loc="upper right",fontsize=9)
    ax2.grid(axis="x",alpha=.16)
    ax2.text(.98,.06,f"22 paired island regions\nNative − introduced median: +{med:.2f} pp\nDescriptive 95%: +{ci[0]:.2f} to +{ci[1]:.2f} pp",
        transform=ax2.transAxes,va="bottom",ha="right",fontsize=10,
        bbox={"facecolor":"white","alpha":.88,"edgecolor":"none","pad":4})
    fig.suptitle("Two observed island plant-geography signatures — not proven transport routes",
         fontweight="bold",fontsize=15)
    fig.supxlabel(
      "A: WPI port records (current) × WCVP alien flora (contemporary), adjusted for area, native richness & bioregion. "
      "B: EuDiS morphology subset (European-biased). Neither panel measures sea drift or ship cargo.",
      fontsize=8.9)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    fn="Fig2_candidate_two_island_plant_supply_signatures_20261117"
    for extension in (".png",".svg"):
        fig.savefig(a.out_dir/(fn+extension),dpi=210,bbox_inches="tight")
    plt.close(fig)
    result={
       "schema":"chocho_science_conservative_20261117_two_panel_island_port_trait_evidence_figure_candidate_v1",
       "status":"SOURCE_VERIFIED_FIG2_CANDIDATE_OBSERVATIONAL_ANNOTATION_NOT_MARINE_SHIPPING_MECHANISM",
       "source_WCVP126_Git_blob":SOURCE_WCVP126_SHA,
       "panel_A_original_registered_port_vs_all_introduced_WCVP_flora_adjusted_beta":slope,
       "panel_A_regions":126,"panel_A_Level1_groups":9,
       "panel_B_original_EuDiS_v5_paired_regions":22,
       "panel_B_native_minus_introduced_median_pp":med,
       "panel_B_descriptive_bootstrap_CI_pp":ci,
       "panel_B_positive_negative_tie":[plus,minus,ties],
       "source_original_EuDiS_CI":43892007427,
       "source_original_WCVP_intro_port_audit":source.get("permutations",{}),
       "output_files":[fn+".png",fn+".svg"],
       "scientific_limits":[
          "This is a candidate evidence-accurate figure, not a replacement of canonical GEB Fig2.",
          "Panel A is a contemporary partial association; cannot establish that ships introduced plants or ports caused local plant arrivals.",
          "Panel B is EuDiS European-biased morphological potential classification within 22 regions, not actual sea-current transport.",
          "No biologically observed natural-current versus artificial ship-mediated supply contrast in the two plotted response variables.",
          "The independent 126-island WCVP native/alien Jaccard cross tests and 2005 NOAA actual directional pilot remain supporting source-audited results not plotted here.",
          "The original single integrated butterfly GEB main positive resource and forecasting outcomes are untouched."
       ]
    }
    (a.out_dir/"fig2_candidate_source_and_inference_receipt.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("FIG2_CANDIDATE_TWO_SIGNALS_SOURCE_LOCKED",json.dumps({
      "status":result["status"],"A_126_port_beta":slope,
      "B_22_native_minus_introduced_pp":med,"B_ci_pp":ci,
      "B_positive_negative_tie":[plus,minus,ties]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
