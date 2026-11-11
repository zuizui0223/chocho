"""Estimate the port-host / later first-index association on the full archived risk set."""
import io
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.optimize import linear_sum_assignment
from scipy.stats import t
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups

root = Path(__file__).resolve().parents[1]
source = "ff642c8:provenance/reviewer_defenses/results/expanded_port_islands/species_islands.csv"
data = pd.read_csv(io.BytesIO(subprocess.check_output(["git", "show", source], cwd=root)))
out = root / "results/port_host_first_index_27"
out.mkdir(parents=True, exist_ok=True)
figdir = root / "figures/integrated_island_port_butterfly_v03"
assert (len(data), data.island.nunique(), data.species.nunique()) == (2994, 27, 132)
assert not data.duplicated(["island", "key"]).any()
assert data[["host", "km", "event"]].notna().all().all() and data.km.gt(0).all()
data["logkm"] = np.log(data.km)
data.groupby("host").agg(n=("event","size"), new=("event","sum"), median_km=("km","median"),
                         islands=("island","nunique"), species=("key","nunique")).to_csv(out/"groups.csv")
# Risk differences. Each risk-set row receives equal weight; intercepts absorb
# average recording differences among destinations and among butterfly species.
fits = []
for name, rhs in [("Crude", "host"), ("Distance", "host + bs(logkm, df=4)"),
                  ("Distance + island + species", "host + bs(logkm, df=4) + C(island) + C(key)")]:
    fit = smf.ols("event ~ " + rhs, data).fit()
    i = fit.model.exog_names.index("host")
    cov, _, _ = cov_cluster_2groups(fit, data.island, data.key)
    se = np.sqrt(cov[i,i])
    b = fit.params["host"]
    q = t.ppf(.975, min(data.island.nunique(), data.key.nunique()) - 1)
    fits.append(dict(model=name, difference_pp=100*b, lower_pp=100*(b-q*se), upper_pp=100*(b+q*se)))
effects = pd.DataFrame(fits)
effects.to_csv(out/"risk_differences.csv",index=False)

# Influence check for concentration of events at any single destination.
influence = []
for island in sorted(data.island.unique()):
    subset = data[data.island.ne(island)]
    model = smf.ols("event ~ host + bs(logkm, df=4) + C(island) + C(key)", subset).fit()
    influence.append(dict(omitted_island=island, difference_pp=100*model.params["host"]))
pd.DataFrame(influence).to_csv(out/"leave_one_island_out.csv",index=False)

# Match without outcomes. Retain maximum cardinality, then minimum total
# log-distance difference. Check outcome extremes among equally optimal ties.
def match(frame, tie=0):
    matched = []
    for _, g in frame.groupby("key"):
        a = g[g.host.eq(1)]
        b = g[g.host.eq(0)]
        if a.empty or b.empty:
            continue
        dist = abs(a.logkm.to_numpy()[:,None] - b.logkm.to_numpy())
        ok = dist <= np.log(1.25) + 1e-12
        delta = a.event.to_numpy()[:,None] - b.event.to_numpy()
        cost = np.where(ok, dist + tie * 1e-9 * delta, 1e6)
        ii,jj = linear_sum_assignment(cost)
        for x,y in zip(ii,jj):
            if ok[x,y]:
                matched.append((a.index[x],b.index[y],dist[x,y]))
    return pd.DataFrame(matched,columns=["a","b","log_gap"])
m = match(data)
comparisons = []
for sign in (0,-1,1):
    z = match(data,sign)
    assert len(z)==len(m) and np.isclose(z.log_gap.sum(),m.log_gap.sum(),rtol=0,atol=1e-8)
    aa,bb=data.loc[z.a],data.loc[z.b]
    comparisons.append(dict(tie=sign,n=len(z),host_events=int(aa.event.sum()),unconfirmed_events=int(bb.event.sum()),
                            species=aa.key.nunique(),destination_islands=len(set(aa.island)|set(bb.island))))
pd.DataFrame(comparisons).to_csv(out/"matching.csv",index=False)
pair_rows = pd.concat([data.loc[m.a,["species","island","km","event"]].reset_index(drop=True).add_prefix("host_"),
                       data.loc[m.b,["species","island","km","event"]].reset_index(drop=True).add_prefix("unconfirmed_")],axis=1)
pair_rows.to_csv(out/"matched_pairs.csv",index=False)
fig,(left,right)=plt.subplots(1,2,figsize=(10.2,3.9),gridspec_kw={"width_ratios":[1,1.22]})
fig.subplots_adjust(left=.12,right=.98,bottom=.23,top=.81,wspace=.78)
rates=[(data[data.host==1].event,data[data.host==0].event),(data.loc[m.a].event,data.loc[m.b].event)]
for j,(label,color) in enumerate([("Prior host record","#146e64"),("Unconfirmed","#75647e")]):
    for k,groups in enumerate(rates):
        vec=groups[j]
        x=k+(j-.5)*.36
        val=100*vec.mean()
        left.bar(x,val,.30,color=color,label=label if k==0 else None)
        left.text(x,val+.08,f"{vec.sum()}/{len(vec)}",ha="center",size=8)
left.set(xticks=[0,1],xticklabels=["Full risk set","Same-species pairs"],ylabel="First indexed in 2024 (%)",ylim=(0,max(3,1.25*max(100*v.mean() for g in rates for v in g))))
left.legend(loc="upper left",bbox_to_anchor=(-.04,-.14),ncol=2,fontsize=8,frameon=False)
left.set_title("A   Observed recording rates",loc="left",fontweight="bold",size=12)
right.axvline(0,color="#84949a",ls="--",lw=1)
for j,r in effects.iterrows():
    right.errorbar(r.difference_pp,j,xerr=[[r.difference_pp-r.lower_pp],[r.upper_pp-r.difference_pp]],fmt="o",color="#22566a",capsize=4)
right.set(yticks=[0,1,2],yticklabels=["Crude","Distance","Distance, island\nand species"],ylim=(2.6,-.6),xlabel="Prior host − unconfirmed (percentage points)")
right.set_title("B   Full risk set · rate differences",loc="left",fontweight="bold",size=12)
for ax in (left,right):
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=9)
for ext in ("png","pdf"):
    fig.savefig(figdir/("FigS2_port_host_first_index."+ext),dpi=200)
print(effects.to_string(index=False))
print(comparisons)
