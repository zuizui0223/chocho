"""Run from choch workspace with original dated port inventories."""
import io
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups

out=Path(__file__).resolve().parent
o=Path("outputs")
w=Path("work/barriers")
rev="3221899dcefa539d72a3751632577d29e4a6a2b7"
f="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
raw=subprocess.check_output(["git","-C","work/remote-sensing-runner","show",f"{rev}:{f}"])
b=pd.read_csv(io.BytesIO(raw))
assert len(b)==225*51 and not b.duplicated(["speciesKey","landmass_id"]).any()
p=pd.read_csv(o/"island_port_early_inventory_status.csv")
assert len(p)==194 and p.complete.all() and set(p.landmass_id)==set(b.landmass_id)
effort=p.groupby("landmass_id").early_plant_records.sum()
h=pd.read_parquet(w/"port_named_hosts_2001_2010.parquet")
h["landmass_id"]=pd.to_numeric(h.landmass_id,errors="raise").astype(int)
assert set(h.landmass_id)<=set(b.landmass_id)
assert h.year_end.eq(2010).all() and h.year_start.eq(2001).all()
diet=pd.read_parquet(w/"all225_named_host_edges_local.parquet").groupby("species").host_name.agg(set).to_dict()
plants=h.groupby("landmass_id").host_name.agg(set).to_dict()
alien=h[h.regional_catalogue_status.eq("contemporary_not_native_in_catalogue")].groupby("landmass_id").host_name.agg(set).to_dict()
d=pd.read_csv(o/"island_host_network_pairs.csv")
dist={(int(r.island_a),int(r.island_b)):r.port_distance_km for r in d.itertuples()}
dist.update({(j,i):v for (i,j),v in list(dist.items())})
assert len(dist)==51*50
rows=[]
for late,window in enumerate(["2011_2017","2018_2022"]):
    prior=b.pre2001_2010_indexed_records.copy()
    if late:prior+=b.future2011_2017_indexed_records
    sources=b[prior.gt(0)].groupby("species").landmass_id.agg(set).to_dict()
    previous_effort=prior.groupby(b.landmass_id).sum()
    risk=b[b["eligible_"+window].eq(1)]
    assert len(risk)==prior.eq(0).sum()
    for r in risk.itertuples():
        island=r.landmass_id
        src=sources.get(r.species,set())
        assert island not in src
        km=min((dist[island,j] for j in src),default=np.nan)
        menu=diet.get(r.species,set())
        recorded=len(menu&plants.get(island,set()))
        rows.append(dict(period=window,island=island,species=r.species,key=r.speciesKey,
                         new=getattr(r,"first_indexed_in_"+window),km=km,known_diet=bool(menu),
                         plant_effort=effort[island],prior_index=previous_effort[island],
                         host=int(recorded>0),alien=int(bool(menu&alien.get(island,set())))))
panel=pd.DataFrame(rows)
summary=[];coverage=[]; fits=[]; distance_checks=[]
for window,full in panel.groupby("period"):
    mask=full.km.notna()&full.known_diet&full.plant_effort.gt(0)&full.prior_index.gt(0)
    q=full[mask].copy()
    coverage.append(dict(period=window,eligible=len(full),events=int(full.new.sum()),source_known=int(full.km.notna().sum()),
                         used=len(q),islands=q.island.nunique(),species=q.key.nunique(),used_events=int(q.new.sum()),
                         zero_plant_islands=int(effort.eq(0).sum())))
    q["near"]=q.km.lt(500).astype(int)
    for host,g in q.groupby("host"):
        distance_checks.append(dict(period=window,host=host,n=len(g),rho=g.km.rank().corr(g.new.rank())))
    for subset,frame in [("main",q),("no_GB_IE",q[~q.island.isin([15,26])])]:
        for predictor in ["host","alien"]:
            for (near,det),g in frame.groupby(["near",predictor]):
                summary.append(dict(period=window,subset=subset,predictor=predictor,near=near,detected=det,
                                    n=len(g),events=int(g.new.sum()),rate=g.new.mean(),islands=g.island.nunique(),species=g.key.nunique()))
    # Island and species fixed effects control all time-invariant differences within each window.
    model=smf.ols("new ~ host * near + C(island) + C(key)",data=q).fit()
    cov=cov_cluster_2groups(model,q.island.to_numpy(),q.key.to_numpy())[0]
    for label,terms in [("host_difference_far",["host"]),("host_difference_near",["host","host:near"]),("interaction",["host:near"])]:
        ix=model.params.index.get_indexer(terms)
        delta=model.params.iloc[ix].sum()
        var=cov[np.ix_(ix,ix)].sum()
        se=np.sqrt(var) if var>=0 else np.nan
        fits.append(dict(period=window,contrast=label,difference=delta,low=delta-1.96*se,high=delta+1.96*se))
m=pd.DataFrame(summary)
m.to_csv(out/"recording_rates.csv",index=False)
pd.DataFrame(coverage).to_csv(out/"coverage.csv",index=False)
pd.DataFrame(fits).to_csv(out/"adjusted.csv",index=False)
pd.DataFrame(distance_checks).to_csv(out/"continuous_distance.csv",index=False)
plt.rcParams["font.family"]="Yu Gothic"
fig,axes=plt.subplots(1,2,figsize=(9,4.8),sharey=True,layout="constrained")
for ax,(window,sub) in zip(axes,m.query("subset=='main' and predictor=='host'").groupby("period")):
    for det,col,name,dx in [(0,"#9099a3","食草の先行記録なし",-.19),(1,"#25806b","食草の先行記録あり",.19)]:
        vals=sub[sub.detected.eq(det)].set_index("near").loc[[1,0]]
        bars=ax.bar(np.arange(2)+dx,vals.rate*100,width=.36,color=col,label=name)
        ax.bar_label(bars,labels=[f"{x.rate:.1%}\n{x.events}/{x.n}" for x in vals.itertuples()],padding=3,fontsize=9)
    ax.set_xticks([0,1],["500 km未満","500 km以上"])
    ax.set_title(window.replace("_","–")+"年")
    ax.set_xlabel("以前の同種記録がある最寄りの島まで")
    ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("島で新しく記録された割合（%）")
axes[0].set_ylim(0,m.query("subset=='main' and predictor=='host'").rate.max()*138)
fig.legend(*axes[0].get_legend_handles_labels(),loc="outside lower center",ncol=2,frameon=False)
fig.suptitle("港の食草記録（2001–2010年）と島間距離による比較")
fig.savefig(out/"host_access.png",dpi=170)
print(pd.DataFrame(coverage).to_string(index=False))
print(m.query("subset=='main'").round(4).to_string(index=False))
print(pd.DataFrame(fits).round(4).to_string(index=False))
