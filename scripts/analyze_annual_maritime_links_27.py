"""Compare annual maritime links within the prior-host stratum (source snapshot ff642c8)."""
import argparse
import io
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("barriers", type=Path, help="Directory containing original species–island annual counts")
args = parser.parse_args()
out = root/"results/annual_maritime_links_27"
out.mkdir(exist_ok=True, parents=True)
def archived(name):
    path = "ff642c8:provenance/reviewer_defenses/results/expanded_port_islands/" + name
    return pd.read_csv(io.BytesIO(subprocess.check_output(["git","show",path],cwd=root)))

universe = set(archived("island_coverage.csv").island)
assert len(universe)==200
edges = archived("links.csv")
nodes = archived("endpoint_coverage.csv")
f = archived("ncfo_segments.csv")
f = pd.concat([f, f.rename(columns={"a":"b","b":"a"})], ignore_index=True)
f = f.assign(source="US BTS NCFO", kind="ferry")
edges = pd.concat([edges, f[edges.columns]]).drop_duplicates()
for col, direction in [("a","out"),("b","in")]:
    nf = f.rename(columns={col:"island"}).assign(direction=direction)
    nodes = pd.concat([nodes,nf[nodes.columns]]).drop_duplicates()
edges=edges[edges.a.isin(universe) & edges.b.isin(universe) & edges.a.ne(edges.b)]
nodes=nodes[nodes.island.isin(universe)]
inbound=nodes[nodes.direction.eq("in")].groupby(["year","island"]).source.agg(set).to_dict()
departure=nodes[nodes.direction.eq("out")].groupby(["year","source"]).island.agg(set).to_dict()
assert set(nodes.year)==set(range(2019,2024))
raw=pd.read_parquet(args.barriers/"global_original_species_island_year.parquet")
raw["island"]=raw.landmass_id.astype(int)
first=raw[raw.records.gt(0)&raw.year.between(2001,2024)&raw.island.isin(universe)].groupby(["speciesKey","island"]).year.min()
panel=pd.read_csv(root/"results/annual_port_host_27/risk_years.csv.gz")
edges_by_year={k:list(g.itertuples()) for k,g in edges.groupby(["year","b"])}
first=first.to_dict()
def classify(lag):
    rr=[]
    for r in panel[panel.year.sub(lag).between(2019,2023)].itertuples():
        service_year=r.year-lag
        registries=inbound.get((service_year,r.island),set())
        eligible=set().union(*(departure.get((service_year,s),set()) for s in registries))-{r.island}
        prior={i for i in eligible if first.get((r.key,i),9999)<service_year}
        e=edges_by_year.get((service_year,r.island),[])
        hits=[edge for edge in e if edge.a in prior]
        rr.append((r.year,r.key,r.island,r.host,r.event,r.km,service_year,bool(registries),bool(prior),
                   int(bool(hits)),"|".join(sorted({edge.source for edge in hits})),
                   "|".join(str(int(x)) for x in sorted({edge.a for edge in hits}))))
    return pd.DataFrame(rr,columns=["year","key","island","host","event","km","service_year","inbound",
                                    "same_source_frame","linked","source","linked_islands"])
# Keep eligible and unresolved observations separately. Port or species combinations
# can recur across years; matching is without replacement within an analysis.
def match(g, strata, tie_target=None, sign=0):
    rows=[]
    for _,s in g.groupby(strata):
        a=s[s.linked.eq(1)]; b=s[s.linked.eq(0)]
        if a.empty or b.empty:
            continue
        dist=abs(np.log(a.km.to_numpy())[:,None]-np.log(b.km.to_numpy()))
        ok=dist<=np.log(1.25)+1e-12
        score={"linked":a.event.to_numpy()[:,None], "unconfirmed":b.event.to_numpy(),
               "difference":a.event.to_numpy()[:,None]-b.event.to_numpy()}.get(tie_target,0)
        ii,jj=linear_sum_assignment(np.where(ok, dist+sign*1e-9*score,1e6))
        for x,y in zip(ii,jj):
            if ok[x,y]:
                rows.append((a.index[x],b.index[y],dist[x,y]))
    return pd.DataFrame(rows,columns=["a","b","log_gap"]).astype({"a":int,"b":int})
summary=[]; matches=[]
for lag in (0,1):
    d=classify(lag)
    d.to_csv(out/f"coverage_lag{lag}.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    h=d[d.host.eq(1)]
    g=h[h.same_source_frame]
    for label,sub in [("host in service years",h),("destination inbound",h[h.inbound]),
                      ("earlier island in same registry",g),("linked",g[g.linked.eq(1)]),
                      ("unconfirmed",g[g.linked.eq(0)])]:
        summary.append(dict(lag=lag,stage=label,n=len(sub),events=int(sub.event.sum()),
                            islands=sub.island.nunique(),species=sub.key.nunique(),island_years=sub.groupby(["year","island"]).ngroups))
    for label,strata in [("Year and distance",["year"]),("Same species",["year","key"]),("Same destination",["year","island"])]:
        z=match(g,strata)
        a,b=g.loc[z.a],g.loc[z.b]
        row=dict(lag=lag,match=label,n=len(z),linked_events=int(a.event.sum()),unconfirmed_events=int(b.event.sum()),
                 linked_islands=a.island.nunique(),unconfirmed_islands=b.island.nunique(),islands=len(set(a.island)|set(b.island)),
                 species=len(set(a.key)|set(b.key)))
        for target in ["linked","unconfirmed","difference"]:
            for sign,side in [(1,"min"),(-1,"max")]:
                alt=match(g,strata,target,sign)
                assert len(alt)==len(z) and abs(alt.log_gap.sum()-z.log_gap.sum())<1e-7
                x=int(g.loc[alt.a].event.sum()); y=int(g.loc[alt.b].event.sum())
                row[f"{target}_{side}"]={"linked":x,"unconfirmed":y,"difference":x-y}[target]
        matches.append(row)
        if lag==0 and label=="Year and distance":
            pairs=pd.concat([a.reset_index(drop=True).add_prefix("linked_"),b.reset_index(drop=True).add_prefix("unconfirmed_")],axis=1)
            assert (pairs.linked_year==pairs.unconfirmed_year).all()
            pairs.to_csv(out/"primary_pairs.csv",index=False)
    if not lag:
        cells=d[d.same_source_frame].groupby(["host","linked"]).agg(n=("event","size"),events=("event","sum"),islands=("island","nunique"))
        cells.to_csv(out/"host_by_link.csv")
        g.groupby(["year","linked"]).agg(n=("event","size"),events=("event","sum"),islands=("island","nunique")).to_csv(out/"annual_groups.csv")
s=pd.DataFrame(summary);m=pd.DataFrame(matches)
s.to_csv(out/"coverage_summary.csv",index=False);m.to_csv(out/"matching.csv",index=False)
print(s.to_string(index=False));print(m.to_string(index=False));print(cells)

# Display observed rates and report tied optimal contrasts in the footnote.
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.4,4.15),gridspec_kw={"width_ratios":[1.0,1.45]})
fig.subplots_adjust(left=.18,right=.97,top=.81,bottom=.28,wspace=.67)
stages=s.query("lag==0").iloc[:3]
for i,r in enumerate(stages.itertuples()):
    ax.barh(i,r.n,color=["#abbdc4","#6696a3","#275e72"][i],height=.45)
    ax.text(r.n+25,i,f"{r.n:,}  ·  {r.islands} islands",va="center",fontsize=9)
ax.set(yticks=range(3),yticklabels=["Prior host record","Inbound coverage","Earlier island in\nsame data source"],ylim=(2.6,-.7),xlim=(0,3200),xlabel="Species–island–years")
ax.set_title("A   Coverage, 2019–2023",loc="left",weight="bold")
groups=s.query("lag==0").set_index("stage")
labels=["All comparable\n15 islands","Year + distance\n24 pairs · 9 islands","Also same species\n5 pairs · 3 islands","Also same island\n5 pairs · 1 island"]
rr=[(groups.loc["linked","events"],groups.loc["linked","n"],groups.loc["unconfirmed","events"],groups.loc["unconfirmed","n"])]
for r in m.query("lag==0").itertuples():
    rr.append((r.linked_events,r.n,r.unconfirmed_events,r.n))
for i,(ea,na,eb,nb) in enumerate(rr):
    for j,(ev,n,color,label) in enumerate([(ea,na,"#166d61","Direct link"),(eb,nb,"#846682","Unconfirmed")]):
        y=i+(j-.5)*.30
        bx.scatter(100*ev/n,y,color=color,s=36,label=label if i==0 else None)
        bx.text(100*ev/n+.6,y,f"{ev}/{n}",va="center",fontsize=9,color=color)
bx.set(yticks=range(4),yticklabels=labels,ylim=(3.56,-.57),xlim=(-1,24),xlabel="First indexed in the service year (%)")
bx.set_title("B   Prior-host comparisons",loc="left",weight="bold")
bx.legend(frameon=False,loc="upper center",bbox_to_anchor=(.52,-.18),ncol=2,fontsize=9)
for axis in (ax,bx):
    axis.spines[["top","right"]].set_visible(False)
    axis.tick_params(labelsize=9)
fig.text(.065,.012,"One optimal pairing shown. Equally close choices in the 24-pair comparison give −1 to +2 excess first records with a link.",size=9)
figdir=root/"figures/integrated_island_port_butterfly_v03"
for ext in ("pdf","png"):
    fig.savefig(figdir/("FigS4_annual_maritime_links."+ext),dpi=200)
