"""Reapply archived same-year distance matching after expanding regional service sources."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def expanded_endpoints(root,edges,nodes):
    add=pd.read_csv(root/"results/regional_maritime_links/positive_routes.csv")
    add["kind"]="regional service"
    legs=pd.concat([add,add[add.two_way].rename(columns={"a":"b","b":"a"})])
    new_nodes=[legs.rename(columns={c:"island"}).assign(direction=d)[nodes.columns] for c,d in [("a","out"),("b","in")]]
    return pd.concat([edges,legs[edges.columns]]).drop_duplicates(),pd.concat([nodes,*new_nodes]).drop_duplicates()

def compare(root,classify,match,archived):
    out=root/"results/regional_maritime_links"
    dk={98,169,276,454,715,893,1577}
    old=pd.read_csv(root/"results/annual_maritime_links_27/coverage_lag0.csv.gz")
    new=classify(0)
    outside=classify(0,excluded=dk)
    assert len(old)==len(new)
    both=old.merge(new,on=["island","year","key"],validate="one_to_one")
    assert (both.linked_x<=both.linked_y).all() and (both.same_source_frame_x<=both.same_source_frame_y).all()
    summary=[]; pairs=[]; groups=[]
    for label,d in [("Previous",old),("Expanded",new),("Outside Denmark",outside)]:
        g=d[d.same_source_frame & d.host.eq(1)]
        for constraint,frame,strata in [("Year and distance",g,["year"]),("Both within 250 km",g[g.km.le(250)],["year"]),
                                         ("Same species",g,["year","key"]),("Same island",g,["year","island"])]:
            z=match(frame,strata)
            a,b=frame.loc[z.a],frame.loc[z.b]
            row=dict(scope=label,restriction=constraint,pairs=len(z),islands=len(set(a.island)|set(b.island)),
                     linked_islands=a.island.nunique(),control_islands=b.island.nunique(),species=len(set(a.key)|set(b.key)),
                     linked_new=int(a.event.sum()),unconfirmed_new=int(b.event.sum()))
            for target in ["linked","unconfirmed","difference"]:
                for sign,side in [(1,"min"),(-1,"max")]:
                    v=match(frame,strata,target,sign)
                    assert len(v)==len(z) and abs(v.log_gap.sum()-z.log_gap.sum())<1e-7
                    x,y=int(frame.loc[v.a].event.sum()),int(frame.loc[v.b].event.sum())
                    row[f"{target}_{side}"]={"linked":x,"unconfirmed":y,"difference":x-y}[target]
            summary.append(row)
            if constraint in ("Year and distance","Both within 250 km"):
                pair=pd.concat([a.reset_index(drop=True).add_prefix("linked_"),b.reset_index(drop=True).add_prefix("unconfirmed_")],axis=1)
                pairs.append(pair.assign(scope=label,restriction=constraint))
        by=g.groupby("linked").agg(n=("event","size"),new=("event","sum"),islands=("island","nunique"),species=("key","nunique")).reset_index()
        groups.append(by.assign(scope=label,all_islands=g.island.nunique(),all_species=g.key.nunique()))
        if label!="Previous":
            g.to_csv(out/("eligible_"+label.lower().replace(" ","_")+".csv"),index=False)
    m=pd.DataFrame(summary);m.to_csv(out/"matched_summary.csv",index=False)
    assert m.iloc[0].pairs==24
    pd.concat(pairs).to_csv(out/"matched_pairs.csv",index=False)
    totals=pd.concat(groups);totals.to_csv(out/"crude_rates.csv",index=False)
    print(m[["scope","restriction","pairs","islands","linked_new","unconfirmed_new","difference_min","difference_max"]].to_string(index=False))
    print(totals.to_string(index=False))
    draw(root,m)

def draw(root,m):
    fig,(ax,bx)=plt.subplots(1,2,figsize=(11.3,4.7),gridspec_kw={"width_ratios":[1,1.3]})
    fig.subplots_adjust(left=.245,right=.96,wspace=.28,top=.84,bottom=.22)
    pick=m[(m.restriction=="Year and distance")|((m.scope!="Previous") & (m.restriction=="Both within 250 km"))].reset_index(drop=True)
    labels=["Previous sources","Expanded sources","Expanded · both ≤250 km","Outside Denmark","Outside Denmark · both ≤250 km"]
    for i,r in pick.iterrows():
        ax.barh(i,r.pairs,height=.50,color="#47818a")
        ax.text(r.pairs+.4,i,f"{r.pairs} pairs · {r.islands} islands",va="center",fontsize=9)
        if not r.pairs:
            bx.text(.2,i,"No comparable pairs",va="center",color="#606c73")
        for j,(ev,col,lab) in enumerate([(r.linked_new,"#156e5c","Direct link"),(r.unconfirmed_new,"#815c7e","Unconfirmed")]):
            if not r.pairs:continue
            val=100*ev/r.pairs
            bx.scatter(val,i+(j-.5)*.34,color=col,s=45,label=lab if i==0 else None)
            bx.annotate(f"{ev}/{r.pairs}",(val,i+(j-.5)*.34),xytext=(6,0),textcoords="offset points",va="center",fontsize=9)
    ax.set(yticks=range(len(pick)),yticklabels=labels,xlim=(0,pick.pairs.max()*1.8),xlabel="Matched pairs")
    bx.set(yticks=[],xlabel="New records per 100 at-risk combinations (%)")
    bx.set_xlim(-.6,max(100*pick.linked_new.div(pick.pairs),default=1)+max(100*pick.unconfirmed_new.div(pick.pairs),default=1)+3)
    for x in (ax,bx):
        x.set_ylim(len(pick)-.45,-.6);x.spines[["top","right"]].set_visible(False)
    ax.set_title("A   Geographic support",loc="left",weight="bold");bx.set_title("B   First records in the service year",loc="left",weight="bold")
    bx.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.17),fontsize=9)
    for ext in ["png","pdf"]:
        fig.savefig(root/"figures/integrated_island_port_butterfly_v03"/("FigS5_regional_routes."+ext),dpi=200)
