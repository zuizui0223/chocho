"""Rebuild two-sided host comparisons and assess single-group influence."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from analyze_native_plus_introduced_27 import estimate
from acquire_route_free_port_flora import OUT as SOURCE, ROOT

OUT=ROOT/"results/route_free_host_influence"
MODEL="bs(np.log(km), df=4) + C(year) + C(island) + C(key)"
# Geographic groupings of destinations; overseas islands follow their location.
GEOGRAPHY={
    "East Asia":{13:"Honshu",27:"Hokkaido",44:"Taiwan"},
    "North Atlantic and Baltic":{15:"Great Britain",24:"Iceland",26:"Ireland",168:"Gotland",232:"Skye",262:"Öland",456:"Orkney Mainland",464:"Isle of Man"},
    "Mediterranean":{85:"Cyprus",89:"Corsica",151:"Euboea",152:"Mallorca"},
    "Macaronesia":{386:"São Miguel",391:"Madeira",549:"Azores (Lajes port)",602:"Terceira"},
    "Pacific North America":{551:"Whidbey",821:"Santa Cruz (CA)",951:"Santa Catalina (CA)"},
    "Hawaii":{214:"Maui"},"Caribbean":{1592:"Saint Martin"},
    "Southwest Indian Ocean":{183:"Réunion"},"Tasmania":{32:"Tasmania"},
}


def overlap(d):
    while len(d):
        keep=(d.groupby("island").both.transform("nunique")==2)&(d.groupby("key").both.transform("nunique")==2)
        if keep.all():break
        d=d[keep]
    return d


def fit(d):
    res=estimate(d,MODEL)
    for h in (0,1):
        group=d[d.both.eq(h)]
        res[f"risk_{h}"]=len(group)
        res[f"events_{h}"]=int(group.event.sum())
    return res


def main(workspace):
    OUT.mkdir(parents=True,exist_ok=True)
    cells=pd.read_csv(SOURCE/"species_island_hosts.csv")
    all_rows=pd.read_csv(SOURCE/"risk_years.csv.gz").merge(cells[cells.native.gt(0)],on=["island","key"],validate="many_to_one").merge(pd.read_csv(SOURCE/"island_registry.csv"),on="island",validate="many_to_one")
    all_rows["both"]=all_rows.introduced.gt(0).astype(int)
    panels={"All overlap":overlap(all_rows),"Additional overlap":overlap(all_rows[~all_rows.old_destination])}
    assert panels["All overlap"].island.nunique()==26 and panels["All overlap"].key.nunique()==81
    linkfile=OUT/"introduced_links.csv"
    if workspace:
        edges=pd.read_parquet(workspace/"work/barriers/all225_named_host_edges_local.parquet").drop_duplicates()
        plants=pd.read_csv(SOURCE/"port_host_status.csv").query("status=='introduced'")[["island","host_name"]].drop_duplicates()
        links=plants.merge(edges).merge(panels["All overlap"][["island","species","key"]].drop_duplicates()).sort_values(["island","key","host_name"])
        links.to_csv(linkfile,index=False)
    links=pd.read_csv(linkfile)
    counts=links.groupby(["island","key"]).size()
    labels={i:(name,area) for area,items in GEOGRAPHY.items() for i,name in items.items()}
    assert set(panels["All overlap"].island)==set(labels)
    baselines=[];changes=[];host_support=[]
    for name,d in panels.items():
        d=d.copy()
        assert np.array_equal(d.introduced,counts.reindex(pd.MultiIndex.from_frame(d[["island","key"]])).fillna(0))
        d["island_name"]=d.island.map(lambda i:labels[i][0])
        d["area"]=d.island.map(lambda i:labels[i][1])
        b=fit(d);baselines.append(dict(panel=name,**b))
        def record(kind,label,remaining,affected,refit=True):
            kept=overlap(remaining) if refit else remaining
            result=fit(kept) if refit else b
            changes.append(dict(panel=name,kind=kind,label=label,affected_years=len(affected),affected_events=int(affected.event.sum()),
                                affected_both_events=int(affected.loc[affected.both.eq(1),"event"].sum()),change_pp=result["diff_pp"]-b["diff_pp"],**result))
        for column,kind in [("area","Region"),("island_name","Island"),("species","Butterfly")]:
            for label,g in d.groupby(column):
                record(kind,label,d[d[column].ne(label)],g)
        # Recode exposure only where this plant was the sole introduced host;
        # a different documented introduced host keeps a combination mixed.
        pool=links.merge(d[["island","key","introduced"]].drop_duplicates())
        pair_stats=d.groupby(["island","key"],as_index=False).agg(at_risk=("event","size"),first_records=("event","sum"))
        for plant,g in pool.groupby("host_name"):
            matched=g.merge(pair_stats,on=["island","key"],validate="one_to_one")
            host_support.append(dict(panel=name,plant=plant,islands=matched.island.nunique(),butterflies=matched.key.nunique(),
                                     plant_pairs=len(matched),associated_first_records=int(matched.first_records.sum()),
                                     sole_host_first_records=int(matched.loc[matched.introduced.eq(1),"first_records"].sum())))
            unique=g[g.introduced.eq(1)]
            mask=pd.MultiIndex.from_frame(d[["island","key"]]).isin(pd.MultiIndex.from_frame(unique[["island","key"]]))
            if mask.any():
                recoded=d.copy()
                recoded.loc[mask,"both"]=0
                record("Introduced host",plant,recoded,d[mask])
            else:
                record("Introduced host",plant,d,d.iloc[0:0],False)
        print(name,b,flush=True)
    pd.DataFrame(host_support).to_csv(OUT/"plant_support.csv",index=False)
    pd.DataFrame(baselines).to_csv(OUT/"panel_estimates.csv",index=False)
    pd.DataFrame(changes).to_csv(OUT/"deletions.csv",index=False)


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--workspace",type=Path)
    with threadpool_limits(limits=2):main(ap.parse_args().workspace)
