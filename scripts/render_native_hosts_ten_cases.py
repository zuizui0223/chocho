"""Summarize the archived GBIF records; one row per preselected butterfly–island event."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/native_hosts_ten_island_cases"
FIG = ROOT / "figures/integrated_island_port_butterfly_v03/FigS7_local_native_hosts"


def summarize():
    cases=pd.read_csv(ROOT/"results/top_host_island_chronology/events.csv").loc[lambda z:z.global_cell.str.startswith("Already")]
    cat=pd.read_csv(OUT/"candidate_hosts.csv")
    records=pd.read_csv(OUT/"occurrences.csv.gz")
    pool=cat.merge(records,on=["island","plant_key"])
    pool=pool[pool.year.ge(2001)&pool.year.lt(pool.first_year)]
    assert not pool.duplicated(["island","butterfly_key","gbif_id"]).any()
    rows=[]
    for case in cases.itertuples():
        evidence=pool[pool.island.eq(case.island)&pool.butterfly_key.eq(case.key)].sort_values(["year","native_host","gbif_id"])
        for name,mask in [("port",evidence.distance_m.le(5000)),("elsewhere",evidence.distance_m.gt(5000))]:
            q=evidence[mask]
            early=q[q.year.le(2010)]
            precise=q[q.uncertainty_m.gt(0)&q.uncertainty_m.le(1000)]
            row=dict(island=case.island_name,butterfly=case.species,first_year=case.year,introduced=case.plant,ports=case.ports,
                     locality=name,records=len(q),native_host_taxa=q.plant_key.nunique(),first_native_year=None,example_host=None,
                     example_gbif_id=None,example_port=None,records_2001_2010=len(early),early_precise_records=int(precise.year.le(2010).sum()),records_uncertainty_le_1km=len(precise))
            if not q.empty:
                r=q.iloc[0]
                row.update(first_native_year=int(r.year),example_host=r.native_host,example_gbif_id=int(r.gbif_id),
                           example_port=r.nearest_focal_port if name=="port" else None)
            rows.append(row)
    result=pd.DataFrame(rows)
    result.to_csv(OUT/"case_comparison.csv",index=False)
    # Plant-specific totals preserve host identities when several taxa support one cell.
    detail=pool.assign(locality=pool.distance_m.le(5000).map({True:"port",False:"elsewhere"})).groupby(
        ["island","butterfly","first_year","locality","native_host","plant_key"],as_index=False).agg(first_record_year=("year","min"),gbif_records=("gbif_id","nunique"))
    detail.to_csv(OUT/"dated_native_hosts.csv",index=False)
    return cases,result


def draw(cases,results):
    fig,ax=plt.subplots(figsize=(11.1,6.7))
    fig.subplots_adjust(left=.035,right=.995,top=.81,bottom=.15)
    ax.set(xlim=(0,1),ylim=(10,0));ax.axis("off")
    for x,label in [(0,"Island"),(.168,"Butterfly · first year"),(.464,"At focal port(s)"),(.740,"Elsewhere on island")]:
        ax.text(x,1.025,label,transform=ax.transAxes,weight="bold",size=11)
    for i,c in enumerate(cases.itertuples()):
        ax.axhline(i,color="#ccd7d8",lw=.7)
        ax.text(0,i+.50,c.island_name,va="center",fontsize=10)
        ax.text(.168,i+.50,f"{c.species}  ·  {c.year}",va="center",fontsize=9.4,style="italic")
        for x,zone in [(.464,"port"),(.740,"elsewhere")]:
            r=results[(results.island==c.island_name)&(results.butterfly==c.species)&(results.locality==zone)].iloc[0]
            confirmed=r.records>0
            ax.add_patch(Rectangle((x,i+.12),.247,.74,facecolor="#e2f0ea" if confirmed else "#f2eee5",edgecolor="none"))
            ax.text(x+.010,i+.49,(f"{r.example_host} · {int(r.first_native_year)}" if confirmed else "No native host recorded"),
                    va="center",size=8.9,color="#125745" if confirmed else "#69573d")
    fig.text(.035,.95,"Were native larval hosts already recorded locally?",size=18,weight="bold")
    at_port=results[results.locality.eq("port")]
    both=int(at_port.records_2001_2010.gt(0).sum())
    fig.text(.035,.893,f"At the same ports in 2001–2010: introduced and native hosts in {both} of {len(cases)} cases; introduced host alone recorded in {len(cases)-both}.",size=11,weight="bold",color="#185c48")
    fig.text(.035,.073,"Port = within 5 km of a port that recorded the focal introduced host. Both columns require coordinates on the same island.",size=9)
    fig.text(.035,.036,"Records span 2001 to the year before the butterfly first record. Native status is regional; a blank search cannot establish absence.",size=9)
    for ext in ("png","pdf"):
        fig.savefig(FIG.with_suffix("."+ext),dpi=200)
    plt.close(fig)


if __name__=="__main__":
    cases,result=summarize()
    draw(cases,result)
    print(result[["island","butterfly","locality","records","first_native_year","example_host","records_2001_2010","records_uncertainty_le_1km"]].to_string(index=False))
