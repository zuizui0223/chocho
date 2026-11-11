"""Rebuild the six-island case supplement from the archived #125 table and annual index."""
import argparse
import io
import subprocess
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("barriers", type=Path, help="Directory with original butterfly annual index and lng_t3 download")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / "manuscript/integrated_butterfly_resource_access_20261109"
figdir = root / "figures/integrated_island_port_butterfly_v03"
source = "ff642c8:provenance/reviewer_defenses/results/expanded_port_islands/"
def archived(name):
    return pd.read_csv(io.BytesIO(subprocess.check_output(["git","show",source+name],cwd=root)))
cases = archived("case_records.csv")
annual = pd.read_parquet(args.barriers / "global_original_species_island_year.parquet")
annual["island"] = annual.landmass_id.astype(int)
first = annual.groupby(["speciesKey","island"]).year.min()
counts = annual.query("year == 2024").groupby(["speciesKey","island"]).records.sum()
terms = archived("through_2023/terminal_join.csv").set_index("name").island
assert terms.index.is_unique
v = pd.read_csv(args.barriers / "lng_t3/LNG_tanker_voyage.csv").drop_duplicates()
v = v[v.start_date.ge("2020-01-01") & v.end_date.le("2023-12-31") & v.end_date.ge(v.start_date)].copy()
v["src"] = v.from_terminal.map(terms)
v = v[v.to_terminal.map(terms).eq(28) & v.src.isin([13,27,43,44])]
assert set(v.voyage)=={"return"} and set(v.to_terminal)=={"Sakhalin II LNG Terminal"}
assert v.groupby("src").size().to_dict()=={13:9,27:1,43:2,44:10}
names={454:"Bornholm",27:"Hokkaido",28:"Sakhalin",98:"Zealand",551:"Whidbey",951:"Santa Catalina",262:"Öland",13:"Honshu",43:"Kyushu",44:"Taiwan",169:"Fyn",50:"Vancouver Island",1149:"San Clemente"}
events=cases[cases.event.eq(1)].sort_values(["island","species"])
assert len(cases)==390 and len(events)==8 and events.island.nunique()==6
order=[454,27,28,98,551,951]
rows=[]
for island in order:
    subset=cases[cases.island.eq(island)]
    for r in events[events.island.eq(island)].itertuples():
        src=[] if pd.isna(r.linked_sources) else [int(i) for i in str(r.linked_sources).split("|")]
        def dated(i):
            yr=first[r.key,i]
            assert yr<=2018
            return f"{names[i]} ({yr})"
        assert first[r.key,island]==2024 and counts[r.key,island]==r.n_records_2024
        rows.append(dict(island=names[island],species=r.species,port_hosts_2001_2010=r.ports_by_host,
                         nearest=dated(r.nearest_island),km=round(r.km,1),
                         linked="; ".join(map(dated,src)) or "Unconfirmed",records=r.n_records_2024,
                         eligible_species=len(subset),linked_eligible=int(subset.route.sum())))
tab=pd.DataFrame(rows)
tab.to_csv(out/"TABLE_S1_six_island_cases.csv",index=False)
fig,ax=plt.subplots(figsize=(10,4.9))
fig.subplots_adjust(left=.22,right=.71,top=.83,bottom=.16)
ax.set_xlim(0,1000); ax.set_ylim(5.65,-.7)
ax.set_xlabel("Distance to nearest earlier island record (km)")
ax.set_yticks(range(6),[names[i] for i in order]); ax.tick_params(length=0,pad=8)
ax.grid(axis="x",color="#d8dee1",lw=.65);ax.set_axisbelow(True)
for spine in ("top","right","left"): ax.spines[spine].set_visible(False)
for y,i in enumerate(order):
    g=tab[tab.island.eq(names[i])];r=g.iloc[0]
    style=("#156b61","o") if i in (27,28) else ("#77519b","s") if i==454 else ("#6a7984","D")
    ax.scatter(r.km,y,color=style[0],marker=style[1],s=68,zorder=3)
    side=-1 if r.km>650 else 1
    ax.annotate(r.nearest.split(" (")[0]+f"  {r.km:.0f}",(r.km,y),xytext=(side*8,0),
                textcoords="offset points",ha="left" if side==1 else "right",va="center",fontsize=9)
    desc={454:"Ferry · other earlier island",27:"Ferry · nearest island",28:"LNG return · nearest + others"}.get(i,"Direct link unconfirmed")
    ax.text(1.03,y,desc,transform=ax.get_yaxis_transform(),va="center",fontsize=9,color=style[0])
    ax.text(-.43,y,f"{len(g)}/{r.eligible_species}",transform=ax.get_yaxis_transform(),va="center",fontsize=10)
ax.text(-.43,-.68,"New / eligible",transform=ax.get_yaxis_transform(),fontsize=9,weight="bold")
ax.text(1.03,-.68,"Service from an earlier island",transform=ax.get_yaxis_transform(),fontsize=9,weight="bold")
fig.suptitle("Earlier host records and maritime connections on six islands",fontsize=15,y=.98)
fig.savefig(figdir/"FigS1_six_island_cases.png",dpi=220)
fig.savefig(figdir/"FigS1_six_island_cases.pdf")
head="""# Table S1. Earlier host records and shipping services on six islands

The table includes every 2024 first-indexed species–island combination with a compatible port host recorded during 2001–2010. Years in parentheses are first entries in the source butterfly series beginning in 2001. Port names delimit 5-km plant searches. The regional flora catalogue classifies each plant as contemporary and outside its native range at botanical-region scale. The year of plant introduction, local feeding and historical occupancy before indexing are unresolved.

| Destination | Butterfly (2024 records) | Host candidates with recording ports, 2001–2010 | Nearest earlier island (first year); km | Earlier islands with a direct service (first year) |
|---|---|---|---|---|"""
lines=[head]
for r in tab.itertuples():
    lines.append(f"| {r.island} | *{r.species}* ({r.records}) | {r.port_hosts_2001_2010} | {r.nearest}; {r.km} | {r.linked} |")
lines+=["","## Services linked to the butterfly cases","",
"Bornholm: Køge (Zealand)–Rønne ferry, positive journeys in 2020–2023. Hokkaido: Ōarai (Honshu)–Tomakomai passenger ferry, positive passenger counts in the three fiscal years ending March 2021–2023. Ishikari Bay New Port, where the host was recorded, is a different Hokkaido port. Both ferry sources report two-way service. For Zealand, Whidbey and Santa Catalina the compiled service records yielded no link from any earlier island for the focal butterfly. This allows other services, including mainland links.","",
"Sakhalin: all qualifying inbound legs listed below were classified as return voyages to Sakhalin II LNG Terminal. Both plant-recording ports (Uglegorsk and Kholmsk) differ from that terminal. Each species uses only the earlier source islands named in its row above. Voyage totals count distinct source records and are shared across species.","",
"| Departure island | Departure terminals | Arrival years (number of return voyages) |","|---|---|---|"]
for i,g in v.groupby("src"):
    yy=g.end_date.str[:4].value_counts().sort_index()
    lines.append(f"| {names[i]} | {', '.join(sorted(set(g.from_terminal)))} | "+", ".join(f"{y} ({n})" for y,n in yy.items())+" |")
lines+=["","Distance is the WGS84 geodesic between island representatives computed from port coordinates. The nearest island is selected from 200 eligible port islands. Positive cases occurred on six of 20 islands with prior host records; 382 of 390 eligible host-positive combinations had no first entry in 2024. Candidate areas, host status, observation histories and route coverage differ among islands. Three Sakhalin species share one destination and year.","",
"[Archived risk set, cases and route compilation](https://github.com/zuizui0223/chocho/tree/ff642c8/provenance/reviewer_defenses/results/expanded_port_islands) · [LNG-T3 voyage source](https://zenodo.org/records/19571058) · [Machine-readable Table S1](TABLE_S1_six_island_cases.csv)"]
(out/"TABLE_S1_six_island_cases.md").write_text("\n".join(lines)+"\n",encoding="utf8")

