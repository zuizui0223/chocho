"""Plot the complete 250-km chronology and distance sensitivity from the archived tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/island_breeding_relay"
FIG = ROOT / "figures/integrated_island_port_butterfly_v03/FigS12_island_breeding_relay"
NAMES = {26: "Ireland", 50: "Vancouver I.", 81: "Hawaiʻi", 214: "Maui",
         402: "Molokaʻi", 464: "Isle of Man", 551: "Whidbey", 778: "Antigua",
         905: "St Croix", 1070: "St Kitts", 3692: "Sint Eustatius"}


def main():
    links = pd.read_csv(DATA / "ordered_links.csv").query("km <= 250").copy()
    plants = pd.read_csv(DATA / "compatible_port_hosts.csv")
    # One illustrative host per island and butterfly, prioritising introduced hosts.
    hosts = plants.sort_values(["status", "host_name", "port_id"]).drop_duplicates(["island", "species"])
    links = links.merge(hosts, validate="many_to_one", on=["island", "species"])
    links["middle_name"] = links.island.map(NAMES)
    links["destination_name"] = links.destination.map(NAMES)
    assert links[["middle_name", "destination_name"]].notna().all().all()
    lays = pd.read_csv(DATA / "ordered_witnesses.csv.gz").query("laying and km <= 250")
    observed = lays.groupby(["key", "island", "destination"]).stage_year.min().rename("laying_year")
    links = links.join(observed, on=["key", "island", "destination"])
    links = links.sort_values(["species", "middle_name", "km"])
    links.to_csv(DATA / "illustrated_sequences.csv", index=False)
    counts = pd.read_csv(DATA / "counts.csv").query("scope == 'larva or explicit laying'")
    fig = plt.figure(figsize=(12, 7), layout=None)
    fig.text(.045, .948, "Where a host record preceded a larva and a neighbouring island record", fontsize=19, weight="bold")
    fig.text(.045, .898, "All nine sequences within 250 km   ·   5 potential relay islands   ·   3 migrant-listed species   ·   8 distinct first records", fontsize=12)
    ax = fig.add_axes([.044, .235, .92, .60])
    ax.set_xlim(2010.7, 2024.8)
    ax.set_ylim(-1, 9.3)
    ax.set_axis_off()
    col = [.002, .19, .39, .55, .81, .975]
    def label(x, y, value, **kw):
        ax.text(x, y, value, transform=ax.get_yaxis_transform(), va="center", **kw)
    for x, t in zip(col, ["Butterfly", "Port host¹, 2001–10", "Island A", "Records by year", "Island B", "km"]):
        label(x, 8.9, t, fontsize=10, weight="bold")
    for i, row in enumerate(links.itertuples()):
        y = 7.95 - i
        if i % 2 == 0:
            ax.axhspan(y-.47, y+.47, color="#f0f4f3", zorder=0)
        label(col[0], y, row.species, fontsize=10, style="italic")
        suffix = "*" if row.status == "introduced" else ""
        label(col[1], y, row.host_name+suffix, fontsize=10, style="italic")
        label(col[2], y, row.middle_name, fontsize=11)
        label(col[4], y, row.destination_name, fontsize=11)
        label(col[5], y, str(round(row.km)), fontsize=10, ha="right")
        def pos(year):
            return .55 + (year - 2012) * .017
        start, end = pos(row.stage_year), pos(row.first_year)
        ax.plot([start,end], [y,y], transform=ax.get_yaxis_transform(), color="#a9b7b7", lw=2)
        if pd.notna(row.laying_year):
            x = pos(row.laying_year)
            ax.plot(x, y, "D", transform=ax.get_yaxis_transform(), color="#51437a", ms=5)
            label(x, y+.23, str(int(row.laying_year)), color="#51437a", fontsize=8, ha="center")
        for x, color, yr, shift in [(start,"#187364",row.stage_year,.23),(end,"#944b29",row.first_year,-.24)]:
            ax.plot(x,y,"o",transform=ax.get_yaxis_transform(),color=color,ms=5)
            label(x,y+shift,str(yr),color=color,fontsize=8,ha="center")
    fig.text(.045,.204,"Green = larva on A.   Purple diamond = observed egg-laying on A.   Brown = first indexed record on B since 2001.",fontsize=11)
    fig.text(.045,.157,"¹ One compatible host per sequence. * Introduced in the botanical region. Feeding on the named plant was not observed in this analysis.",fontsize=10)
    parts = [f"≤{r.km} km: {r.middle_islands} / {r.species} / {r.links}" for r in counts.itertuples()]
    fig.text(.045,.087,"Relay islands / species / sequences     "+"     •     ".join(parts),fontsize=10,weight="bold")
    fig.text(.045,.042,"Distances join representative island port locations. A shared year on B can appear twice when two islands qualify as A.",fontsize=10)
    fig.savefig(FIG.with_suffix(".png"),dpi=185)
    fig.savefig(FIG.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
