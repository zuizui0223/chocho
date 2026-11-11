"""Render the manuscript's three-scale synthesis from the published result tables."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/integrated_island_port_butterfly_v03/Fig1_resource_geography_scales_v03"
PORTS = "provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_validation.json"
ANNUAL = ROOT / "results/annual_port_host_27"


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    groups = json.loads((ROOT / "results/island-mainland-resource-q1/island_mainland_summary.json").read_text())["groups"]
    ports = json.loads((ROOT / PORTS).read_text())
    annual = table(ANNUAL / "annual_counts.csv")
    model = next(row for row in table(ANNUAL / "rate_differences.csv") if row["label"] == "Pooled: + island + species")
    ink, quiet, green, grey = "#172f33", "#42585c", "#007866", "#78929a"
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "DejaVu Sans", "svg.fonttype": "none", "svg.hashsalt": "chocho-fig1"})
    fig = plt.figure(figsize=(8.3, 7.45), facecolor="white")
    def label(x, y, s, size=10, color=ink, bold=False):
        fig.text(x, y, s, fontsize=size, color=color, va="top",
                 weight="bold" if bold else "normal", linespacing=1.35)
    def rule(y):
        fig.add_artist(plt.Line2D([.05, .96], [y, y], transform=fig.transFigure, lw=.65, color="#b7c8c7"))
    def bar(y, value, maximum, tag, tone):
        fig.add_artist(plt.Rectangle((.58, y-.009), .22*value/maximum, .013,
                                    color=tone, transform=fig.transFigure))
        label(.82, y+.007, tag, 12, tone, True)

    label(.05,.97,"Where do introduced plants add butterfly hosts?",17,bold=True)
    label(.05,.910,"Earlier work documented alien host use in California and naturalized plant richness on islands.\nHere, each butterfly keeps its known hosts. A region gains potential resources only when\nan introduced plant supplies its first known host for that butterfly.",10.5)

    rule(.805)
    label(.05,.785,"A   Global resource geography",14,green,True)
    label(.05,.737,"239 butterflies × botanical regions\n126 insular and 237 mainland regions",10.5)
    label(.05,.656,"Comparison: native plant ranges versus\nnative plus introduced ranges.",10,quiet)
    label(.58,.779,"Gain in butterfly–region combinations",10,bold=True)
    for key, y, title, tone in [("island",.698,"Islands",green),("mainland",.625,"Mainland",grey)]:
        native = groups[key]["native_butterfly_region_units"]
        current = groups[key]["contemporary_butterfly_region_units"]
        gain = 100*(current/native-1)
        label(.58,y+.043,f"{title}   {native:,} → {current:,}",9.5)
        bar(y,gain,120,f"+{gain:.1f}%",tone)

    rule(.585)
    label(.05,.565,"B   Vegetation near ports",14,green,True)
    label(.05,.516,f"{ports['complete_ports']:,} ports with plant inventories (5 km)\n{ports['alien_candidate_species']} candidate introduced host species",10.5)
    label(.05,.436,"Regional comparison: port counts and\nintroduced flora in 98 island regions.",10,quiet)
    label(.58,.561,"Ports with a candidate host record",10,bold=True)
    hits, total = ports["ports_with_alien_candidates"], ports["complete_ports"]
    label(.58,.520,f"{hits:,} / {total:,}   ({100*hits/total:.1f}%)",15,green,True)
    label(.58,.459,"Port-rich regions had more host taxa.\nAccounting for all introduced plants\nremoved the positive host association.",9.7)

    rule(.379)
    label(.05,.357,"C   New island butterfly records",14,green,True)
    label(.05,.306,"27 islands × 132 butterflies, 2011–2024\n42,228 at-risk species–island–years",10.5)
    label(.05,.226,"Exposure: a compatible candidate host\nrecorded near a port in 2001–2010.",10,quiet)
    label(.58,.352,"Annual first-record rate on the island",10,bold=True)
    for h, y, desc, tone in [("1",.267,"Earlier host",green),("0",.194,"Unconfirmed",grey)]:
        subset = [r for r in annual if r["host"] == h]
        events = sum(int(r["first_indices"]) for r in subset)
        risk = sum(int(r["at_risk"]) for r in subset)
        rate = 100*events/risk
        label(.58,y+.043,f"{desc}   {events} / {risk:,}",9.5)
        bar(y,rate,1.2,f"{rate:.2f}%",tone)
    rule(.153)
    label(.05,.134,"Adjusted annual difference",10,green,True)
    label(.345,.134,f"+{float(model['diff_pp']):.2f} percentage points (95% interval +{float(model['lower_pp']):.2f} to +{float(model['upper_pp']):.2f})",10,bold=True)
    label(.05,.091,"Adjustment includes year, distance to an earlier same-species island, island and butterfly species.\nThe three samples share host identities. Ship transport, arrival and breeding were unobserved.",9,quiet)
    for ext in ("png","pdf","svg"):
        meta = {"Date": None} if ext == "svg" else {"CreationDate": None} if ext == "pdf" else {}
        fig.savefig(OUT.with_suffix("."+ext), dpi=210, metadata=meta)
    plt.close(fig)


if __name__ == "__main__":
    main()
