"""Draw the paired host-report contrast and held-island prediction contrasts."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/integrated_island_port_butterfly_v03"
data = json.loads((ROOT / "results/island_resource_layers/figure_inputs.json").read_text())
green, blue = "#08775f", "#456680"


def finish(fig, stem):
    for ext in ["pdf", "png"]:
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


fig, ax = plt.subplots(figsize=(8, 3.4))
for offset, (name, group) in zip([.12, -.12], data["hosts"].items()):
    for y, key in enumerate(["difference_family_excess", "difference_family_expectation", "difference_raw"]):
        r = group["paired_comparison"][key]
        value = 100*r["mean_paired_difference"]
        lo, hi = [100*x for x in r["bootstrap95_species_level_mean_descriptive"]]
        label = "All hosts" if offset > 0 else "Excluding listed crop species"
        ax.errorbar(value, y+offset, xerr=[[value-lo],[hi-value]], fmt="o", capsize=3,
                    color=green if offset>0 else blue, label=f"{label} (n={group['shared_species_island_and_mainland']})" if y==0 else None)
ax.axvline(0,color=".65",lw=1)
ax.set(yticks=[0,1,2],yticklabels=["Observed minus expected","Expected from host families","Observed host reports"],xlabel="Island minus mainland introduced-host share (percentage points)")
ax.legend(loc="upper center",bbox_to_anchor=(.48,-.23),ncol=2,frameon=False)
ax.set_title("Introduced plants in the diets reported for the same butterflies",loc="left")
ax.spines[["top","right"]].set_visible(False)
finish(fig,"Fig3_host_reports")


ledger = json.loads((ROOT / "results/port_weed_current_review_20261117/ONE_PAPER_CIRCUIT_RESULTS_LEDGER_20261117.json").read_text())
fig, ax = plt.subplots(figsize=(8, 3))
for y, (entry, interval, color) in enumerate([
    (ledger["floras"]["native"], "circuit_bootstrap95_descriptive", blue),
    (ledger["floras"]["introduced"], "circuit_bootstrap95_descriptive", green),
    (ledger["butterflies"], "bootstrap20k_descriptive95", "#a6602d"),
]):
    value = entry["circuit_pct"]
    lo, hi = entry[interval]
    ax.errorbar(value, y, xerr=[[value-lo],[hi-value]], fmt="o", color=color, capsize=4)
    ax.text(hi+.18,y,f"{value:+.2f}%",va="center")
ax.axvline(0,color=".45",lw=.8)
ax.set(yticks=range(3),yticklabels=["Native flora · 104 regions","Introduced flora · 104 regions","Butterfly records · 51 landmasses"],xlim=(-1,10.3),ylim=(2.5,-.5),xlabel="Change in held-island prediction error (% reduction)")
ax.spines[["top","right"]].set_visible(False)
ax.set_title("Information added by undirected sea resistance",loc="left")
finish(fig,"Fig4_sea_connectivity")
