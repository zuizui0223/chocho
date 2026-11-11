"""Trace the fixed global top ten through the 27-island annual panel."""
import argparse
import io
import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy.stats import binomtest, t
from analyze_butterfly_host_contribution_concentration import load_descriptors, load_pairs, load_units

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/top_host_island_chronology"
FIG = ROOT / "figures/integrated_island_port_butterfly_v03/FigS6_named_host_chronology"
SITES = {13: ("Honshu", "JAP"), 169: ("Fyn", "DEN"), 276: ("Lolland", "DEN"),
         454: ("Bornholm", "DEN"), 551: ("Whidbey", "WAS"), 821: ("Santa Cruz (CA)", "CAL")}


def archived(name):
    ref = "ff642c8:provenance/reviewer_defenses/results/expanded_port_islands/" + name
    return pd.read_csv(io.BytesIO(subprocess.check_output(["git", "show", ref], cwd=ROOT)))


def regional_check(inputs, leaders, events):
    original = inputs / "original"
    diets, meta = load_pairs(original / "insect_host_accepted.csv")
    native = load_units(original / "native_extant_nondoubtful_wgsrpd3.csv")
    present = load_units(original / "extant_nondoubtful_wgsrpd3.csv")
    descriptors = load_descriptors(inputs / "s1_resource_descriptors.csv")
    by_species = {}
    added = lost = 0
    for sp, desc in descriptors.items():
        if not desc["host_family_count"] or not desc["native_resource_units"]:
            continue
        hosts = diets[sp]
        nat = set().union(*(native[h] for h in hosts))
        assert len(nat) == desc["native_resource_units"]
        current = set().union(*(present[h] for h in hosts))
        other = set().union(*(present[h] for h in hosts if meta[h]["accepted_name"] not in leaders))
        added += len(current - nat)
        lost += len(current - (nat | other))
        by_species[sp] = hosts
    assert len(by_species) == 239 and added == 14553
    # Regional status uses the original global reconstruction, independent of port indexing.
    for idx, r in events.iterrows():
        island, region = SITES[r.island]
        events.loc[idx, ["island_name","region"]] = [island, region]
        if r.species not in by_species:
            events.loc[idx,"global_cell"] = "Outside global panel"
            continue
        hosts = by_species[r.species]
        target = {h for h in hosts if meta[h]["accepted_name"] == r.plant}
        if not target:
            events.loc[idx,"global_cell"] = "Focal edge absent from original 239-species catalogue"
        elif not any(region in present[h] - native[h] for h in target):
            events.loc[idx,"global_cell"] = "Focal edge or introduced range unconfirmed"
        elif any(region in native[h] for h in hosts):
            events.loc[idx,"global_cell"] = "Already covered by a native host"
        else:
            events.loc[idx,"global_cell"] = "New regional resource cell"
    return added, lost


def specificity(d):
    # Adjust for other compatible recorded hosts as well as pre-event geography.
    X = np.asarray(dmatrix("other_host + bs(np.log(km), df=4) + C(year) + C(island) + C(key)", d))
    response = d[["top_host","event"]].to_numpy()
    coef, _, rank, _ = np.linalg.lstsq(X, response, rcond=None)
    v, u = (response-X@coef).T
    assert v@v > 0
    b = v@u/(v@v)
    score = pd.Series(v*(u-b*v)/(v@v))
    def meat(group):
        sums = score.groupby(np.asarray(group)).sum()
        return len(sums)/(len(sums)-1)*(sums@sums)
    variance = (meat(d.island)+meat(d.key)-meat(d.groupby(["island","key"]).ngroup()))*(len(d)-1)/(len(d)-rank-1)
    assert variance >= 0
    radius = t.ppf(.975,min(d.island.nunique(),d.key.nunique())-1)*np.sqrt(variance)
    return [100*x for x in (b,b-radius,b+radius)]


def draw(summary, events):
    fig, (a,b) = plt.subplots(1,2,figsize=(12,6.3),gridspec_kw={"width_ratios":[.94,1.46]})
    fig.subplots_adjust(left=.145,right=.99,top=.79,bottom=.19,wspace=.37)
    green = "#006e61"
    a.set_title("A   Global top ten",loc="left",weight="bold")
    a.barh(range(10),100*summary.share,color=[green if k else "#a8babd" for k in summary.events])
    a.set_yticks(range(10),summary.plant,fontstyle="italic",fontsize=9); a.invert_yaxis()
    a.set_xlabel("Share of global added combinations (%)",fontsize=9)
    a.spines[["top","right","left"]].set_visible(False)
    a.tick_params(axis="y",length=0);a.set_xlim(0,5.9)
    for i,r in enumerate(summary.itertuples()):
        a.text(100*r.share+.08,i,f"{r.islands} / {r.events}",va="center",fontsize=8)
    b.set(xlim=(0,1),ylim=(0,1)); b.axis("off")
    b.set_title("B   Earlier host at port → later island first record",loc="left",weight="bold")
    y=.965
    for plant,group in events.groupby("plant",sort=False):
        b.text(0,y,plant,weight="bold",style="italic",color=green,fontsize=11)
        y-=.057
        for island,sub in group.groupby("island_name",sort=False):
            b.text(.01,y,island,fontsize=9,weight="bold")
            b.text(.34,y,", ".join(sub.ports.iloc[0].split(" | ")).title(),fontsize=8)
            y-=.051
            for row in sub.itertuples():
                b.text(.34,y,f"{row.year}    {row.species}",fontsize=9,style="italic")
                y-=.042
            y-=.012
        y-=.013
    fig.text(.035,.95,"Which globally important hosts precede butterfly records on islands?",weight="bold",fontsize=17)
    fig.text(.035,.885,"Plant records within 5 km of ports, 2001–2010   /   Butterfly first records anywhere on the island, 2011–2024",fontsize=10)
    fig.text(.035,.075,"A: labels give islands with a prior host record / later compatible butterfly–island events (27 islands screened).",fontsize=9)
    fig.text(.035,.039,"B: catalogue host matches; port names locate the plant records. Butterfly arrival, local feeding and establishment were unobserved.",fontsize=9)
    for ext in ["png","pdf"]:
        fig.savefig(FIG.with_suffix("."+ext),dpi=195)
    plt.close(fig)


def main(workspace):
    RESULT.mkdir(exist_ok=True, parents=True)
    base = workspace / "work/barriers"
    top = pd.DataFrame(json.loads((ROOT/"provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json").read_text())["plant_species"]["top10"]).rename(columns={"accepted_name":"plant"})
    risk = pd.read_csv(ROOT/"results/annual_port_host_27/risk_years.csv.gz")
    species = archived("species_islands.csv")[["key","species"]].drop_duplicates()
    ports = archived("prior_hosts.csv").rename(columns={"host":"plant"})
    ports = ports[ports.island.isin(risk.island)]
    reg = json.loads((base/"ports_2019.geojson").read_text())["features"]
    ports["port"] = ports.port_id.map({f["properties"]["OBJECTID"]:f["properties"]["PORT_NAME"] for f in reg})
    assert ports.port.notna().all()
    edges = pd.read_parquet(base/"all225_named_host_edges_local.parquet").rename(columns={"host_name":"plant"}).drop_duplicates()
    exposure = ports.merge(edges,on="plant").merge(species,on="species")
    ex = exposure.assign(top_host=exposure.plant.isin(top.plant),other_host=~exposure.plant.isin(top.plant))
    flags = ex.groupby(["island","key"])[["top_host","other_host"]].any()
    d = risk.join(flags,on=["island","key"]); d[flags.columns] = d[flags.columns].eq(True).astype(int)
    assert (d.host == d[flags.columns].max(axis=1)).all() and len(d)==42228
    hits = ex[ex.top_host].merge(d[d.event.eq(1)][["island","key","year"]],on=["island","key"])
    events = hits.groupby(["plant","island","key","species","year"],as_index=False).agg(ports=("port",lambda x:" | ".join(sorted(set(x)))))
    added,lost = regional_check(workspace/"work/remote_sensing/host_source/extracted/inputs",set(top.plant),events)
    events["rank"] = events.plant.map(dict(zip(top.plant,range(1,11))))
    events = events.sort_values(["rank","island","year","species"])
    events.to_csv(RESULT/"events.csv",index=False)
    summary = top.copy()
    for col,series in [("islands",ports.groupby("plant").island.nunique()),("ports",ports.groupby("plant").port_id.nunique()),("events",events.groupby("plant").size())]:
        summary[col] = summary.plant.map(series).fillna(0).astype(int)
    summary.to_csv(RESULT/"ranked_hosts.csv",index=False)
    focal = pd.read_csv(ROOT/"provenance/reviewer_defenses/results/paired_island_port_nonport_coast_20261109/focal_named_butterfly_host_time_order.csv")
    focal = focal[focal.compatible_plant.eq("Medicago sativa")][["landmass_id","port_id","side","compatible_plant_recorded_2001_2010"]].drop_duplicates()
    pairs = focal.pivot(index=["landmass_id","port_id"],columns="side",values="compatible_plant_recorded_2001_2010")
    assert len(pairs)==24 and not pairs.isna().any().any()
    islands = pairs.groupby(level=0).mean()
    delta = islands.port-islands.nonport
    metrics = dict(global_added=added, lost_without_top10_introduced_ranges=lost,
                   paired_sites={str(k):int(v) for k,v in pairs.value_counts().items()},
                   matched_islands=len(delta),port_higher=int((delta>0).sum()),coast_higher=int((delta<0).sum()),
                   island_sign_p=binomtest(int((delta>0).sum()),int((delta!=0).sum())).pvalue if any(delta!=0) else None,
                   top10_diff_pp=specificity(d),
                   top10_diff_pp_other_butterfly_recorded=specificity(d[d.other_record.eq(1)]),
                   at_risk_by_top_host=d.groupby("top_host").event.agg(["size","sum"]).to_dict("index"))
    pairs.to_csv(RESULT/"medicago_matched_coasts.csv")
    (RESULT/"checks.json").write_text(json.dumps(metrics,indent=2)+"\n")
    draw(summary,events)
    print(events.to_string(index=False)); print(json.dumps(metrics,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("workspace",type=Path)
    main(parser.parse_args().workspace)
