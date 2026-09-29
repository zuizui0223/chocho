#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def jitter(label: str, width: float = 0.12):
    import hashlib
    x = int(hashlib.sha256(label.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF
    return (x - 0.5) * 2 * width


def fig1_resource_and_null(anth, matched, hostbias, prominence, concentration, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.7))
    for label, ax in zip(("a", "b", "c"), axes):
        ax.text(
            -0.10, 1.06, label,
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=12,
            va="bottom",
        )

    native = np.asarray([int(r["native_resource_units"]) for r in anth], float)
    contemporary = np.asarray(
        [int(r["contemporary_resource_units"]) for r in anth], float
    )
    axes[0].scatter(native, contemporary, s=16, alpha=0.5)
    lo = max(1, float(min(native.min(), contemporary.min())))
    hi = float(max(native.max(), contemporary.max()))
    axes[0].plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Native resource units")
    axes[0].set_ylabel("Contemporary resource units")
    expanded = sum(int(r["introduced_added_units"]) > 0 for r in anth)
    axes[0].set_title(f"{expanded}/{len(anth)} species expanded")
    axes[0].text(
        0.04, 0.95,
        "26,530 → 41,083\n+54.9%",
        transform=axes[0].transAxes,
        va="top",
        fontsize=10,
    )

    order = [
        ("1_family", "1 family"),
        ("2_families", "2 families"),
        ("3_to_5_families", "3–5 families"),
        ("6plus_families", "6+ families"),
    ]
    groups = []
    labels = []
    ns = []
    for key, label in order:
        vals = np.asarray(
            [
                float(r["log_resource_expansion"])
                for r in anth
                if r["host_breadth_stratum"] == key
            ],
            dtype=float,
        )
        groups.append(vals)
        labels.append(label)
        ns.append(len(vals))

    positions = np.arange(1, len(groups) + 1)
    axes[1].boxplot(
        groups,
        positions=positions,
        widths=0.55,
        showfliers=False,
    )
    for pos, vals in zip(positions, groups):
        offsets = np.asarray(
            [jitter(f"fig1|{pos}|{i}", 0.13) for i in range(len(vals))],
            dtype=float,
        )
        axes[1].scatter(
            np.full(len(vals), pos, dtype=float) + offsets,
            vals,
            s=12,
            alpha=0.45,
        )
    axes[1].axhline(0, linewidth=0.8)
    axes[1].set_xticks(
        positions,
        [f"{lab}\n(n={n})" for lab, n in zip(labels, ns)],
        fontsize=8,
    )
    axes[1].set_ylabel("Log proportional resource expansion")
    axes[1].set_title("Little proportional gradient across diet breadth")
    axes[1].text(
        0.03, 0.96,
        "Spearman ρ = 0.008\n"
        "bootstrap 95% CI: -0.111 to 0.128",
        transform=axes[1].transAxes,
        va="top",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "alpha": 0.85,
        },
    )

    ranks = np.asarray([1, 5, 10, 20, 50, 100], dtype=float)
    plant = concentration["plant_species"]
    genera = concentration["genera"]
    plant_cum = np.asarray(
        [100.0 * float(plant["top_k_share"][str(int(k))]) for k in ranks],
        dtype=float,
    )
    genus_cum = np.asarray(
        [100.0 * float(genera["top_k_share"][str(int(k))]) for k in ranks],
        dtype=float,
    )
    axes[2].plot(ranks, plant_cum, marker="o", label="Host species")
    axes[2].plot(ranks, genus_cum, marker="s", linestyle="--", label="Host genera")
    axes[2].axhline(50, linestyle=":", linewidth=1)
    axes[2].set_xscale("log")
    axes[2].set_xticks(ranks, [str(int(k)) for k in ranks])
    axes[2].set_ylim(0, 90)
    axes[2].set_xlabel("Top-ranked contributors")
    axes[2].set_ylabel("Cumulative share of added units (%)")
    axes[2].set_title("Added opportunity is concentrated among hosts")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].text(
        0.04, 0.95,
        "50% reached by\n38/670 species\n25/431 genera",
        transform=axes[2].transAxes,
        va="top",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "alpha": 0.85,
        },
    )

    fig.suptitle(
        "Human redistribution expands butterfly resource geography across diet-breadth classes",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, outdir, "Figure1_resource_expansion_matched_null")

def fig2_occurrence(occ, overlap, species_robustness, outdir):
    rows=[r for r in occ if int(r["observed_outside_native"])>0]
    rows.sort(key=lambda r:(float(r["fraction_outside_native_explained_by_introduced"] or 0),int(r["observed_outside_native"])))
    fig,axes=plt.subplots(1,2,figsize=(12.2,7.2),gridspec_kw={"width_ratios":[1.6,1]})
    for label,ax in zip(("a","b"),axes):
        ax.text(-0.12,1.03,label,transform=ax.transAxes,fontweight="bold",fontsize=12)
    y=np.arange(len(rows))
    frac=np.asarray([float(r["fraction_outside_native_explained_by_introduced"]) for r in rows])
    axes[0].barh(y,frac)
    axes[0].set_yticks(y,[r["species"] for r in rows],fontsize=8)
    axes[0].set_xlim(0,1.03)
    axes[0].set_xlabel("Fraction of outside-native occurrence units recovered")
    for yi,row,v in zip(y,rows,frac):
        k=int(row["outside_native_units_explained_by_introduced_host_ranges"]); n=int(row["observed_outside_native"])
        axes[0].text(min(v+0.02,0.95),yi,f"{k}/{n}",va="center",fontsize=7)
    axes[0].set_title("Recovery across 23 climate-stratified-panel species")

    loo = species_robustness["leave_one_out"]
    pyrgus = loo["pyrgus_communis_excluded"]
    pieris = loo["pieris_brassicae_excluded"]
    labels=["All 23","− Pyrgus\ncommunis","− Pieris\nbrassicae"]
    observed=[
        overlap["observed_recovered_species_x_units"],
        pyrgus["recovered_units"],
        pieris["recovered_units"],
    ]
    nulls=[
        overlap["level1_composition_matched_null"],
        {"median":pyrgus["null_median"],"q025":pyrgus["null_q025"],"q975":pyrgus["null_q975"]},
        {"median":pieris["null_median"],"q025":pieris["null_q025"],"q975":pieris["null_q975"]},
    ]
    x=np.arange(3)
    med=[z["median"] for z in nulls]; lo=[z["q025"] for z in nulls]; hi=[z["q975"] for z in nulls]
    axes[1].errorbar(x,med,yerr=[np.asarray(med)-np.asarray(lo),np.asarray(hi)-np.asarray(med)],fmt="o",capsize=5,label="Region-matched null 95%")
    axes[1].scatter(x,observed,marker="D",s=42,label="Observed")
    axes[1].set_xticks(x,labels)
    axes[1].set_ylabel("Recovered species × WGSRPD3 units")
    axes[1].set_title("Recovery survives high-leverage species removal")
    axes[1].legend(frameon=False,fontsize=8)
    axes[1].text(
        0.03,0.04,
        "Species-level mean recovery: 0.583\nnull median: 0.385; 0/99,999 exceedances",
        transform=axes[1].transAxes,fontsize=8.5
    )

    fig.suptitle("Introduced host geography aligns with butterfly occurrence beyond structural null expectation",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.95))
    save(fig,outdir,"Figure2_occurrence_validation")


def fig3_robustness(ceiling,regional,outdir):
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.6))
    for label,ax in zip(("a","b"),axes):
        ax.text(0.0,1.06,label,transform=ax.transAxes,fontweight="bold",fontsize=12,va="bottom")

    thresh=ceiling["restricted_native_breadth_thresholds"]
    xs=[100,150,168,200,250]
    ys=[thresh[str(x)]["rho_family_vs_log_expansion"] for x in xs]
    axes[0].plot(xs,ys,marker="o")
    axes[0].axhline(ceiling["rho_host_family_vs_log_expansion"],linestyle="--",linewidth=1,label="Full panel ρ=0.008")
    axes[0].axhline(0,linewidth=0.8)
    axes[0].set_xlabel("Maximum native resource breadth retained")
    axes[0].set_ylabel("Spearman ρ: host families vs log expansion")
    axes[0].set_title("Associations remain small under narrower map support")
    axes[0].legend(frameon=False,fontsize=8)

    region=regional["dominant_native_resource_region"]
    order=["Northern America","Europe","Africa","Temperate Asia","Tropical Asia","Southern America"]
    vals=[region[k]["rho"] if "rho" in region[k] else region[k].get("rho_host_family_vs_log_expansion") for k in order]
    ns=[region[k]["species"] for k in order]
    xx=np.arange(len(order))
    axes[1].bar(xx,vals)
    axes[1].axhline(0,linewidth=0.8)
    axes[1].set_xticks(xx,[f"{name}\n(n={n})" for name,n in zip(order,ns)],rotation=28,ha="right")
    axes[1].set_ylabel("Within-region Spearman ρ")
    axes[1].set_title("Regional diet-breadth associations are heterogeneous")
    fig.suptitle("Ceiling sensitivities remain weak; regional diet-breadth associations are heterogeneous",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig,outdir,"Figure3_expansion_robustness")


def fig4_climate(climate_rows,effect,outdir):
    cols=[
      ("original_score_recomputed","Primary"),
      ("same_level1_score","Same region"),
      ("distance_matched_score_250km","250 km"),
      ("distance_matched_score_500km","500 km"),
      ("distance_matched_score_1000km","1000 km"),
    ]
    values=[np.asarray([float(r[k]) for r in climate_rows if r[k]!=""],float) for k,_ in cols]
    fig,axes=plt.subplots(1,2,figsize=(11.2,4.7),gridspec_kw={"width_ratios":[1.5,1]})
    for label,ax in zip(("a","b"),axes):
        ax.text(-0.08,1.14,label,transform=ax.transAxes,fontweight="bold",fontsize=12,va="bottom")
    pos=np.arange(1,len(cols)+1)
    axes[0].boxplot(values,positions=pos,widths=0.55,showfliers=False)
    for p,vals in zip(pos,values):
        for i,v in enumerate(vals):
            axes[0].scatter(p+jitter(f"{p}|{i}",0.11),v,s=12,alpha=0.55)
        axes[0].text(p,min(1.03,max(vals)+0.035),f"{np.median(vals):.3f}",ha="center",fontsize=8)
    axes[0].axhline(0.5,linestyle="--",linewidth=1)
    axes[0].set_ylim(0.2,1.08)
    axes[0].set_xticks(pos,[label for _,label in cols],rotation=20)
    axes[0].set_ylabel("Climate-filtering score")
    axes[0].set_title("Filtering persists after geographic matching")

    r=effect["observed_partial_spearman"]; ci=effect["bootstrap_percentile_95_ci"]
    axes[1].errorbar([0],[r],yerr=[[r-ci[0]],[ci[1]-r]],fmt="o",capsize=6)
    axes[1].axhline(0,linewidth=0.8)
    mde=effect["approximate_minimum_detectable_absolute_partial_correlation"]["fisher_z_approximation"]
    axes[1].axhline(mde,linestyle=":",linewidth=1)
    axes[1].axhline(-mde,linestyle=":",linewidth=1)
    axes[1].set_xlim(-0.7,0.7); axes[1].set_xticks([0],["Host-family breadth\n(primary partial effect)"])
    axes[1].set_ylabel("Partial Spearman ρ")
    axes[1].set_ylim(-0.75,0.75)
    axes[1].set_title("Predicted negative effect is unsupported and imprecise")
    axes[1].text(0.04,0.04,f"ρ = {r:.3f}\n95% bootstrap: {ci[0]:.3f} to {ci[1]:.3f}\none-sided p = 0.2237\napprox. 80% power at |ρ|≈{mde:.3f}",transform=axes[1].transAxes,fontsize=8.5)
    fig.suptitle("Climate-associated filtering persists; host-family breadth did not detectably explain its strength",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.94))
    save(fig,outdir,"FigureS1_climate_filtering_and_precision")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--anthropogenic-csv",type=Path,required=True)
    ap.add_argument("--matched-null-json",type=Path,required=True)
    ap.add_argument("--hostbias-null-json",type=Path,required=True)
    ap.add_argument("--plant-prominence-json",type=Path,required=True)
    ap.add_argument("--host-concentration-json",type=Path,required=True)
    ap.add_argument("--occurrence-csv",type=Path,required=True)
    ap.add_argument("--occurrence-null-json",type=Path,required=True)
    ap.add_argument("--occurrence-species-robustness-json",type=Path,required=True)
    ap.add_argument("--ceiling-json",type=Path,required=True)
    ap.add_argument("--regional-json",type=Path,required=True)
    ap.add_argument("--climate-csv",type=Path,required=True)
    ap.add_argument("--climate-effect-json",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args()
    fig1_resource_and_null(read_csv(a.anthropogenic_csv),json.loads(a.matched_null_json.read_text()),json.loads(a.hostbias_null_json.read_text()),json.loads(a.plant_prominence_json.read_text()),json.loads(a.host_concentration_json.read_text()),a.output_dir)
    fig2_occurrence(read_csv(a.occurrence_csv),json.loads(a.occurrence_null_json.read_text()),json.loads(a.occurrence_species_robustness_json.read_text()),a.output_dir)
    fig3_robustness(json.loads(a.ceiling_json.read_text()),json.loads(a.regional_json.read_text()),a.output_dir)
    fig4_climate(read_csv(a.climate_csv),json.loads(a.climate_effect_json.read_text()),a.output_dir)
    print(json.dumps({"figures":4,"output_dir":str(a.output_dir)},indent=2))


if __name__=="__main__":
    main()
